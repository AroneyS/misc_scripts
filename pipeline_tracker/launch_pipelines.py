#!/usr/bin/env python3

"""
Author: Samuel Aroney
Launch Snakemake pipelines from a TSV manifest and track them in a shared
ledger (host, PID, log, submitted cluster job IDs), so `status` can report
what's running/failed/complete across nodes, and `restart` only re-runs a
pipeline once its previously submitted job IDs are confirmed gone from qstat.
"""

import argparse
import fcntl
import json
import logging
import os
import re
import shlex
import socket
import subprocess
import sys
import csv
from datetime import datetime, timezone
from pathlib import Path

JOBID_RE = re.compile(r"external jobid '([^']+)'")
REQUIRED_MANIFEST_COLUMNS = {"name", "workdir"}


def default_ledger_path():
    state_home = os.environ.get("XDG_STATE_HOME") or str(Path.home() / ".local" / "state")
    return Path(state_home) / "pipeline-tracker" / "ledger.jsonl"


def parse_manifest(lines):
    content_lines = [line for line in lines if line.strip() and not line.lstrip().startswith("#")]
    if not content_lines:
        raise ValueError("manifest has no entries")

    reader = csv.DictReader(content_lines, delimiter="\t")
    missing = REQUIRED_MANIFEST_COLUMNS - set(reader.fieldnames or [])
    if missing:
        raise ValueError(f"manifest missing required column(s): {', '.join(sorted(missing))}")

    entries = []
    seen_names = set()
    for row in reader:
        name = (row.get("name") or "").strip()
        workdir = (row.get("workdir") or "").strip()
        if not name or not workdir:
            raise ValueError(f"manifest row missing name/workdir: {row}")
        if name in seen_names:
            raise ValueError(f"duplicate pipeline name in manifest: {name}")
        seen_names.add(name)
        entries.append({
            "name": name,
            "workdir": workdir,
            "profile": (row.get("profile") or "").strip(),
            "extra_args": (row.get("extra_args") or "").strip(),
        })
    return entries


def job_ids_from_log(text):
    seen = []
    for match in JOBID_RE.finditer(text):
        job_id = match.group(1)
        if job_id not in seen:
            seen.append(job_id)
    return seen


def pid_alive(pid):
    return Path(f"/proc/{pid}").exists()


def read_exit_code(exit_path):
    if not exit_path.exists():
        return None
    content = exit_path.read_text().strip()
    return int(content) if content else None


def derive_status(exit_code, alive):
    if exit_code is not None:
        return "completed" if exit_code == 0 else f"failed (exit {exit_code})"
    if alive:
        return "running"
    return "unknown (process gone, no exit code recorded)"


def build_snakemake_command(entry, unlock=False):
    cmd = ["snakemake", "--unlock"] if unlock else ["snakemake"]
    if entry["profile"]:
        cmd += ["--profile", entry["profile"]]
    if entry["extra_args"]:
        cmd += shlex.split(entry["extra_args"])
    return cmd


def append_ledger(record, ledger_path):
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(record) + "\n"
    with open(ledger_path, "a") as f:
        # lustre is mounted with flock support; this guards concurrent launches from different nodes
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            f.write(line)
            f.flush()
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def read_ledger(ledger_path):
    if not ledger_path.exists():
        return []
    with open(ledger_path) as f:
        return [json.loads(line) for line in f if line.strip()]


def latest_by_name(records):
    latest = {}
    for record in records:
        latest[record["name"]] = record
    return latest


def check_jobs_cleared(job_ids):
    still_present = []
    for job_id in job_ids:
        result = subprocess.run(["qstat", job_id], capture_output=True)
        if result.returncode == 0:
            still_present.append(job_id)
    return still_present


def launch_one(entry, manifest_path, ledger_path, action="launch"):
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    log_dir = ledger_path.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{entry['name']}.{timestamp}.log"
    exit_path = log_dir / f"{entry['name']}.{timestamp}.exit"

    quoted_cmd = " ".join(shlex.quote(part) for part in build_snakemake_command(entry))
    # detached shell wrapper: our own process exits right after Popen, so the
    # exit code has to be captured by something that outlives it
    shell_command = (
        f"cd {shlex.quote(entry['workdir'])} && "
        f"{{ {quoted_cmd} ; }} > {shlex.quote(str(log_path))} 2>&1; "
        f"echo $? > {shlex.quote(str(exit_path))}"
    )
    process = subprocess.Popen(
        ["bash", "-c", shell_command],
        stdin=subprocess.DEVNULL,
        start_new_session=True,
    )

    record = {
        "name": entry["name"],
        "workdir": entry["workdir"],
        "profile": entry["profile"],
        "extra_args": entry["extra_args"],
        "manifest_path": str(manifest_path) if manifest_path else None,
        "host": socket.gethostname(),
        "pid": process.pid,
        "log_path": str(log_path),
        "exit_path": str(exit_path),
        "launched_at": datetime.now(timezone.utc).isoformat(),
        "action": action,
    }
    append_ledger(record, ledger_path)
    logging.info(f"{action} '{entry['name']}': pid={process.pid} log={log_path}")
    return record


def cmd_launch(args):
    try:
        entries = parse_manifest(args.manifest.read_text().splitlines())
    except ValueError as e:
        logging.error(f"{args.manifest}: {e}")
        return 1

    ledger_latest = latest_by_name(read_ledger(args.ledger))
    for entry in entries:
        previous = ledger_latest.get(entry["name"])
        if previous is not None:
            status = derive_status(read_exit_code(Path(previous["exit_path"])), pid_alive(previous["pid"]))
            if status == "running":
                logging.warning(f"skipping '{entry['name']}': already running (pid={previous['pid']}, host={previous['host']})")
                continue
        launch_one(entry, args.manifest, args.ledger)
    return 0


def cmd_status(args):
    latest = latest_by_name(read_ledger(args.ledger))
    names = args.names or sorted(latest.keys())

    header = ["NAME", "STATUS", "HOST", "PID", "LAUNCHED_AT", "LOG"]
    rows = []
    for name in names:
        record = latest.get(name)
        if record is None:
            rows.append([name, "unknown (not in ledger)", "-", "-", "-", "-"])
            continue
        status = derive_status(read_exit_code(Path(record["exit_path"])), pid_alive(record["pid"]))
        if args.check_jobs and status == "running":
            job_ids = job_ids_from_log(Path(record["log_path"]).read_text())
            still_present = check_jobs_cleared(job_ids)
            status += f" ({len(still_present)}/{len(job_ids)} cluster jobs active)"
        rows.append([name, status, record["host"], str(record["pid"]), record["launched_at"], record["log_path"]])

    all_rows = [header] + rows
    widths = [max(len(r[i]) for r in all_rows) for i in range(len(header))]
    for row in all_rows:
        print("  ".join(c.ljust(w) for c, w in zip(row, widths)))
    return 0


def cmd_restart(args):
    record = latest_by_name(read_ledger(args.ledger)).get(args.name)
    if record is None:
        logging.error(f"'{args.name}' not found in ledger")
        return 1

    status = derive_status(read_exit_code(Path(record["exit_path"])), pid_alive(record["pid"]))
    if status == "running":
        logging.error(f"'{args.name}' is still running (pid={record['pid']}, host={record['host']}); not restarting")
        return 1

    job_ids = job_ids_from_log(Path(record["log_path"]).read_text())
    still_present = check_jobs_cleared(job_ids)
    if still_present:
        logging.error(
            f"'{args.name}' still has cluster jobs queued/running: {', '.join(still_present)}; "
            "wait for them to finish or qdel them before restarting"
        )
        return 1

    entry = {
        "name": record["name"],
        "workdir": record["workdir"],
        "profile": record["profile"],
        "extra_args": record["extra_args"],
    }
    unlock_cmd = build_snakemake_command(entry, unlock=True)
    logging.info(f"unlocking '{args.name}': {' '.join(unlock_cmd)}")
    subprocess.run(unlock_cmd, cwd=entry["workdir"], check=True)

    launch_one(entry, record.get("manifest_path"), args.ledger, action="restart")
    return 0


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")
    parser.add_argument("--ledger", type=Path, default=default_ledger_path(), help="ledger JSONL path (default: %(default)s)")

    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    launch_parser = subparsers.add_parser("launch", help="launch every pipeline in a manifest that isn't already running")
    launch_parser.add_argument("manifest", type=Path, help="TSV manifest: name, workdir, profile, extra_args")

    status_parser = subparsers.add_parser("status", help="show the latest known status of tracked pipelines")
    status_parser.add_argument("names", nargs="*", help="pipeline names to show (default: all tracked)")
    status_parser.add_argument("--check-jobs", action="store_true", help="also query qstat for running pipelines' cluster job IDs")

    restart_parser = subparsers.add_parser("restart", help="restart a tracked pipeline, after confirming its previous cluster jobs are gone")
    restart_parser.add_argument("name", help="pipeline name to restart")

    args = parser.parse_args(arguments)

    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%Y/%m/%d %I:%M:%S %p")

    if args.subcommand == "launch":
        return cmd_launch(args)
    elif args.subcommand == "status":
        return cmd_status(args)
    elif args.subcommand == "restart":
        return cmd_restart(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
