#!/usr/bin/env python3

"""
Author: Samuel Aroney
Wrap an arbitrary long-running command (typically a Snakemake invocation) so
it runs detached, and record it in a shared TSV ledger: host, PID, cwd,
command, log, and exit code. `status` reports running/completed/failed across
every node that launched into the same ledger. `restart` only re-runs a
tracked command once its previously submitted cluster job IDs (scraped from
the log) are confirmed gone from qstat, then unlocks and relaunches.
"""

import argparse
import csv
import fcntl
import logging
import os
import re
import shlex
import socket
import subprocess
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

JOBID_RE = re.compile(r"external jobid '([^']+)'")
LEDGER_FIELDS = ["name", "host", "pid", "cwd", "command", "log_path", "exit_path", "launched_at", "action"]


def default_ledger_path():
    state_home = os.environ.get("XDG_STATE_HOME") or str(Path.home() / ".local" / "state")
    return Path(state_home) / "pipeline-tracker" / "ledger.tsv"


@contextmanager
def locked_ledger(ledger_path):
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    with open(ledger_path, "a+") as f:
        # held across the whole check-then-act below, so two concurrent launches/restarts
        # of the same name (e.g. from different nodes) can't both pass the "not running" check
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            yield f
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def read_records(f):
    f.seek(0)
    return list(csv.DictReader(f, delimiter="\t"))


def write_record(f, record):
    f.seek(0, os.SEEK_END)
    writer = csv.DictWriter(f, fieldnames=LEDGER_FIELDS, delimiter="\t")
    if f.tell() == 0:
        writer.writeheader()
    writer.writerow(record)
    f.flush()


def latest_by_name(records):
    latest = {}
    for record in records:
        latest[record["name"]] = record
    return latest


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


def build_unlock_command(command):
    for i, token in enumerate(command):
        if token == "snakemake" or token.endswith("/snakemake"):
            return command[: i + 1] + ["--unlock"] + command[i + 1 :]
    return None


def check_jobs_cleared(job_ids):
    still_present = []
    for job_id in job_ids:
        result = subprocess.run(["qstat", job_id], capture_output=True)
        if result.returncode == 0:
            still_present.append(job_id)
    return still_present


def start_process(name, command, ledger_path, action, cwd):
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    log_dir = ledger_path.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{name}.{timestamp}.log"
    exit_path = log_dir / f"{name}.{timestamp}.exit"

    quoted_cmd = " ".join(shlex.quote(part) for part in command)
    # detached shell wrapper: our own process exits right after Popen, so the
    # exit code has to be captured by something that outlives it
    shell_command = (
        f"cd {shlex.quote(cwd)} && "
        f"{{ {quoted_cmd} ; }} > {shlex.quote(str(log_path))} 2>&1; "
        f"echo $? > {shlex.quote(str(exit_path))}"
    )
    process = subprocess.Popen(
        ["bash", "-c", shell_command],
        stdin=subprocess.DEVNULL,
        start_new_session=True,
    )

    return {
        "name": name,
        "host": socket.gethostname(),
        "pid": process.pid,
        "cwd": cwd,
        "command": shlex.join(command),
        "log_path": str(log_path),
        "exit_path": str(exit_path),
        "launched_at": datetime.now(timezone.utc).isoformat(),
        "action": action,
    }


def cmd_launch(args):
    with locked_ledger(args.ledger) as f:
        previous = latest_by_name(read_records(f)).get(args.name)
        if previous is not None:
            status = derive_status(read_exit_code(Path(previous["exit_path"])), pid_alive(previous["pid"]))
            if status == "running":
                logging.error(f"'{args.name}' is already running (pid={previous['pid']}, host={previous['host']}); not launching")
                return 1
        record = start_process(args.name, args.command, args.ledger, "launch", os.getcwd())
        write_record(f, record)
    logging.info(f"launch '{args.name}': pid={record['pid']} log={record['log_path']}")
    return 0


def cmd_status(args):
    latest = latest_by_name(read_records_from_path(args.ledger))
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
        rows.append([name, status, record["host"], record["pid"], record["launched_at"], record["log_path"]])

    all_rows = [header] + rows
    widths = [max(len(r[i]) for r in all_rows) for i in range(len(header))]
    for row in all_rows:
        print("  ".join(c.ljust(w) for c, w in zip(row, widths)))
    return 0


def read_records_from_path(ledger_path):
    if not ledger_path.exists():
        return []
    with open(ledger_path) as f:
        fcntl.flock(f, fcntl.LOCK_SH)
        try:
            return read_records(f)
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def cmd_restart(args):
    with locked_ledger(args.ledger) as f:
        record = latest_by_name(read_records(f)).get(args.name)
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

        command = shlex.split(record["command"])
        unlock_command = build_unlock_command(command)
        if unlock_command is not None:
            logging.info(f"unlocking '{args.name}': {' '.join(unlock_command)}")
            subprocess.run(unlock_command, cwd=record["cwd"], check=True)
        else:
            logging.warning(f"'{args.name}' command has no 'snakemake' token; skipping --unlock")

        new_record = start_process(args.name, command, args.ledger, "restart", record["cwd"])
        write_record(f, new_record)
    logging.info(f"restart '{args.name}': pid={new_record['pid']} log={new_record['log_path']}")
    return 0


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")
    parser.add_argument("--ledger", type=Path, default=default_ledger_path(), help="ledger TSV path (default: %(default)s)")

    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    launch_parser = subparsers.add_parser("launch", help="run a command detached and track it, unless it's already running under this name")
    launch_parser.add_argument("--name", required=True, help="unique tracking name for this run")
    launch_parser.add_argument("command", nargs="+", help="command to run, e.g. -- pixi run --as-is snakemake -s foo.smk all --profile aqua")

    status_parser = subparsers.add_parser("status", help="show the latest known status of tracked runs")
    status_parser.add_argument("names", nargs="*", help="names to show (default: all tracked)")
    status_parser.add_argument("--check-jobs", action="store_true", help="also query qstat for running runs' cluster job IDs")

    restart_parser = subparsers.add_parser("restart", help="restart a tracked run, after confirming its previous cluster jobs are gone")
    restart_parser.add_argument("name", help="name to restart")

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
