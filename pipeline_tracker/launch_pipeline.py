#!/usr/bin/env python3

"""
Author: Samuel Aroney
Wrap an arbitrary long-running command (typically a Snakemake invocation) so
it runs detached, and record it in a shared TSV ledger: host, PID, cwd,
command, log, and exit code. `status` reports running/completed/failed (plus
% steps done, parsed from the log) across every node that launched into the
same ledger -- for a run on another host, liveness is checked over SSH, since
PIDs are only meaningful on the host that owns them. `check` reports whether
a tracked run still has cluster jobs queued/running (e.g. after it was
cancelled or crashed), without touching anything. `cancel` qdels any such
active cluster jobs and stops the local or remote process. `restart` does the
same check `cancel` would report and, only once it's clear, unlocks and
relaunches the tracked command.
"""

import argparse
import csv
import fcntl
import logging
import os
import re
import shlex
import signal
import socket
import subprocess
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

JOBID_RE = re.compile(r"external jobid '([^']+)'")
PROGRESS_RE = re.compile(r"(\d+) of (\d+) steps \(([\d.]+)%\) done")
LEDGER_FIELDS = ["name", "host", "pid", "cwd", "command", "log_path", "exit_path", "launched_at", "action"]
SSH_OPTS = ["-o", "BatchMode=yes", "-o", "ConnectTimeout=5", "-o", "StrictHostKeyChecking=accept-new"]
# a run is only safe to relaunch/unlock/treat-as-not-running once we're sure it isn't still going;
# "unreachable" (couldn't SSH to its host) must NOT count as safe, or a flaky connection could
# cause a double-launch of a pipeline that's actually still running elsewhere
SAFE_TO_ACT_STATES = {"completed", "failed", "gone"}


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


def remote_pids_alive(host, pids):
    """One SSH round trip per host, checking /proc for every pid at once.
    Returns {pid: True/False/None}; None means the host couldn't be reached."""
    pids = list(dict.fromkeys(str(p) for p in pids))
    if not pids:
        return {}
    probe = "; ".join(f"[ -e /proc/{pid} ] && echo ALIVE:{pid}" for pid in pids) + "; true"
    try:
        result = subprocess.run(["ssh", *SSH_OPTS, host, probe], capture_output=True, text=True, timeout=15)
    except subprocess.TimeoutExpired:
        logging.warning(f"ssh to {host} timed out after 15s")
        return {pid: None for pid in pids}
    except OSError as e:
        logging.warning(f"ssh to {host} failed to start: {e}")
        return {pid: None for pid in pids}
    if result.returncode != 0:
        logging.warning(f"ssh to {host} failed (exit {result.returncode}): {result.stderr.strip()}")
        return {pid: None for pid in pids}
    alive = {line.split(":", 1)[1] for line in result.stdout.splitlines() if line.startswith("ALIVE:")}
    return {pid: (pid in alive) for pid in pids}


def is_alive(host, pid):
    if host == socket.gethostname():
        return pid_alive(pid)
    return remote_pids_alive(host, [pid]).get(str(pid))


def terminate_process(host, pid):
    """Sends SIGTERM to the whole process group (see start_process), locally or over SSH."""
    if host == socket.gethostname():
        try:
            os.killpg(int(pid), signal.SIGTERM)
        except ProcessLookupError:
            pass
    else:
        subprocess.run(["ssh", *SSH_OPTS, host, f"kill -TERM -{pid}"], capture_output=True)


def read_exit_code(exit_path):
    if not exit_path.exists():
        return None
    content = exit_path.read_text().strip()
    return int(content) if content else None


def run_state(exit_code, alive):
    """One of: completed, failed, running, gone (confirmed dead locally/remotely,
    but no exit code -- e.g. killed), unreachable (host couldn't be checked)."""
    if exit_code is not None:
        return "completed" if exit_code == 0 else "failed"
    if alive is True:
        return "running"
    if alive is False:
        return "gone"
    return "unreachable"


def describe_state(state, exit_code=None):
    return {
        "completed": "completed",
        "failed": f"failed (exit {exit_code})",
        "running": "running",
        "gone": "unknown (process gone, no exit code recorded)",
        "unreachable": "unknown",
    }[state]


def read_tail(path, size=65536):
    with open(path, "rb") as f:
        f.seek(0, os.SEEK_END)
        f.seek(max(0, f.tell() - size))
        return f.read().decode(errors="replace")


def progress_from_log(text):
    matches = PROGRESS_RE.findall(text)
    if not matches:
        return None
    done, total, pct = matches[-1]
    return f"{pct}% ({done}/{total})"


def format_launched_at(launched_at):
    return datetime.fromisoformat(launched_at).astimezone().strftime("%m-%d %H:%M")


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
    echo_cmd = f"echo {shlex.quote('+ ' + quoted_cmd)}"
    # detached shell wrapper: our own process exits right after Popen, so the
    # exit code has to be captured by something that outlives it
    shell_command = (
        f"cd {shlex.quote(cwd)} && "
        f"{{ {echo_cmd}; {quoted_cmd} ; }} > {shlex.quote(str(log_path))} 2>&1; "
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
            exit_code = read_exit_code(Path(previous["exit_path"]))
            alive = None if exit_code is not None else is_alive(previous["host"], previous["pid"])
            state = run_state(exit_code, alive)
            if state not in SAFE_TO_ACT_STATES:
                logging.error(
                    f"'{args.name}' is {describe_state(state, exit_code)} "
                    f"(pid={previous['pid']}, host={previous['host']}); not launching"
                )
                return 1
        record = start_process(args.name, args.command, args.ledger, "launch", os.getcwd())
        write_record(f, record)
    logging.info(f"launch '{args.name}': pid={record['pid']} log={record['log_path']}")
    return 0


def cmd_status(args):
    latest = latest_by_name(read_records_from_path(args.ledger))
    names = args.names or sorted(latest.keys())
    current_host = socket.gethostname()

    # batch the liveness checks one SSH round trip per remote host, not one per run
    pending_by_host = {}
    for name in names:
        record = latest.get(name)
        if record is not None and read_exit_code(Path(record["exit_path"])) is None:
            pending_by_host.setdefault(record["host"], []).append(record["pid"])

    alive_lookup = {}
    for host, pids in pending_by_host.items():
        if host == current_host:
            for pid in pids:
                alive_lookup[(host, pid)] = pid_alive(pid)
        elif not args.no_remote_check:
            for pid, alive in remote_pids_alive(host, pids).items():
                alive_lookup[(host, pid)] = alive

    header = ["NAME", "STATUS", "PROGRESS", "HOST", "PID", "LAUNCHED_AT", "LOG"]
    rows = []
    for name in names:
        record = latest.get(name)
        if record is None:
            rows.append([name, "unknown (not in ledger)", "-", "-", "-", "-", "-"])
            continue

        exit_code = read_exit_code(Path(record["exit_path"]))
        alive = None if exit_code is not None else alive_lookup.get((record["host"], record["pid"]))
        state = run_state(exit_code, alive)
        status = describe_state(state, exit_code)

        progress = "-"
        log_path = Path(record["log_path"])
        if log_path.exists():
            progress = progress_from_log(read_tail(log_path)) or "-"

        if args.check_jobs and state == "running":
            job_ids = job_ids_from_log(log_path.read_text())
            still_present = check_jobs_cleared(job_ids)
            status += f" ({len(still_present)}/{len(job_ids)} cluster jobs active)"

        rows.append([name, status, progress, record["host"], record["pid"], format_launched_at(record["launched_at"]), record["log_path"]])

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


def restart_safety(record):
    """State plus cluster-job-clearance check shared by `check` and `restart`."""
    exit_code = read_exit_code(Path(record["exit_path"]))
    alive = None if exit_code is not None else is_alive(record["host"], record["pid"])
    state = run_state(exit_code, alive)
    job_ids = job_ids_from_log(Path(record["log_path"]).read_text())
    still_present = check_jobs_cleared(job_ids)
    return state, exit_code, job_ids, still_present


def cmd_check(args):
    record = latest_by_name(read_records_from_path(args.ledger)).get(args.name)
    if record is None:
        logging.error(f"'{args.name}' not found in ledger")
        return 1

    state, exit_code, job_ids, still_present = restart_safety(record)
    print(f"name:        {record['name']}")
    print(f"status:      {describe_state(state, exit_code)}")
    print(f"host:        {record['host']}")
    print(f"pid:         {record['pid']}")
    print(f"log:         {record['log_path']}")
    print(f"cluster jobs: {len(job_ids)} submitted, {len(still_present)} still queued/running")
    if still_present:
        print(f"still active: {', '.join(still_present)}")

    if state not in SAFE_TO_ACT_STATES:
        print(f"verdict: {describe_state(state, exit_code)}; not safe to restart")
        return 1
    if still_present:
        print("verdict: cluster jobs still queued/running; wait for them to finish or qdel them before restarting")
        return 1
    print("verdict: safe to restart")
    return 0


def cmd_cancel(args):
    record = latest_by_name(read_records_from_path(args.ledger)).get(args.name)
    if record is None:
        logging.error(f"'{args.name}' not found in ledger")
        return 1

    job_ids = job_ids_from_log(Path(record["log_path"]).read_text())
    still_present = check_jobs_cleared(job_ids)

    failed = []
    for job_id in still_present:
        result = subprocess.run(["qdel", job_id], capture_output=True, text=True)
        if result.returncode == 0:
            logging.info(f"qdel {job_id}: cancelled")
        else:
            failed.append(job_id)
            logging.error(f"qdel {job_id} failed: {(result.stderr or result.stdout).strip()}")
    if not still_present:
        logging.info(f"'{args.name}': no active cluster jobs to cancel")

    exit_code = read_exit_code(Path(record["exit_path"]))
    if exit_code is None:
        alive = is_alive(record["host"], record["pid"])
        if alive:
            logging.info(f"'{args.name}': terminating process group (pid={record['pid']}, host={record['host']})")
            terminate_process(record["host"], record["pid"])
        elif alive is None:
            logging.warning(f"'{args.name}': could not reach host {record['host']} to check/terminate its local process")

    return 1 if failed else 0


def cmd_restart(args):
    with locked_ledger(args.ledger) as f:
        record = latest_by_name(read_records(f)).get(args.name)
        if record is None:
            logging.error(f"'{args.name}' not found in ledger")
            return 1

        state, exit_code, job_ids, still_present = restart_safety(record)
        if state not in SAFE_TO_ACT_STATES:
            logging.error(
                f"'{args.name}' is {describe_state(state, exit_code)} "
                f"(pid={record['pid']}, host={record['host']}); not restarting"
            )
            return 1
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
    status_parser.add_argument("--no-remote-check", action="store_true", help="skip SSH liveness checks for runs on other hosts (faster; they show as unreachable)")

    check_parser = subparsers.add_parser("check", help="report whether a tracked run still has cluster jobs queued/running, without restarting it")
    check_parser.add_argument("name", help="name to check")

    cancel_parser = subparsers.add_parser("cancel", help="qdel all cluster jobs still active for a tracked run, and stop its local process if still running")
    cancel_parser.add_argument("name", help="name to cancel")

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
    elif args.subcommand == "check":
        return cmd_check(args)
    elif args.subcommand == "cancel":
        return cmd_cancel(args)
    elif args.subcommand == "restart":
        return cmd_restart(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
