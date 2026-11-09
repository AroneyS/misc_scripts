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
active cluster jobs and stops the local or remote process. `forget` removes a
run from the ledger so it stops showing up in `status`. `restart` does the
same check `cancel` would report and, only once it's clear, unlocks and
relaunches the tracked command. `watch` is an interactive, auto-refreshing
version of `status`, with keys to view/tail a run's log and to check/cancel/
forget/restart it.
"""

import argparse
import csv
import curses
import fcntl
import logging
import os
import re
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import textwrap
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

JOBID_RE = re.compile(r"external jobid '([^']+)'")
PROGRESS_RE = re.compile(r"(\d+) of (\d+) steps \(([\d.]+)%\) done")
ERROR_RE = re.compile(r"Error in rule \S+:")
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


def read_recent_log(path, start_size=65536, max_size=8 * 1024 * 1024):
    """Tail that grows until it finds a progress line, instead of a fixed window --
    a burst of job-submission lines with no completions yet can push the most recent
    "N of M steps" line out of a small fixed tail even though it's still in the log."""
    file_size = path.stat().st_size
    size = start_size
    while True:
        text = read_tail(path, size)
        if PROGRESS_RE.search(text) or size >= file_size or size >= max_size:
            return text
        size *= 4


def progress_from_log(text):
    matches = PROGRESS_RE.findall(text)
    if not matches:
        return None
    done, total, pct = matches[-1]
    return f"{pct}% ({done}/{total})"


def rule_errors_from_log(text):
    return len(ERROR_RE.findall(text))


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


STATUS_HEADER = ["NAME", "PROJECT", "STATUS", "PROGRESS", "ERRS", "HOST", "LAUNCHED_AT", "LOG"]


def project_from_cwd(cwd):
    """Short project label for a run's working directory.

    ~/projects/<project>/... -> <project>; any other path -> its basename."""
    parts = Path(cwd).parts
    if "projects" in parts:
        idx = parts.index("projects")
        if idx + 1 < len(parts):
            return parts[idx + 1]
    return Path(cwd).name or "-"


def gather_rows(latest, names, check_jobs=False, no_remote_check=False):
    """Compute one status row per name, from the ledger's latest records.

    Shared by `status` (plain table) and `watch` (curses table, plus the raw
    `record` dict each row carries so key presses can act on it)."""
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
        elif not no_remote_check:
            for pid, alive in remote_pids_alive(host, pids).items():
                alive_lookup[(host, pid)] = alive

    rows = []
    for name in names:
        record = latest.get(name)
        if record is None:
            rows.append({
                "name": name, "record": None, "state": None, "exit_code": None,
                "status": "unknown (not in ledger)", "progress": "-", "errors": "-",
                "host": "-", "launched_at": "-", "log_path": "-", "project": "-",
            })
            continue

        exit_code = read_exit_code(Path(record["exit_path"]))
        alive = None if exit_code is not None else alive_lookup.get((record["host"], record["pid"]))
        state = run_state(exit_code, alive)
        status = describe_state(state, exit_code)

        progress, errors = "-", "0"
        log_path = Path(record["log_path"])
        if log_path.exists():
            tail_text = read_recent_log(log_path)
            progress = progress_from_log(tail_text) or "-"
            errors = str(rule_errors_from_log(tail_text))

        if check_jobs and state == "running":
            job_ids = job_ids_from_log(log_path.read_text())
            still_present = check_jobs_cleared(job_ids)
            status += f" ({len(still_present)}/{len(job_ids)} cluster jobs active)"

        rows.append({
            "name": name, "record": record, "state": state, "exit_code": exit_code,
            "status": status, "progress": progress, "errors": errors,
            "host": record["host"], "launched_at": format_launched_at(record["launched_at"]),
            "log_path": record["log_path"], "project": project_from_cwd(record["cwd"]),
        })
    return rows


def render_status_table(rows, term_width=None):
    """Format `gather_rows()` output into aligned header+row strings.

    LOG is last and lowest-priority, so it's the one that gets truncated
    (with an ellipsis) to fit `term_width` -- but only when actually too
    narrow. Pass `term_width=None` (e.g. output is piped) to never truncate,
    so the full path stays copy/paste-able."""
    table_rows = [
        [r["name"], r["project"], r["status"], r["progress"], r["errors"], r["host"], r["launched_at"], r["log_path"]]
        for r in rows
    ]
    all_rows = [STATUS_HEADER] + table_rows
    log_col = len(STATUS_HEADER) - 1
    widths = [max(len(r[i]) for r in all_rows) for i in range(len(STATUS_HEADER))]

    if term_width is not None:
        available = term_width - sum(widths[:log_col]) - 2 * log_col
        if 10 <= available < widths[log_col]:
            widths[log_col] = available
            for row in all_rows[1:]:
                if len(row[log_col]) > available:
                    row[log_col] = row[log_col][: available - 1] + "…"

    return ["  ".join(c.ljust(w) for c, w in zip(row, widths)) for row in all_rows]


def cmd_status(args):
    latest = latest_by_name(read_records_from_path(args.ledger))
    names = args.names or sorted(latest.keys())
    rows = gather_rows(latest, names, check_jobs=args.check_jobs, no_remote_check=args.no_remote_check)
    term_width = shutil.get_terminal_size().columns if sys.stdout.isatty() else None
    for line in render_status_table(rows, term_width):
        print(line)
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


def cmd_forget(args):
    if not args.ledger.exists():
        logging.error(f"ledger not found: {args.ledger}")
        return 1

    # the ledger is append-only elsewhere (see locked_ledger/"a+"), but removing a row needs an
    # in-place rewrite, so this opens its own "r+" handle and takes the flock directly
    with open(args.ledger, "r+") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            records = read_records(f)
            matching = [r for r in records if r["name"] == args.name]
            if not matching:
                logging.error(f"'{args.name}' not found in ledger")
                return 1

            latest = matching[-1]
            exit_code = read_exit_code(Path(latest["exit_path"]))
            alive = None if exit_code is not None else is_alive(latest["host"], latest["pid"])
            state = run_state(exit_code, alive)
            if state not in SAFE_TO_ACT_STATES:
                logging.warning(
                    f"'{args.name}' is {describe_state(state, exit_code)}; forgetting it anyway, "
                    "but you'll lose status/restart/cancel tracking for it"
                )

            remaining = [r for r in records if r["name"] != args.name]
            f.seek(0)
            f.truncate()
            writer = csv.DictWriter(f, fieldnames=LEDGER_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(remaining)
            f.flush()
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)

    logging.info(f"'{args.name}' removed from ledger")
    return 0


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


DEFAULT_WATCH_INTERVAL = 15.0

WATCH_HELP_TEXT = textwrap.dedent(
    """
    watch key bindings:
      q           quit
      h           show this help
      up/down     move selection
      PgUp/PgDn   page
      Home        jump to top
      /           search by name
      r           refresh now (otherwise auto-refreshes every %(interval)gs)
      o / l       view full log in less
      t           tail -f the log (Ctrl-C to stop)
      c           check cluster-job status (read-only)
      k           cancel: qdel active cluster jobs, stop the process
      f           forget: remove from the ledger (loses status/restart/cancel tracking)
      R           restart: unlock + relaunch, only once previous jobs are confirmed clear

    Row colours:
      green:  completed (exit 0)
      red:    failed (non-zero exit)
      yellow: gone/unreachable (process missing with no exit code, or host unreachable)
    """
) % {"interval": DEFAULT_WATCH_INTERVAL}
WATCH_HELP_LINES = WATCH_HELP_TEXT.splitlines()


def _row_attr(row):
    state = row.get("state")
    if state == "completed":
        return curses.color_pair(1)
    if state == "failed":
        return curses.color_pair(2)
    if state in ("gone", "unreachable"):
        return curses.color_pair(3)
    return 0


def _addstr_safe(stdscr, y, x, text, maxx, attr=0):
    text = text[: max(0, maxx - x - 1)]
    if text:
        try:
            stdscr.addstr(y, x, text, attr)
        except curses.error:
            pass


def _show_watch_help(stdscr):
    stdscr.nodelay(False)
    top = 0
    while True:
        maxy, maxx = stdscr.getmaxyx()
        stdscr.erase()
        for idx, line in enumerate(WATCH_HELP_LINES[top : top + maxy - 1]):
            _addstr_safe(stdscr, idx, 0, line, maxx)
        _addstr_safe(stdscr, maxy - 1, 0, "q to exit help", maxx, curses.A_REVERSE)
        stdscr.refresh()
        key = stdscr.getch()
        if key in (ord("q"), ord("h")):
            break
        elif key == curses.KEY_DOWN:
            top = min(max(0, len(WATCH_HELP_LINES) - (maxy - 1)), top + 1)
        elif key == curses.KEY_UP:
            top = max(0, top - 1)
    stdscr.nodelay(True)
    stdscr.clear()
    curses.curs_set(0)


def cmd_watch(args):
    names_filter = args.names
    ledger = args.ledger

    def fetch_rows():
        latest = latest_by_name(read_records_from_path(ledger))
        names = names_filter or sorted(latest.keys())
        # Warnings (e.g. ssh timeouts) would be written over the curses screen on
        # every refresh; unreachable hosts already show as such in the table.
        previous = logging.root.manager.disable
        logging.disable(logging.WARNING)
        try:
            return gather_rows(latest, names, no_remote_check=args.no_remote_check)
        finally:
            logging.disable(previous)

    def _draw(stdscr):
        curses.curs_set(0)
        curses.start_color()
        curses.use_default_colors()
        curses.init_pair(1, curses.COLOR_GREEN, -1)
        curses.init_pair(2, curses.COLOR_RED, -1)
        curses.init_pair(3, curses.COLOR_YELLOW, -1)
        stdscr.nodelay(True)
        stdscr.keypad(True)

        rows, lines = [], []
        selected = 1
        offset = 0
        last_refresh = 0.0
        message, message_until = None, 0.0

        def refresh():
            nonlocal rows, lines, last_refresh
            rows = fetch_rows()
            _, maxx = stdscr.getmaxyx()
            lines = render_status_table(rows, term_width=maxx)
            last_refresh = time.time()

        def flash(text, seconds=2.5):
            nonlocal message, message_until
            message, message_until = text, time.time() + seconds

        def confirm(prompt):
            maxy, maxx = stdscr.getmaxyx()
            stdscr.move(maxy - 1, 0)
            stdscr.clrtoeol()
            _addstr_safe(stdscr, maxy - 1, 0, prompt + " [y/N] ", maxx, curses.A_REVERSE)
            stdscr.refresh()
            stdscr.nodelay(False)
            key = stdscr.getch()
            stdscr.nodelay(True)
            return key in (ord("y"), ord("Y"))

        def run_action(action_fn):
            """Leave curses so the action's normal logging/printing is
            visible, let the user read it, then resume the table."""
            curses.endwin()
            try:
                action_fn()
            finally:
                input("\nPress Enter to return...")
                stdscr.clear()
                curses.curs_set(0)
            refresh()

        refresh()

        while True:
            maxy, maxx = stdscr.getmaxyx()
            key = stdscr.getch()

            if message and time.time() > message_until:
                message = None

            selected_row = rows[selected - 1] if rows and 0 < selected <= len(rows) else None

            if key == -1:
                if time.time() - last_refresh >= args.interval:
                    refresh()
                time.sleep(0.1)
            elif key == ord("q"):
                break
            elif key == curses.KEY_UP:
                selected = max(1, selected - 1)
            elif key == curses.KEY_DOWN:
                selected = min(len(lines) - 1, selected + 1) if lines else 1
            elif key == curses.KEY_NPAGE:
                if lines:
                    selected = min(len(lines) - 1, selected + (maxy - 2))
            elif key == curses.KEY_PPAGE:
                selected = max(1, selected - (maxy - 2))
            elif key == curses.KEY_HOME:
                selected, offset = 1, 0
            elif key == ord("r"):
                refresh()
            elif key == ord("h"):
                _show_watch_help(stdscr)
            elif key == ord("/"):
                curses.echo()
                stdscr.nodelay(False)
                stdscr.move(maxy - 1, 0)
                stdscr.clrtoeol()
                stdscr.addstr(maxy - 1, 0, "/")
                try:
                    query = stdscr.getstr(maxy - 1, 1).decode()
                except Exception:
                    query = ""
                stdscr.nodelay(True)
                curses.noecho()
                if query:
                    ql = query.lower()
                    for idx, row in enumerate(rows, start=1):
                        if ql in row["name"].lower():
                            selected = idx
                            break
            elif key in (ord("o"), ord("l")):
                if selected_row is None:
                    flash("no run selected")
                elif selected_row["log_path"] != "-" and os.path.exists(selected_row["log_path"]):
                    curses.endwin()
                    subprocess.call(["less", selected_row["log_path"]])
                    stdscr.clear()
                    curses.curs_set(0)
                else:
                    flash(f"no log for '{selected_row['name']}'")
            elif key == ord("t"):
                if selected_row is None:
                    flash("no run selected")
                elif selected_row["log_path"] != "-" and os.path.exists(selected_row["log_path"]):
                    curses.endwin()
                    try:
                        subprocess.call(["tail", "-f", selected_row["log_path"]])
                    except KeyboardInterrupt:
                        pass
                    stdscr.clear()
                    curses.curs_set(0)
                else:
                    flash(f"no log for '{selected_row['name']}'")
            elif key == ord("c"):
                if selected_row is None or selected_row["record"] is None:
                    flash("no run selected")
                else:
                    name = selected_row["name"]
                    run_action(lambda: cmd_check(argparse.Namespace(name=name, ledger=ledger)))
            elif key == ord("k"):
                if selected_row is None or selected_row["record"] is None:
                    flash("no run selected")
                elif confirm(f"cancel '{selected_row['name']}'? qdels its cluster jobs and stops its process"):
                    name = selected_row["name"]
                    run_action(lambda: cmd_cancel(argparse.Namespace(name=name, ledger=ledger)))
            elif key == ord("f"):
                if selected_row is None or selected_row["record"] is None:
                    flash("no run selected")
                elif confirm(f"forget '{selected_row['name']}'? removes it from status/restart/cancel tracking"):
                    name = selected_row["name"]
                    run_action(lambda: cmd_forget(argparse.Namespace(name=name, ledger=ledger)))
            elif key == ord("R"):
                if selected_row is None or selected_row["record"] is None:
                    flash("no run selected")
                elif confirm(f"restart '{selected_row['name']}'? only runs if its previous jobs are confirmed clear"):
                    name = selected_row["name"]
                    run_action(lambda: cmd_restart(argparse.Namespace(name=name, ledger=ledger)))
            elif key not in (-1, curses.KEY_RESIZE):
                flash(f"unknown key '{chr(key)}'" if 32 <= key <= 126 else "unknown key")

            stdscr.erase()
            if lines:
                max_sel = len(lines) - 1
                selected = min(max(1, selected), max_sel)
                page = maxy - 2
                offset = max(0, min(offset, max(0, max_sel - page)))
                if selected < offset + 1:
                    offset = selected - 1
                elif selected > offset + page:
                    offset = selected - page
                _addstr_safe(stdscr, 0, 0, lines[0], maxx, curses.A_BOLD)
                visible = lines[1 + offset : 1 + offset + page]
                for i, line in enumerate(visible, start=1):
                    row_idx = offset + i - 1
                    attr = curses.A_REVERSE if (offset + i == selected) else _row_attr(rows[row_idx])
                    _addstr_safe(stdscr, i, 0, line, maxx, attr)
            else:
                offset, selected = 0, 1
            if message:
                _addstr_safe(stdscr, maxy - 2, 0, message, maxx, curses.color_pair(2))
            help_text = "q quit  h help  /  search  o log  t tail  c check  k cancel  f forget  R restart  r refresh"
            _addstr_safe(stdscr, maxy - 1, 0, help_text, maxx, curses.A_REVERSE)
            stdscr.refresh()

    curses.wrapper(_draw)
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

    forget_parser = subparsers.add_parser("forget", help="remove a tracked run from the ledger so it no longer appears in status")
    forget_parser.add_argument("name", help="name to forget")

    restart_parser = subparsers.add_parser("restart", help="restart a tracked run, after confirming its previous cluster jobs are gone")
    restart_parser.add_argument("name", help="name to restart")

    watch_parser = subparsers.add_parser("watch", help="interactively view tracked runs, auto-refreshing, with keys to view/tail a log and check/cancel/forget/restart a run")
    watch_parser.add_argument("names", nargs="*", help="names to show (default: all tracked)")
    watch_parser.add_argument("--no-remote-check", action="store_true", help="skip SSH liveness checks for runs on other hosts (faster; they show as unreachable)")
    watch_parser.add_argument("--interval", type=float, default=DEFAULT_WATCH_INTERVAL, help="auto-refresh interval in seconds (default: %(default)s)")

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
    elif args.subcommand == "forget":
        return cmd_forget(args)
    elif args.subcommand == "restart":
        return cmd_restart(args)
    elif args.subcommand == "watch":
        return cmd_watch(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
