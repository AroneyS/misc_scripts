#!/usr/bin/env python3

import subprocess
import sys
import time
import unittest
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch, Mock
import launch_pipeline
from launch_pipeline import (
    job_ids_from_log,
    pid_alive,
    read_exit_code,
    run_state,
    describe_state,
    remote_pids_alive,
    is_alive,
    progress_from_log,
    read_recent_log,
    rule_errors_from_log,
    format_launched_at,
    build_unlock_command,
    latest_by_name,
    locked_ledger,
    read_records,
    read_records_from_path,
    write_record,
    restart_safety,
    start_process,
    gather_rows,
    render_status_table,
    STATUS_HEADER,
    project_from_cwd,
    _row_attr,
)

NONEXISTENT_PID = 999999999


class Tests(unittest.TestCase):
    def test_job_ids_from_log(self):
        text = (
            "Submitted job 1 with external jobid '111'.\n"
            "Submitted job 2 with external jobid '222'.\n"
            "Submitted job 3 with external jobid '111'.\n"
        )
        self.assertEqual(job_ids_from_log(text), ["111", "222"])

    def test_job_ids_from_log_none(self):
        self.assertEqual(job_ids_from_log("nothing interesting here\n"), [])

    def test_progress_from_log(self):
        text = "3 of 10 steps (30%) done\nsome other line\n7 of 10 steps (70%) done\n"
        self.assertEqual(progress_from_log(text), "70% (7/10)")

    def test_progress_from_log_none(self):
        self.assertIsNone(progress_from_log("nothing here\n"))

    def test_progress_from_log_fractional_percent(self):
        # snakemake prints fractional percentages once step counts are large
        text = "Finished jobid: 6531 (Rule: genome_prodigal)\n6 of 6531 steps (0.1%) done\n"
        self.assertEqual(progress_from_log(text), "0.1% (6/6531)")

    def test_read_recent_log_grows_past_a_small_tail(self):
        with TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            # progress line followed by enough submission chatter to push it out of a 64KB tail
            log_path.write_text("3 of 6531 steps (0.0%) done\n" + "Submitted job with external jobid 'x'.\n" * 5000)
            text = read_recent_log(log_path, start_size=64)
        self.assertEqual(progress_from_log(text), "0.0% (3/6531)")

    def test_rule_errors_from_log(self):
        text = "Error in rule failing_rule:\n    jobid: 5\nError in rule another_rule:\n"
        self.assertEqual(rule_errors_from_log(text), 2)

    def test_rule_errors_from_log_none(self):
        self.assertEqual(rule_errors_from_log("3 of 10 steps (30%) done\n"), 0)

    def test_format_launched_at(self):
        self.assertEqual(format_launched_at("2026-11-02T21:48:54.576829+00:00"), datetime.fromisoformat("2026-11-02T21:48:54.576829+00:00").astimezone().strftime("%m-%d %H:%M"))

    def test_start_process_echoes_command_to_log(self):
        with TemporaryDirectory() as tmp:
            ledger_path = Path(tmp) / "ledger.tsv"
            record = start_process("echotest", ["echo", "hi there"], ledger_path, "launch", tmp)
            for _ in range(50):
                if Path(record["exit_path"]).exists():
                    break
                time.sleep(0.1)
            log_text = Path(record["log_path"]).read_text()
        self.assertEqual(log_text.splitlines()[0], "+ echo 'hi there'")
        self.assertIn("hi there", log_text.splitlines()[1])

    def test_pid_alive(self):
        self.assertTrue(pid_alive(1))
        self.assertFalse(pid_alive(NONEXISTENT_PID))

    def test_read_exit_code(self):
        with TemporaryDirectory() as tmp:
            exit_path = Path(tmp) / "run.exit"
            self.assertIsNone(read_exit_code(exit_path))

            exit_path.write_text("0\n")
            self.assertEqual(read_exit_code(exit_path), 0)

            exit_path.write_text("1\n")
            self.assertEqual(read_exit_code(exit_path), 1)

    def test_run_state_and_describe_state(self):
        self.assertEqual(run_state(0, alive=None), "completed")
        self.assertEqual(run_state(1, alive=None), "failed")
        self.assertEqual(run_state(None, alive=True), "running")
        self.assertEqual(run_state(None, alive=False), "gone")
        self.assertEqual(run_state(None, alive=None), "unreachable")

        self.assertEqual(describe_state("completed"), "completed")
        self.assertEqual(describe_state("failed", exit_code=1), "failed (exit 1)")
        self.assertEqual(describe_state("running"), "running")
        self.assertEqual(describe_state("gone"), "unknown (process gone, no exit code recorded)")
        self.assertEqual(describe_state("unreachable"), "unknown")

    def test_remote_pids_alive_parses_output(self):
        mock_result = Mock(returncode=0, stdout="ALIVE:111\n")
        with patch("launch_pipeline.subprocess.run", return_value=mock_result) as mock_run:
            result = remote_pids_alive("otherhost", ["111", "222"])
        self.assertEqual(result, {"111": True, "222": False})
        ssh_cmd = mock_run.call_args[0][0]
        self.assertEqual(ssh_cmd[0], "ssh")
        self.assertIn("otherhost", ssh_cmd)

    def test_remote_pids_alive_unreachable(self):
        mock_result = Mock(returncode=255, stdout="", stderr="")
        with patch("launch_pipeline.subprocess.run", return_value=mock_result):
            result = remote_pids_alive("otherhost", ["111"])
        self.assertEqual(result, {"111": None})

    def test_remote_pids_alive_empty(self):
        self.assertEqual(remote_pids_alive("otherhost", []), {})

    def test_is_alive_dispatches_local_vs_remote(self):
        with patch("launch_pipeline.socket.gethostname", return_value="thishost"):
            with patch("launch_pipeline.pid_alive", return_value=True) as mock_local:
                self.assertTrue(is_alive("thishost", 123))
                mock_local.assert_called_once_with(123)

            with patch("launch_pipeline.remote_pids_alive", return_value={"123": False}) as mock_remote:
                self.assertFalse(is_alive("otherhost", 123))
                mock_remote.assert_called_once_with("otherhost", [123])

    def test_build_unlock_command_bare(self):
        command = ["snakemake", "-s", "foo.smk", "all", "--profile", "aqua"]
        self.assertEqual(
            build_unlock_command(command),
            ["snakemake", "--unlock", "-s", "foo.smk", "all", "--profile", "aqua"],
        )

    def test_build_unlock_command_wrapped(self):
        command = ["pixi", "run", "--as-is", "snakemake", "-s", "foo.smk", "all", "--profile", "aqua"]
        self.assertEqual(
            build_unlock_command(command),
            ["pixi", "run", "--as-is", "snakemake", "--unlock", "-s", "foo.smk", "all", "--profile", "aqua"],
        )

    def test_build_unlock_command_no_snakemake(self):
        self.assertIsNone(build_unlock_command(["echo", "hello"]))

    def test_latest_by_name(self):
        records = [
            {"name": "run1", "pid": "1"},
            {"name": "run2", "pid": "2"},
            {"name": "run1", "pid": "3"},
        ]
        self.assertEqual(latest_by_name(records), {"run1": {"name": "run1", "pid": "3"}, "run2": {"name": "run2", "pid": "2"}})

    def test_restart_safety_still_queued(self):
        with TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            log_path.write_text("Submitted job 1 with external jobid '111'.\n")
            record = {"host": "thishost", "exit_path": str(Path(tmp) / "run.exit"), "pid": NONEXISTENT_PID, "log_path": str(log_path)}
            with patch("launch_pipeline.socket.gethostname", return_value="thishost"), \
                 patch("launch_pipeline.check_jobs_cleared", return_value=["111"]):
                state, exit_code, job_ids, still_present = restart_safety(record)
        self.assertEqual(state, "gone")
        self.assertIsNone(exit_code)
        self.assertEqual(job_ids, ["111"])
        self.assertEqual(still_present, ["111"])

    def test_restart_safety_clear(self):
        with TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            log_path.write_text("Submitted job 1 with external jobid '111'.\n")
            exit_path = Path(tmp) / "run.exit"
            exit_path.write_text("1\n")
            record = {"host": "thishost", "exit_path": str(exit_path), "pid": NONEXISTENT_PID, "log_path": str(log_path)}
            with patch("launch_pipeline.check_jobs_cleared", return_value=[]):
                state, exit_code, job_ids, still_present = restart_safety(record)
        self.assertEqual(state, "failed")
        self.assertEqual(exit_code, 1)
        self.assertEqual(job_ids, ["111"])
        self.assertEqual(still_present, [])

    def test_restart_safety_unreachable_host_is_not_safe(self):
        with TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            log_path.write_text("")
            record = {"host": "otherhost", "exit_path": str(Path(tmp) / "run.exit"), "pid": "123", "log_path": str(log_path)}
            with patch("launch_pipeline.socket.gethostname", return_value="thishost"), \
                 patch("launch_pipeline.remote_pids_alive", return_value={"123": None}), \
                 patch("launch_pipeline.check_jobs_cleared", return_value=[]):
                state, exit_code, job_ids, still_present = restart_safety(record)
        self.assertEqual(state, "unreachable")
        self.assertNotIn(state, launch_pipeline.SAFE_TO_ACT_STATES)

    def test_cmd_cancel_qdels_still_active_jobs(self):
        with TemporaryDirectory() as tmp:
            ledger_path = Path(tmp) / "ledger.tsv"
            log_path = Path(tmp) / "run.log"
            log_path.write_text("Submitted job 1 with external jobid '111'.\n")
            record = {
                "name": "run1", "host": "node1", "pid": str(NONEXISTENT_PID), "cwd": "/tmp",
                "command": "snakemake --profile aqua", "log_path": str(log_path),
                "exit_path": str(Path(tmp) / "run.exit"), "launched_at": "2026-01-01T00:00:00+00:00", "action": "launch",
            }
            with locked_ledger(ledger_path) as f:
                write_record(f, record)

            class Args:
                name = "run1"
                ledger = ledger_path

            with patch("launch_pipeline.check_jobs_cleared", return_value=["111"]), \
                 patch("launch_pipeline.is_alive", return_value=False), \
                 patch("launch_pipeline.subprocess.run") as mock_run:
                mock_run.return_value.returncode = 0
                rc = launch_pipeline.cmd_cancel(Args())

        mock_run.assert_called_once_with(["qdel", "111"], capture_output=True, text=True)
        self.assertEqual(rc, 0)

    def test_cmd_forget_removes_only_matching_name(self):
        with TemporaryDirectory() as tmp:
            ledger_path = Path(tmp) / "ledger.tsv"
            for name in ("run1", "run2"):
                record = {
                    "name": name, "host": "node1", "pid": str(NONEXISTENT_PID), "cwd": "/tmp",
                    "command": "snakemake --profile aqua", "log_path": "/tmp/none.log",
                    "exit_path": "/tmp/none.exit", "launched_at": "2026-01-01T00:00:00+00:00", "action": "launch",
                }
                with locked_ledger(ledger_path) as f:
                    write_record(f, record)

            class Args:
                name = "run1"
                ledger = ledger_path

            rc = launch_pipeline.cmd_forget(Args())
            remaining = latest_by_name(read_records_from_path(ledger_path))

        self.assertEqual(rc, 0)
        self.assertNotIn("run1", remaining)
        self.assertIn("run2", remaining)

    def test_cmd_forget_not_found(self):
        with TemporaryDirectory() as tmp:
            ledger_path = Path(tmp) / "ledger.tsv"
            record = {
                "name": "run1", "host": "node1", "pid": "1", "cwd": "/tmp",
                "command": "snakemake --profile aqua", "log_path": "/tmp/none.log",
                "exit_path": "/tmp/none.exit", "launched_at": "2026-01-01T00:00:00+00:00", "action": "launch",
            }
            with locked_ledger(ledger_path) as f:
                write_record(f, record)

            class Args:
                name = "nosuchname"
                ledger = ledger_path

            rc = launch_pipeline.cmd_forget(Args())

        self.assertEqual(rc, 1)

    def test_ledger_round_trip_and_header_written_once(self):
        with TemporaryDirectory() as tmp:
            ledger_path = Path(tmp) / "ledger.tsv"
            record = {
                "name": "run1", "host": "node1", "pid": "123", "cwd": "/tmp",
                "command": "snakemake --profile aqua", "log_path": "/tmp/run1.log",
                "exit_path": "/tmp/run1.exit", "launched_at": "2026-01-01T00:00:00+00:00", "action": "launch",
            }
            with locked_ledger(ledger_path) as f:
                write_record(f, record)
            with locked_ledger(ledger_path) as f:
                records = read_records(f)
            self.assertEqual(records, [record])
            self.assertEqual(ledger_path.read_text().splitlines()[0].split("\t"), [
                "name", "host", "pid", "cwd", "command", "log_path", "exit_path", "launched_at", "action",
            ])


    def test_gather_rows_not_in_ledger(self):
        rows = gather_rows({}, ["nosuchname"])
        self.assertEqual(rows, [{
            "name": "nosuchname", "record": None, "state": None, "exit_code": None,
            "project": "-", "status": "unknown (not in ledger)", "progress": "-", "errors": "-",
            "host": "-", "launched_at": "-", "log_path": "-",
        }])

    def test_gather_rows_completed(self):
        with TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            log_path.write_text("3 of 10 steps (30%) done\n")
            exit_path = Path(tmp) / "run.exit"
            exit_path.write_text("0\n")
            record = {
                "name": "run1", "host": "node1", "pid": "123", "cwd": "/tmp",
                "command": "snakemake --profile aqua", "log_path": str(log_path),
                "exit_path": str(exit_path), "launched_at": "2026-01-01T00:00:00+00:00", "action": "launch",
            }
            rows = gather_rows({"run1": record}, ["run1"])
        self.assertEqual(rows[0]["name"], "run1")
        self.assertEqual(rows[0]["state"], "completed")
        self.assertEqual(rows[0]["exit_code"], 0)
        self.assertEqual(rows[0]["status"], "completed")
        self.assertEqual(rows[0]["progress"], "30% (3/10)")
        self.assertEqual(rows[0]["errors"], "0")
        self.assertIs(rows[0]["record"], record)

    def test_gather_rows_check_jobs_appends_cluster_job_count(self):
        with TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            log_path.write_text("Submitted job 1 with external jobid '111'.\n")
            record = {
                "name": "run1", "host": "thishost", "pid": str(NONEXISTENT_PID), "cwd": "/tmp",
                "command": "snakemake --profile aqua", "log_path": str(log_path),
                "exit_path": str(Path(tmp) / "run.exit"), "launched_at": "2026-01-01T00:00:00+00:00", "action": "launch",
            }
            with patch("launch_pipeline.socket.gethostname", return_value="thishost"), \
                 patch("launch_pipeline.pid_alive", return_value=True), \
                 patch("launch_pipeline.check_jobs_cleared", return_value=["111"]):
                rows = gather_rows({"run1": record}, ["run1"], check_jobs=True)
        self.assertEqual(rows[0]["state"], "running")
        self.assertEqual(rows[0]["status"], "running (1/1 cluster jobs active)")

    def test_render_status_table_header(self):
        lines = render_status_table([])
        self.assertEqual(lines, ["  ".join(STATUS_HEADER)])

    def test_render_status_table_truncates_log_when_narrow(self):
        rows = [{
            "name": "run1", "project": "proj", "status": "running", "progress": "-", "errors": "0",
            "host": "node1", "launched_at": "01-01 00:00", "log_path": "/" + "a" * 300 + "/run1.log",
        }]
        lines = render_status_table(rows, term_width=80)
        log_field = lines[1].split("  ")[-1]
        self.assertTrue(log_field.endswith("…"))
        self.assertLess(len(log_field), len(rows[0]["log_path"]))

    def test_render_status_table_no_truncation_when_term_width_none(self):
        rows = [{
            "name": "run1", "project": "proj", "status": "running", "progress": "-", "errors": "0",
            "host": "node1", "launched_at": "01-01 00:00", "log_path": "/a/very/long/path/to/a/run1.log",
        }]
        lines = render_status_table(rows, term_width=None)
        self.assertIn(rows[0]["log_path"], lines[1])

    def test_row_attr(self):
        # curses.color_pair() needs an initscr'd terminal, so stub it rather
        # than actually initialising curses just for this mapping check
        with patch("launch_pipeline.curses.color_pair", side_effect=lambda n: n):
            self.assertEqual(_row_attr({"state": "completed"}), 1)
            self.assertEqual(_row_attr({"state": "failed"}), 2)
            self.assertEqual(_row_attr({"state": "gone"}), 3)
            self.assertEqual(_row_attr({"state": "unreachable"}), 3)
        self.assertEqual(_row_attr({"state": "running"}), 0)

    def test_watch_subcommand_help(self):
        script = Path(__file__).resolve().parent / "launch_pipeline.py"
        result = subprocess.run([sys.executable, str(script), "watch", "--help"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("--interval", result.stdout)



class ProjectFromCwdTests(unittest.TestCase):
    def test_projects_subdir(self):
        self.assertEqual(project_from_cwd("/home/u/projects/foo/results/x"), "foo")

    def test_other_path_uses_basename(self):
        self.assertEqual(project_from_cwd("/tmp/bar"), "bar")

if __name__ == "__main__":
    unittest.main()
