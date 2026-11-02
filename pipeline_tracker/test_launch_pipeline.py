#!/usr/bin/env python3

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from launch_pipeline import (
    job_ids_from_log,
    pid_alive,
    read_exit_code,
    derive_status,
    build_unlock_command,
    latest_by_name,
    locked_ledger,
    read_records,
    write_record,
    restart_safety,
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

    def test_derive_status(self):
        self.assertEqual(derive_status(0, alive=False), "completed")
        self.assertEqual(derive_status(1, alive=False), "failed (exit 1)")
        self.assertEqual(derive_status(None, alive=True), "running")
        self.assertEqual(derive_status(None, alive=False), "unknown (process gone, no exit code recorded)")

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
            record = {"exit_path": str(Path(tmp) / "run.exit"), "pid": NONEXISTENT_PID, "log_path": str(log_path)}
            with patch("launch_pipeline.check_jobs_cleared", return_value=["111"]):
                status, job_ids, still_present = restart_safety(record)
        self.assertEqual(status, "unknown (process gone, no exit code recorded)")
        self.assertEqual(job_ids, ["111"])
        self.assertEqual(still_present, ["111"])

    def test_restart_safety_clear(self):
        with TemporaryDirectory() as tmp:
            log_path = Path(tmp) / "run.log"
            log_path.write_text("Submitted job 1 with external jobid '111'.\n")
            exit_path = Path(tmp) / "run.exit"
            exit_path.write_text("1\n")
            record = {"exit_path": str(exit_path), "pid": NONEXISTENT_PID, "log_path": str(log_path)}
            with patch("launch_pipeline.check_jobs_cleared", return_value=[]):
                status, job_ids, still_present = restart_safety(record)
        self.assertEqual(status, "failed (exit 1)")
        self.assertEqual(job_ids, ["111"])
        self.assertEqual(still_present, [])

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


if __name__ == "__main__":
    unittest.main()
