#!/usr/bin/env python3

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from launch_pipelines import (
    parse_manifest,
    job_ids_from_log,
    pid_alive,
    read_exit_code,
    derive_status,
    build_snakemake_command,
    latest_by_name,
)

NONEXISTENT_PID = 999999999


class Tests(unittest.TestCase):
    def test_parse_manifest(self):
        lines = [
            "# comment, ignored",
            "name\tworkdir\tprofile\textra_args",
            "",
            "run1\t/some/dir\taqua\t--configfile x.yaml",
            "run2\t/other/dir\t\t",
        ]
        observed = parse_manifest(lines)
        self.assertEqual(observed, [
            {"name": "run1", "workdir": "/some/dir", "profile": "aqua", "extra_args": "--configfile x.yaml"},
            {"name": "run2", "workdir": "/other/dir", "profile": "", "extra_args": ""},
        ])

    def test_parse_manifest_missing_column(self):
        with self.assertRaises(ValueError):
            parse_manifest(["name\tprofile", "run1\taqua"])

    def test_parse_manifest_duplicate_name(self):
        lines = ["name\tworkdir", "run1\t/a", "run1\t/b"]
        with self.assertRaises(ValueError):
            parse_manifest(lines)

    def test_parse_manifest_empty(self):
        with self.assertRaises(ValueError):
            parse_manifest(["# only comments", ""])

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

    def test_build_snakemake_command(self):
        entry = {"profile": "aqua", "extra_args": "--configfile x.yaml -n"}
        self.assertEqual(
            build_snakemake_command(entry),
            ["snakemake", "--profile", "aqua", "--configfile", "x.yaml", "-n"],
        )
        self.assertEqual(
            build_snakemake_command(entry, unlock=True),
            ["snakemake", "--unlock", "--profile", "aqua", "--configfile", "x.yaml", "-n"],
        )

    def test_build_snakemake_command_bare(self):
        entry = {"profile": "", "extra_args": ""}
        self.assertEqual(build_snakemake_command(entry), ["snakemake"])

    def test_latest_by_name(self):
        records = [
            {"name": "run1", "pid": 1},
            {"name": "run2", "pid": 2},
            {"name": "run1", "pid": 3},
        ]
        self.assertEqual(latest_by_name(records), {"run1": {"name": "run1", "pid": 3}, "run2": {"name": "run2", "pid": 2}})


if __name__ == "__main__":
    unittest.main()
