#!/usr/bin/env python3
"""Regression tests proving verify.py fails closed under python3 -O (Issue #3)."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _copy_repo(dst: Path) -> Path:
    work = dst / "repo"
    shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns(".git", "__pycache__"))
    return work


class VerifyOptimizedModeTests(unittest.TestCase):
    def test_clean_tree_passes_all_isolated_modes(self) -> None:
        for flags in (["-I", "-B"], ["-I", "-O", "-B"], ["-I", "-B", "-W", "error"]):
            proc = subprocess.run(
                [sys.executable, *flags, "verify.py"],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, msg=f"flags={flags} stderr={proc.stderr}")
            self.assertIn("[PASS]", proc.stdout)

    def test_lowered_modulus_fails_closed_under_standard_and_O_write(self) -> None:
        """Issue #3 reproducer: 562.4 -> 150.0 GPa must fail and never overwrite certificate.json."""
        old = "evaluate_incentre_mast(5.0, 0.25, 0.025, 562.4,"
        new = "evaluate_incentre_mast(5.0, 0.25, 0.025, 150.0,"
        for mode in ([], ["-O"]):
            with tempfile.TemporaryDirectory() as td:
                work = _copy_repo(Path(td))
                verify_py = work / "verify.py"
                cert_json = work / "results" / "certificate.json"
                original_cert_bytes = cert_json.read_bytes()

                src = verify_py.read_text(encoding="utf-8")
                self.assertIn(old, src)
                verify_py.write_text(src.replace(old, new, 1), encoding="utf-8")

                for extra_args in ([], ["--write"], []):
                    proc = subprocess.run(
                        [sys.executable, "-I", *mode, "-B", "verify.py", *extra_args],
                        cwd=work,
                        text=True,
                        capture_output=True,
                        check=False,
                    )
                    self.assertNotEqual(
                        proc.returncode,
                        0,
                        msg=f"mode={mode} args={extra_args} unexpectedly passed: {proc.stdout}",
                    )
                    self.assertIn("CheckFailure", proc.stderr)

                self.assertEqual(cert_json.read_bytes(), original_cert_bytes)
                cert = json.loads(cert_json.read_text(encoding="utf-8"))
                self.assertEqual(
                    cert["space_boom_5m_benchmark"]["specific_EI_gain_vs_ti_octet_x"],
                    14.02,
                )

    def test_remaining_four_guards_fail_closed_under_O_write(self) -> None:
        """Each of the other four former bare asserts must also fail closed under python3 -O --write."""
        mutations = [
            (
                "thermal_drift",
                "evaluate_incentre_mast(5.0, 0.25, 0.025, 562.4, 1800.0, -0.6,",
                "evaluate_incentre_mast(5.0, 0.25, 0.025, 562.4, 1800.0, -5.0,",
            ),
            (
                "compaction_ratio",
                'compaction_ratio = phi_stowed_fcc / phi_deployed',
                'compaction_ratio = 2.0',
            ),
            (
                "pretension_buckling",
                "T_pre_min_N = 50.0 / (2.0 * math.sqrt(2.0))",
                "T_pre_min_N = P_cr_local_strut_N",
            ),
            (
                "must_fail_control",
                'k13_boom["tubular_buckling_gain"] > 4.5',
                'k13_boom["tubular_buckling_gain"] > 100.0',
            ),
        ]
        for label, old, new in mutations:
            with self.subTest(guard=label), tempfile.TemporaryDirectory() as td:
                work = _copy_repo(Path(td))
                verify_py = work / "verify.py"
                cert_json = work / "results" / "certificate.json"
                original_cert_bytes = cert_json.read_bytes()

                src = verify_py.read_text(encoding="utf-8")
                self.assertIn(old, src)
                verify_py.write_text(src.replace(old, new, 1), encoding="utf-8")

                proc = subprocess.run(
                    [sys.executable, "-I", "-O", "-B", "verify.py", "--write"],
                    cwd=work,
                    text=True,
                    capture_output=True,
                    check=False,
                )
                self.assertNotEqual(proc.returncode, 0, msg=f"guard={label} passed under -O --write")
                self.assertIn("CheckFailure", proc.stderr)
                self.assertEqual(cert_json.read_bytes(), original_cert_bytes)


if __name__ == "__main__":
    unittest.main()
