"""Physical compression-limit regression tests for aska-physical/cfrp-masts.

Covers:
  - Pin-jointed (K=1) three-mode tubular CFRP capacity P_allow = min(P_Euler, P_kink, P_shell_SP8007)
  - Honest rejection of unphysical 22.85 GPa Euler-only compressive stress on K13D2U tubes
    (governs at 350 MPa kink-band = 42.61 kN, SF = 803.54x at 3x pretension rather than 52,453x)
  - NASA SP-8007 thin-wall shell-buckling dominance at eta=0.98 (R_mid/t=99), capping effective
    tube gain below 1.0 even though elastic (1+eta^2)/(1-eta^2) = 49.5x
  - Fail-closed validation on NaN/Inf/bool/out-of-domain inputs and zero bare asserts across -I modes.
"""

from __future__ import annotations

import ast
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import verify


class TestCfrpMastsPhysicalHardening(unittest.TestCase):
    def test_k13d2u_strut_governed_by_kink_band_not_22p8GPa_euler(self) -> None:
        cert = verify.run_verification()
        win = cert["space_boom_5m_benchmark"]["cable_pretension_window"]
        self.assertEqual(win["governing_compressive_mode"], "kink_band_compression")
        self.assertGreater(win["local_tube_euler_stress_GPa"], 22.0)
        self.assertEqual(win["local_tube_allowable_stress_MPa"], 350.0)
        self.assertEqual(win["local_tube_buckling_N"], win["local_tube_kink_band_N"])
        self.assertLess(win["local_tube_buckling_N"], 0.02 * win["local_tube_euler_elastic_N"])
        self.assertLess(win["safety_factor_at_3x_pretension"], 1000.0)
        self.assertGreater(win["safety_factor_at_3x_pretension"], 500.0)
        self.assertLess(win["sp8007_knockdown_gamma"], 0.95)
        self.assertGreater(win["sp8007_knockdown_gamma"], 0.80)

    def test_ultra_thin_shell_eta099_switches_to_sp8007_shell_buckling(self) -> None:
        A_m2 = 1.21755e-4
        L_m = 0.25 * (math.sqrt(2.0) - 1.0)
        thin = verify.evaluate_local_strut_compression(A_m2, L_m, 562.4, 0.99, sigma_kink_MPa=350.0)
        self.assertEqual(thin["governing_mode"], "sp8007_local_shell_buckling")
        self.assertLess(thin["P_allow_N"], thin["P_kink_N"])
        self.assertGreater(thin["elastic_euler_tube_gain"], 99.0)
        self.assertLess(thin["effective_compressive_tube_gain"], 1.0)

        # Slender strut where Euler genuinely governs (L = 2.0 m, small area)
        slender = verify.evaluate_local_strut_compression(1.0e-5, 2.0, 195.8, 0.50, sigma_kink_MPa=1500.0)
        self.assertEqual(slender["governing_mode"], "euler_column_buckling")
        self.assertAlmostEqual(slender["P_allow_N"], slender["P_euler_N"], places=6)
        self.assertGreater(slender["effective_compressive_tube_gain"], 1.6)

    def test_invalid_inputs_fail_closed(self) -> None:
        for bad in (0.0, -1.0, float("nan"), float("inf"), True):
            with self.assertRaises(verify.CheckFailure):
                verify.evaluate_local_strut_compression(bad, 0.1, 562.4, 0.8)
            with self.assertRaises(verify.CheckFailure):
                verify.evaluate_local_strut_compression(1e-4, bad, 562.4, 0.8)
            with self.assertRaises(verify.CheckFailure):
                verify.evaluate_local_strut_compression(1e-4, 0.1, 562.4, 0.8, sigma_kink_MPa=bad)
        for bad_eta in (-0.1, 1.0, 1.2, float("nan"), True):
            with self.assertRaises(verify.CheckFailure):
                verify.evaluate_local_strut_compression(1e-4, 0.1, 562.4, bad_eta)

    def test_zero_bare_asserts(self) -> None:
        for rel in ("verify.py", "tests/test_physical_hardening.py"):
            tree = ast.parse((ROOT / rel).read_text(encoding="utf-8"))
            self.assertEqual([n.lineno for n in ast.walk(tree) if isinstance(n, ast.Assert)], [])


if __name__ == "__main__":
    unittest.main()
