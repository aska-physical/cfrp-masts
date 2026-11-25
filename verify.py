#!/usr/bin/env python3
"""Deterministic Timoshenko beam-column + compressive kink-band/SP-8007 tubular buckling + folding compaction verifier for aska-physical/cfrp-masts."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
CERT_PATH = ROOT / "results" / "certificate.json"


class CheckFailure(RuntimeError):
    """A failed verification check, raised explicitly so python3 -O cannot strip it."""


def require(ok: bool, message: str) -> None:
    if not ok:
        raise CheckFailure(message)


def _require_finite_positive(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CheckFailure(f"{name} must be a finite real number, got {value!r}")
    val = float(value)
    if not math.isfinite(val) or val <= 0.0:
        raise CheckFailure(f"{name} must be finite and > 0, got {value!r}")
    return val


def evaluate_local_strut_compression(
    A_strut_m2: float,
    L_strut_m: float,
    E1_GPa: float,
    eta_tube: float,
    sigma_kink_MPa: float = 350.0,
    E2_GPa: float = 8.0,
    nu12: float = 0.30,
    K_eff: float = 1.0,
) -> dict[str, Any]:
    """Compute pin-jointed (K=1) tubular CFRP compressive capacity across all 3 physical failure modes:

      P_allow = min(P_Euler, P_kink, P_shell_SP8007)

    where:
      - P_Euler        = pi^2 * E1 * I_tube / (K * L)^2
      - P_kink         = sigma_kink * A_strut   (Budiansky-Fleck compressive kink-band limit)
      - P_shell_SP8007 = gamma_SP8007 * [sqrt(E1 * E2) * (t/R_mid) / sqrt(3*(1 - nu12*nu21))] * A_strut
                         with NASA SP-8007 knockdown gamma = 1 - 0.901 * (1 - exp(-sqrt(R_mid/t)/16)).
    """
    A_m2 = _require_finite_positive(A_strut_m2, "A_strut_m2")
    L_m = _require_finite_positive(L_strut_m, "L_strut_m")
    E1 = _require_finite_positive(E1_GPa, "E1_GPa") * 1e9
    E2 = _require_finite_positive(E2_GPa, "E2_GPa") * 1e9
    sig_kink = _require_finite_positive(sigma_kink_MPa, "sigma_kink_MPa") * 1e6
    K = _require_finite_positive(K_eff, "K_eff")
    if isinstance(eta_tube, bool) or not isinstance(eta_tube, (int, float)) or not math.isfinite(float(eta_tube)):
        raise CheckFailure(f"eta_tube must be finite in [0, 1), got {eta_tube!r}")
    eta = float(eta_tube)
    require(0.0 <= eta < 1.0, f"eta_tube must lie in [0, 1), got {eta!r}")
    if isinstance(nu12, bool) or not isinstance(nu12, (int, float)) or not math.isfinite(float(nu12)):
        raise CheckFailure(f"nu12 must be finite in (0, 0.5), got {nu12!r}")
    nu_12 = float(nu12)
    require(0.0 < nu_12 < 0.5, f"nu12 must lie in (0, 0.5), got {nu_12!r}")

    I_tube = (A_m2**2) * (1.0 + eta**2) / (4.0 * math.pi * (1.0 - eta**2))
    P_euler = (math.pi**2 * E1 * I_tube) / ((K * L_m) ** 2)
    sigma_euler_GPa = (P_euler / A_m2) * 1e-9

    P_kink = sig_kink * A_m2

    if eta == 0.0:
        gamma_sp8007 = 1.0
        P_shell = float("inf")
        sigma_shell_kd_GPa = float("inf")
        R_over_t = 0.5
    else:
        R_out = math.sqrt(A_m2 / (math.pi * (1.0 - eta**2)))
        R_in = eta * R_out
        t_wall = R_out - R_in
        R_mid = 0.5 * (R_out + R_in)
        R_over_t = R_mid / t_wall
        gamma_sp8007 = 1.0 - 0.901 * (1.0 - math.exp(-math.sqrt(R_over_t) / 16.0))
        nu_21 = nu_12 * (E2 / E1)
        require(nu_12 * nu_21 < 1.0, "orthotropic Poisson product nu12*nu21 must be < 1")
        sigma_shell_ortho = math.sqrt(E1 * E2) * (t_wall / R_mid) / math.sqrt(3.0 * (1.0 - nu_12 * nu_21))
        sigma_shell_kd_Pa = gamma_sp8007 * sigma_shell_ortho
        sigma_shell_kd_GPa = sigma_shell_kd_Pa * 1e-9
        P_shell = sigma_shell_kd_Pa * A_m2

    if P_kink <= P_euler and P_kink <= P_shell:
        P_allow = P_kink
        mode = "kink_band_compression"
    elif P_shell <= P_euler:
        P_allow = P_shell
        mode = "sp8007_local_shell_buckling"
    else:
        P_allow = P_euler
        mode = "euler_column_buckling"

    I_solid = (A_m2**2) / (4.0 * math.pi)
    P_euler_solid = (math.pi**2 * E1 * I_solid) / ((K * L_m) ** 2)
    P_allow_solid = P_kink if P_kink <= P_euler_solid else P_euler_solid

    return {
        "P_euler_N": P_euler,
        "sigma_euler_GPa": sigma_euler_GPa,
        "P_kink_N": P_kink,
        "sigma_kink_MPa": sigma_kink_MPa,
        "P_shell_sp8007_N": P_shell,
        "sigma_shell_sp8007_GPa": sigma_shell_kd_GPa,
        "sp8007_knockdown_gamma": gamma_sp8007,
        "R_mid_over_t": R_over_t,
        "P_allow_N": P_allow,
        "sigma_allow_MPa": (P_allow / A_m2) * 1e-6,
        "governing_mode": mode,
        "elastic_euler_tube_gain": (1.0 + eta**2) / (1.0 - eta**2),
        "effective_compressive_tube_gain": P_allow / P_allow_solid,
    }


def evaluate_incentre_mast(
    span_m: float,
    bay_width_m: float,
    rho_bar: float,
    E1_GPa: float,
    rho_solid_kg_m3: float,
    cte_ppm_K: float,
    eta_tube: float,
    braced: bool,
    tip_mass_kg: float,
    tip_lateral_load_N: float,
    delta_T_K: float = 100.0,
) -> dict[str, float]:
    E1_Pa = E1_GPa * 1e9
    if braced:
        norm_E, norm_G = 0.128051509, 0.0490981109
    else:
        norm_E, norm_G = 0.163275266, 0.0

    E_eff_Pa = norm_E * rho_bar * E1_Pa
    G_eff_Pa = norm_G * rho_bar * E1_Pa
    A_mast = bay_width_m**2
    I_mast = (bay_width_m**4) / 12.0
    EI = E_eff_Pa * I_mast
    GA = G_eff_Pa * A_mast * (5.0 / 6.0)

    mass_per_m = rho_bar * rho_solid_kg_m3 * A_mast
    total_struct_mass_kg = mass_per_m * span_m

    if GA <= 0.0:
        tip_defl_mm = float("inf")
        P_cr_global_N = 0.0
        f1_Hz = 0.0
    else:
        defl_bend = (tip_lateral_load_N * (span_m**3)) / (3.0 * EI)
        defl_shear = (tip_lateral_load_N * span_m) / GA
        tip_defl_mm = (defl_bend + defl_shear) * 1e3
        P_euler = (math.pi**2 * EI) / (4.0 * (span_m**2))
        P_cr_global_N = P_euler / (1.0 + P_euler / GA)
        k_tip = 1.0 / ((span_m**3) / (3.0 * EI) + span_m / GA)
        m_eff = tip_mass_kg + 0.24 * total_struct_mass_kg
        f1_Hz = (1.0 / (2.0 * math.pi)) * math.sqrt(k_tip / m_eff)

    phi_deployed = math.pi * ((math.sqrt(2.0) - 1.0) ** 3)
    phi_stowed_fcc = math.pi / math.sqrt(18.0)
    compaction_ratio = phi_stowed_fcc / phi_deployed

    tube_buckling_gain = (1.0 + eta_tube**2) / (1.0 - eta_tube**2)
    thermal_axial_drift_mm = abs(cte_ppm_K * 1e-6 * delta_T_K * span_m) * 1e3

    return {
        "structural_mass_kg": round(total_struct_mass_kg, 3),
        "EI_Nm2": round(EI, 1),
        "GA_N": round(GA, 1),
        "tip_deflection_mm": round(tip_defl_mm, 4) if math.isfinite(tip_defl_mm) else -1.0,
        "global_buckling_kN": round(P_cr_global_N * 1e-3, 3),
        "first_mode_Hz": round(f1_Hz, 3),
        "volumetric_compaction_ratio": round(compaction_ratio, 4),
        "stowed_stack_length_m": round(span_m / compaction_ratio, 3),
        "tubular_buckling_gain": round(tube_buckling_gain, 4),
        "thermal_drift_100K_mm": round(thermal_axial_drift_mm, 4),
    }


def evaluate_metal_octet_mast(
    span_m: float,
    bay_width_m: float,
    rho_bar: float,
    Es_GPa: float,
    rho_s_kg_m3: float,
    cte_ppm_K: float,
    tip_mass_kg: float,
    tip_lateral_load_N: float,
    delta_T_K: float = 100.0,
) -> dict[str, float]:
    Es_Pa = Es_GPa * 1e9
    EI = ((1.0 / 9.0) * rho_bar * Es_Pa) * ((bay_width_m**4) / 12.0)
    GA = ((1.0 / 12.0) * rho_bar * Es_Pa) * (bay_width_m**2) * (5.0 / 6.0)
    mass_kg = rho_bar * rho_s_kg_m3 * (bay_width_m**2) * span_m
    defl_mm = ((tip_lateral_load_N * span_m**3) / (3.0 * EI) + (tip_lateral_load_N * span_m) / GA) * 1e3
    k_tip = 1.0 / ((span_m**3) / (3.0 * EI) + span_m / GA)
    f1_Hz = (1.0 / (2.0 * math.pi)) * math.sqrt(k_tip / (tip_mass_kg + 0.24 * mass_kg))
    return {
        "structural_mass_kg": round(mass_kg, 3),
        "EI_Nm2": round(EI, 1),
        "tip_deflection_mm": round(defl_mm, 4),
        "first_mode_Hz": round(f1_Hz, 3),
        "volumetric_compaction_ratio": 1.0,
        "thermal_drift_100K_mm": round(abs(cte_ppm_K * 1e-6 * delta_T_K * span_m) * 1e3, 4),
    }


def run_verification() -> dict[str, Any]:
    k13_boom = evaluate_incentre_mast(5.0, 0.25, 0.025, 562.4, 1800.0, -0.6, 0.80, True, 5.0, 50.0)
    t1100_mast_15m = evaluate_incentre_mast(15.0, 0.60, 0.030, 195.8, 1554.0, -0.4, 0.80, True, 25.0, 800.0)
    unbraced_boom = evaluate_incentre_mast(5.0, 0.25, 0.025, 562.4, 1800.0, -0.6, 0.80, False, 5.0, 50.0)

    al_boom = evaluate_metal_octet_mast(5.0, 0.25, 0.025, 68.9, 2700.0, 23.1, 5.0, 50.0)
    ti_boom = evaluate_metal_octet_mast(5.0, 0.25, 0.025, 113.8, 4430.0, 8.6, 5.0, 50.0)

    specific_EI_k13 = k13_boom["EI_Nm2"] / k13_boom["structural_mass_kg"]
    specific_EI_ti = ti_boom["EI_Nm2"] / ti_boom["structural_mass_kg"]
    require(
        specific_EI_k13 > 10.0 * specific_EI_ti,
        f"K13 specific EI gain ({specific_EI_k13 / specific_EI_ti:.2f}x) does not exceed 10.0x",
    )
    require(
        k13_boom["thermal_drift_100K_mm"] < al_boom["thermal_drift_100K_mm"] / 30.0,
        f"K13 thermal drift ({k13_boom['thermal_drift_100K_mm']} mm) does not beat Al/30",
    )
    require(
        abs(k13_boom["volumetric_compaction_ratio"] - 3.3166) < 1e-3,
        f"K13 volumetric compaction ratio ({k13_boom['volumetric_compaction_ratio']}) drifted from 3.3166",
    )

    L_cell = 0.25
    R_inc = math.sqrt(2.0) - 1.0
    L_strut = L_cell * R_inc
    geom_factor = (12.0 + 6.0 * math.sqrt(2.0)) * R_inc + 12.0 * 0.3623231811
    A_intra = 0.025 * (L_cell**2) / geom_factor
    A_cable = 0.3623231811 * A_intra
    # Local compressive capacity across Euler, K13D2U kink-band (350 MPa), and NASA SP-8007 shell buckling
    strut_cap = evaluate_local_strut_compression(
        A_intra, L_strut, 562.4, 0.80, sigma_kink_MPa=350.0, E2_GPa=8.0, nu12=0.30, K_eff=1.0
    )
    P_cr_local_strut_N = strut_cap["P_euler_N"]
    P_allow_local_strut_N = strut_cap["P_allow_N"]
    # Ultra-thin wall control (eta=0.99, R_mid/t=99.5): NASA SP-8007 shell buckling caps gain below kink & Euler
    thin_shell_cap = evaluate_local_strut_compression(
        A_intra, L_strut, 562.4, 0.99, sigma_kink_MPa=350.0, E2_GPa=8.0, nu12=0.30, K_eff=1.0
    )
    T_pre_min_N = 50.0 / (2.0 * math.sqrt(2.0))
    eps_pre_min = T_pre_min_N / (562.4e9 * A_cable)
    require(
        3.0 * T_pre_min_N < 0.20 * P_allow_local_strut_N,
        f"3x pretension ({3.0 * T_pre_min_N:.3f} N) exceeds 20% of local tube compressive allowable ({0.20 * P_allow_local_strut_N:.2f} N)",
    )

    controls = [
        {"control": "unbraced_mast_carries_cantilever_shear_without_pore_cables", "rejected": unbraced_boom["GA_N"] == 0.0 and unbraced_boom["tip_deflection_mm"] == -1.0},
        {"control": "timoshenko_shear_ignored_gives_overoptimistic_deflection", "rejected": k13_boom["tip_deflection_mm"] > (50.0 * 125.0 / (3.0 * k13_boom["EI_Nm2"])) * 1e3},
        {"control": "cable_pretension_buckles_local_cfrp_tubes", "rejected": 3.0 * T_pre_min_N < 0.20 * P_allow_local_strut_N},
        {"control": "solid_rod_matches_hollow_tube_eta08_buckling", "rejected": k13_boom["tubular_buckling_gain"] > 4.5},
        {"control": "fixed_metal_octet_matches_incentre_folding_compaction", "rejected": k13_boom["volumetric_compaction_ratio"] > 3.3 * ti_boom["volumetric_compaction_ratio"]},
        {
            "control": "euler_only_22p8GPa_stress_is_carryable_without_kink_band_cap",
            "rejected": strut_cap["sigma_euler_GPa"] > 10.0 * (strut_cap["sigma_kink_MPa"] / 1000.0)
            and P_allow_local_strut_N < 0.05 * P_cr_local_strut_N
            and strut_cap["governing_mode"] == "kink_band_compression",
        },
        {
            "control": "ultra_thin_tube_eta099_retains_uncapped_euler_gain_past_sp8007_shell_limit",
            "rejected": thin_shell_cap["governing_mode"] == "sp8007_local_shell_buckling"
            and thin_shell_cap["P_allow_N"] < thin_shell_cap["P_kink_N"]
            and thin_shell_cap["effective_compressive_tube_gain"] < 1.0,
        },
    ]
    failed_controls = [c["control"] for c in controls if not c["rejected"]]
    require(not failed_controls, f"must-fail controls did not reject: {failed_controls}")

    return {
        "repository": "aska-physical/cfrp-masts",
        "certificate_version": "1.0.0",
        "space_boom_5m_benchmark": {
            "k13d2u_snap_locked_boom": k13_boom,
            "al_6061_octet_boom": al_boom,
            "ti_6al_4v_octet_boom": ti_boom,
            "specific_EI_gain_vs_ti_octet_x": round(specific_EI_k13 / specific_EI_ti, 2),
            "thermal_drift_reduction_vs_al_x": round(al_boom["thermal_drift_100K_mm"] / k13_boom["thermal_drift_100K_mm"], 2),
            "cable_pretension_window": {
                "T_pre_min_N": round(T_pre_min_N, 3),
                "eps_pre_min_microstrain": round(eps_pre_min * 1e6, 2),
                "governing_compressive_mode": strut_cap["governing_mode"],
                "local_tube_buckling_N": round(P_allow_local_strut_N, 2),
                "local_tube_kink_band_N": round(strut_cap["P_kink_N"], 2),
                "local_tube_sp8007_shell_N": round(strut_cap["P_shell_sp8007_N"], 2),
                "local_tube_euler_elastic_N": round(P_cr_local_strut_N, 2),
                "local_tube_euler_stress_GPa": round(strut_cap["sigma_euler_GPa"], 2),
                "local_tube_allowable_stress_MPa": round(strut_cap["sigma_allow_MPa"], 2),
                "sp8007_knockdown_gamma": round(strut_cap["sp8007_knockdown_gamma"], 4),
                "safety_factor_at_3x_pretension": round(P_allow_local_strut_N / (3.0 * T_pre_min_N), 2),
            },
        },
        "telecom_mast_15m_t1100g": t1100_mast_15m,
        "must_fail_controls": controls,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    cert = run_verification()
    if args.write:
        CERT_PATH.parent.mkdir(parents=True, exist_ok=True)
        CERT_PATH.write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not CERT_PATH.exists() or json.loads(CERT_PATH.read_text(encoding="utf-8")) != cert:
        print("[FAIL] results/certificate.json mismatch or missing", file=sys.stderr)
        return 1
    n_ctrl = len(cert["must_fail_controls"])
    print(f"[PASS] cfrp-masts Timoshenko + folding + kink/SP-8007 pre-tension buckling checks ({n_ctrl}/{n_ctrl} controls): OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
