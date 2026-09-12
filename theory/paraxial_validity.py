"""Quantifies where the paraxial/thin-lens approximation (theory.pdf Sec.
4.5, Eq. paraxial_condition) is actually being pushed for THIS system.

D_LS(t) >> 2*R_Sch(M_L) is required for the Fresnel/paraxial diffraction
integral (the derivation behind F(w,y)) to be valid. D_LS(t) shrinks to
zero twice per outer orbit here (unlike textbook cosmological lensing),
so this is checked explicitly rather than assumed. Produces the numbers in
theory.tex's Table 4.1 -- run this script to regenerate them, not by hand.

Run: `python3 theory/paraxial_validity.py`
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "src"))

from gwlens import system, units, geometry as geo

NUMBERS = {}


def main():
    R_sch_m = 2.0 * units.msun_to_meters(system.M_LENS_MSUN)
    R_sch_AU = R_sch_m / units.AU_SI
    NUMBERS["R_Sch_AU"] = R_sch_AU

    # D_LS, y, theta_E at the actual y_min found by scanning the orbit
    # (same search used in cases/case_A_chirp/run.py -- not re-derived here,
    # imported logic would be nicer, but this is a one-off validity check,
    # not part of the physics pipeline, so a fresh scan is fine and simpler).
    t_scan = np.linspace(0.0, system.P_OUT_S, 4000, endpoint=False)
    r, nu = geo.relative_separation(t_scan, system.A_OUT_M / units.PC_SI,
                                     system.E_OUT, system.P_OUT_S, system.T_PERI_OUT)
    x_sky, y_sky, z_los = geo.orbital_plane_to_sky(
        r, nu, system.LITTLE_OMEGA_OUT, np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT)
    y_t, lensed_mask, theta_E = geo.impact_parameter_of_time(
        t_scan, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
        system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)

    idx_min = np.argmin(np.where(lensed_mask, y_t, np.inf))
    D_LS_at_ymin_AU = z_los[idx_min] * units.PC_SI / units.AU_SI
    NUMBERS["D_LS_at_ymin_AU"] = float(D_LS_at_ymin_AU)
    NUMBERS["D_LS_over_2Rsch_at_ymin"] = float(D_LS_at_ymin_AU / (2 * R_sch_AU))

    xi0_at_ymin_AU = np.sqrt(2 * R_sch_AU * D_LS_at_ymin_AU)
    NUMBERS["xi0_over_DLS_at_ymin"] = float(xi0_at_ymin_AU / D_LS_at_ymin_AU)

    # time per crossing spent with D_LS < 2*R_Sch (circular orbit, linearize
    # z_los near each crossing: z_los(t) ~ [2*pi/P * a_out * sin(i)] * t)
    P = system.P_OUT_S
    a_out_AU = system.A_OUT_M / units.AU_SI
    slope_AU_per_s = 2 * np.pi / P * a_out_AU * np.sin(np.deg2rad(system.I_OUT_DEG))
    t_marginal_s = 2 * R_sch_AU / slope_AU_per_s
    NUMBERS["marginal_window_seconds"] = float(t_marginal_s)
    NUMBERS["marginal_window_fraction_of_Pout"] = float(t_marginal_s / P)

    out = HERE / "paraxial_validity_numbers.json"
    out.write_text(json.dumps(NUMBERS, indent=2))
    print(json.dumps(NUMBERS, indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
