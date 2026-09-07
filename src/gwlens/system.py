"""The one hierarchical triple used throughout this project (wiki/conventions.md).

Single source of truth for every physical parameter -- case scripts import
from here, never hand-type a mass or distance again. Run this file directly
to print the derived quantities and the hierarchy/regime checks that justify
the whole setup (also exercised in tests/test_system.py).
"""
import numpy as np

from . import units, chirp, geometry, waveoptics

# --- inner binary (the GW source) ---------------------------------------
M1_MSUN = 20.0
M2_MSUN = 15.0
MTOT_MSUN = M1_MSUN + M2_MSUN
MCHIRP_MSUN = chirp.chirp_mass(M1_MSUN, M2_MSUN)
F_ISCO_HZ = chirp.isco_frequency(MTOT_MSUN)

# --- outer body (the lens) -----------------------------------------------
# "Massive IMBH" (~5e4 Msun): chosen, not textbook-typical, so that the
# wave-optics parameter w=8*pi*G*M_L*f/c^3 is O(0.1-10) at BOTH epochs used
# below -- logged as a choice in wiki/log.md.
M_LENS_MSUN = 5.0e4
D_L_PC = 5000.0  # distance to the whole triple (lens ~ source, see geometry.py)

# --- outer orbit (source around the lens) ---------------------------------
P_OUT_S = 4.0 * units.YEAR_SI / 365.25  # 4 days, in seconds
E_OUT = 0.0          # circular outer orbit (kept simple; see wiki/todo.md)
I_OUT_DEG = 87.0     # near edge-on: makes "repeated lensing" pulses possible
                     # without guaranteeing an exact, over-fine-tuned transit
OMEGA_OUT = 0.0
LITTLE_OMEGA_OUT = 0.0
T_PERI_OUT = 0.0


def a_out_from_kepler(P_out_s, M_lens_msun, M_binary_msun):
    """Kepler's third law: a^3 = G(M_lens+M_binary) P^2 / (4 pi^2). Returns
    a_out in metres; case scripts convert to pc for geometry.py."""
    GM_total_SI = units.G_SI * (M_lens_msun + M_binary_msun) * units.M_SUN_SI
    a3 = GM_total_SI * P_out_s ** 2 / (4.0 * np.pi ** 2)
    return a3 ** (1.0 / 3.0)


A_OUT_M = a_out_from_kepler(P_OUT_S, M_LENS_MSUN, MTOT_MSUN)
A_OUT_PC = A_OUT_M / units.PC_SI
A_OUT_AU = A_OUT_M / units.AU_SI

# --- the two epochs on the SAME inspiral track ----------------------------
F_A_START_HZ = 10.0   # Case A ("chirp"): observation starts here, ends near merger
F_B_HZ = 0.05          # Case B ("quasi-monochromatic"): fixed observing frequency
T_OBS_B_S = 6.0 * P_OUT_S  # Case B observing baseline: 6 outer periods


def inner_separation_m(f_gw_hz, m1_msun, m2_msun):
    """Inner-binary physical separation at GW frequency f_gw (Kepler, GW freq
    = 2x orbital freq): a_in^3 = G Mtot / (pi f_gw)^2 / 4."""
    Mtot_SI = units.G_SI * (m1_msun + m2_msun) * units.M_SUN_SI
    f_orb = f_gw_hz / 2.0
    return (Mtot_SI / (2.0 * np.pi * f_orb) ** 2) ** (1.0 / 3.0)


def summary():
    lines = []
    lines.append(f"Inner binary: m1={M1_MSUN}, m2={M2_MSUN} Msun, "
                  f"Mtot={MTOT_MSUN} Msun, Mchirp={MCHIRP_MSUN:.4f} Msun")
    lines.append(f"f_isco = {F_ISCO_HZ:.2f} Hz")
    lines.append(f"Lens: M_L={M_LENS_MSUN:.3e} Msun at D_L={D_L_PC} pc")
    lines.append(f"Outer orbit: P_out={P_OUT_S/86400:.3f} d -> "
                  f"a_out={A_OUT_AU:.4f} AU = {A_OUT_PC:.3e} pc "
                  f"(Kepler's third law, i_out={I_OUT_DEG} deg, e_out={E_OUT})")

    a_in_A = inner_separation_m(F_A_START_HZ, M1_MSUN, M2_MSUN)
    a_in_B = inner_separation_m(F_B_HZ, M1_MSUN, M2_MSUN)
    lines.append(f"Inner separation at f_A_start={F_A_START_HZ} Hz: "
                  f"{a_in_A/units.AU_SI:.5f} AU "
                  f"(hierarchy a_in/a_out = {a_in_A/A_OUT_M:.2e})")
    lines.append(f"Inner separation at f_B={F_B_HZ} Hz: "
                  f"{a_in_B/units.AU_SI:.5f} AU "
                  f"(hierarchy a_in/a_out = {a_in_B/A_OUT_M:.2e})")

    tau_A = chirp.time_to_merger(F_A_START_HZ, MCHIRP_MSUN)
    tau_B = chirp.time_to_merger(F_B_HZ, MCHIRP_MSUN)
    lines.append(f"Time to merger from f_A_start: {tau_A:.3f} s "
                  f"({tau_A/P_OUT_S:.2e} outer periods -- static-lens regime)")
    lines.append(f"Time to merger from f_B: {tau_B/86400:.2f} d "
                  f"({tau_B/P_OUT_S:.1f} outer periods -- quasi-monochromatic regime)")

    w_A_start = waveoptics.w_of_frequency(F_A_START_HZ, M_LENS_MSUN)
    w_A_isco = waveoptics.w_of_frequency(F_ISCO_HZ, M_LENS_MSUN)
    w_B = waveoptics.w_of_frequency(F_B_HZ, M_LENS_MSUN)
    lines.append(f"w at f_A_start={F_A_START_HZ}Hz: {w_A_start:.3f}; "
                 f"w at f_isco={F_ISCO_HZ:.1f}Hz: {w_A_isco:.3f}")
    lines.append(f"w at f_B={F_B_HZ}Hz: {w_B:.4f}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
