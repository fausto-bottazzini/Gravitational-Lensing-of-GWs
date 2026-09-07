"""Physical constants and unit conversions.

Pinned once, here, and nowhere else in the repo (see wiki/conventions.md).

Sources:
    G, c        : CODATA 2018 (exact/defined SI values where applicable;
                   c is exact by definition of the metre, G is CODATA 2018
                   recommended value, the dominant source of uncertainty in
                   GM_sun below).
    GM_sun      : IAU 2015 Resolution B3 nominal solar mass parameter
                   (defined exactly, does not depend on G's uncertainty).
"""
import numpy as np

# --- SI base constants -------------------------------------------------
C_SI = 299_792_458.0                 # m/s, exact (SI definition of the metre)
G_SI = 6.674_30e-11                  # m^3 kg^-1 s^-2, CODATA 2018

# --- IAU 2015 nominal solar parameters (exact by definition) -----------
GM_SUN_SI = 1.327_124_400_18e20      # m^3/s^2 (IAU 2015 nominal GM_sun)
M_SUN_SI = GM_SUN_SI / G_SI          # kg, derived (inherits G's uncertainty)

# --- derived geometrized units of one solar mass ------------------------
# length: G Msun / c^2 ; time: G Msun / c^3
MSUN_IN_METERS = GM_SUN_SI / C_SI**2      # m
MSUN_IN_SECONDS = GM_SUN_SI / C_SI**3     # s

# --- distance units ------------------------------------------------------
PC_SI = 3.085_677_581_49e16          # m, IAU 2015 definition
MPC_SI = 1.0e6 * PC_SI
AU_SI = 1.495_978_707_0e11           # m, IAU 2012 definition (exact)
YEAR_SI = 365.25 * 86400.0           # Julian year, s


def msun_to_seconds(m_msun):
    """Mass in solar masses -> geometrized time unit (seconds)."""
    return np.asarray(m_msun) * MSUN_IN_SECONDS


def msun_to_meters(m_msun):
    """Mass in solar masses -> geometrized length unit (metres)."""
    return np.asarray(m_msun) * MSUN_IN_METERS


if __name__ == "__main__":
    # Self-check: print the numbers this module is trusted for, so a reader
    # (or claims.yaml) can point straight at reproducible output.
    print(f"G_SI            = {G_SI:.6e} m^3 kg^-1 s^-2  (CODATA 2018)")
    print(f"c_SI            = {C_SI:.0f} m/s (exact)")
    print(f"GM_sun_SI       = {GM_SUN_SI:.8e} m^3/s^2 (IAU 2015 nominal)")
    print(f"M_sun_SI        = {M_SUN_SI:.6e} kg")
    print(f"GMsun/c^2       = {MSUN_IN_METERS/1e3:.6f} km")
    print(f"GMsun/c^3       = {MSUN_IN_SECONDS*1e6:.6f} us")
