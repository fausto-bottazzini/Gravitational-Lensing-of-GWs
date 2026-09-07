"""Checks for src/gwlens/chirp.py. Run with `python3 tests/test_chirp.py`."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from gwlens import chirp


def check_freq_of_time_matches_f0():
    """freq_of_time(t=0) must reproduce the f0 that time_to_merger(f0) used
    to define t_c -- i.e. the algebraic inversion is self-consistent."""
    Mc = chirp.chirp_mass(30.0, 30.0)
    f0 = 20.0
    t_c = chirp.time_to_merger(f0, Mc)
    f_at_0 = chirp.freq_of_time(np.array([0.0]), t_c, Mc)[0]
    relerr = abs(f_at_0 - f0) / f0
    ok = relerr < 1e-8
    return ok, f"f(t=0)={f_at_0:.6f} Hz vs input f0={f0} Hz, relerr={relerr:.2e}"


def check_closed_form_matches_numeric_ode():
    """The closed-form f(t) must agree with direct numerical integration of
    df/dt = (96/5) pi^(8/3) (GMc/c^3)^(5/3) f^(11/3) (solve_ivp, RK45,
    tight tolerances) -- an independent method, not just an algebra check."""
    Mc = chirp.chirp_mass(30.0, 30.0)
    f0 = 20.0
    t_c = chirp.time_to_merger(f0, Mc)
    t = np.linspace(0.0, 0.999 * t_c, 500)
    f_closed = chirp.freq_of_time(t, t_c, Mc)
    f_numeric = chirp.freq_of_time_numeric(t, f0, Mc)
    relerr = np.max(np.abs(f_closed - f_numeric) / f_numeric)
    ok = relerr < 1e-6
    return ok, f"max relative error closed-form vs solve_ivp = {relerr:.2e} (tol 1e-6)"


def check_isco_frequency_reasonable():
    """f_isco(M) = c^3/(pi*6^1.5*GM) matches the well-known approximate
    engineering formula f_isco[Hz] ~ 4400 Hz * (Msun/M) quoted throughout the
    PN/merger literature -- an independent numerical cross-check, not a
    tautology (4400 is not used anywhere in `isco_frequency`'s own code)."""
    f_isco = chirp.isco_frequency(60.0)
    approx = 4400.0 / 60.0  # Hz, the commonly quoted round-number formula
    relerr = abs(f_isco - approx) / approx
    ok = relerr < 0.01
    return ok, f"f_isco(M=60 Msun) = {f_isco:.2f} Hz vs ~4400/M approx = {approx:.2f} Hz, relerr={relerr:.2e}"


def check_chirp_mass_symmetry():
    """Mchirp(m1,m2) must be symmetric, and at equal mass m1=m2=m:
    Mchirp = m/2^(1/5) (from (m*m)^(3/5)/(2m)^(1/5) = m/2^(1/5)), so
    Mchirp/Mtotal = 1/(2*2^(1/5)) = 2^(-6/5) -- re-derived here algebraically
    rather than assumed, after a first attempt at this check used the wrong
    exponent (2^-1/5) and failed against the code, which was correct (see
    wiki/log.md)."""
    m1, m2 = 20.0, 40.0
    a = chirp.chirp_mass(m1, m2)
    b = chirp.chirp_mass(m2, m1)
    eq_mass_ratio = chirp.chirp_mass(30.0, 30.0) / 60.0
    expected = 2 ** (-6.0 / 5.0)
    ok = abs(a - b) < 1e-12 and abs(eq_mass_ratio - expected) < 1e-12
    return ok, f"Mc(20,40)={a:.6f}=Mc(40,20)={b:.6f}; Mc/Mtot at equal mass={eq_mass_ratio:.6f} vs 2^-6/5={expected:.6f}"


CHECKS = [
    ("freq_of_time_matches_f0", check_freq_of_time_matches_f0),
    ("closed_form_matches_numeric_ode", check_closed_form_matches_numeric_ode),
    ("isco_frequency_reasonable", check_isco_frequency_reasonable),
    ("chirp_mass_symmetry", check_chirp_mass_symmetry),
]


def main():
    results = {}
    all_ok = True
    for name, fn in CHECKS:
        try:
            ok, msg = fn()
        except Exception as exc:  # noqa: BLE001
            ok, msg = False, f"EXCEPTION: {exc!r}"
        all_ok &= ok
        results[name] = {"passed": bool(ok), "message": msg}
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")
    out = Path(__file__).parent / "CHECKS_chirp.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
