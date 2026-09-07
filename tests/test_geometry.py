"""Checks for src/gwlens/geometry.py. Run with `python3 tests/test_geometry.py`."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from gwlens import geometry as geo


def check_kepler_equation_inverts():
    worst = 0.0
    for e in [0.0, 0.3, 0.7, 0.95]:
        M = np.linspace(0, 2 * np.pi, 37)
        E = geo.solve_kepler_equation(M, e)
        M_back = E - e * np.sin(E)
        err = np.max(np.abs(np.mod(M_back - M + np.pi, 2 * np.pi) - np.pi))
        worst = max(worst, err)
    ok = worst < 1e-10
    return ok, f"max |M - (E - e sinE)| over e in [0,0.3,0.7,0.95] = {worst:.2e}"


def check_circular_orbit_constant_radius():
    r, nu = geo.relative_separation(np.linspace(0, 20, 9), a=1.0, e=0.0, period=10.0)
    ok = np.allclose(r, 1.0, atol=1e-12)
    return ok, f"circular-orbit r range [{r.min():.12f}, {r.max():.12f}] (expect exactly 1.0)"

def check_kepler_third_law_period():
    """A circular orbit's true anomaly must advance by exactly 2*pi per
    period P -- checks relative_separation's mean-motion n=2*pi/P is wired
    correctly, independent of the Newton solver itself."""
    a, e, P = 2.0, 0.0, 5.0
    _, nu0 = geo.relative_separation(np.array([0.0]), a, e, P)
    _, nu1 = geo.relative_separation(np.array([P]), a, e, P)
    ok = abs(np.mod(nu1[0] - nu0[0], 2 * np.pi)) < 1e-10 or abs(np.mod(nu1[0] - nu0[0], 2*np.pi) - 2*np.pi) < 1e-10
    return ok, f"nu(t=0)={nu0[0]:.10f}, nu(t=P)={nu1[0]:.10f} (should differ by 2*pi mod 2*pi)"


def check_edge_on_orbit_projects_to_a_line():
    """i=pi/2, Omega=omega=0: the sky-plane y-coordinate must vanish
    identically (orbit projects to a 1D line -- the eclipsing/self-lensing
    geometry)."""
    r, nu = geo.relative_separation(np.linspace(0, 10, 25), 1.0, 0.3, 10.0)
    x, y, z = geo.orbital_plane_to_sky(r, nu, omega=0.0, i=np.pi / 2, Omega=0.0)
    ok = np.max(np.abs(y)) < 1e-10
    return ok, f"max|y_sky| for edge-on orbit = {np.max(np.abs(y)):.2e} (expect ~0)"


def check_face_on_orbit_projects_full_ellipse():
    """i=0: the sky-plane projection must equal the orbital-plane position
    exactly (no foreshortening), i.e. x_sky^2+y_sky^2 = r^2."""
    r, nu = geo.relative_separation(np.linspace(0, 10, 25), 1.0, 0.3, 10.0)
    x, y, z = geo.orbital_plane_to_sky(r, nu, omega=0.3, i=0.0, Omega=0.7)
    rho = np.hypot(x, y)
    ok = np.allclose(rho, r, atol=1e-10) and np.allclose(z, 0.0, atol=1e-10)
    return ok, f"face-on: max|rho-r|={np.max(np.abs(rho-r)):.2e}, max|z|={np.max(np.abs(z)):.2e}"


def check_einstein_radius_scaling():
    """theta_E ~ sqrt(M_L): doubling the lens mass at fixed geometry must
    scale theta_E by exactly sqrt(2)."""
    t1 = geo.einstein_radius_angle(1000.0, 1000.0, 8000.0)
    t2 = geo.einstein_radius_angle(2000.0, 1000.0, 8000.0)
    ratio = t2 / t1
    ok = abs(ratio - np.sqrt(2.0)) < 1e-10
    return ok, f"theta_E(2M)/theta_E(M) = {ratio:.10f} (expect sqrt(2)={np.sqrt(2):.10f})"


CHECKS = [
    ("kepler_equation_inverts", check_kepler_equation_inverts),
    ("circular_orbit_constant_radius", check_circular_orbit_constant_radius),
    ("kepler_third_law_period", check_kepler_third_law_period),
    ("edge_on_orbit_projects_to_a_line", check_edge_on_orbit_projects_to_a_line),
    ("face_on_orbit_projects_full_ellipse", check_face_on_orbit_projects_full_ellipse),
    ("einstein_radius_scaling", check_einstein_radius_scaling),
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
    out = Path(__file__).parent / "CHECKS_geometry.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
