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
    scale theta_E by exactly sqrt(2). Exercised through
    `impact_parameter_of_time`, which is the code path the cases actually
    use and the only place theta_E is computed."""
    a_out_pc = 1.0 * 4.84814e-6  # 1 AU in pc
    t = np.array([0.25])          # quarter period: source behind the lens
    args = (a_out_pc, 0.0, 1.0, np.deg2rad(87.0), 0.0, 0.0, 0.0)
    _, _, tE1 = geo.impact_parameter_of_time(t, *args, 1000.0, 1000.0)
    _, _, tE2 = geo.impact_parameter_of_time(t, *args, 2000.0, 1000.0)
    ratio = tE2[0] / tE1[0]
    ok = abs(ratio - np.sqrt(2.0)) < 1e-10
    return ok, f"theta_E(2M)/theta_E(M) = {ratio:.10f} (expect sqrt(2)={np.sqrt(2):.10f})"


def check_one_lensing_pulse_per_outer_period():
    """The structural claim Case B rests on, which the previous version of
    this check did not test: over one outer period, y(t) must have exactly
    ONE minimum, it must sit at conjunction, and its value must be the
    closed form obtained by putting theta = nu+omega = pi/2 into y(theta)
    (theory.pdf Eq. 4.7),

        y_min = sqrt(a_out c^2 / (4 G M_L)) * cos(i) / sqrt(sin i).

    That is what makes the signal a train of separated pulses, one per
    period, rather than a continuous modulation. Counting how much of the
    orbit is lensed -- the old assertion, kept below as a third condition --
    would pass for a y(t) with any number of minima, anywhere."""
    from gwlens import units
    a_out_pc = 1.0 * 4.84814e-6  # 1 AU in pc
    a_out_m = a_out_pc * units.PC_SI
    P_out, i_out, M_L = 4 * 86400.0, np.deg2rad(87.0), 5.0e4
    t = np.linspace(0.0, P_out, 4001, endpoint=False)
    y, lensed, _ = geo.impact_parameter_of_time(
        t, a_out_pc, 0.0, P_out, i_out, 0.0, 0.0, 0.0, M_L, 5000.0)

    # local minima of y on the lensed half
    interior = np.where(lensed)[0]
    interior = interior[(interior > 0) & (interior < len(t) - 1)]
    minima = [k for k in interior if y[k] < y[k - 1] and y[k] < y[k + 1]]
    n_min = len(minima)

    # e=0 and omega=0, so theta = 2*pi*t/P and conjunction is t = P/4
    t_min = t[minima[0]] if n_min == 1 else float("nan")
    phase_err = abs(t_min / P_out - 0.25) if n_min == 1 else np.inf

    y_closed = (np.sqrt(a_out_m / (4.0 * units.msun_to_meters(M_L)))
                * np.cos(i_out) / np.sqrt(np.sin(i_out)))
    y_err = abs(y[minima[0]] - y_closed) / y_closed if n_min == 1 else np.inf

    frac = float(np.mean(lensed))
    ok = (n_min == 1 and phase_err < 1e-3 and y_err < 1e-4
          and abs(frac - 0.5) < 0.05
          and np.all(np.isfinite(y[lensed])) and np.all(y[lensed] > 0))
    return ok, (f"{n_min} minimum of y per outer period (expect exactly 1) at "
                f"t/P={t_min/P_out:.4f} (expect 0.25); y_min={y[minima[0]]:.6f} vs "
                f"closed form {y_closed:.6f}, relative error {y_err:.1e}; "
                f"lensed fraction {frac:.3f}")


def check_impact_parameter_closed_form_and_D_L_independence():
    """`impact_parameter_of_time` is the one place this project's geometry
    genuinely departs from textbook lensing, and until now the only thing
    checked about it was that ~half the orbit comes out lensed -- which a
    formula off by any overall factor would still satisfy. Two assertions
    that would not:

    1. Substituting theta = rho/D_L and theta_E^2 = 4GM_L D_LS/(c^2 D_L^2)
       into y = theta/theta_E gives, exactly,

           y = rho / sqrt(4 G M_L D_LS / c^2)

       with every distance in metres and D_L cancelling identically. Checked
       against the function's own output to machine precision.
    2. Therefore y must not depend on D_L AT ALL. This is the physically
       non-obvious content -- in cosmological lensing y very much does depend
       on the distances -- and it is what a stray D_L, a D_L/D_S mix-up, or a
       pc/metre slip in either factor would break. Checked by moving the
       whole triple from 5 kpc to 500 pc and requiring y to be unchanged.
    """
    from gwlens import units
    a_out_pc = 1.8172 * units.AU_SI / units.PC_SI
    P_out, M_L = 4 * 86400.0, 5.0e4
    args = (a_out_pc, 0.0, P_out, np.deg2rad(87.0), 0.0, 0.0, 0.0, M_L)
    t = np.linspace(0.0, P_out, 401)

    y, lensed, _ = geo.impact_parameter_of_time(t, *args, 5000.0)
    r, nu = geo.relative_separation(t, a_out_pc, 0.0, P_out, 0.0)
    x_sky, y_sky, z_los = geo.orbital_plane_to_sky(r, nu, 0.0, np.deg2rad(87.0), 0.0)
    rho_m = np.hypot(x_sky, y_sky) * units.PC_SI
    D_LS_m = z_los * units.PC_SI
    y_closed = rho_m[lensed] / np.sqrt(
        4.0 * units.msun_to_meters(M_L) * D_LS_m[lensed])
    err_closed = float(np.max(np.abs(y[lensed] - y_closed) / y_closed))

    y_near, lensed_near, _ = geo.impact_parameter_of_time(t, *args, 500.0)
    err_D_L = float(np.max(np.abs(y_near[lensed] - y[lensed]) / y[lensed]))

    ok = err_closed < 1e-13 and err_D_L < 1e-13 and np.array_equal(lensed, lensed_near)
    return ok, (f"max relative error vs y = rho/sqrt(4GM_L D_LS/c^2): {err_closed:.2e}; "
                f"max relative change when D_L goes 5000 pc -> 500 pc: {err_D_L:.2e} "
                f"(y must be independent of D_L; both tol 1e-13)")


CHECKS = [
    ("kepler_equation_inverts", check_kepler_equation_inverts),
    ("circular_orbit_constant_radius", check_circular_orbit_constant_radius),
    ("kepler_third_law_period", check_kepler_third_law_period),
    ("edge_on_orbit_projects_to_a_line", check_edge_on_orbit_projects_to_a_line),
    ("face_on_orbit_projects_full_ellipse", check_face_on_orbit_projects_full_ellipse),
    ("einstein_radius_scaling", check_einstein_radius_scaling),
    ("one_lensing_pulse_per_outer_period", check_one_lensing_pulse_per_outer_period),
    ("impact_parameter_closed_form_and_D_L_independence",
     check_impact_parameter_closed_form_and_D_L_independence),
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
