"""Checks for src/gwlens/waveoptics.py.

No pytest dependency (kept out on purpose, see wiki/conventions.md): this is
a plain script, run with `python3 tests/test_waveoptics.py`, following the
course's day5 convention of runnable, self-reporting check scripts. Each
`check_*` function returns (passed: bool, message: str); `main()` runs all of
them, prints a PASS/FAIL table, writes tests/CHECKS.json, and exits non-zero
if anything failed -- so this can also gate `reproduce.sh`.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from gwlens import waveoptics as wo


def check_paczynski_magnification():
    """Image-sum magnification (from `image_positions`+`magnification`, our
    own derivation of the point-lens Jacobian) must equal the independent,
    textbook Paczynski (1986) total-magnification formula, for several y."""
    worst = 0.0
    for y in [0.05, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 4.0]:
        xp, xm = wo.image_positions(y)
        assert abs(xp * xm + 1.0) < 1e-12, f"x_+ x_- != -1 at y={y}"
        mu_sum = abs(wo.magnification(xp)) + abs(wo.magnification(xm))
        pac = wo.total_magnification_paczynski(y)
        worst = max(worst, abs(mu_sum - pac) / pac)
    ok = worst < 1e-10
    return ok, f"max relative error vs Paczynski(1986) = {worst:.2e} (tol 1e-10)"


def check_low_w_limit():
    """F(w,y) -> 1 as w -> 0 (long-wavelength waves are not lensed)."""
    worst = 0.0
    for y in [0.3, 1.0, 2.0]:
        F = wo.F_point_lens(1e-4, y)
        worst = max(worst, abs(F - 1.0))
    ok = worst < 2e-3
    return ok, f"max |F(w=1e-4,y)-1| = {worst:.2e} (tol 2e-3)"


def check_geometric_optics_limit():
    """F_point_lens (exact) vs F_geometric_optics (independent stationary-
    phase derivation, including the exp(i*w*phi_+) reference-phase term --
    see wiki/log.md for the bug this caught) must agree, and the relative
    error must SHRINK as w grows (it's an asymptotic expansion)."""
    y = 1.0
    ws = [10, 50, 200, 1000]
    errs = []
    for w in ws:
        Fw = wo.F_point_lens(w, y)
        Fg = wo.F_geometric_optics(w, y)
        errs.append(abs(Fw - Fg) / abs(Fw))
    # An asymptotic (1/w) expansion: not pairwise monotonic (the sub-leading
    # correction is itself oscillatory in w), but bounded overall and smaller
    # at the largest w than at the smallest.
    all_small = all(e < 5e-3 for e in errs)
    overall_shrinks = errs[-1] < errs[0]
    ok = all_small and overall_shrinks
    msg = f"relative errors at w={ws} -> {['%.1e'%e for e in errs]}, all<5e-3: {all_small}, err(w=1000)<err(w=10): {overall_shrinks}"
    return ok, msg


def check_radial_1d_converges_to_closed_form():
    """F_radial_1d (from-scratch Bessel-reduced integral, i-epsilon
    regularized) must converge to F_point_lens as eta -> 0, at the expected
    O(eta) rate (checked, not assumed)."""
    w, y = 1.0, 1.0
    F_exact = wo.F_point_lens(w, y)
    etas = [0.02, 0.01, 0.005]
    errs = [abs(wo.F_radial_1d(w, y, eta=e) - F_exact) / abs(F_exact) for e in etas]
    decreasing = errs[0] > errs[1] > errs[2]
    # O(eta): halving eta should roughly halve the error
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    linear_ok = all(1.3 < r < 3.0 for r in ratios)
    ok = decreasing and errs[-1] < 0.01 and linear_ok
    msg = f"eta={etas} -> relerr={['%.2e'%e for e in errs]}, error ratios={['%.2f'%r for r in ratios]} (expect ~2, O(eta))"
    return ok, msg


def check_bruteforce_2d_agrees():
    """F_bruteforce_2d (2D quadrature, no Bessel-identity shortcut at all) at
    modest resolution agrees with F_point_lens at the few-% level, the same
    O(eta) regularization error seen in the 1D check above."""
    w, y = 1.0, 1.0
    F_exact = wo.F_point_lens(w, y)
    F_2d = wo.F_bruteforce_2d(w, y, x_max=25, n=1200, eta=0.005)
    relerr = abs(F_2d - F_exact) / abs(F_exact)
    ok = relerr < 0.03
    return ok, f"2D brute-force relative error = {relerr:.2e} (tol 3e-2, x_max=25,n=1200,eta=0.005)"


def check_einstein_ring_regularization():
    """At y=0 (perfect alignment, an Einstein ring) geometric optics
    DIVERGES (mu(x_+=1) is singular) but the exact wave-optics F(w,0) must
    stay finite -- this is the textbook statement that diffraction
    regularizes the point-source caustic; check it is actually true here."""
    xp, _ = wo.image_positions(1e-6)  # y->0 image position -> 1
    mu_blowup = abs(wo.magnification(xp))
    F0 = [abs(wo.F_point_lens(w, 0.0)) for w in [1, 5, 20]]
    finite = all(np.isfinite(v) and v < 100 for v in F0)
    ok = finite and mu_blowup > 1e4
    return ok, f"geometric mu(y->0) = {mu_blowup:.2e} (diverges, as expected); |F(w,0)| for w=1,5,20 -> {['%.3f'%v for v in F0]} (stays finite)"


CHECKS = [
    ("paczynski_magnification", check_paczynski_magnification),
    ("low_w_limit", check_low_w_limit),
    ("geometric_optics_limit", check_geometric_optics_limit),
    ("radial_1d_converges_to_closed_form", check_radial_1d_converges_to_closed_form),
    ("bruteforce_2d_agrees", check_bruteforce_2d_agrees),
    ("einstein_ring_regularization", check_einstein_ring_regularization),
]


def main():
    results = {}
    all_ok = True
    for name, fn in CHECKS:
        try:
            ok, msg = fn()
        except Exception as exc:  # noqa: BLE001 - report, don't hide
            ok, msg = False, f"EXCEPTION: {exc!r}"
        all_ok &= ok
        results[name] = {"passed": bool(ok), "message": msg}
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")

    out = Path(__file__).parent / "CHECKS.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
