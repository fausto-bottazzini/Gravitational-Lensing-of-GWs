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
    phase derivation) must agree, and the error at the largest w must be
    smaller than at the smallest, as an asymptotic expansion requires. It is
    not monotonic in between: |F| oscillates and so does the relative error,
    so what this asserts is the two endpoints and a ceiling, not
    monotonicity. Both evaluators are referenced to the strong image's own
    arrival (see the module docstring of waveoptics.py)."""
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
    regularized) must converge to F_point_lens as eta -> 0, at the linear
    rate, checked rather than assumed. This is the ONLY independent
    validation `F_point_lens` has: theory.pdf Sec. 3.4 is explicit that
    `F_geometric_optics` is not a second opinion but a second evaluation of
    the same integral, so agreement with it proves nothing about the
    closed form.

    Run at TWO points, and the second one is the reason: the regulator's
    small parameter is eta*w*x_+^2/2, not eta, and w=1 is the single place
    in the plane where those two readings agree numerically. A check there
    alone therefore cannot tell them apart, and says nothing about the rest
    of the range. `F_hybrid` calls `F_point_lens` for every w < 30 and
    switches to the asymptotic form above, so the demanding end of its
    operating range is w -> 30; the second point sits just under it, at
    w=25 with Case A's y=1.589. There eta must be ~40x smaller for the same
    accuracy -- at eta=1e-3 the error is 6.1e-2, sixty times what a bare
    "O(eta)" reading would predict. Case B's own w_B = 0.31 is easier than
    either point, so it is covered a fortiori rather than spending test
    time on it."""
    puntos = [(1.0, 1.0, [0.02, 0.01, 0.005], 1.0e-2),
              (25.0, 1.589, [1.0e-3, 2.5e-4], 2.0e-2)]
    ok, partes = True, []
    for w, y, etas, tol in puntos:
        F_exact = wo.F_point_lens(w, y)
        errs = [abs(wo.F_radial_1d(w, y, eta=e) - F_exact) / abs(F_exact)
                for e in etas]
        ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
        esperado = [etas[i] / etas[i + 1] for i in range(len(etas) - 1)]
        # linear in eta: the error ratio must track the eta ratio
        lineal = all(0.7 * s < r < 1.4 * s for r, s in zip(ratios, esperado))
        baja = all(errs[i] > errs[i + 1] for i in range(len(errs) - 1))
        ok = ok and baja and lineal and errs[-1] < tol
        partes.append(
            f"(w={w:.2f}, y={y:.3f}) eta={etas} -> "
            f"relerr={['%.1e' % e for e in errs]}, "
            f"ratios={['%.2f' % r for r in ratios]} "
            f"(expect {['%.1f' % s for s in esperado]}, linear in eta)")
    return ok, "; ".join(partes)


def check_bruteforce_2d_agrees():
    """F_bruteforce_2d (2D quadrature, no Bessel-identity shortcut at all) at
    modest resolution agrees with F_point_lens at the few-% level, the same
    regularization error seen in the 1D check above. Deliberately left at
    w=y=1: the fixed grid (x_max, n) cannot follow the integrand as w
    grows, and with these defaults the error is 7.5e-1 by w=10. Raising
    the resolution enough to make this meaningful at Case B's w costs
    minutes, and `F_radial_1d` already provides that check there at a
    thousandth of the cost. This one is a shape check on the 2D form."""
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


def check_hybrid_matches_at_threshold():
    """At w=w_geo_threshold exactly, F_point_lens (exact) and
    F_geometric_optics (what F_hybrid switches to just above threshold)
    must already agree well, at the y values Case A/B actually use --
    isolates the approximation error at the switch point, rather than
    conflating it with real F(w) variation between two different w (an
    earlier version of this check compared w=19.9 to w=20.1 and mixed the
    two effects together)."""
    w = 30.0
    ys = [0.3, 1.0, 1.589, 1.6, 2.5]
    errs = {}
    for y in ys:
        exact = wo.F_point_lens(w, y)
        errs[y] = abs(exact - wo.F_geometric_optics(w, y)) / abs(exact)
    worst = max(errs.values())
    ok = worst < 0.04
    # Report per-y, not just the max: the max is dominated by the SMALLEST y,
    # where the two images are close together and highly magnified and so the
    # stationary-phase approximation is at its worst. Quoting only that number
    # (as this check used to) makes the switch look ~100x cruder than it is
    # anywhere the project actually evaluates F -- at Case A's y_A=1.589 the
    # discontinuity at the threshold is 3.1e-4, and Case A's whole band sits
    # above w=30 in any case. y=0.3 is kept in the list precisely because it
    # is the hard case, not because anything here uses it.
    detail = ", ".join(f"y={y}: {errs[y]:.1e}" for y in ys)
    return ok, (f"relative error |F_point_lens-F_geometric_optics| at w=30 -- {detail} "
                f"(max {worst:.2e}, tol 4e-2). The max is set by the smallest y; "
                f"at Case A's y_A=1.589 the switch is smooth to {errs[1.589]:.1e}.")


def _lensed_test_pulse(y, M_L_msun, f0, t_ref, sigma, w_evaluator, fs=4096.0, T=16.0):
    """Shared machinery for the two causality checks below: build a short
    Gaussian-enveloped test pulse, lens it via FFT/multiply-by-F/IFFT exactly
    as cases/case_A_chirp/run.py does, using `w_evaluator` (a callable
    w,y -> complex) for F at every positive frequency. Returns (t, lensed,
    dT_seconds, mu_plus, mu_minus)."""
    x_plus, x_minus = wo.image_positions(y)
    mu_plus = abs(wo.magnification(x_plus))
    mu_minus = abs(wo.magnification(x_minus))
    dT_dimless = wo.time_delay_difference(y)

    from gwlens import units  # already on sys.path, see top of this file

    t_char = 4.0 * units.msun_to_seconds(M_L_msun)  # seconds, = 4GM/c^3
    dT_seconds = dT_dimless * t_char

    n = int(T * fs)
    n = 1 << (n - 1).bit_length()
    t = np.arange(n) / fs
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)

    pulse = np.exp(-0.5 * ((t - t_ref) / sigma) ** 2) * np.cos(2 * np.pi * f0 * (t - t_ref))
    H = np.fft.rfft(pulse)
    w_arr = wo.w_of_frequency(np.abs(freqs), M_L_msun)
    F = np.ones_like(H, dtype=complex)
    pos = freqs > 0
    F[pos] = np.array([w_evaluator(w, y) for w in w_arr[pos]])
    lensed = np.fft.irfft(H * F, n=n)
    return t, lensed, dT_seconds, mu_plus, mu_minus


def _local_peak(t, lensed, center, width=0.3):
    """(time, |amplitude|) of the largest sample of `lensed` within `width`
    of `center`."""
    mask = (t > center - width) & (t < center + width)
    idx = np.where(mask)[0]
    j = idx[np.argmax(np.abs(lensed[idx]))]
    return t[j], np.abs(lensed[j])


def check_causality_of_lensed_pulse():
    """THE check that catches the sign-convention bug documented at the top
    of waveoptics.py -- every other check here compares two evaluators of F
    to each other, and |F| (all any of them test) is conjugation-invariant,
    so none of them could have caught a global sign error. This one
    actually applies F(f) to a short test pulse via FFT/IFFT (exactly as
    cases/case_A_chirp/run.py does) and checks the physical requirement
    that the weaker (saddle-point) image comes out AFTER the stronger
    (minimum) image, not before -- a genuine causality/reconstruction test,
    not a self-consistency one. Also checks the reference-phase fix
    (2026-09-07, see wiki/log.md and the module docstring): with it, the
    strong image sits at EXACTLY t_ref (not offset by an arbitrary
    convention-dependent amount) and at amplitude sqrt(mu_plus) (not 1),
    and the weak image at amplitude sqrt(mu_minus) (not sqrt(mu_minus/mu_plus)
    -- the earlier version of this check used the wrong normalization here
    too, see wiki/log.md). y=1.0 gives w>>30 across the whole pulse band, so
    this exercises F_geometric_optics only -- see
    check_causality_of_lensed_pulse_low_w below for F_point_lens."""
    y = 1.0
    M_L_msun = 5.0e4
    f0, t_ref, sigma = 40.0, 6.0, 0.05
    t, lensed, dT_seconds, mu_plus, mu_minus = _lensed_test_pulse(
        y, M_L_msun, f0, t_ref, sigma, wo.F_hybrid)

    t_strong, amp_strong = _local_peak(t, lensed, t_ref, 0.3)
    t_after, amp_after = _local_peak(t, lensed, t_ref + dT_seconds, 0.3)
    _, amp_before = _local_peak(t, lensed, t_ref - dT_seconds, 0.3)
    expected_strong = np.sqrt(mu_plus)
    expected_weak = np.sqrt(mu_minus)

    strong_at_t_ref = abs(t_strong - t_ref) < 2.0 / 4096.0  # within a couple of samples
    strong_amp_ok = abs(amp_strong - expected_strong) < 0.05 * expected_strong
    # the weak image must show up clearly AFTER, not before
    # two-sided, and tight: a one-sided `> 0.2 * expected_weak` passed with
    # the sqrt(mu_-/mu_+) normalization this check exists to rule out.
    weak_after_ok = (abs(amp_after - expected_weak) < 0.05 * expected_weak
                     and amp_before < 0.05 * expected_weak)
    ok = strong_at_t_ref and strong_amp_ok and weak_after_ok
    return ok, (f"strong image at t_ref+{t_strong - t_ref:.4f}s, amplitude {amp_strong:.4f} "
                f"(expect sqrt(mu_+)={expected_strong:.4f}); weak image AFTER = {amp_after:.4f} "
                f"(expect sqrt(mu_-)={expected_weak:.4f}), BEFORE = {amp_before:.4f} (expect ~0) "
                f"-- dT={dT_seconds:.3f}s")


def check_causality_of_lensed_pulse_low_w():
    """Same as check_causality_of_lensed_pulse but at low w (~6, well below
    w_geo_threshold=30), so F_hybrid dispatches to F_point_lens -- the
    evaluator all of Case B and the sub-threshold part of Case A actually
    use, and which the check above never exercises (an independent review
    found this: with y=1.0/f0=40Hz/M_L=5e4Msun, 1.6e-34 of the pulse's power
    sits below w_geo_threshold, see wiki/log.md). A lower carrier frequency
    keeps the same M_L but puts the whole pulse's band below threshold.

    What it does NOT assert, and why (fixed 2026-09-07, see wiki/log.md):
    an earlier version of this check required the trailing response to reach
    0.2*sqrt(mu_-), and passed with a measured 0.254 against an "expected"
    0.413 -- a 38% miss inside a 5x-wide tolerance window, i.e. a check that
    looked like an amplitude validation while validating nothing. The
    expectation was the wrong one: sqrt(mu_-) is the GEOMETRIC-optics weak-
    image amplitude, and at w~6 there are no separated images to have an
    amplitude -- the second image is smeared over a time comparable to its
    own delay, which is the entire content of being in the wave-optics
    regime. Replaced by the assertion that is exact at any w: CAUSALITY.
    All of the lensing-induced response must arrive after the strong image,
    so the energy in the window preceding it must sit at the numerical floor.
    That is what a conjugation error breaks, and it breaks it by orders of
    magnitude, not by tens of percent. The geometric-optics comparison is
    still reported in the message, as information rather than as a
    pass/fail criterion.
    """
    y = 1.0
    M_L_msun = 5.0e4
    f0, t_ref, sigma = 1.0, 6.0, 0.05
    t, lensed, dT_seconds, mu_plus, mu_minus = _lensed_test_pulse(
        y, M_L_msun, f0, t_ref, sigma, wo.F_hybrid)

    w_at_f0 = wo.w_of_frequency(f0, M_L_msun)
    t_strong, amp_strong = _local_peak(t, lensed, t_ref, 0.3)
    expected_strong = np.sqrt(mu_plus)

    # energy strictly before / strictly after the strong image, excluding the
    # pulse itself (+-5 sigma) so this measures the LENSING response only
    pad = 5.0 * sigma
    before = (t > t_ref - dT_seconds - pad) & (t < t_ref - pad)
    after = (t > t_ref + pad) & (t < t_ref + dT_seconds + pad)
    e_before = float(np.sum(lensed[before] ** 2))
    e_after = float(np.sum(lensed[after] ** 2))

    strong_at_t_ref = abs(t_strong - t_ref) < 2.0 / 4096.0
    strong_amp_ok = abs(amp_strong - expected_strong) < 0.02 * expected_strong
    causal = e_after > 1e4 * e_before
    ok = (w_at_f0 < 30.0) and strong_at_t_ref and strong_amp_ok and causal
    _, amp_after = _local_peak(t, lensed, t_ref + dT_seconds, 0.3)
    return ok, (f"w(f0)={w_at_f0:.2f} (<30, so F_point_lens is exercised); strong image at "
                f"t_ref+{t_strong - t_ref:.4f}s, amplitude {amp_strong:.4f} "
                f"(expect sqrt(mu_+)={expected_strong:.4f}, tol 2%); trailing/leading "
                f"energy ratio = {e_after / e_before:.1e} (require > 1e4). For information "
                f"only: peak of the trailing response is {amp_after:.4f}, vs the "
                f"geometric-optics sqrt(mu_-)={np.sqrt(mu_minus):.4f} it is NOT expected to "
                f"reach at this w -- dT={dT_seconds:.3f}s")


CHECKS = [
    ("causality_of_lensed_pulse", check_causality_of_lensed_pulse),
    ("causality_of_lensed_pulse_low_w", check_causality_of_lensed_pulse_low_w),
    ("paczynski_magnification", check_paczynski_magnification),
    ("low_w_limit", check_low_w_limit),
    ("geometric_optics_limit", check_geometric_optics_limit),
    ("radial_1d_converges_to_closed_form", check_radial_1d_converges_to_closed_form),
    ("bruteforce_2d_agrees", check_bruteforce_2d_agrees),
    ("einstein_ring_regularization", check_einstein_ring_regularization),
    ("hybrid_matches_at_threshold", check_hybrid_matches_at_threshold),
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
