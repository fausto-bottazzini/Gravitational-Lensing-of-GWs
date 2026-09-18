"""Checks for src/gwlens/taylorf2.py. Run with `python3 tests/test_taylorf2.py`."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from gwlens import taylorf2, chirp, units


def check_leading_order_group_delay():
    """The SPA stationary-phase relation, in THIS project's FT convention
    (see spa_phase docstring), gives group delay (1/2pi)dPsi/df =
    -(t_c - tau(f)) for the LEADING (v^0) term alone, tau from chirp.py's
    independently-derived (and separately checked, tests/test_chirp.py)
    time_to_merger. NOTE what this does and does not touch: `psi_leading`
    below is a local hand-written copy of the leading term, so what is
    validated here is that ALGEBRA against chirp.py, plus `_v_of_f`, which
    is the only thing it imports from the module. It does NOT exercise
    `spa_phase` -- check_spa_phase_series_matches_literature does that. It
    does NOT by itself catch an overall sign flip relative to numpy's
    ifft (a constant-derivative-sign error), which is what
    check_ifft_reconstructs_chirp_at_tc below catches instead."""
    m1, m2 = 20.0, 15.0
    Mc = chirp.chirp_mass(m1, m2)
    f0 = 10.0
    t_c = chirp.time_to_merger(f0, Mc)

    def psi_leading(f):
        mtot = m1 + m2
        eta = m1 * m2 / mtot ** 2
        v = taylorf2._v_of_f(f, mtot)
        return -(2 * np.pi * f * t_c - np.pi / 4.0 + (3.0 / (128 * eta * v ** 5)))

    f = np.linspace(10.0, 100.0, 500)
    psi = psi_leading(f)
    group_delay = np.gradient(psi, f) / (2 * np.pi)
    predicted = -(t_c - chirp.time_to_merger(f, Mc))

    # trim edges (np.gradient is less accurate there)
    relerr = np.abs(group_delay - predicted)[10:-10] / np.abs(predicted)[10:-10]
    ok = np.nanmax(relerr) < 2e-3
    return ok, f"max relative error, leading-order group delay vs -(t_c-tau(f)) = {np.nanmax(relerr):.2e} (tol 2e-3)"


def check_ifft_reconstructs_chirp_at_tc():
    """The real, physical test: numpy.fft.irfft of htilde(f) (banded to
    [f0,f_isco], both edges tapered) must reconstruct a time-domain chirp
    peaking AT t_c, and that peak location must NOT move when the amount
    of zero-padding changes (a peak that drifts with the padding is the
    signature of a sign error aliasing power to the wrong end of the
    periodic IFFT window -- exactly the bug this check replaced; see
    spa_phase's docstring and wiki/log.md)."""
    m1, m2 = 20.0, 15.0
    Mc = chirp.chirp_mass(m1, m2)
    f0, f_isco = 10.0, chirp.isco_frequency(m1 + m2)
    t_c = chirp.time_to_merger(f0, Mc)
    fs = 4.0 * f_isco
    peak_times = []
    for pad_mult in [2, 4, 8]:
        n = int(np.ceil((t_c - chirp.time_to_merger(f_isco, Mc)) * fs))
        n = 1 << (n - 1).bit_length()  # round up to a power of two
        n_pad = n * pad_mult  # n already a power of two, so this equals 1 << (n*pad_mult-1).bit_length()
        freqs = np.fft.rfftfreq(n_pad, d=1.0 / fs)
        inband = (freqs >= f0) & (freqs <= f_isco)
        H = np.zeros_like(freqs, dtype=complex)
        H[inband] = taylorf2.htilde(freqs[inband], m1, m2, t_c, 5.0e-3)
        h = np.fft.irfft(H, n=n_pad)
        t = np.arange(n_pad) / fs
        peak_times.append(t[np.argmax(np.abs(h))])
    spread = max(peak_times) - min(peak_times)
    close_to_tc = all(abs(pt - t_c) / t_c < 0.05 for pt in peak_times)
    ok = spread < 0.05 * t_c and close_to_tc
    # float(), not the bare np.float64 list: numpy >=2.0 reprs those as
    # "np.float64(...)", numpy <2.0 as a plain number, so this message (and
    # so this committed CHECKS_taylorf2.json) was not byte-reproducible
    # across the repo's own two supported environments; caught by an
    # independent review, see wiki/log.md.
    peak_times_f = [float(pt) for pt in peak_times]
    return ok, f"peak times at pad_mult=[2,4,8]: {peak_times_f} (t_c={t_c:.3f}); spread={spread:.4f}s"


def check_pn_terms_are_small_corrections():
    """A precondition, not a verification: the PN series is only meaningful
    while its terms are corrections. At f=20 Hz for a 20+15 Msun binary, v
    must be small enough that the 1PN, 1.5PN and 2PN terms stay at the tens
    of percent level rather than order unity. If this failed, nothing else
    in taylorf2.py would mean anything for this system at this frequency,
    however correctly it were implemented."""
    m1, m2 = 20.0, 15.0
    mtot = m1 + m2
    eta = m1 * m2 / mtot ** 2
    v = taylorf2._v_of_f(20.0, mtot)
    phi2 = 3715.0 / 756.0 + 55.0 * eta / 9.0
    phi3 = -16.0 * np.pi
    phi4 = 15293365.0 / 508032.0 + 27145.0 * eta / 504.0 + 3085.0 * eta ** 2 / 72.0
    terms = [float(phi2 * v ** 2), float(phi3 * v ** 3), float(phi4 * v ** 4)]
    ok = all(abs(t) < 1.0 for t in terms)
    return ok, f"v={v:.4f}; PN correction terms at f=20Hz: {terms} (all should be < 1 in magnitude)"


def check_spa_phase_series_matches_literature():
    """Exercises `spa_phase` itself, which the check above does not.

    Inverting the definition in its docstring,

        -Psi = 2*pi*f*t_c - phi_c - pi/4 + (3/(128*eta*v^5)) * P(v),

    recovers the PN polynomial P(v) FROM the function's own output; it is then
    compared against Buonanno et al. (2009), Eq. 3.18, written out here from
    the paper. One comparison therefore pins the overall sign, the
    3/(128*eta*v^5) prefactor, the v = (pi*M*f)^(1/3) definition and all three
    PN coefficients at once -- none of which any check reached before: with the
    leading-order check comparing a local copy against chirp.py, mutating
    phi2, phi3, the v^5 exponent or the amplitude exponent all passed, and so
    did doubling the entire leading-order phase.

    v is rebuilt here from units rather than taken from `_v_of_f`, so the two
    checks do not share a way to be wrong."""
    m1, m2 = 20.0, 15.0
    mtot = m1 + m2
    eta = m1 * m2 / mtot ** 2
    t_c = 15.0
    f = np.array([12.0, 20.0, 35.0, 60.0, 100.0, 160.0])
    v = (np.pi * units.msun_to_seconds(mtot) * f) ** (1.0 / 3.0)

    psi = taylorf2.spa_phase(f, m1, m2, t_c)
    recovered = (-psi - 2.0 * np.pi * f * t_c + np.pi / 4.0) * (128.0 * eta * v ** 5 / 3.0)

    phi2 = 3715.0 / 756.0 + 55.0 * eta / 9.0
    phi3 = -16.0 * np.pi
    phi4 = 15293365.0 / 508032.0 + 27145.0 * eta / 504.0 + 3085.0 * eta ** 2 / 72.0
    literature = 1.0 + phi2 * v ** 2 + phi3 * v ** 3 + phi4 * v ** 4

    relerr = float(np.max(np.abs(recovered - literature) / np.abs(literature)))
    ok = relerr < 1e-12
    return ok, (f"PN series recovered from spa_phase vs Buonanno et al. (2009) "
                f"Eq. 3.18 over f=12..160 Hz: max relative error {relerr:.2e} "
                f"(tol 1e-12)")


def check_spa_amplitude_scaling():
    """|h(f)| ~ f^(-7/6), the exponent spa_amplitude's own docstring states
    and that nothing checked: mutating it to f^(-1/6) passed every check in
    this file. A ratio needs no absolute normalization, so this tests the
    exponent alone."""
    m1, m2 = 20.0, 15.0
    f = np.array([10.0, 20.0, 40.0, 80.0])
    a = taylorf2.spa_amplitude(f, m1, m2, 100.0)
    ratios = a[1:] / a[:-1]
    expected = 2.0 ** (-7.0 / 6.0)
    relerr = float(np.max(np.abs(ratios - expected) / expected))
    ok = relerr < 1e-12
    return ok, (f"spa_amplitude doubling ratio A(2f)/A(f) = {ratios[0]:.9f} vs "
                f"2^(-7/6) = {expected:.9f} (max relative error {relerr:.2e}, "
                f"tol 1e-12)")


CHECKS = [
    ("leading_order_group_delay", check_leading_order_group_delay),
    ("ifft_reconstructs_chirp_at_tc", check_ifft_reconstructs_chirp_at_tc),
    ("pn_terms_are_small_corrections", check_pn_terms_are_small_corrections),
    ("spa_phase_series_matches_literature", check_spa_phase_series_matches_literature),
    ("spa_amplitude_scaling", check_spa_amplitude_scaling),
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
    out = Path(__file__).parent / "CHECKS_taylorf2.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
