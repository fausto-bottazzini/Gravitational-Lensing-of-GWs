"""Checks for src/gwlens/imr_waveform.py. Run with `python3 tests/test_imr_waveform.py`.

Some checks only run if `pycbc` is importable in the current environment
(WSL2/Linux here; not native Windows -- see README.md and wiki/log.md).
Skipped checks are reported as such, not silently passed.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from gwlens import imr_waveform, chirp


def _pycbc_available():
    try:
        import pycbc  # noqa: F401
        return True
    except ImportError:
        return False


def check_fallback_matches_taylorf2_directly():
    """_taylorf2_banded (the fallback path) must reproduce taylorf2.htilde
    exactly in-band -- it is a thin wrapper (band + taper) around it, not a
    reimplementation."""
    from gwlens import taylorf2
    freqs = np.linspace(0, 200, 4001)
    m1, m2 = 20.0, 15.0
    Mc = chirp.chirp_mass(m1, m2)
    t_c = chirp.time_to_merger(10.0, Mc)
    H = imr_waveform._taylorf2_banded(freqs, m1, m2, t_c, 5.0e-3, 10.0)
    f_isco = chirp.isco_frequency(m1 + m2)
    mid = (freqs > 30) & (freqs < f_isco - 20)  # away from the taper edges
    H_direct = taylorf2.htilde(freqs[mid], m1, m2, t_c, 5.0e-3)
    relerr = np.max(np.abs(H[mid] - H_direct) / np.abs(H_direct))
    ok = relerr < 1e-10
    return ok, f"max relative error, banded fallback vs direct taylorf2.htilde (mid-band, untapered) = {relerr:.2e}"


def check_fallback_used_without_pycbc():
    """get_unlensed_htilde_fd must fall back cleanly (no exception) and
    label itself correctly when pycbc genuinely is not importable -- forced
    here by a bogus module name substitution rather than assuming this
    environment's actual pycbc availability either way."""
    import builtins
    real_import = builtins.__import__

    def blocked_import(name, *a, **kw):
        if name == "pycbc" or name.startswith("pycbc."):
            raise ImportError("blocked for this check")
        return real_import(name, *a, **kw)

    freqs = np.linspace(0, 200, 4001)
    m1, m2 = 20.0, 15.0
    Mc = chirp.chirp_mass(m1, m2)
    t_c = chirp.time_to_merger(10.0, Mc)
    builtins.__import__ = blocked_import
    try:
        H, label = imr_waveform.get_unlensed_htilde_fd(freqs, m1, m2, t_c, 5.0e-3, 10.0)
    finally:
        builtins.__import__ = real_import
    ok = "TaylorF2" in label and "fallback" in label and np.any(np.abs(H) > 0)
    return ok, f"label='{label}', nonzero bins={np.sum(np.abs(H)>0)}"


def check_pycbc_merger_time_alignment():
    """When pycbc IS available: the time-shifted IMRPhenomD waveform's
    peak must land within a tiny fraction of a cycle of the requested
    merger_time, and -- the actual regression this guards -- must NOT sit
    at the edge of the FFT window (the unshifted-alignment bug)."""
    if not _pycbc_available():
        return True, "SKIPPED: pycbc not importable in this environment"
    from gwlens import system
    fs = 4.0 * system.F_ISCO_HZ
    n_pad = 16384
    T = n_pad / fs
    freqs = np.fft.rfftfreq(n_pad, d=1.0 / fs)
    target = 0.85 * T
    H, label = imr_waveform.get_unlensed_htilde_fd(
        freqs, system.M1_MSUN, system.M2_MSUN, None, system.D_L_PC / 1e6,
        system.F_A_START_HZ, merger_time=target)
    h = np.fft.irfft(H, n=n_pad)
    t = np.arange(n_pad) / fs
    peak_t = t[np.argmax(np.abs(h))]
    err = abs(peak_t - target)
    not_at_edge = peak_t < 0.95 * T
    ok = err < 0.05 and not_at_edge and "IMRPhenomD" in label
    return ok, f"target={target:.3f}s, peak={peak_t:.3f}s, |diff|={err:.4f}s, not_at_window_edge={not_at_edge}"


CHECKS = [
    ("fallback_matches_taylorf2_directly", check_fallback_matches_taylorf2_directly),
    ("fallback_used_without_pycbc", check_fallback_used_without_pycbc),
    ("pycbc_merger_time_alignment", check_pycbc_merger_time_alignment),
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
    out = Path(__file__).parent / "CHECKS_imr_waveform.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
