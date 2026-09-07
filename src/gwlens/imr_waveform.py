"""The unlensed inner-binary waveform, full inspiral-merger-ringdown (IMR)
where possible.

Tries `pycbc` (IMRPhenomD, Khan et al. 2016 -- an NR-calibrated, published
approximant, exactly what LIGO/Virgo parameter estimation and search
pipelines use) first; if `pycbc` is not importable in the current
environment, falls back to `src/gwlens/taylorf2.py` (the hand-built,
restricted 2PN inspiral-ONLY SPA waveform already validated in
tests/test_taylorf2.py) -- no merger, no ringdown, cut off at f_isco.

Why both exist, not just one: `pycbc`/`lalsimulation` has real, historically
poor native-Windows support (see wiki/log.md -- a first `pip install pycbc`
attempt on native Windows hung indefinitely resolving the lalsuite
dependency chain). It installs cleanly on Linux/macOS, INCLUDING WSL2,
which is where this repository's own committed IMRPhenomD results were
actually generated (see README.md for exactly how). A plain
`pip install -r requirements.txt` on native Windows will NOT have `pycbc`
-- that is expected, not an error, and this module's automatic fallback
keeps every script runnable there too, just with the inspiral-only
waveform instead of the full IMR one.

Which path actually ran is always recorded (`waveform_source` in the
calling script's `numbers.json`) -- never silently assumed.
"""
import numpy as np

from . import taylorf2


def get_unlensed_htilde_fd(freqs, m1_msun, m2_msun, t_c, d_eff_mpc, f_lower,
                            merger_time=None):
    """H(f) on EXACTLY the input `freqs` grid (must be freqs = k*delta_f,
    k=0,1,2,..., numpy.fft.rfftfreq's convention). `merger_time`, if given,
    places the merger at that time (seconds) after an `irfft` of the
    returned array on this same grid -- needed because pycbc's FD
    waveforms are NOT delivered pre-aligned to a convenient positive time
    (see wiki/log.md for how this was diagnosed, the same symptom as the
    earlier TaylorF2 sign bug: an unshifted `irfft` put the merger right at
    the edge of the array, wrapping most of the inspiral to the "wrong"
    end). Returns (H, source_label).
    """
    try:
        H, label = _pycbc_imrphenomd(freqs, m1_msun, m2_msun, d_eff_mpc, f_lower)
    except ImportError:
        H = np.zeros_like(freqs, dtype=complex)
        H[:] = _taylorf2_banded(freqs, m1_msun, m2_msun, t_c, d_eff_mpc, f_lower)
        label = ("TaylorF2 2PN inspiral-only (pycbc not available in this "
                 "environment -- fallback, see wiki/log.md; no merger/ringdown)")

    if merger_time is not None:
        n_pad = 2 * (len(freqs) - 1)
        fs = freqs[-1] * 2.0
        t = np.arange(n_pad) / fs
        h0 = np.fft.irfft(H, n=n_pad)
        current_peak_t = t[np.argmax(np.abs(h0))]
        delta = merger_time - current_peak_t
        H = H * np.exp(-1j * 2.0 * np.pi * freqs * delta)

    return H, label


def _pycbc_imrphenomd(freqs, m1_msun, m2_msun, d_eff_mpc, f_lower):
    from pycbc.waveform import get_fd_waveform  # noqa: local import, optional dependency

    delta_f = float(freqs[1] - freqs[0])
    n = len(freqs)
    hp, _ = get_fd_waveform(
        approximant="IMRPhenomD",
        mass1=float(m1_msun), mass2=float(m2_msun),
        delta_f=delta_f, f_lower=float(f_lower),
        f_final=float(freqs[-1]),
        distance=float(d_eff_mpc),
    )
    src = hp.numpy()
    H = np.zeros(n, dtype=complex)
    m = min(n, len(src))
    H[:m] = src[:m]
    return H, "IMRPhenomD (pycbc/LALSimulation, Khan et al. 2016) -- full inspiral-merger-ringdown"


def _taylorf2_banded(freqs, m1_msun, m2_msun, t_c, d_eff_mpc, f_lower):
    """TaylorF2, banded [f_lower, f_isco] with a half-cosine edge taper --
    same construction Case A used before the pycbc integration."""
    from . import chirp
    f_isco = chirp.isco_frequency(m1_msun + m2_msun)
    taper_frac = 0.03
    band_width = f_isco - f_lower
    taper_hz = taper_frac * band_width
    w = np.zeros_like(freqs)
    inband = (freqs >= f_lower) & (freqs <= f_isco)
    w[inband] = 1.0
    lo_ramp = (freqs >= f_lower) & (freqs < f_lower + taper_hz)
    hi_ramp = (freqs > f_isco - taper_hz) & (freqs <= f_isco)
    w[lo_ramp] = 0.5 * (1 - np.cos(np.pi * (freqs[lo_ramp] - f_lower) / taper_hz))
    w[hi_ramp] = 0.5 * (1 - np.cos(np.pi * (f_isco - freqs[hi_ramp]) / taper_hz))
    H = np.zeros_like(freqs, dtype=complex)
    H[inband] = w[inband] * taylorf2.htilde(freqs[inband], m1_msun, m2_msun, t_c, d_eff_mpc)
    return H
