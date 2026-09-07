"""Case A -- the chirp, static lens.

The inner binary (src/gwlens/system.py) inspirals from f=F_A_START_HZ to
f_isco in ~15 s: 4.4e-5 of the outer orbital period, so the lens is frozen
at a single impact parameter y_A for the whole observation (checked in
tests/test_system.py::check_static_lens_regime).

Produces (in this directory, all from THIS script, nothing hand-edited):
    numbers -> provenance/numbers.json
    figures -> *.png
    this docstring's claims -> provenance/claims.yaml (written by hand,
        pointing at the numbers/figures this script produces)

Run: `python3 cases/case_A_chirp/run.py`
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from gwlens import system, chirp, waveoptics as wo, geometry as geo, units

OUT = HERE
NUMBERS = {}


def log(name, value):
    NUMBERS[name] = value
    return value


def build_waveform():
    """Time-domain restricted-quadrupole chirp from f=F_A_START_HZ to
    f_isco, sampled fast enough to Nyquist-resolve f_isco, windowed with a
    Tukey taper to control FFT edge artefacts."""
    f0 = system.F_A_START_HZ
    Mc = system.MCHIRP_MSUN
    t_c = chirp.time_to_merger(f0, Mc)
    tau_isco = chirp.time_to_merger(system.F_ISCO_HZ, Mc)
    t_end = t_c - tau_isco  # time at which f(t) = f_isco
    log("t_c_seconds", float(t_c))
    log("t_end_seconds", float(t_end))
    log("duration_seconds", float(t_end))

    fs = 4.0 * system.F_ISCO_HZ  # >2x Nyquist margin
    n = int(np.ceil(t_end * fs))
    n = 1 << (n - 1).bit_length()  # next power of 2, for fft speed
    t = np.arange(n) / fs
    t = t[t < t_end]

    f_t = chirp.freq_of_time(t, t_c, Mc)
    _, phase = chirp.phase_of_time(t, lambda tt: chirp.freq_of_time(tt, t_c, Mc), 0.0, t[-1])
    # phase_of_time returns its own (finer) time grid; interpolate onto t
    tt_fine, _ = chirp.phase_of_time(t, lambda tt: chirp.freq_of_time(tt, t_c, Mc), 0.0, t[-1])
    phase_interp = interp1d(tt_fine, phase, kind="cubic")(t)

    D_eff_mpc = system.D_L_PC / 1.0e6  # same physical distance as the lens (~5 kpc)
    amp = chirp.restricted_pn_amplitude_td(f_t, Mc, D_eff_mpc)

    window = np.ones_like(t)
    taper_n = int(0.05 * len(t))
    if taper_n > 1:
        ramp = 0.5 * (1 - np.cos(np.linspace(0, np.pi, taper_n)))
        window[:taper_n] = ramp
    h_t = amp * window * np.cos(phase_interp)
    return t, f_t, h_t, fs, D_eff_mpc


def apply_lensing(t, h_t, fs, y_A):
    n = len(h_t)
    n_pad = 1 << (2 * n - 1).bit_length()  # zero-pad to reduce circular-conv wraparound
    h_padded = np.zeros(n_pad)
    h_padded[:n] = h_t
    H = np.fft.rfft(h_padded)
    freqs = np.fft.rfftfreq(n_pad, d=1.0 / fs)

    w_arr = wo.w_of_frequency(np.abs(freqs), system.M_LENS_MSUN)
    # F at f=0 is ill-defined (w=0 handled by low-w limit -> 1); skip it
    F = np.ones_like(H, dtype=complex)
    positive = freqs > 0
    t0 = time.time()
    F_vals = np.array([wo.F_hybrid(w, y_A) for w in w_arr[positive]])
    log("F_evaluation_seconds", time.time() - t0)
    F[positive] = F_vals

    H_lensed = H * F
    h_lensed_padded = np.fft.irfft(H_lensed, n=n_pad)
    h_lensed = h_lensed_padded[:n]
    return h_lensed, freqs, F, H


def make_ring_pattern(w, y_max=3.0, n_grid=481, n_radial=4000,
                       center=(0.0, 0.0), half_width=None):
    """|F(w,y)|^2 on a 2D (y1,y2) grid, built from a 1D radial evaluation
    (point lens is axisymmetric) and interpolated onto the grid -- avoids
    ~n_grid^2 individual mpmath/geometric-optics calls.

    Two resolutions matter, both must resolve the finest fringe spacing
    ~2*pi/(w*y) present in the region shown: `n_radial` (the underlying 1D
    F(w,r) sampling) AND the DISPLAY grid `n_grid` over the plotted window
    -- a first version of the zoomed high-w panel had a fine `n_radial` but
    still looked aliased because `n_grid` points were spread over the full
    [-y_max,y_max] axis rather than the actual (small) zoomed window; fixed
    by letting the caller pass `center`/`half_width` so the display grid is
    only ever as coarse as the window actually shown (see wiki/log.md).
    """
    cx, cy = center
    if half_width is None:
        half_width = y_max
        x1, x2 = -y_max, y_max
        y1, y2 = -y_max, y_max
    else:
        x1, x2 = cx - half_width, cx + half_width
        y1, y2 = cy - half_width, cy + half_width
    r_max_needed = np.hypot(max(abs(x1), abs(x2)), max(abs(y1), abs(y2)))

    r_grid = np.linspace(1e-3, r_max_needed, n_radial)
    F_r = np.array([wo.F_hybrid(w, r) for r in r_grid])
    mag2 = np.abs(F_r) ** 2
    interp = interp1d(r_grid, mag2, bounds_error=False, fill_value=(mag2[0], mag2[-1]))

    xs = np.linspace(x1, x2, n_grid)
    ys = np.linspace(y1, y2, n_grid)
    Y1, Y2 = np.meshgrid(xs, ys)
    R = np.hypot(Y1 - 0.0, Y2 - 0.0)  # radius from the LENS at the origin
    return xs, ys, interp(R)


def main():
    print(system.summary())

    # y_A: the impact parameter frozen at the moment of merger. The merger
    # can in principle happen at any outer-orbital phase; we illustrate the
    # phase of CLOSEST approach (deepest, most visible lensing) rather than
    # an arbitrary t=0 -- found by scanning one full outer period and taking
    # the minimum y among the lensed (D_LS>0) half. Explicit choice, not a
    # fit; a different merger phase would just move y_A within the pattern
    # plotted in Figure 3 below.
    t_scan = np.linspace(0.0, system.P_OUT_S, 2000, endpoint=False)
    y_scan, lensed_scan, _ = geo.impact_parameter_of_time(
        t_scan, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
        system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)
    idx_min = np.argmin(np.where(lensed_scan, y_scan, np.inf))
    y_A = float(y_scan[idx_min])
    log("y_A", y_A)
    log("y_A_outer_orbit_phase_fraction", float(t_scan[idx_min] / system.P_OUT_S))
    log("lensed_at_merger_phase", bool(lensed_scan[idx_min]))

    t, f_t, h_unlensed, fs, D_eff_mpc = build_waveform()
    log("D_eff_Mpc", D_eff_mpc)
    log("sample_rate_Hz", fs)
    log("n_samples", len(t))

    h_lensed, freqs, F, H = apply_lensing(t, h_unlensed, fs, y_A)

    w_start = wo.w_of_frequency(system.F_A_START_HZ, system.M_LENS_MSUN)
    w_isco = wo.w_of_frequency(system.F_ISCO_HZ, system.M_LENS_MSUN)
    log("w_at_f_start", float(w_start))
    log("w_at_f_isco", float(w_isco))

    F_at_start = wo.F_hybrid(w_start, y_A)
    F_at_isco = wo.F_hybrid(w_isco, y_A)
    log("abs_F_at_f_start", float(abs(F_at_start)))
    log("abs_F_at_f_isco", float(abs(F_at_isco)))

    mag_lensed = np.max(np.abs(h_lensed))
    mag_unlensed = np.max(np.abs(h_unlensed))
    log("peak_strain_amplification", float(mag_lensed / mag_unlensed))
    # Peak-sample amplification is a fragile statistic here: F(f) modulates
    # PHASE as well as amplitude across the band (arg F sweeps through many
    # cycles, Figure 1), so the lensed waveform is not simply a rescaled
    # copy of the unlensed one -- individual peaks can come out lower even
    # with |F|>1 everywhere sampled, because the coherent buildup at any one
    # instant is now built from dephased harmonics. RMS power, integrated
    # over the whole chirp, is the more robust summary:
    rms_lensed = np.sqrt(np.mean(h_lensed ** 2))
    rms_unlensed = np.sqrt(np.mean(h_unlensed ** 2))
    log("rms_strain_amplification", float(rms_lensed / rms_unlensed))

    # ---- Figure 1: F(f) across the chirp (interference fringes) ---------
    # The fringes are far too fine to resolve over the full 10-125 Hz band
    # at print resolution (w up to 778 means thousands of oscillations); a
    # first version of this figure just looked like a filled band. Fixed
    # with a zoomed inset showing a few Hz where the fringes are resolved.
    mask = (freqs >= system.F_A_START_HZ) & (freqs <= system.F_ISCO_HZ)
    f_zoom_lo, f_zoom_hi = system.F_A_START_HZ, system.F_A_START_HZ + 3.0
    zmask = (freqs >= f_zoom_lo) & (freqs <= f_zoom_hi)

    fig, axes = plt.subplots(2, 1, figsize=(7, 6.5))
    axes[0].plot(freqs[mask], np.abs(F[mask]), lw=0.4, color="#2b6cb0")
    axes[0].axvspan(f_zoom_lo, f_zoom_hi, color="gold", alpha=0.3)
    axes[0].set_ylabel(r"$|F(f)|$")
    axes[0].set_xlabel("f [Hz]")
    axes[0].set_title("Case A: amplification factor across the whole chirp\n"
                       f"(static lens, $y={y_A:.3f}$; shaded = zoom below)")
    axes[1].plot(freqs[zmask], np.abs(F[zmask]), lw=1.0, color="#2b6cb0")
    axes[1].set_ylabel(r"$|F(f)|$")
    axes[1].set_xlabel("f [Hz]")
    axes[1].set_title(f"Zoom: {f_zoom_lo:.0f}-{f_zoom_hi:.0f} Hz -- individual "
                       "interference fringes resolved")
    fig.tight_layout()
    fig.savefig(OUT / "caseA_F_of_f.png", dpi=170)
    plt.close(fig)

    # ---- Figure 2: lensed vs unlensed time-domain strain -----------------
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.plot(t, h_unlensed, lw=0.6, color="#888888", label="unlensed")
    ax.plot(t, h_lensed, lw=0.6, color="#2b6cb0", label="lensed", alpha=0.85)
    ax.set_xlabel(f"t [s] (merger at t={NUMBERS['t_end_seconds']:.2f} s)")
    ax.set_ylabel("h(t) [arb. units]")
    ax.set_title("Case A: lensed vs. unlensed chirp")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "caseA_strain_time.png", dpi=170)
    plt.close(fig)

    # ---- Figure 3: extended diffraction/interference pattern in y-plane -
    from matplotlib.colors import LogNorm
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4))

    axis0, axis0y, mag2_0 = make_ring_pattern(w_start, y_max=3.0, n_grid=481)
    im0 = axes[0].pcolormesh(axis0, axis0y, mag2_0, shading="auto", cmap="inferno",
                              norm=LogNorm(vmin=max(mag2_0.min(), 1e-3), vmax=mag2_0.max()))
    axes[0].set_title(f"start of band, f={system.F_A_START_HZ} Hz, w={w_start:.1f}", fontsize=9)
    fig.colorbar(im0, ax=axes[0], shrink=0.8, label=r"$|F|^2$ (log scale)")

    axis1, axis1y, mag2_1 = make_ring_pattern(
        w_isco, n_grid=500, center=(y_A, 0.0), half_width=0.15, n_radial=6000)
    im1 = axes[1].pcolormesh(axis1, axis1y, mag2_1, shading="auto", cmap="inferno",
                              norm=LogNorm(vmin=max(mag2_1.min(), 1e-3), vmax=mag2_1.max()))
    axes[1].set_title(f"near merger, f=f_isco, w={w_isco:.1f}\n(zoomed to |Δy|<0.15 near the source)", fontsize=9)
    fig.colorbar(im1, ax=axes[1], shrink=0.8, label=r"$|F|^2$ (log scale)")

    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel(r"$y_1$")
        ax.plot(y_A, 0, "w+", ms=12, mew=2)
    axes[0].set_ylabel(r"$y_2$")
    fig.suptitle("Case A: extended amplification pattern on the source plane\n"
                 "(point lens: axisymmetric rings; central bright spot = "
                 "Einstein-ring diffraction spike; source at fixed "
                 fr"$y_A={y_A:.3f}$, marked)")
    fig.tight_layout()
    fig.savefig(OUT / "caseA_diffraction_pattern.png", dpi=170)
    plt.close(fig)

    # ---- Figure 4: simple idealized detector view -------------------------
    fig, ax = plt.subplots(figsize=(8, 3.0))
    from scipy.signal import hilbert
    env_u = np.abs(hilbert(h_unlensed))
    env_l = np.abs(hilbert(h_lensed))
    ax.plot(t, env_u, color="#888888", lw=1.0, label="unlensed envelope")
    ax.plot(t, env_l, color="#2b6cb0", lw=1.0, label="lensed envelope")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("strain envelope [arb. units]")
    ax.set_title("Case A: idealized detector view (strain envelope, no noise/antenna pattern)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "caseA_detector_envelope.png", dpi=170)
    plt.close(fig)

    with open(OUT / "provenance" / "numbers.json", "w") as fh:
        json.dump(NUMBERS, fh, indent=2)
    print(json.dumps(NUMBERS, indent=2))
    print("\nFigures written:",
          "caseA_F_of_f.png, caseA_strain_time.png, "
          "caseA_diffraction_pattern.png, caseA_detector_envelope.png")


if __name__ == "__main__":
    main()
