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

from gwlens import system, chirp, imr_waveform, waveoptics as wo, geometry as geo, units

OUT = HERE
NUMBERS = {}


def log(name, value):
    NUMBERS[name] = value
    return value


def build_waveform_fd():
    """The unlensed waveform in the frequency domain, via
    src/gwlens/imr_waveform.py: real inspiral-merger-ringdown (IMRPhenomD,
    pycbc) where pycbc is importable, else the inspiral-only TaylorF2
    fallback (src/gwlens/taylorf2.py) -- see README.md for which
    environments give which, and wiki/log.md for why both exist. Lensing
    is a single complex multiplication by F(f) on this same grid -- the
    natural representation for it, since F(f) IS a frequency-domain object
    (Eq. Fdimless, theory.tex)."""
    f0 = system.F_A_START_HZ
    Mc = system.MCHIRP_MSUN
    t_c = chirp.time_to_merger(f0, Mc)
    tau_isco = chirp.time_to_merger(system.F_ISCO_HZ, Mc)
    t_end = t_c - tau_isco  # time at which the LEADING-ORDER f(t) = f_isco
    log("t_c_seconds", float(t_c))
    log("t_end_seconds", float(t_end))
    log("duration_seconds", float(t_end))

    # NOT 4*f_isco: a real IMR waveform's ringdown has significant power
    # well above f_isco (checked directly against pycbc's own IMRPhenomD
    # output for this system: >1% of peak amplitude out to ~500 Hz, >0.1%
    # to ~640 Hz) -- 4*f_isco=502 Hz undersamples that badly (caught by
    # eye: the lensed merger/ringdown zoom looked suspiciously smooth/flat
    # compared to the unlensed one, a classic aliasing symptom, not a
    # lensing effect -- see wiki/log.md). 2048 Hz gives >3x margin above
    # the observed high-frequency content.
    fs = 2048.0
    n = int(np.ceil(t_end * fs))
    n = 1 << (n - 1).bit_length()  # round up to a power of two
    n_pad = 2 * n  # extra padding: frequency resolution for the fringes (n already a power of two, so this equals 1 << (2*n-1).bit_length())
    freqs = np.fft.rfftfreq(n_pad, d=1.0 / fs)
    log("sample_rate_Hz", fs)
    log("n_samples", n_pad)
    log("frequency_resolution_Hz", float(freqs[1] - freqs[0]))

    D_eff_mpc = system.D_L_PC / 1.0e6  # same physical distance as the lens (~5 kpc)
    merger_time = 0.85 * (n_pad / fs)  # leaves inspiral before and ringdown after in-window
    H_unlensed, source_label = imr_waveform.get_unlensed_htilde_fd(
        freqs, system.M1_MSUN, system.M2_MSUN, t_c, D_eff_mpc, f0,
        merger_time=merger_time)
    log("waveform_source", source_label)
    log("merger_time_target_s", float(merger_time))

    return freqs, H_unlensed, fs, D_eff_mpc, n_pad


def apply_lensing(freqs, H_unlensed, y_A):
    w_arr = wo.w_of_frequency(np.abs(freqs), system.M_LENS_MSUN)
    F = np.ones_like(freqs, dtype=complex)
    positive = freqs > 0
    t0 = time.time()
    F_vals = np.array([wo.F_hybrid(w, y_A) for w in w_arr[positive]])
    log("F_evaluation_seconds", time.time() - t0)
    F[positive] = F_vals
    H_lensed = H_unlensed * F
    return H_lensed, F


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
    # F_hybrid switches evaluator on w alone, not on r -- at high w (as used
    # for the LEFT/full-pattern panel here) it always returns
    # F_geometric_optics, which formally diverges as r->0 (mu_+ -> infinity).
    # An independent review caught this: the plotted "Einstein-ring
    # diffraction spike" was actually that divergence, cut off only by the
    # innermost grid sample, not the true (finite) wave-optics value -- the
    # opposite of what the figure's caption claims. Force the exact closed
    # form near the axis regardless of w; see wiki/log.md.
    F_r = np.array([wo.F_point_lens(w, r) if r < 0.5 else wo.F_hybrid(w, r) for r in r_grid])
    mag2 = np.abs(F_r) ** 2
    interp = interp1d(r_grid, mag2, bounds_error=False, fill_value=(mag2[0], mag2[-1]))

    xs = np.linspace(x1, x2, n_grid)
    ys = np.linspace(y1, y2, n_grid)
    Y1, Y2 = np.meshgrid(xs, ys)
    R = np.hypot(Y1, Y2)  # radius from the LENS at the origin
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

    freqs, H_unlensed, fs, D_eff_mpc, n_pad = build_waveform_fd()
    log("D_eff_Mpc", D_eff_mpc)
    H_lensed, F = apply_lensing(freqs, H_unlensed, y_A)

    t = np.arange(n_pad) / fs
    h_unlensed = np.fft.irfft(H_unlensed, n=n_pad)
    h_lensed = np.fft.irfft(H_lensed, n=n_pad)
    log("peak_time_unlensed_s", float(t[np.argmax(np.abs(h_unlensed))]))

    w_start = wo.w_of_frequency(system.F_A_START_HZ, system.M_LENS_MSUN)
    w_isco = wo.w_of_frequency(system.F_ISCO_HZ, system.M_LENS_MSUN)
    log("w_at_f_start", float(w_start))
    log("w_at_f_isco", float(w_isco))

    F_at_start = wo.F_hybrid(w_start, y_A)
    F_at_isco = wo.F_hybrid(w_isco, y_A)
    log("abs_F_at_f_start", float(abs(F_at_start)))
    log("abs_F_at_f_isco", float(abs(F_at_isco)))

    band_mask = (freqs >= system.F_A_START_HZ) & (freqs <= system.F_ISCO_HZ)
    abs_F_band = np.abs(F[band_mask])
    log("min_abs_F_in_band", float(abs_F_band.min()))
    log("max_abs_F_in_band", float(abs_F_band.max()))

    i_peak_lensed = np.argmax(np.abs(h_lensed))
    mag_lensed = np.abs(h_lensed[i_peak_lensed])
    mag_unlensed = np.max(np.abs(h_unlensed))
    log("peak_strain_amplification", float(mag_lensed / mag_unlensed))
    log("peak_time_lensed_s", float(t[i_peak_lensed]))
    log("peak_time_lensed_minus_unlensed_s", float(t[i_peak_lensed] - NUMBERS["peak_time_unlensed_s"]))

    dT_dimensionless = wo.time_delay_difference(y_A)
    t_char_seconds = 4.0 * units.msun_to_seconds(system.M_LENS_MSUN)  # 4GM/c^3
    log("image_time_delay_seconds", float(dT_dimensionless * t_char_seconds))
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

    # ---- Genuine second-image echo: locate it empirically, don't assume --
    # F(w,y) is referenced to the strong image's own arrival (waveoptics.py
    # module docstring, "REFERENCE-PHASE NORMALIZATION FIX"), so the naive
    # prediction "echo at t_peak_unlensed + image_time_delay" should now be
    # exact -- search near it and report how well it lines up, rather than
    # assume: this is the more honest, and more easily reproduced, check.
    # (An earlier version of this code, before that normalization fix, found
    # this naive prediction off by ~0.6s and mis-attributed the offset to
    # interference between the two images; it was actually the un-subtracted
    # reference phase this fix removes -- see wiki/log.md.)
    from scipy.signal import hilbert
    env_u = np.abs(hilbert(h_unlensed))
    env_l = np.abs(hilbert(h_lensed))
    t_peak_unlensed = NUMBERS["peak_time_unlensed_s"]
    predicted_echo_t = t_peak_unlensed + NUMBERS["image_time_delay_seconds"]
    search = (t > predicted_echo_t - 1.0) & (t < predicted_echo_t + 1.0)
    echo_idx_local = np.argmax(env_l[search])
    echo_t = float(t[search][echo_idx_local])
    log("echo_peak_time_s", echo_t)
    log("echo_peak_time_minus_prediction_s", float(echo_t - predicted_echo_t))
    log("echo_peak_env_lensed", float(env_l[search][echo_idx_local]))
    # same instant on the unlensed envelope, as a noise-floor / sanity check
    # (the unlensed signal has already ended by here -- this should be small)
    i_at_echo_unlensed = int(np.argmin(np.abs(t - echo_t)))
    log("echo_time_env_unlensed", float(env_u[i_at_echo_unlensed]))

    # How much fainter than the MERGER PEAK itself (not the noise floor
    # above) is the echo -- the number that actually answers "how faint is
    # it", and the cleanest independent validation of F available: in the
    # geometric-optics limit (which the whole in-band F is, here -- see
    # RESULTS.md) the echo/peak envelope ratio should equal sqrt(|mu_-|),
    # the weak image's own magnification, exactly.
    unlensed_peak_env = float(np.max(env_u))
    log("unlensed_peak_env", unlensed_peak_env)
    echo_to_unlensed_peak_ratio = float(NUMBERS["echo_peak_env_lensed"] / unlensed_peak_env)
    log("echo_to_unlensed_peak_ratio", echo_to_unlensed_peak_ratio)
    x_plus_A, x_minus_A = wo.image_positions(y_A)
    sqrt_mu_minus = float(np.sqrt(np.abs(wo.magnification(x_minus_A))))
    sqrt_mu_plus = float(np.sqrt(wo.magnification(x_plus_A)))
    log("sqrt_mu_minus_geometric", sqrt_mu_minus)
    log("echo_ratio_vs_sqrt_mu_minus_relerr",
        float(abs(echo_to_unlensed_peak_ratio - sqrt_mu_minus) / sqrt_mu_minus))

    # ---- Figure 1: F(f) across the chirp (interference fringes) ---------
    # The fringes are far too fine to resolve over the full 10-125 Hz band
    # at print resolution (w up to 778 means thousands of oscillations) --
    # plotting the raw |F(f)| curve there just renders as an unresolved,
    # visually noisy blob with no information content beyond its envelope
    # (an earlier version of this figure did exactly that). The envelope
    # itself, though, is EXACTLY two constant lines here, not something
    # that needs a plot to discover: F is in the geometric-optics regime
    # throughout the band (min/max_abs_F_in_band above), so
    # |F(f)| = |sqrt(mu_+) + i*sqrt(|mu_-|)*exp(-2*pi*i*f*DeltaT_seconds)|
    # oscillates strictly between sqrt(mu_+)-sqrt(|mu_-|) and
    # sqrt(mu_+)+sqrt(|mu_-|) for every f -- so the top panel shows that
    # constant band directly (shaded), which is the actually-informative
    # summary of the whole-band behavior. The oscillation ITSELF has a
    # period in frequency of exactly 1/Delta_T_seconds
    # (image_time_delay_seconds), independent of f (Delta_T is fixed for a
    # static lens, and phase = 2*pi*f*Delta_T is linear in f) -- so the
    # zoomed inset below, though drawn at the start of the band for
    # concreteness, is representative of any equal-width window anywhere
    # in it, not particular to that location.
    mask = (freqs >= system.F_A_START_HZ) & (freqs <= system.F_ISCO_HZ)
    f_zoom_lo, f_zoom_hi = system.F_A_START_HZ, system.F_A_START_HZ + 3.0
    zmask = (freqs >= f_zoom_lo) & (freqs <= f_zoom_hi)
    fringe_period_hz = 1.0 / NUMBERS["image_time_delay_seconds"]
    log("fringe_period_hz", float(fringe_period_hz))

    # Three panels: constant-envelope overview, |F(f)| zoomed (fringes
    # resolved), and arg F(f) over that same zoom -- an earlier version of
    # this figure had only the two |F| panels while RESULTS.md/claims.yaml
    # already described an "arg F" sub-panel that did not exist; an
    # independent review caught the mismatch, see wiki/log.md.
    fig, axes = plt.subplots(3, 1, figsize=(7, 9.2))
    axes[0].axhspan(sqrt_mu_plus - sqrt_mu_minus, sqrt_mu_plus + sqrt_mu_minus,
                     color="#2b6cb0", alpha=0.25,
                     label=r"$\sqrt{\mu_+}\pm\sqrt{|\mu_-|}$ (exact envelope)")
    axes[0].axhline(sqrt_mu_plus, color="#2b6cb0", lw=1.0, ls="--",
                     label=r"$\sqrt{\mu_+}$ (strong-image-only value)")
    axes[0].axvspan(f_zoom_lo, f_zoom_hi, color="gold", alpha=0.3)
    axes[0].set_ylabel(r"$|F(f)|$")
    axes[0].set_xlabel("f [Hz]")
    axes[0].set_xlim(system.F_A_START_HZ, system.F_ISCO_HZ)
    axes[0].legend(fontsize=7, loc="lower right")
    axes[0].set_title("Case A: amplification factor across the whole chirp\n"
                       f"(static lens, $y={y_A:.3f}$; band is exactly this width at every f; zoom below)",
                       fontsize=10)
    axes[1].plot(freqs[zmask], np.abs(F[zmask]), lw=1.0, color="#2b6cb0")
    axes[1].set_ylabel(r"$|F(f)|$")
    axes[1].set_xlabel("f [Hz]")
    axes[1].set_title(f"Zoom: {f_zoom_lo:.0f}-{f_zoom_hi:.0f} Hz (fringe period "
                       f"{fringe_period_hz:.2f} Hz, the same anywhere in the band)")
    axes[2].plot(freqs[zmask], np.angle(F[zmask]), lw=1.0, color="#c05621")
    axes[2].set_ylabel(r"$\arg F(f)$ [rad]")
    axes[2].set_xlabel("f [Hz]")
    # F traces a circle of radius sqrt(mu_-) centered on sqrt(mu_+) in the
    # complex plane as w*DeltaT sweeps 2*pi per fringe (the strong image is
    # the fixed reference, weak image the rotating arm) -- since mu_- < mu_+
    # here, that circle does not enclose the origin, so arg F oscillates
    # (bounded, not a full 2*pi wrap) once per fringe, in phase with |F|
    # above.
    axes[2].set_title("Same zoom: phase -- bounded oscillation in phase with each fringe above")
    fig.tight_layout()
    fig.savefig(OUT / "caseA_F_of_f.png", dpi=170)
    plt.close(fig)

    # ---- Figure 2: lensed vs unlensed time-domain strain -----------------
    # Three panels: overview (extended past merger to include the
    # second-image echo), a merger/ringdown zoom, and a dedicated echo zoom
    # -- the echo is the unambiguous causality signature (a second, weaker,
    # LATER copy of the signal). With F referenced to the strong image's own
    # arrival (waveoptics.py, "REFERENCE-PHASE NORMALIZATION FIX"), the
    # lensed merger now peaks at (to within a sample of) the SAME time as
    # the unlensed one -- see wiki/log.md for the earlier, incorrect version
    # of this story, where an un-subtracted reference phase looked like a
    # genuine 0.6s early shift and was mis-attributed to interference.
    t_peak = t_peak_unlensed
    fig, axes = plt.subplots(3, 1, figsize=(8, 9.0))
    t_hi = min(t[-1], echo_t + 1.0)
    view = (t >= t_peak - 1.15 * NUMBERS["t_end_seconds"]) & (t <= t_hi)
    axes[0].plot(t[view], h_unlensed[view], lw=0.5, color="#888888", label="unlensed")
    axes[0].plot(t[view], h_lensed[view], lw=0.5, color="#2b6cb0", label="lensed", alpha=0.85)
    axes[0].set_xlabel(f"t [s] (merger at t={t_peak:.2f} s)")
    axes[0].set_ylabel("h(t) [arb. units]")
    axes[0].axvspan(t_peak - 0.05, t_peak + 0.3, color="gold", alpha=0.3)
    axes[0].axvspan(echo_t - 0.3, echo_t + 0.3, color="mediumseagreen", alpha=0.3)
    axes[0].set_title(f"Case A: lensed vs. unlensed waveform\n({NUMBERS['waveform_source']})", fontsize=10)
    axes[0].legend(loc="upper left", fontsize=8)

    # window wide enough to comfortably contain BOTH the unlensed peak (at
    # t_peak, by definition) and the lensed curve's own peak, which is not
    # guaranteed to sit at t_peak too -- F(f) dephases the lensed merger,
    # so its extremum can be offset from the unlensed one (an earlier,
    # narrower window here missed it; see wiki/log.md)
    zoom = (t > t_peak - 0.05) & (t < t_peak + 0.3)
    axes[1].plot(t[zoom] - t_peak, h_unlensed[zoom], lw=1.0, color="#888888", label="unlensed")
    axes[1].plot(t[zoom] - t_peak, h_lensed[zoom], lw=1.0, color="#2b6cb0", alpha=0.85, label="lensed")
    axes[1].set_xlabel("t - t_merger [s]")
    axes[1].set_ylabel("h(t) [arb. units]")
    axes[1].set_title("Zoom: merger and ringdown (gold band above)", fontsize=10)
    axes[1].legend(loc="upper right", fontsize=8)

    echo_zoom = (t > echo_t - 0.3) & (t < echo_t + 0.3)
    axes[2].plot(t[echo_zoom] - t_peak, h_unlensed[echo_zoom], lw=1.0, color="#888888", label="unlensed (should be ~0: signal already ended)")
    axes[2].plot(t[echo_zoom] - t_peak, h_lensed[echo_zoom], lw=1.0, color="#2ca02c", alpha=0.9, label="lensed (the echo)")
    axes[2].set_xlabel("t - t_merger [s]")
    axes[2].set_ylabel("h(t) [arb. units]")
    axes[2].set_title(
        f"Zoom: the second-image echo, observed at t_merger+{echo_t - t_peak:.2f} s "
        f"(geometric-optics prediction: +{NUMBERS['image_time_delay_seconds']:.2f} s; "
        f"green band above) -- the causality check made visible", fontsize=10)
    axes[2].legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "caseA_strain_time.png", dpi=170)
    plt.close(fig)

    # ---- Figure 3: extended diffraction/interference pattern in y-plane -
    from matplotlib.colors import LogNorm
    import matplotlib.ticker as mticker

    def _log_tick(value, _pos):
        # "1" instead of matplotlib's default "10^0" for the unit-value
        # tick -- correct either way, but "10^0" reads as unnecessary
        # decoration for a tick that just means 1.
        if abs(value - 1.0) < 1e-9:
            return "1"
        exp = round(np.log10(value))
        return fr"$10^{{{exp}}}$"

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4))

    axis0, axis0y, mag2_0 = make_ring_pattern(w_start, y_max=3.0, n_grid=481)
    im0 = axes[0].pcolormesh(axis0, axis0y, mag2_0, shading="auto", cmap="inferno",
                              norm=LogNorm(vmin=max(mag2_0.min(), 1e-3), vmax=mag2_0.max()))
    axes[0].set_title(f"start of band, f={system.F_A_START_HZ} Hz, w={w_start:.1f}", fontsize=9)
    cb0 = fig.colorbar(im0, ax=axes[0], shrink=0.8, label=r"$|F|^2$ (log scale)")
    cb0.ax.yaxis.set_major_formatter(mticker.FuncFormatter(_log_tick))

    axis1, axis1y, mag2_1 = make_ring_pattern(
        w_isco, n_grid=500, center=(y_A, 0.0), half_width=0.15, n_radial=6000)
    im1 = axes[1].pcolormesh(axis1, axis1y, mag2_1, shading="auto", cmap="inferno",
                              norm=LogNorm(vmin=max(mag2_1.min(), 1e-3), vmax=mag2_1.max()))
    axes[1].set_title(f"near merger, f=f_isco, w={w_isco:.1f}\n(zoomed to |Δy|<0.15 near the source)", fontsize=9)
    cb1 = fig.colorbar(im1, ax=axes[1], shrink=0.8, label=r"$|F|^2$ (log scale)")
    cb1.ax.yaxis.set_major_formatter(mticker.FuncFormatter(_log_tick))

    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel(r"$y_1$")
        ax.plot(y_A, 0, "w+", ms=12, mew=2)
    axes[0].set_ylabel(r"$y_2$")
    # A wider, single-line-per-row suptitle was getting clipped at both
    # edges of this fairly narrow figure; shortened and set explicitly
    # inside the axes bounding box (independent review, see wiki/log.md).
    fig.suptitle("Case A: amplification pattern on the source plane\n"
                 "(axisymmetric rings; central spot = Einstein-ring "
                 fr"diffraction spike; source at $y_A={y_A:.3f}$, marked)",
                 fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90))
    fig.savefig(OUT / "caseA_diffraction_pattern.png", dpi=170)
    plt.close(fig)

    # ---- Figure 4: simple idealized detector view -------------------------
    # Log-scale y-axis and a window extended past the echo: the echo is
    # ~4x (0.6 dex) below the merger peak envelope
    # (echo_to_unlensed_peak_ratio, matching sqrt(|mu_-|) to <1%, see above)
    # -- clear on a log axis, easy to under-sell on a linear one.
    fig, ax = plt.subplots(figsize=(8, 3.6))
    # same window as Figure 2's overview panel above -- identical expression, reused rather than rebuilt
    ax.plot(t[view], env_u[view], color="#888888", lw=1.0, label="unlensed envelope")
    ax.plot(t[view], env_l[view], color="#2b6cb0", lw=1.0, label="lensed envelope")
    ax.axvspan(echo_t - 0.3, echo_t + 0.3, color="mediumseagreen", alpha=0.25,
               label=f"observed echo (t_merger+{echo_t - t_peak:.2f} s)")
    ax.set_yscale("log")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("strain envelope [arb. units, log scale]")
    ax.set_title("Case A: idealized detector view (strain envelope, no noise/antenna pattern)\n"
                 "note the second, weaker peak: the second-image echo", fontsize=10)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "caseA_detector_envelope.png", dpi=170)
    plt.close(fig)

    # ---- Animation data for report/report.html (JSON, reduced resolution) -
    win = (t >= t_peak - 1.15 * NUMBERS["t_end_seconds"]) & (t <= t_peak + 0.3)
    idx = np.where(win)[0]
    stride = max(1, len(idx) // 1400)
    idx = idx[::stride]
    f_isco = system.F_ISCO_HZ
    # (instantaneous frequency during inspiral, leading-order estimate, for
    # syncing the animation marker only -- not a physics claim; simpler to
    # recompute directly than thread through build_waveform_fd's internals)
    tau_isco = chirp.time_to_merger(system.F_ISCO_HZ, system.MCHIRP_MSUN)
    t_c_local = t_peak + tau_isco
    with np.errstate(invalid="ignore"):
        f_inspiral = chirp.freq_of_time(t[idx], t_c_local, system.MCHIRP_MSUN)
    f_of_t = np.where(t[idx] < t_peak, np.minimum(f_inspiral, f_isco), f_isco)

    zmask_idx = np.where(zmask)[0]
    anim = {
        "t_peak": float(t_peak),
        "t_end": float(NUMBERS["t_end_seconds"]),
        "f_isco": float(f_isco),
        "f_zoom_lo": float(f_zoom_lo),
        "f_zoom_hi": float(f_zoom_hi),
        "waveform": {
            "t": t[idx].tolist(),
            "h_unlensed": h_unlensed[idx].tolist(),
            "h_lensed": h_lensed[idx].tolist(),
            "env_unlensed": env_u[idx].tolist(),
            "env_lensed": env_l[idx].tolist(),
            "f_of_t": f_of_t.tolist(),
        },
        "fringes": {
            "freqs": freqs[zmask_idx].tolist(),
            "F_abs": np.abs(F[zmask_idx]).tolist(),
        },
    }
    with open(OUT / "caseA_animation_data.json", "w") as fh:
        json.dump(anim, fh)

    with open(OUT / "provenance" / "numbers.json", "w") as fh:
        json.dump(NUMBERS, fh, indent=2)
    print(json.dumps(NUMBERS, indent=2))
    print("\nFigures written:",
          "caseA_F_of_f.png, caseA_strain_time.png, "
          "caseA_diffraction_pattern.png, caseA_detector_envelope.png")


if __name__ == "__main__":
    main()
