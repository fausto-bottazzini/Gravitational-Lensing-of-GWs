"""Case A -- the chirp, static lens.

The inner binary (src/gwlens/system.py) inspirals from f=F_A_START_HZ to
f_isco in ~15 s: 4.4e-5 of the outer orbital period, so the lens is frozen
at a single impact parameter y_A for the whole observation (checked in
tests/test_system.py::check_static_lens_regime).

The source is frozen too, kinematically -- but that is a separate statement
and it is checked separately, not folded into the one above. The outer orbit
also Doppler-shifts the source, and in Case B that effect dominates
everything the lens does; here it does not, for the quantitative reasons
computed in `main()` below (mean_beta_los_over_chirp,
roemer_residual_phase_at_fisco_rad, orbital_redshift_factor_minus_1 in
provenance/numbers.json). See src/gwlens/doppler.py and wiki/log.md.

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
from gwlens import doppler as dp

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

    # --- why this case can treat the source as kinematically static --------
    # The outer orbit moves the source along the line of sight at
    # beta_los ~ 1.6e-2, which in Case B (src/gwlens/doppler.py) dominates
    # everything the lens does. Here it does not, and the reason is
    # quantitative, not hand-waved -- so it is computed, not asserted:
    #
    #  * The merger is placed at the orbital phase of closest approach
    #    (t/P_out = 0.25, chosen just above for the LENSING, not for this),
    #    which is exactly where the line-of-sight velocity crosses zero.
    #  * Over a 15 s window out of a 4-day orbit, a constant delay and a
    #    constant slope in the delay are anyway EXACTLY degenerate with t_c
    #    and with the chirp mass -- unmeasurable in a single event. Only the
    #    curvature (the line-of-sight acceleration) is non-degenerate.
    #  * That residual is `roemer_residual_ptp_s` below: ~9e-6 s, i.e.
    #    ~0.007 rad of phase at the highest in-band frequency.
    #
    # It is therefore NOT applied to the waveform. Applying it would mean
    # resampling h(t) onto a retarded time grid, and the cubic-interpolation
    # error of that resampling would itself exceed the 0.007 rad it was meant
    # to model -- machinery that adds numerical noise in place of a physical
    # effect smaller than the noise. Reported instead; the check behind it is
    # tests/test_doppler.py::check_case_A_roemer_is_negligible.
    orb = dict(a_out_m=system.A_OUT_M, e_out=system.E_OUT,
               period_out=system.P_OUT_S, i_out=np.deg2rad(system.I_OUT_DEG),
               Omega_out=system.OMEGA_OUT, omega_out=system.LITTLE_OMEGA_OUT,
               t_peri=system.T_PERI_OUT, M_lens_msun=system.M_LENS_MSUN,
               M_binary_msun=system.MTOT_MSUN)
    tau_chirp = chirp.time_to_merger(system.F_A_START_HZ, system.MCHIRP_MSUN)
    t_window = np.linspace(t_scan[idx_min] - tau_chirp, t_scan[idx_min], 4001)
    roemer_res, roemer_slope, _ = dp.roemer_residual(t_window, **orb)
    # the best-fit slope over the window, i.e. the MEAN beta_los across the
    # chirp -- not the instantaneous value, which is exactly zero at t/P=0.25
    log("mean_beta_los_over_chirp", float(roemer_slope))
    log("roemer_residual_ptp_s", float(np.ptp(roemer_res)))
    log("roemer_residual_phase_at_fisco_rad",
        float(2.0 * np.pi * system.F_ISCO_HZ * np.ptp(roemer_res)))
    # The one kinematic effect that is NOT small -- and is also not a waveform
    # feature: constant for a circular orbit, hence exactly degenerate with
    # the chirp mass. A 4.1e-4 bias on any M_chirp inferred from this signal,
    # not something visible in any figure here.
    log("orbital_redshift_factor_minus_1",
        float(dp.orbital_redshift_factor(system.A_OUT_M, system.M_LENS_MSUN) - 1.0))

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

    # How good the geometric-optics evaluator actually is AT THIS y -- the
    # approximation the whole case rests on, since F_hybrid is above
    # w_geo_threshold across the entire band. Quoted here at y_A rather than
    # as the worst case over y, because `check_hybrid_matches_at_threshold`'s
    # 3.1% is set entirely by y=0.3, where the two images sit close together
    # and highly magnified and stationary phase is at its weakest -- a y
    # nothing in this project evaluates.
    for label, w_val in (("f_start", w_start), ("f_isco", w_isco)):
        exact = wo.F_point_lens(w_val, y_A)
        log(f"geo_vs_exact_relerr_at_{label}",
            float(abs(exact - wo.F_geometric_optics(w_val, y_A)) / abs(exact)))
    log("geo_vs_exact_relerr_at_threshold_w30",
        float(abs(wo.F_point_lens(30.0, y_A) - wo.F_geometric_optics(30.0, y_A))
              / abs(wo.F_point_lens(30.0, y_A))))

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

    # ...and that rms is not a free number: it has a closed-form prediction,
    # which makes it the one end-to-end validation of this entire pipeline
    # (waveform generator, F evaluator, the FFT convention, the reference
    # phase, the window) against an independent textbook result.
    #   By Parseval, (rms_lensed/rms_unlensed)^2 is the |H(f)|^2-weighted mean
    #   of |F(f)|^2 across the band. In the geometric-optics regime -- which
    #   the whole in-band F is here, w = 62 to 778 --
    #       |F|^2 = mu_+ + |mu_-| + 2 sqrt(mu_+ |mu_-|) sin(2 pi f Delta_T),
    #   and the oscillating term averages away over the ~400 fringes in the
    #   band, leaving exactly mu_+ + |mu_-|: the Paczynski (1986) TOTAL
    #   magnification A(y), already implemented independently in
    #   waveoptics.total_magnification_paczynski and separately checked
    #   against the image sum in tests/test_waveoptics.py. So
    #       rms_strain_amplification = sqrt(A(y_A)),
    #   with no free parameter and nothing fitted. The residual below is the
    #   finite-fringe-count sampling of that average, not a modelling error.
    A_paczynski = float(wo.total_magnification_paczynski(y_A))
    log("total_magnification_paczynski", A_paczynski)
    log("rms_amplification_predicted", float(np.sqrt(A_paczynski)))
    log("rms_amplification_relerr",
        float(abs(NUMBERS["rms_strain_amplification"] - np.sqrt(A_paczynski))
              / np.sqrt(A_paczynski)))

    # ---- Genuine second-image echo: locate it empirically, don't assume --
    # F(w,y) is referenced to the strong image's own arrival (waveoptics.py
    # module docstring, under "Reference phase"), so the naive
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

    # Two windows of equal width at OPPOSITE ends of the band, not one window
    # plus a full-band overview. The old top panel drew the constant envelope
    # as a shaded span across 10-125 Hz with no curve in it, which rendered
    # as a featureless blue rectangle: it asserted the envelope is constant
    # instead of showing it. Two resolved windows at 10 Hz and at f_isco show
    # the same envelope and the same fringe period at w=62 and at w=778,
    # which is the claim, made visible.
    #
    # F is evaluated on its own dense grid here, NOT on the FFT grid: the
    # fringe period is 1/Delta_T = 0.29 Hz and the FFT bin spacing is
    # ~0.04 Hz, i.e. 7 points per fringe, which drew visibly polygonal
    # sinusoids. F(f) is a closed form, so the display grid costs nothing
    # and need not inherit the waveform's resolution.
    f_hi_lo = system.F_ISCO_HZ - (f_zoom_hi - f_zoom_lo)

    def _F_dense(f_lo, f_hi, n=2000):
        ff = np.linspace(f_lo, f_hi, n)
        ww = wo.w_of_frequency(ff, system.M_LENS_MSUN)
        return ff, np.array([wo.F_hybrid(w, y_A) for w in ww])

    f_lo_grid, F_lo = _F_dense(f_zoom_lo, f_zoom_hi)
    f_hi_grid, F_hi = _F_dense(f_hi_lo, system.F_ISCO_HZ)

    fig, axes = plt.subplots(3, 1, figsize=(7, 8.4))
    for axk, (ff, FF, w_here, etiqueta) in zip(
            axes[:2],
            [(f_lo_grid, F_lo, w_start, "inicio de banda"),
             (f_hi_grid, F_hi, w_isco, "cerca del merger")]):
        axk.plot(ff, np.abs(FF), lw=1.0, color="#2b6cb0")
        axk.axhline(sqrt_mu_plus + sqrt_mu_minus, color="#c05621", lw=0.9, ls="--")
        axk.axhline(sqrt_mu_plus - sqrt_mu_minus, color="#c05621", lw=0.9, ls="--")
        axk.set_ylabel(r"$|F(f)|$")
        axk.set_xlabel("$f$ [Hz]")
        axk.set_title(f"{etiqueta}: {ff[0]:.0f}–{ff[-1]:.0f} Hz, $w={w_here:.0f}$",
                      fontsize=9.5)
    # The overall title goes on the first panel rather than in a suptitle.
    # A suptitle has to be positioned by hand relative to a layout that
    # tight_layout computes afterwards, and every value either leaves a band
    # of white or overlaps the first panel's own title; folding it in removes
    # the guess entirely.
    axes[0].set_title(
        f"Caso A: franjas de interferencia en $F(f)$ (lente estático, "
        f"$y={y_A:.3f}$)\n"
        f"inicio de banda: {f_lo_grid[0]:.0f}–{f_lo_grid[-1]:.0f} Hz, "
        f"$w={w_start:.0f}$", fontsize=10)
    axes[1].set_title(
        f"cerca del merger: {f_hi_grid[0]:.0f}–{f_hi_grid[-1]:.0f} Hz, "
        f"$w={w_isco:.0f}$", fontsize=9.5)
    axes[0].text(0.985, 0.90, r"$\sqrt{\mu_+}\pm\sqrt{|\mu_-|}$",
                 transform=axes[0].transAxes, ha="right", va="top",
                 fontsize=8, color="#c05621")
    axes[2].plot(f_lo_grid, np.angle(F_lo), lw=1.0, color="#c05621")
    axes[2].set_ylabel(r"$\arg F(f)$ [rad]")
    axes[2].set_xlabel("$f$ [Hz]")
    # F traces a circle of radius sqrt(mu_-) centered on sqrt(mu_+) in the
    # complex plane as w*DeltaT sweeps 2*pi per fringe (the strong image is
    # the fixed reference, weak image the rotating arm) -- since mu_- < mu_+
    # here, that circle does not enclose the origin, so arg F oscillates
    # (bounded, not a full 2*pi wrap) once per fringe, in phase with |F|
    # above.
    axes[2].set_title("la fase, sobre la misma ventana: oscilación acotada, "
                      "una vuelta por franja", fontsize=9.5)
    # tight_layout FIRST, suptitle after: tight_layout reserves vertical space
    # for a suptitle even when that suptitle is placed outside the canvas, so
    # calling it second leaves a band of white between the title and the top
    # panel that bbox_inches="tight" does not crop. Placed afterwards, at
    # y just above 1, the tight bounding box simply expands to include it.
    fig.tight_layout()
    fig.savefig(OUT / "caseA_F_of_f.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 2: lensed vs unlensed time-domain strain -----------------
    # Three panels: overview (extended past merger to include the
    # second-image echo), a merger/ringdown zoom, and a dedicated echo zoom
    # -- the echo is the unambiguous causality signature (a second, weaker,
    # LATER copy of the signal). With F referenced to the strong image's own
    # arrival (waveoptics.py, "Reference phase"), the
    # lensed merger now peaks at (to within a sample of) the SAME time as
    # the unlensed one -- see wiki/log.md for the earlier, incorrect version
    # of this story, where an un-subtracted reference phase looked like a
    # genuine 0.6s early shift and was mis-attributed to interference.
    t_peak = t_peak_unlensed
    t_hi = min(t[-1], echo_t + 1.0)
    view = (t >= t_peak - 1.15 * NUMBERS["t_end_seconds"]) & (t <= t_hi)
    # Four panels, in the order a reader needs them: the signal without the
    # lens, the same signal with it, a zoom centred on the merger, and the
    # DIFFERENCE. The first two were previously drawn on top of each other in
    # one panel, where the lensed curve simply hid the unlensed one and
    # nothing could be compared -- at this time resolution the two are
    # indistinguishable by eye anyway, which is itself the point: what the
    # lens does is not visible in the raw strain.
    #
    # The difference panel is where it becomes visible, and it is worth its
    # own panel for a reason. In geometric optics F = sqrt(mu_+) - i
    # sqrt(|mu_-|) exp(i w Delta_T), so h_lensed - h_unlensed is
    # (sqrt(mu_+)-1) times the signal, plus a full sqrt(|mu_-|) copy of it
    # delayed by Delta_T. The first term is 3% and rides under the chirp; the
    # second is 24% and sits where the chirp has already ended, so the
    # subtraction isolates the second image against nothing at all.
    h_diff = h_lensed - h_unlensed
    fig, axes = plt.subplots(4, 1, figsize=(8, 10.6))
    y_span = 1.08 * float(np.max(np.abs(h_lensed[view])))
    axes[0].plot(t[view], h_unlensed[view], lw=0.5, color="#888888")
    axes[0].set_ylabel("$h(t)$ [u. arb.]")
    axes[0].set_ylim(-y_span, y_span)
    axes[0].set_title("sin lente", fontsize=9.5)
    axes[1].plot(t[view], h_lensed[view], lw=0.5, color="#2b6cb0")
    axes[1].set_xlabel(f"$t$ [s] (merger en $t={t_peak:.2f}$ s)")
    axes[1].set_ylabel("$h(t)$ [u. arb.]")
    axes[1].set_ylim(-y_span, y_span)
    axes[1].set_title("con lente, misma escala", fontsize=9.5)
    # Only the model identity, not the full provenance string: everything
    # after the " -- " is an English descriptor that belongs in numbers.json
    # and RESULTS.md, not in a figure inside a Spanish document. Splitting
    # rather than hardcoding keeps the TaylorF2 fallback distinguishable
    # here too, without maintaining a second, translated copy of a label
    # whose authoritative version is `waveform_source`.

    # window wide enough to comfortably contain BOTH the unlensed peak (at
    # t_peak, by definition) and the lensed curve's own peak, which is not
    # guaranteed to sit at t_peak too -- F(f) dephases the lensed merger,
    # so its extremum can be offset from the unlensed one (an earlier,
    # narrower window here missed it; see wiki/log.md)
    # Centred on the merger, not starting at it: a window running from -0.05
    # to +0.30 s put the merger hard against the left edge with a quarter
    # second of ringdown taking up the rest of the frame.
    zoom = (t > t_peak - 0.15) & (t < t_peak + 0.15)
    axes[2].plot(t[zoom] - t_peak, h_unlensed[zoom], lw=1.1, color="#888888", label="sin lente")
    axes[2].plot(t[zoom] - t_peak, h_lensed[zoom], lw=1.1, color="#2b6cb0", alpha=0.9, label="con lente")
    axes[2].set_xlabel(r"$t - t_\mathrm{merger}$ [s]")
    axes[2].set_ylabel("$h(t)$ [u. arb.]")
    axes[2].set_title("ampliación: la primera imagen", fontsize=9.5)
    axes[2].legend(loc="lower left", fontsize=8)

    # The RATIO of the envelopes, not their difference. The difference isolates
    # the second image -- a delayed copy of the whole signal -- which is a
    # statement about arrival times. The ratio is a statement about
    # interference, which is what the lens actually does to the waveform here:
    # dividing the lensed analytic signal by the unlensed one leaves |F|, and
    # because the chirp sweeps frequency monotonically, plotting it against
    # time sweeps out F's interference fringes. It oscillates between
    # sqrt(mu_+) -/+ sqrt(|mu_-|) throughout, and the fringes crowd together
    # towards the merger because df/dt does: measured, 1.0 fringes per second
    # at 10 Hz against 93 per second at 36 Hz. The fringe period is constant
    # in FREQUENCY (0.29 Hz, caseA_F_of_f.png); this is that same curve seen
    # through f(t).
    #
    # Cut before the merger: past it the unlensed signal decays into the
    # ringdown while the lensed one still carries the second image, so the
    # ratio stops measuring interference and starts diverging.
    rview = (t >= t_peak - NUMBERS["t_end_seconds"] + 1.0) & (t <= t_peak - 0.25)
    ratio = np.abs(env_l[rview]) / np.abs(env_u[rview])
    axes[3].plot(t[rview] - t_peak, ratio, lw=0.6, color="#6b46c1")
    for lvl in (sqrt_mu_plus - sqrt_mu_minus, sqrt_mu_plus + sqrt_mu_minus):
        axes[3].axhline(lvl, color="#c05621", lw=0.9, ls="--")
    axes[3].set_xlabel(r"$t - t_\mathrm{merger}$ [s]")
    axes[3].set_ylabel("con lente / sin lente")
    axes[3].set_title("el cociente: la interferencia entre las dos imágenes",
                      fontsize=9.5)
    fig.tight_layout()   # before the suptitle; see caseA_F_of_f above for why
    fig.suptitle("Caso A: qué le hace el lente a la forma de onda",
                 fontsize=11, y=1.012)
    # Model identity as a footnote, not in the title: the string carries the
    # approximant name and, on the fallback path, a parenthetical about why,
    # which together are wider than the figure.
    fig.text(0.5, -0.006, NUMBERS["waveform_source"].split(" -- ")[0],
             ha="center", va="top", fontsize=8, color="0.35")
    fig.savefig(OUT / "caseA_strain_time.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 3: extended diffraction/interference pattern in y-plane -
    from matplotlib.colors import LogNorm
    import matplotlib.ticker as mticker

    def _label_colorbar(cb, im):
        """Readable ticks on a log colour bar, whichever range it spans.

        Two cases, and the default formatter is wrong for both. Over many
        decades it writes 10^0 for the tick that just means 1; inside a
        single decade the log locator places exactly one major tick and
        leaves the rest to the minor formatter, which renders them as
        "6 x 10^-1" -- scientific notation whose exponent carries nothing,
        since every tick shares it. So: plain decimals and explicit ticks
        when the span is narrow, plain powers of ten when it is wide.
        """
        lo, hi = im.get_clim()
        if hi / max(lo, 1e-12) < 20.0:
            step = 0.1 if hi - lo < 1.2 else 0.25
            ticks = [round(v, 2) for v in np.arange(0.0, hi + step, step)
                     if lo <= v <= hi]
            cb.set_ticks(ticks)
            cb.ax.yaxis.set_major_formatter(
                mticker.FuncFormatter(lambda v, _: f"{v:.1f}"))
            cb.ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        else:
            cb.ax.yaxis.set_major_formatter(mticker.FuncFormatter(
                lambda v, _: ("1" if abs(v - 1.0) < 1e-9
                              else fr"$10^{{{round(np.log10(v))}}}$")))

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4))

    # n_grid=1201, not 481: the fringe period in y shrinks with radius (it is
    # 2*pi/(w * dDeltaT/dy), = 0.050 at y=0.3 but 0.028 at y=3), so 481 points
    # across the 6-wide axis gave dy=0.0125 -- only 2.2 samples per fringe on
    # the OUTER rings, below what is needed to render them without moire.
    # 1201 gives dy=0.005, ~5.6 samples per fringe there. `n_radial` was
    # already fine (dr=0.0011, ~26 per fringe); this was the display grid
    # only, the same distinction the docstring above records.
    # y_max=1.6 rather than 3.0, and a finer display grid: the fringe period
    # in y is ~0.03-0.05 here, so 3.0 across 1201 pixels left ~5 px per ring
    # on the outer half and the rings read as moire rather than as rings.
    # Nothing is lost: the source sits at y_A=1.589 and the pattern is
    # axisymmetric, so the outer rings are more of the same.
    # y_max=2.0: at 1.6 the source marker at y_A=1.589 sat exactly on the
    # frame edge and was drawn half outside it, while the caption said it
    # was marked. 2.0 still leaves ~11 px per fringe at 1601 points.
    axis0, axis0y, mag2_0 = make_ring_pattern(w_start, y_max=2.0, n_grid=1601)
    im0 = axes[0].pcolormesh(axis0, axis0y, mag2_0, shading="auto", cmap="inferno",
                              norm=LogNorm(vmin=max(mag2_0.min(), 1e-3), vmax=mag2_0.max()))
    axes[0].set_title(f"inicio de banda, $f={system.F_A_START_HZ}$ Hz, $w={w_start:.1f}$", fontsize=9)
    cb0 = fig.colorbar(im0, ax=axes[0], shrink=0.8, label=r"$|F|^2$ (escala logarítmica)")
    _label_colorbar(cb0, im0)

    # half_width 0.02, not 0.15: the fringe period in y scales as 1/w, so at
    # w=778 it is ~0.003 and a 0.30-wide window packed ~100 rings into 500
    # pixels -- the moire the whole figure was accused of. 0.04 wide leaves
    # ~13 fringes, which is what a fringe pattern has to look like to read as
    # one.
    axis1, axis1y, mag2_1 = make_ring_pattern(
        w_isco, n_grid=500, center=(y_A, 0.0), half_width=0.02, n_radial=12000)
    im1 = axes[1].pcolormesh(axis1, axis1y, mag2_1, shading="auto", cmap="inferno",
                              norm=LogNorm(vmin=max(mag2_1.min(), 1e-3), vmax=mag2_1.max()))
    axes[1].set_title(f"cerca del merger, $f=f_\\mathrm{{isco}}$, $w={w_isco:.0f}$\n"
                       r"(ampliado a $|\Delta y|<0.02$ en torno a la fuente)", fontsize=9)
    cb1 = fig.colorbar(im1, ax=axes[1], shrink=0.8, label=r"$|F|^2$ (escala logarítmica)")
    _label_colorbar(cb1, im1)

    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel(r"$y_1$")
        ax.plot(y_A, 0, "w+", ms=12, mew=2)
    axes[0].set_ylabel(r"$y_2$")
    # A wider, single-line-per-row suptitle was getting clipped at both
    # edges of this fairly narrow figure; shortened and set explicitly
    # inside the axes bounding box (independent review, see wiki/log.md).
    fig.tight_layout()   # before the suptitle; see caseA_F_of_f above for why
    fig.suptitle("Caso A: patrón de amplificación en el plano de la fuente\n"
                 fr"(anillos axisimétricos; la fuente, en $y_A={y_A:.3f}$, está marcada)",
                 fontsize=10, y=1.04)
    # bbox_inches="tight": with set_aspect("equal") the axes are resized AFTER
    # tight_layout has computed positions, so the x-axis label ended up drawn
    # below the figure canvas and was cut off in the committed PNG (both
    # panels lost their "$y_1$"). Cropping to the actual artist extents at
    # save time is the fix that cannot silently regress the next time this
    # figure's geometry changes.
    fig.savefig(OUT / "caseA_diffraction_pattern.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 4: simple idealized detector view -------------------------
    # Log-scale y-axis and a window extended past the echo: the echo is
    # ~4x (0.6 dex) below the merger peak envelope
    # (echo_to_unlensed_peak_ratio, matching sqrt(|mu_-|) to <1%, see above)
    # -- clear on a log axis, easy to under-sell on a linear one.
    fig, ax = plt.subplots(figsize=(8, 3.6))
    # same window as Figure 2's overview panel above -- identical expression, reused rather than rebuilt
    ax.plot(t[view], env_u[view], color="#888888", lw=1.0, label="sin lente")
    ax.plot(t[view], env_l[view], color="#2b6cb0", lw=1.0, label="con lente")
    ax.annotate("1ra imagen", xy=(t_peak, float(np.max(env_l))),
                xytext=(-46, -14), textcoords="offset points", fontsize=8,
                color="#2b6cb0")
    ax.annotate("2da imagen\n($+%.2f$ s)" % (echo_t - t_peak),
                xy=(echo_t, NUMBERS["echo_peak_env_lensed"]),
                xytext=(-8, 12), textcoords="offset points", fontsize=8,
                color="#2b6cb0")
    ax.set_yscale("log")
    ax.set_xlabel("$t$ [s]")
    ax.set_ylabel("$|h(t)|$ [u. arb.]")
    # Retitled: "vista de detector idealizada (envolvente)" said what the
    # figure was made of, not what it shows. What it shows is that the lens
    # turns one event into two arrivals.
    ax.set_title("Caso A: con lente llegan dos pulsos, no uno\n"
                 "(amplitud instantánea $|h(t)|$, escala logarítmica)", fontsize=10)
    # Lower left: the top right of this axes is where the merger peak and the
    # echo both are, and a legend there covered them. The green band that
    # used to mark the echo is gone for the same reason -- it shaded the
    # feature it was pointing at.
    ax.legend(fontsize=8, loc="lower left")
    # Headroom above the merger peak (the default limits clipped it, so the
    # tallest thing in the figure ran off the top of the frame) and a floor
    # three decades down. Without the floor the axis auto-scales to whatever
    # the post-merger tail decays to, which on the inspiral-only fallback
    # path is ~1e-27 and stretches the useful part of the plot into a band a
    # few pixels tall.
    env_top = float(np.max(env_l[view]))
    ax.set_ylim(env_top / 1.0e3, 3.0 * env_top)
    fig.tight_layout()
    fig.savefig(OUT / "caseA_detector_envelope.png", dpi=170, bbox_inches="tight")
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
