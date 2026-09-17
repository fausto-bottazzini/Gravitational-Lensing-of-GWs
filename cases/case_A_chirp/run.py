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

    # Dos paneles, no tres. La version anterior eran tres ventanas de |F(f)|
    # y arg F(f) que, medidas, son la MISMA curva: el periodo de franja vale
    # 0.29120 Hz a 10 Hz, a 60 y a 123, y las cotas de |F| y de arg F
    # coinciden en las tres a cuatro decimales. Repetir una curva tres veces
    # para demostrar que es la misma es un mal uso de tres paneles.
    #
    # Lo que si hay que mostrar es de donde salen esas cotas. En optica
    # geometrica F = sqrt(mu_+) - i sqrt(|mu_-|) exp(i w Delta_T), asi que al
    # barrer f el punto recorre una CIRCUNFERENCIA de radio sqrt(|mu_-|)
    # centrada en sqrt(mu_+), una vuelta por franja. De ahi se leen las tres
    # cosas a la vez: |F| esta entre los puntos mas cercano y mas lejano de
    # esa circunferencia al origen; arg F no da nunca la vuelta porque la
    # circunferencia no encierra al origen, y las tangentes desde el origen
    # dan la cota arcsin(sqrt(|mu_-|)/sqrt(mu_+)); y el periodo es 1/Delta_T
    # porque w Delta_T es lineal en f.
    un_giro = 1.0 / NUMBERS["image_time_delay_seconds"]        # una franja, en Hz
    f_giro, F_giro = _F_dense(f_zoom_lo, f_zoom_lo + un_giro, n=800)
    arg_max = float(np.arcsin(sqrt_mu_minus / sqrt_mu_plus))

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.3),
                             gridspec_kw={"width_ratios": [1, 1.25]})

    ax = axes[0]
    ax.axhline(0, color="0.6", lw=0.8)
    ax.axvline(0, color="0.6", lw=0.8)
    for signo in (+1, -1):
        ax.plot([0, 1.45 * np.cos(signo * arg_max)],
                [0, 1.45 * np.sin(signo * arg_max)],
                color="#c05621", lw=0.9, ls="--")
    ax.plot(F_giro.real, F_giro.imag, color="#2b6cb0", lw=2.0)
    ax.plot([sqrt_mu_plus - sqrt_mu_minus, sqrt_mu_plus + sqrt_mu_minus], [0, 0],
            "o", color="#c05621", ms=5)
    ax.set_aspect("equal")
    ax.set_xlim(-0.12, 1.45)
    ax.set_ylim(-0.62, 0.62)
    ax.set_xlabel(r"$\mathrm{Re}\,F$")
    ax.set_ylabel(r"$\mathrm{Im}\,F$")
    ax.set_title("$F$ en el plano complejo: una vuelta por franja", fontsize=9.5)

    axes[1].plot(f_lo_grid, np.angle(F_lo), lw=1.1, color="#c05621")
    for signo in (+1, -1):
        axes[1].axhline(signo * arg_max, color="#c05621", lw=0.9, ls="--")
    axes[1].set_ylabel(r"$\arg F(f)$ [rad]")
    axes[1].set_xlabel("$f$ [Hz]")
    axes[1].set_title("la fase, acotada: nunca da la vuelta", fontsize=9.5)
    # La cota va dentro de los ejes, pero el hueco hay que FABRICARLO: con
    # los límites automáticos la curva llena toda la banda entre las dos
    # punteadas, y el rótulo quedaba escrito encima de los mínimos de la
    # fase (y con su extremo izquierdo fuera del eje, alineado a la derecha).
    # Se baja el límite inferior hasta abrir una franja vacía, el rótulo va
    # centrado en ella, y encima lleva recuadro opaco: eso es lo que lo
    # mantiene legible aunque la geometría de la figura vuelva a cambiar.
    axes[1].set_ylim(-2.10 * arg_max, 1.22 * arg_max)
    axes[1].text(0.5, 0.045,
                 r"$|\arg F|\leq\arcsin(\sqrt{|\mu_-|}/\sqrt{\mu_+})="
                 + f"{arg_max:.4f}$",
                 transform=axes[1].transAxes, ha="center", va="bottom",
                 fontsize=8.5, color="#c05621",
                 bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                           edgecolor="#c05621", linewidth=0.7, alpha=0.95))
    fig.tight_layout()
    fig.suptitle(
        f"Caso A: dos imágenes interfieren (lente estático, $y={y_A:.3f}$)\n"
        f"radio $\\sqrt{{|\\mu_-|}}={sqrt_mu_minus:.4f}$, centro "
        f"$\\sqrt{{\\mu_+}}={sqrt_mu_plus:.4f}$, "
        f"período de franja $1/\\Delta T={un_giro:.4f}$ Hz",
        fontsize=10.5, y=1.02)
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
    # The window is set by the envelope, not by taste. env_l and env_u are
    # analytic-signal magnitudes, and an analytic-signal envelope can only
    # follow a modulation SLOWER than its own carrier. The fringe rate here
    # is Delta_T * df/dt, so the ratio Delta_T*(df/dt)/f is the figure of
    # merit: it runs 0.10 at 13 s before the merger, 0.43 at 3 s, 0.86 at
    # 1.5 s and 2.17 at 0.6 s. Past ~0.25 the drawn fringes are no longer the
    # real ones -- they shrink, and an earlier version of this panel showed
    # exactly that, an amplitude decaying into the merger that is a property
    # of the Hilbert transform and not of the lens. Cut where the criterion
    # says to.
    f_of_t = chirp.freq_of_time(t - (t_peak - NUMBERS["t_c_seconds"]),
                                NUMBERS["t_c_seconds"], system.MCHIRP_MSUN)
    dfdt = np.gradient(f_of_t, t)
    slow = NUMBERS["image_time_delay_seconds"] * np.abs(dfdt) < 0.40 * f_of_t
    rview = ((t >= t_peak - NUMBERS["t_end_seconds"] + 1.0)
             & (t <= t_peak) & slow)
    t_cut = float(t[rview][-1] - t_peak)
    log("ratio_panel_cut_before_merger_s", -t_cut)
    ratio = np.abs(env_l[rview]) / np.abs(env_u[rview])
    axes[3].plot(t[rview] - t_peak, ratio, lw=0.6, color="#6b46c1")
    for lvl, lab in ((sqrt_mu_plus + sqrt_mu_minus,
                      r"$\sqrt{\mu_+}\pm\sqrt{|\mu_-|}$"),
                     (sqrt_mu_plus - sqrt_mu_minus, None)):
        axes[3].axhline(lvl, color="#c05621", lw=0.9, ls="--", label=lab)
    axes[3].legend(fontsize=8, loc="lower left")
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
                 fr"(anillos axisimétricos; la fuente, en $y_A={y_A:.3f}$, marcada)",
                 fontsize=10, y=1.04)
    # bbox_inches="tight": with set_aspect("equal") the axes are resized AFTER
    # tight_layout has computed positions, so the x-axis label ended up drawn
    # below the figure canvas and was cut off in the committed PNG (both
    # panels lost their "$y_1$"). Cropping to the actual artist extents at
    # save time is the fix that cannot silently regress the next time this
    # figure's geometry changes.
    fig.savefig(OUT / "caseA_diffraction_pattern.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 4: time-frequency, the way a chirp is usually shown -------
    # Replaces an amplitude-vs-time plot of the same two waveforms, which
    # said nothing the strain figure above did not already say. A
    # spectrogram does: the inspiral is a track sweeping up in frequency,
    # and lensing puts a SECOND track on the plane, the same sweep displaced
    # by Delta_T. Two images become two tracks, which is the whole of
    # geometric-optics lensing in one picture.
    from scipy.signal import spectrogram

    # nperseg=1024 at fs=2048 is a 0.5 s window: five cycles at the 10 Hz
    # band start, which is the least it can be and still localise the low
    # end, and short enough that the 0.29 Hz fringes do not resolve into the
    # picture (they are caseA_F_of_f.png's subject, not this one's).
    nperseg, noverlap = 1024, 1024 - 64
    tf_view = (t >= t_peak - NUMBERS["t_end_seconds"]) & (t <= echo_t + 0.6)
    fig, axes = plt.subplots(2, 1, figsize=(8, 6.2), sharex=True, sharey=True)
    espectros = []
    for h_sig in (h_unlensed[tf_view], h_lensed[tf_view]):
        # Hann, not the default Tukey(0.25): its sidelobes are ~30 dB
        # lower, and with three decades of colour range the Tukey ones
        # drew a fan of spurious arcs above the real track.
        f_sp, t_sp, S = spectrogram(h_sig, fs=fs, window="hann",
                                    nperseg=nperseg, noverlap=noverlap,
                                    scaling="spectrum", mode="magnitude")
        espectros.append((f_sp, t_sp, S))
    vmax = max(S.max() for _, _, S in espectros)
    for axk, (f_sp, t_sp, S), etiqueta in zip(
            axes, espectros, ("sin lente", "con lente")):
        band = (f_sp >= 8.0) & (f_sp <= 400.0)
        im = axk.pcolormesh(t_sp + t[tf_view][0] - t_peak, f_sp[band],
                            S[band], shading="auto", cmap="viridis",
                            norm=LogNorm(vmin=vmax / 2.0e2, vmax=vmax))
        axk.set_yscale("log")
        axk.set_ylabel("$f$ [Hz]")
        axk.set_title(etiqueta, fontsize=9.5)
    axes[1].set_xlabel(r"$t - t_\mathrm{merger}$ [s]")
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label("amplitud espectral [u. arb.]")
    fig.suptitle("Caso A: el chirp en tiempo-frecuencia — el lente lo duplica",
                 fontsize=11, y=0.975)
    fig.savefig(OUT / "caseA_spectrogram.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # ---- Animation data for report/report.html (JSON, reduced resolution) -
    # La ventana llega hasta pasado el eco: la envolvente de la segunda imagen
    # es parte de lo que las diapositivas muestran, y cortarla en t_peak+0.3
    # dejaba afuera justamente eso.
    win = ((t >= t_peak - 1.15 * NUMBERS["t_end_seconds"])
           & (t <= t_peak + NUMBERS["image_time_delay_seconds"] + 0.9))
    idx = np.where(win)[0]
    stride = max(1, len(idx) // 2000)
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

    # Ventana densa, para dibujar la SEÑAL y no su envolvente. La grilla de
    # arriba tiene 1455 puntos en 17.7 s: 82 muestras por segundo contra una
    # portadora que va de 10 a 126 Hz. Con eso solo se puede dibujar una
    # envolvente; la onda sale aliaseada. Aca va h(t) a 1024 Hz sobre una
    # ventana que llega hasta pasado el eco, normalizada al pico sin lente
    # (la diapositiva dibuja formas, no unidades) y redondeada a cinco
    # decimales para no inflar el HTML.
    hi_lo = t_peak - 6.0
    hi_hi = t_peak + NUMBERS["image_time_delay_seconds"] + 0.9
    paso_hi = max(1, int(round(fs / 1024.0)))
    hidx = np.where((t >= hi_lo) & (t <= hi_hi))[0][::paso_hi]
    escala_hi = float(np.abs(h_unlensed[hidx]).max())
    anim_wave_hi = {
        "t0": float(t[hidx[0]] - t_peak),
        "dt": float(paso_hi / fs),
        "scale": escala_hi,
        "h_unlensed": np.round(h_unlensed[hidx] / escala_hi, 5).tolist(),
        "h_lensed": np.round(h_lensed[hidx] / escala_hi, 5).tolist(),
    }

    zmask_idx = np.where(zmask)[0]
    anim = {
        "t_peak": float(t_peak),
        "t_end": float(NUMBERS["t_end_seconds"]),
        "f_isco": float(f_isco),
        "f_start": 10.0,
        "sqrt_mu_plus": sqrt_mu_plus,
        "sqrt_mu_minus": sqrt_mu_minus,
        # a_in va como f^(-2/3); alcanza con el valor en 1 Hz para que el
        # esquema de la coalescencia lo reconstruya sin repetir constantes
        "a_in_km_at_1hz": float(
            system.inner_separation_m(1.0, system.M1_MSUN, system.M2_MSUN) / 1e3),
        "wave_hi": anim_wave_hi,
        "f_zoom_lo": float(f_zoom_lo),
        "f_zoom_hi": float(f_zoom_hi),
        "image_time_delay_s": float(NUMBERS["image_time_delay_seconds"]),
        "waveform": {
            "t": t[idx].tolist(),
            "h_unlensed": h_unlensed[idx].tolist(),
            "h_lensed": h_lensed[idx].tolist(),
            "env_unlensed": env_u[idx].tolist(),
            "env_lensed": env_l[idx].tolist(),
            "f_of_t": f_of_t.tolist(),
        },
        # Grilla propia, no la de la FFT. La de la FFT tiene 0.031 Hz de
        # paso contra un periodo de franja de 0.291: nueve puntos por
        # franja, con los que la fase sale poligonal y el lugar geometrico
        # en el plano complejo se dibuja a cuerdas que cortan la
        # circunferencia. F(f) es forma cerrada, asi que la grilla de
        # dibujo no tiene por que heredar la resolucion de la onda.
        "fringes": {
            "freqs": f_lo_grid.tolist(),
            "F_abs": np.abs(F_lo).tolist(),
            "F_arg": np.angle(F_lo).tolist(),
            # una sola vuelta, para dibujar la circunferencia sin repetirla
            "giro_re": F_giro.real.tolist(),
            "giro_im": F_giro.imag.tolist(),
            "arg_max": arg_max,
        },
        # Parameters only, no grid: across this whole band w > 30, so the
        # pattern is the two-image geometric-optics form and report.html
        # evaluates it directly rather than being shipped a sampled one --
        # |F|^2 = mu_+ + |mu_-| + 2 sqrt(mu_+|mu_-|) sin(w Delta_T(y)), with
        # mu_+-(y) and Delta_T(y) elementary. Same expression as
        # waveoptics.F_geometric_optics; agreement checked below.
        "pattern": {
            "y_A": float(y_A),
            "w_start": float(w_start),
            "w_isco": float(w_isco),
            "f_start": float(system.F_A_START_HZ),
            "f_isco": float(f_isco),
            "y_max": 2.0,
            "zoom_half_width": 0.05,
        },
    }
    # The closed form report.html will evaluate, checked here against this
    # repo's own evaluator so the page cannot quietly drift from it.
    _r = np.linspace(0.05, 2.0, 400)
    _worst = 0.0
    for _w in (w_start, 0.5 * (w_start + w_isco), w_isco):
        _xp, _xm = wo.image_positions(_r)
        _mp, _mm = wo.magnification(_xp), np.abs(wo.magnification(_xm))
        _dT = wo.time_delay_difference(_r)
        _closed = _mp + _mm + 2.0 * np.sqrt(_mp * _mm) * np.sin(_w * _dT)
        _ref = np.array([abs(complex(wo.F_geometric_optics(_w, rr))) ** 2 for rr in _r])
        _worst = max(_worst, float(np.max(np.abs(_closed - _ref) / _ref)))
    log("pattern_closed_form_vs_evaluator_relerr", _worst)
    with open(OUT / "caseA_animation_data.json", "w") as fh:
        json.dump(anim, fh)

    with open(OUT / "provenance" / "numbers.json", "w") as fh:
        json.dump(NUMBERS, fh, indent=2)
    print(json.dumps(NUMBERS, indent=2))
    print("\nFigures written:",
          "caseA_F_of_f.png, caseA_strain_time.png, "
          "caseA_diffraction_pattern.png, caseA_spectrogram.png")


if __name__ == "__main__":
    main()
