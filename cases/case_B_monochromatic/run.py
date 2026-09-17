"""Case B -- quasi-monochromatic source, moving lens.

Same triple as Case A (src/gwlens/system.py), earlier in the SAME inner
binary's inspiral: f=F_B_HZ=0.05 Hz, changing by only ~2% over the
T_OBS_B_S=3-outer-period (12 d) observation window used here (checked in
tests/test_system.py::check_quasi_monochromatic_regime). Over that window
the outer orbit is NOT frozen: the impact parameter y(t) sweeps through the
diffraction pattern derived in Case A, producing periodic amplification
pulses once per outer period -- "repeated lensing" (D'Orazio & Loeb 2020),
here from wave optics rather than their geometric-optics treatment.

THE OUTER ORBIT DOES TWO THINGS, NOT ONE (added 2026-09-07, see
src/gwlens/doppler.py and wiki/log.md). It sweeps y(t) through the
diffraction pattern -- the lensing pulses this case was built for -- and it
also moves the source along the line of sight at beta_los~1.6e-2, imprinting
a Roemer/Doppler modulation of +-905 s, i.e. ~91 GW cycles peak-to-peak at
f_B. That is ~527x the 3.43 s image time delay at the pulse. Until this was
added, this script modelled only the smaller of the two and described it as
the outer orbit's signature. Both are now included: the Doppler enters as a
retarded emission time (nothing multiplied anywhere), it changes the phase
by ~90 cycles and the strain ENVELOPE by ~1e-5, and so every amplitude
result in this case -- |F|^2 pulses, rms amplification, every diffraction
figure -- is numerically unchanged by it. Figure 5
(caseB_doppler_vs_lensing.png) puts the two side by side.

KEY DIFFERENCE FROM CASE A: the lens itself moves on a timescale (hours-days,
via y(t)) that is enormously slower than the GW oscillation period
(1/f_B=20 s), so instead of one FFT/IFFT multiplication by a FIXED F(f), we
multiply the (complex, analytic-signal) waveform by the ADIABATIC/
instantaneous F(w_B, y(t)) at each time sample -- valid because y(t) is, to
outstanding accuracy, constant over any one GW cycle (theory/theory.tex
Sec. 7 derives the timescale hierarchy that justifies this).

Run: `python3 cases/case_B_monochromatic/run.py`
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from gwlens import system, chirp, waveoptics as wo, geometry as geo, units
from gwlens import doppler as dp

OUT = HERE
NUMBERS = {}


def log(name, value):
    NUMBERS[name] = value
    return value


def main():
    print(system.summary())

    w_B = wo.w_of_frequency(system.F_B_HZ, system.M_LENS_MSUN)
    log("w_B", float(w_B))

    # --- outer-orbit impact parameter over the observation ---------------
    n_per_period = 600
    n_periods = system.T_OBS_B_S / system.P_OUT_S
    n_samples = int(n_per_period * n_periods)
    t = np.linspace(0.0, system.T_OBS_B_S, n_samples)
    y_t, lensed_mask, theta_E = geo.impact_parameter_of_time(
        t, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
        system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)
    log("fraction_of_time_lensed", float(np.mean(lensed_mask)))
    log("min_y_during_observation", float(np.min(y_t[lensed_mask])))
    log("n_outer_periods_observed", float(n_periods))
    log("n_samples", n_samples)

    # --- F(w_B, y(t)): exact evaluator throughout (w_B << w_geo_threshold) -
    t0 = time.time()
    F_t = np.array([wo.F_hybrid(w_B, y) if lm else 1.0 + 0.0j
                     for y, lm in zip(y_t, lensed_mask)])
    log("F_evaluation_seconds", time.time() - t0)
    log("max_abs_F", float(np.max(np.abs(F_t))))
    log("mean_abs_F2_over_observation", float(np.mean(np.abs(F_t) ** 2)))

    # --- the SAME outer orbit's kinematic imprint: Roemer / Doppler --------
    # The outer orbit does not only move the source through the diffraction
    # pattern; it also moves the source along the line of sight, at
    # beta_los ~ 1.6e-2 with a light-crossing time a_out/c ~ 900 s. That is
    # the first-order Doppler effect, it is imprinted by exactly the same
    # orbit, and it is far LARGER than the lensing delay this case exists to
    # show -- so leaving it out (as this script did until 2026-09-07) made
    # the case's headline claim about what the outer orbit does to the
    # waveform simply wrong by two orders of magnitude. See
    # src/gwlens/doppler.py and wiki/log.md. It is applied the only way that
    # needs no separate factor anywhere: evaluate the emitted waveform at the
    # retarded emission time.
    orb = dict(a_out_m=system.A_OUT_M, e_out=system.E_OUT,
               period_out=system.P_OUT_S, i_out=np.deg2rad(system.I_OUT_DEG),
               Omega_out=system.OMEGA_OUT, omega_out=system.LITTLE_OMEGA_OUT,
               t_peri=system.T_PERI_OUT, M_lens_msun=system.M_LENS_MSUN,
               M_binary_msun=system.MTOT_MSUN)
    orb_novec = {k: v for k, v in orb.items() if k != "Omega_out"}
    roemer_t = dp.roemer_delay(t, **orb)
    t_em, n_iter_em = dp.emission_time(t, **orb)
    beta_los_t = dp.los_velocity(t, **orb_novec) / units.C_SI
    log("roemer_delay_ptp_s", float(np.ptp(roemer_t)))
    log("roemer_phase_ptp_cycles", float(system.F_B_HZ * np.ptp(roemer_t)))
    log("roemer_phase_ptp_rad", float(2.0 * np.pi * system.F_B_HZ * np.ptp(roemer_t)))
    log("beta_los_max", float(np.max(np.abs(beta_los_t))))
    log("emission_time_iterations", int(n_iter_em))
    log("orbital_redshift_factor_minus_1",
        float(dp.orbital_redshift_factor(system.A_OUT_M, system.M_LENS_MSUN) - 1.0))

    # --- underlying quasi-monochromatic waveform (complex analytic signal) -
    # Everything below is evaluated at the EMISSION time t_em, not the
    # observation time t: that single substitution is the complete
    # first-order Doppler treatment (doppler.emission_time's docstring; the
    # identity d t_em/d t_obs = 1/(1+beta) is checked in
    # tests/test_doppler.py). The remaining, second-order + gravitational
    # piece is exactly constant for this circular orbit and therefore
    # degenerate with M_chirp -- logged above as
    # orbital_redshift_factor_minus_1, deliberately not applied.
    t_c = chirp.time_to_merger(system.F_B_HZ, system.MCHIRP_MSUN)
    f_t = chirp.freq_of_time(t_em, t_c, system.MCHIRP_MSUN)
    log("f_start_Hz", float(f_t[0]))
    log("f_end_Hz", float(f_t[-1]))
    log("fractional_freq_drift", float((f_t[-1] - f_t[0]) / f_t[0]))
    # observed frequency = emitted / (1+beta_los): the Doppler swing is
    # itself larger than the whole intrinsic drift across the observation
    f_obs_t = f_t * dp.classical_doppler_factor(t_em, **orb_novec)
    log("doppler_freq_swing_Hz", float(np.ptp(f_obs_t)))
    log("intrinsic_freq_drift_Hz", float(f_t[-1] - f_t[0]))

    # phase grid padded by more than the Roemer amplitude on each side, so
    # t_em never falls outside it (the interpolation below would extrapolate
    # silently otherwise)
    pad = 2.0 * np.ptp(roemer_t)
    tt_fine, phase = chirp.phase_of_time(
        lambda tt: chirp.freq_of_time(tt, t_c, system.MCHIRP_MSUN),
        t[0] - pad, t[-1] + pad)
    phase_interp = interp1d(tt_fine, phase, kind="cubic")
    phase_t = phase_interp(t_em)

    D_eff_mpc = system.D_L_PC / 1.0e6
    amp_t = chirp.restricted_pn_amplitude_td(f_t, system.MCHIRP_MSUN, D_eff_mpc)

    z_unlensed = amp_t * np.exp(1j * phase_t)          # analytic signal
    z_lensed = F_t * z_unlensed                          # adiabatic lensing
    h_unlensed = z_unlensed.real
    h_lensed = z_lensed.real

    # The same waveform WITHOUT the Roemer substitution, kept only to
    # measure how much the Doppler term actually changes: the phase, by
    # ~90 cycles; the envelope, by ~1e-5 (doppler.py's docstring, and
    # tests/test_doppler.py::check_envelope_unaffected_by_roemer). The second
    # number is why every amplitude result in this case is unchanged by the
    # 2026-09-07 addition of Doppler, and the first is why the waveform
    # itself is not.
    amp_static = chirp.restricted_pn_amplitude_td(
        chirp.freq_of_time(t, t_c, system.MCHIRP_MSUN), system.MCHIRP_MSUN, D_eff_mpc)
    log("doppler_envelope_max_fractional_change",
        float(np.max(np.abs(amp_t - amp_static) / amp_static)))
    log("doppler_phase_ptp_cycles_vs_static",
        float(np.ptp(phase_t - phase_interp(t)) / (2.0 * np.pi)))

    rms_amp = np.sqrt(np.mean(np.abs(z_lensed) ** 2) / np.mean(np.abs(z_unlensed) ** 2))
    log("rms_amplification_over_observation", float(rms_amp))

    # The lensing's own timing, for the like-for-like comparison against the
    # Roemer delay above (Figure 5). Two different, both legitimate, ways to
    # state it, kept separate rather than conflated:
    #   * Delta_T(y(t)): the physical, convention-independent delay BETWEEN
    #     the two images -- the quantity this whole case is about, and the
    #     one that does not depend on how F is phase-referenced.
    #   * arg F/(2 pi f_B): the phase the lens actually writes onto the
    #     observed waveform at this frequency, expressed as a time. Small
    #     here because F is referenced to the strong image's own arrival
    #     (waveoptics.py), i.e. this is the RESIDUAL timing after the strong
    #     image, not the image separation.
    # Two different, both legitimate, statements of "how much timing does the
    # lens contribute", kept separate rather than conflated:
    #   * Delta_T(y) at the pulse -- the physical, convention-independent
    #     delay BETWEEN the two images. Reported as a number only. It is NOT
    #     plotted: Delta_T increases with y, so it formally diverges toward
    #     each D_LS -> 0 crossing, where the second image's magnification
    #     goes to zero -- that divergence is the second image ceasing to
    #     exist, not a large observable delay, and a raw Delta_T(t) curve
    #     puts it on screen as if it were the dominant effect, which is the
    #     exact opposite of the finding.
    #   * arg F/(2 pi) -- the phase the lens actually writes onto the
    #     observed waveform. Bounded, goes to zero where the lensing does,
    #     and directly comparable to the Roemer phase. This is what Figure 5
    #     plots.
    t_char = 4.0 * units.msun_to_seconds(system.M_LENS_MSUN)  # 4GM/c^3
    i_ymin = int(np.argmin(np.where(lensed_mask, y_t, np.inf)))
    image_dT_at_pulse = float(wo.time_delay_difference(y_t[i_ymin]) * t_char)
    lens_phase_cycles = np.where(lensed_mask, np.angle(F_t) / (2.0 * np.pi), 0.0)
    log("image_delay_at_pulse_s", image_dT_at_pulse)
    log("lens_phase_ptp_cycles", float(np.ptp(lens_phase_cycles)))
    log("roemer_over_image_delay_ratio", float(np.ptp(roemer_t) / image_dT_at_pulse))
    log("roemer_over_lens_phase_ratio",
        float(system.F_B_HZ * np.ptp(roemer_t) / np.ptp(lens_phase_cycles)))

    # --- one approximation this case makes, quantified rather than assumed --
    # F is evaluated at the fixed w_B = w(F_B_HZ), i.e. at the nominal
    # emitted frequency: neither the 4% intrinsic drift nor the Doppler shift
    # is fed back into the amplification factor. D'Orazio & Loeb (2020) make
    # the same omission and give the reason -- lensing only happens where the
    # line-of-sight velocity crosses zero, which for this near-edge-on
    # circular orbit is exactly the orbital phase of the lensing pulse, so
    # the Doppler shift is at its minimum precisely where F matters most.
    # That is an argument, not a number, so here is the number: |F| at the
    # true instantaneous OBSERVED frequency versus |F| at fixed w_B, over
    # every lensed sample.
    # Two numbers, because they differ by 50x and only reporting the smaller
    # one would be the same overclaim the D'Orazio-Loeb argument invites: AT
    # the pulse the error really is tiny (their argument holds), but the
    # worst case over the whole lensed half is larger, because far from the
    # pulse |F| sits on the steep flank of a diffraction fringe where a small
    # shift in w moves it comparatively much more -- on a value that is close
    # to 1 and carries no weight in any result quoted here.
    w_obs_t = wo.w_of_frequency(f_obs_t, system.M_LENS_MSUN)
    abs_F_true = np.array([abs(wo.F_hybrid(wj, yj))
                           for wj, yj in zip(w_obs_t[lensed_mask], y_t[lensed_mask])])
    abs_F_fixed = np.abs(F_t[lensed_mask])
    relerr_w = np.abs(abs_F_true - abs_F_fixed) / abs_F_fixed
    log("abs_F_relerr_from_fixed_w_max", float(np.max(relerr_w)))
    log("abs_F_relerr_from_fixed_w_at_pulse",
        float(relerr_w[int(np.argmax(abs_F_fixed))]))

    # ================= Figures ============================================
    # (A y(t)-over-time panel was tried here first -- the divergence at each
    # D_LS=0 crossing dominates the plot and left little room to actually
    # see anything else, so it added no real information beyond what
    # fraction_of_time_lensed already shows;
    # dropped rather than kept for its own sake.)
    one_p = t <= 1.5 * system.P_OUT_S  # used below for the orbit-track figure

    # Fig 1 (was Fig 2): |F(t)|^2 over the whole observation -- repeated lensing pulses,
    # plus a one-period zoom: the diffraction ringing flanking each pulse is
    # genuine (not aliasing -- re-checked on a 40x finer grid, agrees
    # point-for-point) but is compressed into a few pixels at the full
    # 12-day width, so it needs its own panel to be visible at all.
    abs_F2 = np.abs(F_t) ** 2
    fig, axes = plt.subplots(2, 1, figsize=(8, 6.0))
    axes[0].plot(t / 86400.0, abs_F2, color="#c0392b", lw=0.8)
    axes[0].set_xlabel("$t$ [días]")
    axes[0].set_ylabel(r"$|F(w_B,y(t))|^2$")
    axes[0].set_title(f"Caso B: pulsos de lensing repetido a lo largo de {n_periods:.0f} períodos externos "
                       f"($w_B={w_B:.2f}$)")

    i_peak = int(np.argmax(abs_F2))
    t_peak_days = t[i_peak] / 86400.0
    # Dedicated fine grid for the zoom, not a slice of the coarse one. The
    # observation grid is 600 samples per outer period (576 s apart), which is
    # right for three pulses across 12 days but leaves only ~300 points across
    # this +-1 day window -- and the ringing oscillates fast enough there that
    # the curve came out visibly polygonal. 6000 points is 29 s apart and costs
    # about three seconds of F evaluation.
    t_zoom = np.linspace(t[i_peak] - 86400.0, t[i_peak] + 86400.0, 6000)
    y_zoom, lensed_zoom, _ = geo.impact_parameter_of_time(
        t_zoom, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
        system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)
    abs_F2_zoom = np.ones_like(t_zoom)
    near = lensed_zoom & (y_zoom < 60.0)
    abs_F2_zoom[near] = np.abs(
        [complex(wo.F_hybrid(w_B, yy)) for yy in y_zoom[near]]) ** 2
    axes[1].plot((t_zoom - t[i_peak]) / 86400.0, abs_F2_zoom,
                 color="#c0392b", lw=1.0)
    axes[1].set_xlabel(r"$t - t_\mathrm{pico}$ [días]")
    axes[1].set_ylabel(r"$|F(w_B,y(t))|^2$")
    axes[1].set_title("Ampliación: un pulso — el anillado de difracción que lo flanquea")
    fig.tight_layout()
    fig.savefig(OUT / "caseB_repeated_pulses.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # Fig 3: extended diffraction pattern with the orbit trajectory overlaid.
    # y_max from min_y, NOT max_y: max(y_t) is dominated by the near-
    # divergent points close to the D_LS=0 crossing (Figure 1), which would
    # zoom the whole interesting pattern down to an invisible dot (caught by
    # eye on the first version -- see wiki/log.md). 8*min_y (~12.7) was
    # still too tight, though: y(t) only sits that close to y_min for a
    # narrow sliver of the lensed half -- the MEDIAN y among lensed samples
    # is ~25.5, so the interactive marker (report.html) built on that
    # y_max sat clamped at the frame edge 78% of the orbit (caught by an
    # independent review; see wiki/log.md). y truly is unbounded toward
    # each crossing, so no finite y_max keeps the marker on-frame the whole
    # orbit -- this just widens the window enough that it is on-frame for
    # roughly the lensed half's median, not just its sharpest dip.
    y_max = 25.0
    r_grid = np.linspace(1e-3, y_max, 2000)
    F_r = np.array([wo.F_hybrid(w_B, r) for r in r_grid])
    mag2_r = np.abs(F_r) ** 2
    interp = interp1d(r_grid, mag2_r, bounds_error=False, fill_value=(mag2_r[0], mag2_r[-1]))
    axis = np.linspace(-y_max, y_max, 361)
    Y1, Y2 = np.meshgrid(axis, axis)
    R = np.hypot(Y1, Y2)
    pattern = interp(R)

    r, nu = geo.relative_separation(t[one_p], system.A_OUT_M / units.PC_SI, system.E_OUT,
                                     system.P_OUT_S, system.T_PERI_OUT)
    x_sky, y_sky, z_los = geo.orbital_plane_to_sky(r, nu, system.LITTLE_OMEGA_OUT,
                                                    np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT)
    # theta_E(t), NOT a single fixed/median value: it scales as sqrt(D_LS(t))
    # (geometry.py, impact_parameter_of_time docstring) and varies
    # substantially over the orbit -- a fixed median theta_E, applied to
    # every point of the trajectory, distorted the plotted track by ~19% at
    # the point that matters most (peak lensing, where theta_E is largest).
    # Caught by an independent review; see wiki/log.md. Recompute it exactly
    # on this same t[one_p] grid (same formula as geometry.py, using
    # |z_los| rather than geometry.py's own lensed-only masking: this plot
    # is illustrative and shows BOTH halves of the orbit continuously, and
    # |D_LS| gives a smooth, well-defined theta_E(t) for the unlensed half
    # too, purely for a consistent common length scale to plot in -- no
    # lensing is claimed to happen there, per the D_LS<=0 discussion in
    # theory.tex Sec. 4.2).
    # Floored away from exactly zero: at the two front/back crossings each
    # orbit D_LS_traj_pc=0 exactly (a measure-zero set of samples on this
    # grid), which would divide by zero and emit non-finite NaN/Infinity
    # into the exported JSON below (not valid JSON, though it happens to
    # work as inlined JavaScript, where those are valid identifiers -- an
    # independent review caught this; see wiki/log.md). The floor only ever
    # affects those isolated samples; theta_E(t) is still continuous and
    # correct everywhere else.
    D_LS_traj_pc = np.maximum(np.abs(z_los), 1e-12)
    theta_E_traj = np.sqrt(
        4.0 * units.msun_to_meters(system.M_LENS_MSUN) * (D_LS_traj_pc * units.PC_SI)
        / (system.D_L_PC * units.PC_SI) ** 2
    )
    y1_traj = x_sky / system.D_L_PC / theta_E_traj
    y2_traj = y_sky / system.D_L_PC / theta_E_traj

    lensed_traj = z_los > 0   # still needed by the animation export below

    # Fig 3: the pattern, with axes, colour bar and caption -- but WITHOUT
    # the orbit track that used to be drawn over it. The track is a
    # near-straight line crossing rings whose spacing is the entire content
    # of the figure, and it hid what it was meant to locate. Where the source
    # sits on the pattern is said in time instead, in
    # caseB_repeated_pulses.png, and the animation in report/report.html
    # draws a moving marker over the bare render below.
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.pcolormesh(axis, axis, pattern, shading="auto", cmap="inferno",
                        norm=LogNorm(vmin=max(pattern.min(), 1e-2), vmax=pattern.max()))
    ax.set_aspect("equal")
    ax.set_xlim(-y_max, y_max)
    ax.set_ylim(-y_max, y_max)
    # One marker at the closest approach, not the whole track. This is the
    # single piece of information the removed orbit line carried that the
    # pattern itself does not: where on it the source actually gets to.
    y_min_obs = float(np.min(y_t[lensed_mask]))
    ax.plot(y_min_obs, 0.0, "w+", ms=11, mew=2)
    # Label placed right beside the marker, with no leader line. Anchored up
    # in the corner it needed an arrow, and that arrow was a white diagonal
    # drawn straight across the rings -- doing to the figure exactly what the
    # removed orbit track did.
    ax.annotate(f"$y_\mathrm{{min}}={y_min_obs:.2f}$",
                xy=(y_min_obs, 0.0), xytext=(8, -14),
                textcoords="offset points", color="white", fontsize=8)
    ax.set_xlabel(r"$y_1$")
    ax.set_ylabel(r"$y_2$")
    ax.set_title(f"Caso B: patrón de difracción en el plano de la fuente "
                 f"($w_B={w_B:.2f}$)\n"
                 "anillos concéntricos, con el pico del anillo de Einstein al centro",
                 fontsize=10)
    # Plain decimals on the colour bar. The whole range here is 1.0 to 1.55,
    # inside a single decade, and matplotlib's default log formatter renders
    # that as "1.1 x 10^0", "1.2 x 10^0" -- scientific notation carrying no
    # information, since the exponent is zero on every tick.
    import matplotlib.ticker as mticker
    cb = fig.colorbar(im, ax=ax, shrink=0.85,
                      label=r"$|F|^2$ (escala logarítmica)")
    # Explicit ticks as well as an explicit format: inside a single decade the
    # log locator puts exactly one major tick (at 1) and leaves the rest to the
    # minor formatter, so silencing the minors without setting the majors
    # leaves a colour bar with one label on it.
    lo_cb, hi_cb = im.get_clim()
    ticks = [v for v in np.arange(0.8, 1.61, 0.1) if lo_cb <= v <= hi_cb]
    cb.set_ticks(ticks)
    cb.ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.1f}"))
    cb.ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    fig.tight_layout()
    fig.savefig(OUT / "caseB_pattern.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # Fig 4: idealized detector view -- lensed vs unlensed strain envelope
    fig, axes = plt.subplots(2, 1, figsize=(8, 6))
    axes[0].plot(t / 86400.0, np.abs(z_unlensed), color="#888888", lw=0.8, label="envolvente sin lente")
    axes[0].plot(t / 86400.0, np.abs(z_lensed), color="#2b6cb0", lw=0.8, label="envolvente con lente")
    axes[0].set_xlabel("$t$ [días]")
    axes[0].set_ylabel("envolvente de $h(t)$ [u. arb.]")
    axes[0].set_title("Caso B: vista de detector idealizada sobre toda la observación")
    # Headroom above the data so the legend box has empty space to sit in,
    # rather than covering the pulses/carrier ripple beneath it (it did, at
    # a fixed 'upper right'/'best' location, in an earlier version of this
    # figure -- caught by eye, not by any check, since nothing here tests
    # figure legibility; see wiki/log.md).
    # A tight margin around the data instead of 28% of dead canvas above it:
    # the pulses are the whole point and they were being drawn into the
    # bottom half of the panel. The legend may overlap the curve; that is
    # cheaper than shrinking the feature to avoid it.
    lo = float(np.min(np.abs(z_lensed)))
    hi = float(np.max(np.abs(z_lensed)))
    pad = 0.10 * (hi - lo)
    axes[0].set_ylim(lo - pad, hi + pad)
    axes[0].legend(fontsize=8, loc="upper right")

    # The coarse t-grid above (1800 samples over 12 d, ~576 s spacing) is
    # fine for the slowly-varying envelope but wildly undersamples the
    # 1/f_B=20 s carrier -- plotting "raw waveform" from it just aliases
    # (an earlier version of this figure did exactly that; see wiki/log.md).
    # Build a dedicated fine grid around the first pulse instead.
    i_pulse = int(np.argmax(np.abs(F_t)))
    t_pulse = t[i_pulse]
    # +-2500 s, not +-200 s. The pulse is 6653 s wide (FWHM), so over 200 s
    # |F| changes by 0.2% of its peak height: the panel resolved carrier
    # cycles beautifully and showed an amplification that was, correctly but
    # uselessly, constant across the frame. At 2500 s the envelope moves
    # through 28% of the pulse and the amplification visibly grows and
    # fades, at the cost of the individual cycles becoming dense -- which is
    # unavoidable, since the pulse is 333 carrier cycles wide and no single
    # window can resolve both. The envelopes are drawn over the carrier so
    # the pulse shape reads even where the cycles merge.
    t_fine = np.linspace(t_pulse - 50.0, t_pulse + 50.0, 4000)   # dt=0.025 s
    y_fine, lensed_fine, _ = geo.impact_parameter_of_time(
        t_fine, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
        system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)
    F_fine = np.array([wo.F_hybrid(w_B, y) if lm else 1.0 + 0.0j
                        for y, lm in zip(y_fine, lensed_fine)])
    # Retarded emission time here too, for the same reason as the coarse grid
    # above -- and it matters more here, not less: this panel resolves
    # individual 20 s carrier cycles, and the Roemer delay slews by
    # beta_los*600 s ~ 9.9 s across this window, half a cycle.
    # Using the observation time directly would draw a carrier at the wrong
    # instantaneous frequency.
    t_em_fine, _ = dp.emission_time(t_fine, **orb)
    f_fine = chirp.freq_of_time(t_em_fine, t_c, system.MCHIRP_MSUN)
    phase_fine = phase_interp(t_em_fine)
    amp_fine = chirp.restricted_pn_amplitude_td(f_fine, system.MCHIRP_MSUN, D_eff_mpc)
    z_u_fine = amp_fine * np.exp(1j * phase_fine)
    z_l_fine = F_fine * z_u_fine

    tf = t_fine - t_pulse
    # Carriers only, no envelopes. A pulse is 333 carrier cycles wide, so any
    # window that resolves cycles has a flat |F| and the envelopes over it are
    # just horizontal rules, and any window where they move has the cycles
    # merged into a band. The pulse SHAPE is caseB_repeated_pulses.png's job;
    # this panel is the carrier, and what it shows is the crest height.
    # Lensed first, unlensed over it, so the smaller stays visible inside.
    axes[1].plot(tf, z_l_fine.real, color="#2b6cb0", lw=0.9, label="con lente")
    axes[1].plot(tf, z_u_fine.real, color="#777777", lw=0.9, label="sin lente")
    axes[1].set_xlabel("tiempo respecto del pico del pulso [s]")
    axes[1].set_ylabel("$h(t)$ [u. arb.]")
    axes[1].legend(fontsize=8, loc="upper right")
    # What this panel shows is the AMPLITUDE difference: max_abs_F = 1.1766,
    # i.e. 17.7% taller crests, and that is the whole of it.
    #
    # It does not show a dephasing, and -- worth being exact about, since an
    # earlier title claimed it did -- it could not show one even if there were
    # a bigger one to show. Both curves are evaluated at the same retarded
    # time, so the Roemer delay is common to them and cancels; what is left is
    # F, whose phase runs over lens_phase_ptp_cycles = 0.0174 cycles across
    # the WHOLE observation. At this scale that is sub-pixel by construction,
    # so crests that line up here are not evidence of anything. The evidence
    # that the lens writes amplitude and not phase is the number itself, 0.017
    # cycles against Roemer's 90.5, and its figure is
    # caseB_doppler_vs_lensing.png.
    # Nor is it a zoom on the pulse SHAPE: that (the sinc-like diffraction
    # ringing) lives on a ~day timescale and is caseB_repeated_pulses.png.
    # Why +-50 s, i.e. FIVE cycles. The window has to be narrow enough that a
    # 17.7% difference in crest height is legible, and that is a stronger
    # constraint than it looks. At +-300 s (30 cycles, the previous choice)
    # the two curves cross 60 times across the frame and read as one sinusoid
    # drawn twice -- the figure looked like it was UNDERSTATING the
    # amplification, which is how this got caught. At +-50 s each cycle gets
    # ~270 px and the gap between the crests is unmistakable.
    #
    # Nothing is given up by narrowing, because |F| was already flat: it moves
    # 0.064% of its peak across +-300 s and 0.002% across +-50 s. Going the
    # other way, to where |F| really does move (+-2500 s, 4.4%), drops the
    # cycles to 5 px and the carrier fills in as a solid block. The pulse is
    # 333 cycles wide, so no window shows both; its shape belongs to
    # caseB_repeated_pulses.png.
    # Titles name the panel and stop there; the rest is cases/FIGURES.md.
    axes[1].set_title("Ampliación en el pico: las crestas con lente son un "
                       "$%.1f\\%%$ más altas"
                       % (100 * (abs(F_t[i_pulse]) - 1.0)), fontsize=10)
    # Centred on zero, with a symmetric margin: the curves are a symmetric
    # oscillation, so an asymmetric frame reads as an offset that is not
    # there. No legend box to displace them (labels are on the curves).
    span = 1.18 * float(np.max(np.abs(z_l_fine.real)))
    axes[1].set_ylim(-span, span)
    fig.tight_layout()
    fig.savefig(OUT / "caseB_detector_view.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # Fig 5: the two things the SAME outer orbit does to the waveform, on the
    # same axes and in the same units -- the honest version of "what does the
    # outer orbit imprint". Top: timing, in seconds. Bottom: the phase each
    # one actually writes at f_B, in GW cycles. The point of the figure is
    # the vertical scale, which is why the top panel is symlog: the Roemer
    # curve is off the chart of the lensing one by ~2.5 orders of magnitude,
    # and a linear axis would show the lensing curve as a flat line at zero
    # (which, at this scale, is very nearly the honest answer).
    # Deliberately NOT a log or symlog y-axis: the Roemer delay is a pure
    # sinusoid, and symlog renders it as a flat-topped square wave -- a plot
    # artefact that misstates the shape of the very curve the figure is
    # about. One linear panel at the Roemer scale, where the lensing delay
    # correctly sits on the zero line (that IS the finding), plus a small
    # inset at ~150x so the lensing curve is still shown resolved rather than
    # only shown to be invisible.
    # Both curves in GW CYCLES, the unit the waveform is actually written in,
    # so the comparison is like-for-like and both are bounded.
    # One curve per panel, sharing the time axis, instead of two curves on one
    # axis where the second is necessarily a flat line at zero. The point of
    # the figure is the RATIO of two phases, and a ratio of 5206 cannot be
    # drawn on a common axis: either the small curve is invisible or the large
    # one is off-screen. Stacked panels give each its own scale and keep the
    # comparison honest, with the ratio stated rather than drawn.
    fig, axes = plt.subplots(2, 1, figsize=(8, 5.4), sharex=True)
    tp = t[one_p] / 86400.0
    roemer_cycles = system.F_B_HZ * roemer_t

    # Conjunction: where y(t) is smallest, i.e. where the pulses are. It is
    # also the extremum of z_src, hence of the Roemer delay -- which is why
    # the lensing lines up with the Roemer PEAK and not with its zero: the
    # delay is extremal exactly where the line-of-sight VELOCITY vanishes.
    # That coincidence is the reason the fixed-w approximation is safe here,
    # so the figure marks it instead of leaving it to look like an accident.
    t_conj = float(tp[int(np.argmax(np.abs(F_t[one_p])))])
    for k, axk in enumerate(axes):
        axk.axvline(t_conj, color='0.55', lw=0.8, ls=':',
                    label='conjunción' if k == 0 else None)
        axk.axhline(0.0, color="0.85", lw=0.6)

    axes[0].plot(tp, roemer_cycles[one_p], color="#6b46c1", lw=1.6)
    axes[0].set_ylabel("Roemer / Doppler\n[ciclos]")
    axes[0].set_title("Caso B: las dos huellas de la misma órbita, a escala propia\n"
                      f"{NUMBERS['roemer_phase_ptp_cycles']:.0f} ciclos contra "
                      f"{NUMBERS['lens_phase_ptp_cycles']:.3f} — "
                      f"{NUMBERS['roemer_over_lens_phase_ratio']:.0f}$\\times$",
                      fontsize=10)

    axes[1].plot(tp, lens_phase_cycles[one_p], color="#c0392b", lw=1.4)
    axes[1].set_ylabel("lensing, $\\arg F/2\\pi$\n[ciclos]")
    axes[1].set_xlabel("$t$ [días]")
    # The dotted line marks conjunction, and it is in a legend rather than a
    # floating caption: the pulses land on the Roemer MAXIMUM, not on its
    # zero, because the delay is extremal exactly where the line-of-sight
    # velocity vanishes. That is said in cases/FIGURES.md, not written across
    # the plot.
    axes[0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "caseB_doppler_vs_lensing.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # Fig 6: pattern only (no track), fixed square canvas -- background for
    # the animated visualization in report/report.html. A clean render (no
    # overlaid line) so the JS animation draws its own marker on top without
    # duplicating a static track underneath it.
    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.pcolormesh(axis, axis, pattern, shading="auto", cmap="inferno",
                        norm=LogNorm(vmin=max(pattern.min(), 1e-2), vmax=pattern.max()))
    ax.set_aspect("equal")
    ax.set_xlim(-y_max, y_max)
    ax.set_ylim(-y_max, y_max)
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    # Deliberately NOT bbox_inches="tight" here, unlike every other figure in
    # this repo: this PNG is the background the report.html animation draws
    # its marker on, positioned by percentage of the image. It must stay an
    # exact square with zero margin (subplots_adjust above), or the marker
    # silently drifts off the pattern it is supposed to be tracking.
    fig.savefig(OUT / "caseB_pattern_only.png", dpi=200)
    plt.close(fig)

    # ---- Animation data for report/report.html (JSON, reduced resolution) -
    P_days = system.P_OUT_S / 86400.0
    anim = {
        "y_max": float(y_max),
        "P_days": float(P_days),
        "w_B": float(w_B),
        "orbit": {
            "t_days": (t[one_p] / 86400.0).tolist(),
            "y1": y1_traj.tolist(),
            "y2": y2_traj.tolist(),
            "lensed": lensed_traj.tolist(),
        },
    }
    # envelope + pulses: subsample the full 3-period arrays to ~900 points
    stride = max(1, n_samples // 900)
    anim["envelope"] = {
        "t_days": (t[::stride] / 86400.0).tolist(),
        "unlensed": np.abs(z_unlensed[::stride]).tolist(),
        "lensed": np.abs(z_lensed[::stride]).tolist(),
        "F2": (np.abs(F_t[::stride]) ** 2).tolist(),
    }
    # The kinematic (Roemer/Doppler) imprint of the same orbit, in seconds and
    # in GW cycles, on the same subsampled grid. An additional key, not a
    # changed one: report.html's existing JS reads only the keys above and is
    # unaffected by its presence.
    anim["kinematics"] = {
        "roemer_s": roemer_t[::stride].tolist(),
        "roemer_cycles": (system.F_B_HZ * roemer_t[::stride]).tolist(),
        "lens_phase_cycles": lens_phase_cycles[::stride].tolist(),
        "beta_los": beta_los_t[::stride].tolist(),
    }
    with open(OUT / "caseB_animation_data.json", "w") as fh:
        json.dump(anim, fh)

    with open(OUT / "provenance" / "numbers.json", "w") as fh:
        json.dump(NUMBERS, fh, indent=2)
    print(json.dumps(NUMBERS, indent=2))
    print("\nFigures written: caseB_repeated_pulses.png, "
          "caseB_pattern.png, caseB_pattern_only.png, "
          "caseB_detector_view.png, caseB_doppler_vs_lensing.png, "
          "caseB_animation_data.json")


if __name__ == "__main__":
    main()
