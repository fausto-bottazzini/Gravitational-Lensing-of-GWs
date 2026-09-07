"""Case B -- quasi-monochromatic source, moving lens.

Same triple as Case A (src/gwlens/system.py), earlier in the SAME inner
binary's inspiral: f=F_B_HZ=0.05 Hz, changing by only ~4% over the
T_OBS_B_S=6-outer-period (24 d) observation window used here (checked in
tests/test_system.py::check_quasi_monochromatic_regime). Over that window
the outer orbit is NOT frozen: the impact parameter y(t) sweeps through the
diffraction pattern derived in Case A, producing periodic amplification
pulses once per outer period -- "repeated lensing" (D'Orazio & Loeb 2020),
here from wave optics rather than their geometric-optics treatment.

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

    # --- underlying quasi-monochromatic waveform (complex analytic signal) -
    t_c = chirp.time_to_merger(system.F_B_HZ, system.MCHIRP_MSUN)
    f_t = chirp.freq_of_time(t, t_c, system.MCHIRP_MSUN)
    log("f_start_Hz", float(f_t[0]))
    log("f_end_Hz", float(f_t[-1]))
    log("fractional_freq_drift", float((f_t[-1] - f_t[0]) / f_t[0]))
    tt_fine, phase = chirp.phase_of_time(
        lambda tt: chirp.freq_of_time(tt, t_c, system.MCHIRP_MSUN), 0.0, t[-1])
    phase_t = interp1d(tt_fine, phase, kind="cubic")(t)

    D_eff_mpc = system.D_L_PC / 1.0e6
    amp_t = chirp.restricted_pn_amplitude_td(f_t, system.MCHIRP_MSUN, D_eff_mpc)

    z_unlensed = amp_t * np.exp(1j * phase_t)          # analytic signal
    z_lensed = F_t * z_unlensed                          # adiabatic lensing
    h_unlensed = z_unlensed.real
    h_lensed = z_lensed.real

    rms_amp = np.sqrt(np.mean(np.abs(z_lensed) ** 2) / np.mean(np.abs(z_unlensed) ** 2))
    log("rms_amplification_over_observation", float(rms_amp))

    # ================= Figures ============================================

    # Fig 1: y(t) over ~1.5 outer periods, lensed/unlensed shading
    fig, ax = plt.subplots(figsize=(8, 3.2))
    one_p = t <= 1.5 * system.P_OUT_S
    days = t[one_p] / 86400.0
    y_show = np.where(lensed_mask[one_p], y_t[one_p], np.nan)
    ax.plot(days, y_show, color="#2b6cb0", lw=1.2, label="y(t), lensed half")
    ax.fill_between(days, 0, 1, where=~lensed_mask[one_p], transform=ax.get_xaxis_transform(),
                     color="grey", alpha=0.15, label="unlensed half (source in front)")
    ax.set_xlabel("t [days]")
    ax.set_ylabel("impact parameter y(t)")
    ax.set_ylim(0, 10 * NUMBERS["min_y_during_observation"])  # y->inf as D_LS->0+; crop to the interesting dip
    ax.set_title("Case B: impact parameter over 1.5 outer periods\n"
                 "(y formally diverges at the D_LS=0 crossing, cropped here)")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "caseB_y_of_t.png", dpi=170)
    plt.close(fig)

    # Fig 2: |F(t)|^2 over the whole observation -- repeated lensing pulses,
    # plus a one-period zoom: the diffraction ringing flanking each pulse is
    # genuine (not aliasing -- re-checked on a 40x finer grid, agrees
    # point-for-point) but is compressed into a few pixels at the full
    # 24-day width and was not actually visible in a single-panel version of
    # this figure despite RESULTS.md/report.tex describing it as such; an
    # independent review caught this, see wiki/log.md.
    abs_F2 = np.abs(F_t) ** 2
    fig, axes = plt.subplots(2, 1, figsize=(8, 6.0))
    axes[0].plot(t / 86400.0, abs_F2, color="#c0392b", lw=0.8)
    axes[0].set_xlabel("t [days]")
    axes[0].set_ylabel(r"$|F(w_B,y(t))|^2$")
    axes[0].set_title(f"Case B: repeated lensing pulses over {n_periods:.0f} outer periods "
                       f"(w_B={w_B:.2f})")

    i_peak = int(np.argmax(abs_F2))
    t_peak_days = t[i_peak] / 86400.0
    zoom = (np.abs(t / 86400.0 - t_peak_days) < 1.0)
    axes[1].plot(t[zoom] / 86400.0 - t_peak_days, abs_F2[zoom], color="#c0392b", lw=1.0)
    axes[1].set_xlabel("t - t_peak [days]")
    axes[1].set_ylabel(r"$|F(w_B,y(t))|^2$")
    axes[1].set_title("Zoom: one pulse -- the diffraction ringing flanking it")
    fig.tight_layout()
    fig.savefig(OUT / "caseB_repeated_pulses.png", dpi=170)
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

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.pcolormesh(axis, axis, pattern, shading="auto", cmap="inferno",
                        norm=LogNorm(vmin=max(pattern.min(), 1e-2), vmax=pattern.max()))
    lensed_traj = z_los > 0
    ax.plot(y1_traj[lensed_traj], y2_traj[lensed_traj], "-", color="cyan", lw=1.5,
            label="source track (lensed half)")
    ax.plot(y1_traj[~lensed_traj], y2_traj[~lensed_traj], "--", color="lime", lw=1.0,
            label="source track (unlensed half, in front of lens)")
    ax.set_aspect("equal")
    ax.set_xlim(-y_max, y_max)
    ax.set_ylim(-y_max, y_max)
    ax.set_xlabel(r"$y_1$")
    ax.set_ylabel(r"$y_2$")
    ax.set_title(f"Case B: fixed diffraction pattern (w_B={w_B:.2f}) with the\n"
                 "source's outer-orbit track over 1.5 periods")
    fig.colorbar(im, ax=ax, shrink=0.85, label=r"$|F|^2$ (log scale)")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "caseB_pattern_with_orbit.png", dpi=170)
    plt.close(fig)

    # Fig 4: idealized detector view -- lensed vs unlensed strain envelope
    fig, axes = plt.subplots(2, 1, figsize=(8, 6))
    axes[0].plot(t / 86400.0, np.abs(z_unlensed), color="#888888", lw=0.8, label="unlensed envelope")
    axes[0].plot(t / 86400.0, np.abs(z_lensed), color="#2b6cb0", lw=0.8, label="lensed envelope")
    axes[0].set_xlabel("t [days]")
    axes[0].set_ylabel("strain envelope [arb. units]")
    axes[0].set_title("Case B: idealized detector view over the full observation")
    axes[0].legend(fontsize=8)

    # The coarse t-grid above (3600 samples over 24 d, ~576 s spacing) is
    # fine for the slowly-varying envelope but wildly undersamples the
    # 1/f_B=20 s carrier -- plotting "raw waveform" from it just aliases
    # (an earlier version of this figure did exactly that; see wiki/log.md).
    # Build a dedicated fine grid around the first pulse instead.
    i_pulse = int(np.argmax(np.abs(F_t)))
    t_pulse = t[i_pulse]
    t_fine = np.linspace(t_pulse - 200.0, t_pulse + 200.0, 4000)  # dt=0.1s, >>Nyquist for f_B=0.05Hz
    y_fine, lensed_fine, _ = geo.impact_parameter_of_time(
        t_fine, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
        system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)
    F_fine = np.array([wo.F_hybrid(w_B, y) if lm else 1.0 + 0.0j
                        for y, lm in zip(y_fine, lensed_fine)])
    f_fine = chirp.freq_of_time(t_fine, t_c, system.MCHIRP_MSUN)
    phase_fine = 2.0 * np.pi * f_fine * (t_fine - t_fine[0]) + phase_t[i_pulse]
    amp_fine = chirp.restricted_pn_amplitude_td(f_fine, system.MCHIRP_MSUN, D_eff_mpc)
    z_u_fine = amp_fine * np.exp(1j * phase_fine)
    z_l_fine = F_fine * z_u_fine

    axes[1].plot(t_fine - t_pulse, z_u_fine.real, color="#888888", lw=0.7, label="unlensed")
    axes[1].plot(t_fine - t_pulse, z_l_fine.real, color="#2b6cb0", lw=0.7, alpha=0.85, label="lensed")
    axes[1].set_xlabel(f"t - {t_pulse/86400:.3f} d  [s]")
    axes[1].set_ylabel("h(t) [arb. units]")
    axes[1].set_title("Zoom: raw waveform through the peak of one lensing pulse (fine time grid)")
    axes[1].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "caseB_detector_view.png", dpi=170)
    plt.close(fig)

    # Fig 5: pattern only (no track), fixed square canvas -- background for
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
    # envelope + pulses: subsample the full 6-period arrays to ~900 points
    stride = max(1, n_samples // 900)
    anim["envelope"] = {
        "t_days": (t[::stride] / 86400.0).tolist(),
        "unlensed": np.abs(z_unlensed[::stride]).tolist(),
        "lensed": np.abs(z_lensed[::stride]).tolist(),
        "F2": (np.abs(F_t[::stride]) ** 2).tolist(),
    }
    with open(OUT / "caseB_animation_data.json", "w") as fh:
        json.dump(anim, fh)

    with open(OUT / "provenance" / "numbers.json", "w") as fh:
        json.dump(NUMBERS, fh, indent=2)
    print(json.dumps(NUMBERS, indent=2))
    print("\nFigures written: caseB_y_of_t.png, caseB_repeated_pulses.png, "
          "caseB_pattern_with_orbit.png, caseB_pattern_only.png, "
          "caseB_detector_view.png, caseB_animation_data.json")


if __name__ == "__main__":
    main()
