"""Checks that the chosen system parameters (src/gwlens/system.py) actually
satisfy the regime assumptions stated in wiki/conventions.md. If these fail,
the whole "static lens" / "quasi-monochromatic" story is not justified for
the numbers actually used in cases/, so this gates them.
Run with `python3 tests/test_system.py`.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from gwlens import system, chirp, geometry as geo, units
from gwlens import waveoptics as wo


def check_hierarchy():
    """The triple has to actually BE a hierarchical triple -- tested against
    the real criteria, not against the `a_in/a_out < 1e-3` proxy this check
    used until 2026-09-07. That proxy is not the stability condition, and a
    system could satisfy it while being dynamically unstable or having its
    inner binary tidally unbound by the tertiary; both would invalidate the
    whole setup, not just a detail of it. Three things, at both epochs:

      1. Dynamical stability, Mardling & Aarseth (2001):
             a_out/a_in > 2.8 [(1+q_out)(1+e_out)/sqrt(1-e_out)]^(2/5),
         q_out = m3/(m1+m2) = M_lens/M_binary, the tertiary over the
         inner pair -- the lens IS the tertiary here. Margin 93x at
         Case B (the inverted ratio would read ~1700x, which is how to
         recognise it if this ever regresses).
      2. The inner binary sits well inside its own Hill radius about the
         tertiary, a_in << a_out (M_b/3M_L)^(1/3) -- i.e. it is bound to
         itself, not tidally stripped. Margin ~290x at Case B.
      3. The inner binary is a POINT SOURCE to the lens, in both of the
         two inequivalent senses. theory.pdf Sec. 1.3 derives a_in << b
         (no differential deflection across the binary) and then says in
         as many words that this is "the weaker of the two": in wave
         optics what erases the diffraction pattern is not the images
         separating but different points of the source producing patterns
         offset from one another, so the size to beat is the EINSTEIN
         RADIUS, not the impact parameter. theory.pdf asserts that
         stronger condition holds "with room to spare" and gives no
         number; this is the number. Projected into the source plane the
         Einstein radius is eta_0 = theta_E * D_S, evaluated at
         conjunction, where the lensing actually happens. a_in/eta_0 is
         6.4e-3 at Case B against a_in/a_out of 2.1e-4 -- the governing
         condition is ~30x tighter than the one that used to stand in for
         it, and it is the one that rules the regime this whole project
         works in.
    """
    # q_out = m3/(m1+m2), the definition in Mardling & Aarseth (2001) Eq. 90.
    q_out = system.M_LENS_MSUN / system.MTOT_MSUN
    e = system.E_OUT
    crit = 2.8 * ((1.0 + q_out) * (1.0 + e) / np.sqrt(1.0 - e)) ** 0.4
    r_hill = system.A_OUT_M * (system.MTOT_MSUN / (3.0 * system.M_LENS_MSUN)) ** (1.0 / 3.0)

    # Einstein radius in the SOURCE plane at conjunction: eta_0 = theta_E * D_S,
    # with D_S = D_L + D_LS and D_LS the line-of-sight projection geometry.py uses.
    t = np.linspace(0.0, system.P_OUT_S, 4001, endpoint=False)
    y_t, lensed, theta_E = geo.impact_parameter_of_time(
        t, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
        system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)
    k = int(np.argmin(np.where(lensed, y_t, np.inf)))  # conjunction
    r_orb, nu = geo.relative_separation(t, system.A_OUT_M, system.E_OUT,
                                       system.P_OUT_S, system.T_PERI_OUT)
    _, _, z_los = geo.orbital_plane_to_sky(r_orb, nu, system.LITTLE_OMEGA_OUT,
                                          np.deg2rad(system.I_OUT_DEG),
                                          system.OMEGA_OUT)
    eta_0 = theta_E[k] * (system.D_L_PC * units.PC_SI + z_los[k])

    worst_stab, worst_hill, worst_point, worst_wave = np.inf, np.inf, 0.0, 0.0
    for f in (system.F_A_START_HZ, system.F_B_HZ):
        a_in = system.inner_separation_m(f, system.M1_MSUN, system.M2_MSUN)
        worst_stab = min(worst_stab, (system.A_OUT_M / a_in) / crit)
        worst_hill = min(worst_hill, r_hill / a_in)
        worst_point = max(worst_point, a_in / system.A_OUT_M)
        worst_wave = max(worst_wave, a_in / eta_0)

    ok = (worst_stab > 10.0 and worst_hill > 10.0
          and worst_point < 1e-3 and worst_wave < 1e-2)
    return ok, (f"Mardling-Aarseth stability: {worst_stab:.0f}x the critical "
                f"a_out/a_in={crit:.2f} at the tighter epoch; inner binary "
                f"{worst_hill:.0f}x inside its Hill radius about the lens; "
                f"point source, geometric a_in/a_out = {worst_point:.2e} and "
                f"wave-optics a_in/eta_0 = {worst_wave:.2e} at conjunction "
                f"(eta_0 = {eta_0:.3e} m), the latter {worst_wave/worst_point:.0f}x "
                f"tighter and the one that governs. Requires 10x margin on the "
                f"first two, <1e-3 and <1e-2 on the last two.")


def check_outer_orbit_is_static_over_the_observation():
    """Everything here models the outer orbit as a FIXED Keplerian ellipse
    and the inner binary as circular. Two things could break that on the
    timescales involved, and neither was checked anywhere before 2026-09-07:

      1. The outer orbit radiates too. Its own GW inspiral time
         tau = (5/256) c^5 a^4 / (G^3 M_L M_b (M_L+M_b)) must be enormously
         longer than anything observed here, or a_out is not fixed.
      2. Kozai-Lidov oscillations would pump the inner binary's eccentricity
         on t_KL ~ (8/15pi)((M_b+M_L)/M_L)(P_out^2/P_in)(1-e_out^2)^(3/2),
         which would invalidate the circular-inspiral waveform. Here that is
         16 yr -- not negligible against the 240 d to merger on its own (4%),
         which is exactly why the second half of the argument matters: KL is
         quenched when the inner binary's own GR periastron precession is
         faster, t_GR = P_in a_in c^2 (1-e^2)/(3 G M_b) << t_KL, and here it
         is faster by ~3e4. Quoting only the first ratio would have made this
         look marginal; quoting only the second would have skipped why it is
         needed.
    """
    G, c, Msun = units.G_SI, units.C_SI, units.M_SUN_SI
    M_L, M_b = system.M_LENS_MSUN * Msun, system.MTOT_MSUN * Msun
    a_out, P_out = system.A_OUT_M, system.P_OUT_S

    tau_out = (5.0 / 256.0) * c ** 5 * a_out ** 4 / (G ** 3 * M_L * M_b * (M_L + M_b))
    tau_inner_merger = chirp.time_to_merger(system.F_B_HZ, system.MCHIRP_MSUN)
    decay_frac = tau_inner_merger / tau_out

    P_in = 2.0 / system.F_B_HZ            # GW frequency is twice the orbital
    t_kl = (8.0 / (15.0 * np.pi)) * ((M_b + M_L) / M_L) * P_out ** 2 / P_in
    a_in = system.inner_separation_m(system.F_B_HZ, system.M1_MSUN, system.M2_MSUN)
    t_gr = P_in * a_in * c ** 2 / (3.0 * G * M_b)

    ok = decay_frac < 1e-4 and (t_kl / t_gr) > 100.0
    return ok, (f"outer orbit's own GW decay time = {tau_out/units.YEAR_SI:.2e} yr, of which "
                f"{decay_frac:.1e} elapses before the inner binary merges (require <1e-4, so "
                f"a_out is fixed); Kozai-Lidov t_KL = {t_kl/units.YEAR_SI:.1f} yr is "
                f"{t_kl/tau_inner_merger:.0f}x the time to merger, and further quenched "
                f"by inner GR precession, t_KL/t_GR = {t_kl/t_gr:.1e} "
                f"(require >100)")


def check_static_lens_regime():
    """Not a check on any code: an assertion that the pinned system sits in
    the regime Case A claims. The merger duration, from f_A_start to formal
    coalescence, must be a tiny fraction of the outer orbital period, so
    that the lens does not move during the observed chirp. Recorded here
    because the report states it; it would fail if someone changed the
    system parameters out from under that claim."""
    tau_A = chirp.time_to_merger(system.F_A_START_HZ, system.MCHIRP_MSUN)
    frac = tau_A / system.P_OUT_S
    ok = frac < 1e-3
    return ok, f"merger duration / outer period = {frac:.2e} (require < 1e-3)"


def check_quasi_monochromatic_regime():
    """The same kind of assertion as the previous check, for Case B: the
    time to merger from f_B must span several outer periods, so the orbital
    modulation is observable at all, while f itself barely changes over the
    observing baseline."""
    tau_B = chirp.time_to_merger(system.F_B_HZ, system.MCHIRP_MSUN)
    n_periods_to_merger = tau_B / system.P_OUT_S
    # fractional frequency drift over the ACTUAL Case B observing baseline
    # (system.T_OBS_B_S = 3 outer periods), not over the full remaining
    # time to merger: f formally diverges as t -> t_c whatever the window,
    # so a drift measured to merger would say nothing about the window.
    f_start = chirp.freq_of_time(np.array([0.0]), tau_B, system.MCHIRP_MSUN)[0]
    f_end = chirp.freq_of_time(np.array([system.T_OBS_B_S]), tau_B, system.MCHIRP_MSUN)[0]
    drift = abs(f_end - f_start) / f_start
    n_obs_periods = system.T_OBS_B_S / system.P_OUT_S
    ok = n_periods_to_merger > 3 * n_obs_periods and drift < 0.05
    return ok, (f"{n_periods_to_merger:.1f} outer periods before merger "
                f"(observing {n_obs_periods:.0f} of them); fractional freq "
                f"drift over the observing baseline = {drift:.2e} "
                f"(require merger>>observing window, drift<5%)")


def check_wave_optics_regime_nontrivial():
    """The two epochs must land on OPPOSITE sides of the optics transition,
    which is the whole reason for studying both. Case B must be in the
    diffractive regime, w_B < 1, where neither geometric optics nor the
    long-wavelength limit applies and the exact F is needed; Case A must be
    above w=30, the threshold at which F_hybrid switches to the asymptotic
    form, for its whole band. The earlier version of this check only asked
    that w be somewhere between 1e-2 and 1e4 at each epoch, which is no
    constraint at all."""
    from gwlens import waveoptics as wo
    w_B = wo.w_of_frequency(system.F_B_HZ, system.M_LENS_MSUN)
    w_A_start = wo.w_of_frequency(system.F_A_START_HZ, system.M_LENS_MSUN)
    w_A_isco = wo.w_of_frequency(system.F_ISCO_HZ, system.M_LENS_MSUN)
    ok = w_B < 1.0 and w_A_start > 30.0 and w_A_isco > w_A_start
    return ok, (f"w_B={w_B:.3f} (require <1, diffractive); Case A band "
                f"w={w_A_start:.1f} to {w_A_isco:.1f} (require >30 throughout, "
                f"geometric optics)")


def check_adiabatic_approximation_where_it_is_applied():
    """Case B evaluates F(w_B, y(t)) instant by instant. The licence for that
    (theory.pdf, the adiabatic condition) is that the configuration changes
    little while the lens RESPONDS, on t_resp ~ 4 G M_L / c^3:

        t_resp << |d ln y / dt|^{-1}.

    theory.pdf writes the right-hand side as "~P_out". That is the scale in
    the middle of the lensed window but it is NOT a bound, and the difference
    matters: at both ends of the window D_LS -> 0, so theta_E -> 0 and
    y -> infinity, and |d ln y/dt|^{-1} -> 0 there for any system whatsoever.
    Minimised over the whole lensed half the quantity is therefore not a
    property of this system at all -- it is a property of the sampling grid,
    and reads 43 s at 2e4 samples per period, 0.26 s at 2e5, 0.03 s at 2e6.
    Reporting any of those as "the adiabatic margin" would be reporting the
    grid spacing.

    Where the condition has to hold is where the prescription does something:
    the adiabatic factor is F=1 to within a part in 1e2 over almost the whole
    orbit, and only near conjunction does the lens modulate anything. So the
    minimum is taken over |F| > 1.01, which is grid-independent because y is
    finite and smooth there, and the margin is ~5800x."""
    from gwlens import waveoptics as wo

    t_resp = 4.0 * units.msun_to_seconds(system.M_LENS_MSUN)
    w_B = wo.w_of_frequency(system.F_B_HZ, system.M_LENS_MSUN)

    margenes = []
    for n in (8001, 20001):
        t = np.linspace(0.0, system.P_OUT_S, n, endpoint=False)
        y_t, lensed, _ = geo.impact_parameter_of_time(
            t, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
            np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
            system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)
        F_mag = np.ones_like(y_t)
        F_mag[lensed] = np.abs([complex(wo.F_hybrid(w_B, yy)) for yy in y_t[lensed]])
        activo = F_mag > 1.01
        ln_y = np.where(lensed, np.log(np.where(lensed, y_t, 1.0)), np.nan)
        tau = np.abs(1.0 / np.gradient(ln_y, t))
        margenes.append(float(np.nanmin(tau[activo]) / t_resp))

    grid_independiente = abs(margenes[1] - margenes[0]) / margenes[0] < 0.05
    ok = min(margenes) > 100.0 and grid_independiente
    return ok, (f"t_resp = 4GM_L/c^3 = {t_resp:.3f} s; where |F|>1.01, "
                f"min |dln y/dt|^-1 / t_resp = {margenes[0]:.0f}x at 8e3 samples and "
                f"{margenes[1]:.0f}x at 2e4 (require >100x, and agreement to 5% "
                f"between the two so the number is the system's and not the grid's)")


def check_deflection_is_weak_where_the_rays_pass():
    """theory.pdf's optics chapter is built to first order in U/c^2 -- the
    Helmholtz equation, the refractive index, and the deflection
    alpha = 4 G M_L/(c^2 b) all drop O(U^2). The text says the expansion
    parameter is small, which is the right thing for a textbook to say;
    whether it is small for THESE numbers is a statement about this repo's
    system, and belongs here (and in report/, not in theory/).

    The subtlety worth recording: the ray does NOT pass at the source offset
    y*xi_0. It passes at the IMAGE position, x_+- * xi_0, and as y -> 0 the
    two images go to the Einstein ring |x| = 1. So the closest approach of
    any ray is xi_0 itself, for any y whatsoever -- small y does not push
    the problem into the strong field. What is left is a fixed fractional
    correction, the second-order Schwarzschild deflection:

        alpha_2/alpha_1 = (15 pi/16) * (G M_L / c^2 b).

    At the Einstein ring (the y -> 0 worst case) that is 2.4%; at the y this
    project actually uses it is 1.2%, since the strong image sits out at
    x_+ = 2.07. It scales as M_L^(-1/3) at fixed outer period, so a heavier
    lens would be worse."""
    r_g = units.msun_to_meters(system.M_LENS_MSUN)          # G M_L / c^2
    i_out = np.deg2rad(system.I_OUT_DEG)
    D_LS = system.A_OUT_M * np.sin(i_out)                   # at conjunction
    xi_0 = np.sqrt(4.0 * r_g * D_LS)          # Einstein radius, LENS plane

    y_min = np.sqrt(system.A_OUT_M / (4.0 * r_g)) * np.cos(i_out) / np.sqrt(np.sin(i_out))
    x_plus, _ = wo.image_positions(y_min)

    worst = (15.0 * np.pi / 16.0) * (r_g / xi_0)            # y -> 0
    used = (15.0 * np.pi / 16.0) * (r_g / (x_plus * xi_0))  # y = y_min

    ok = worst < 0.05 and xi_0 / r_g > 50.0
    return ok, (f"closest any ray comes is the Einstein ring, xi_0 = {xi_0/r_g:.0f} r_g "
                f"(require >50); second-order deflection is {100*worst:.2f}% there "
                f"(the y->0 bound, require <5%) and {100*used:.2f}% at the y actually "
                f"used, y_min={y_min:.3f} -> x_+={x_plus:.3f}")


CHECKS = [
    ("hierarchy", check_hierarchy),
    ("outer_orbit_is_static_over_the_observation", check_outer_orbit_is_static_over_the_observation),
    ("static_lens_regime", check_static_lens_regime),
    ("quasi_monochromatic_regime", check_quasi_monochromatic_regime),
    ("wave_optics_regime_nontrivial", check_wave_optics_regime_nontrivial),
    ("adiabatic_approximation_where_it_is_applied",
     check_adiabatic_approximation_where_it_is_applied),
    ("deflection_is_weak_where_the_rays_pass",
     check_deflection_is_weak_where_the_rays_pass),
]


def main():
    print(system.summary())
    print()
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
    out = Path(__file__).parent / "CHECKS_system.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
