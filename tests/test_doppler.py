"""Checks for src/gwlens/doppler.py -- the line-of-sight kinematics of the
outer orbit (Roemer delay, Doppler, constant orbital redshift).

Run with `python3 tests/test_doppler.py`. Same plain-script, PASS/FAIL,
writes-its-own-JSON convention as the rest of tests/ (no pytest, see
wiki/conventions.md).

Six checks, deliberately. Each one, if it failed, would invalidate a number
that cases/ or report/ actually states; nothing here re-checks algebra for
its own sake (2 a sin(i)/c is 2 a sin(i)/c, a face-on orbit has no
line-of-sight motion, and a fixed-point iteration converges -- none of those
would ever have been wrong in a way that changed a result, so none of them
are here). What IS here: the tie between the closed-form kinematics and the
orbit code the lensing uses, the identity that makes the retarded-time
substitution a complete first-order treatment, the coefficient in the
constant orbital redshift (which no other check would catch, since a
constant redshift makes no waveform feature), and the three quantitative
claims Case A and Case B rest on.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from gwlens import doppler as dp, geometry as geo, system, units, waveoptics as wo, chirp


def _orbit_kwargs():
    """The pinned system's outer orbit, in the argument form doppler.py takes
    (a_out in METRES here, unlike geometry.py's parsecs -- see
    doppler.source_los_offset)."""
    return dict(a_out_m=system.A_OUT_M, e_out=system.E_OUT,
                period_out=system.P_OUT_S, i_out=np.deg2rad(system.I_OUT_DEG),
                Omega_out=system.OMEGA_OUT, omega_out=system.LITTLE_OMEGA_OUT,
                t_peri=system.T_PERI_OUT, M_lens_msun=system.M_LENS_MSUN,
                M_binary_msun=system.MTOT_MSUN)


def _no_Omega(kw):
    """los_velocity/doppler_factor take no Omega -- z_los is Omega-independent
    (geometry.orbital_plane_to_sky: z = r sin(nu+omega) sin i)."""
    return {k: v for k, v in kw.items() if k != "Omega_out"}


def check_los_velocity_matches_the_orbit_code():
    """The closed-form radial-velocity curve v_z = K[cos(nu+omega)+e cos omega]
    must equal d/dt of `source_los_offset`, which is built from
    geometry.py -- the same orbit code the lensing side uses.

    This is the check that ties the textbook formula to THIS repo's orbit
    rather than to an assumed one: a wrong mass fraction, a wrong sin(i), a
    wrong sign of the line-of-sight axis, or a K with the wrong period would
    all show up here and nowhere else. Run at e=0 (the orbit actually used)
    and e=0.4 (so the eccentric branch of the formula is exercised too,
    since wiki/todo.md leaves eccentric outer orbits open)."""
    worst = 0.0
    for e in [0.0, 0.4]:
        kw = _orbit_kwargs()
        kw["e_out"] = e
        t = np.linspace(0.0, kw["period_out"], 400001)
        v_fd = np.gradient(dp.source_los_offset(t, **kw), t)
        v_closed = dp.los_velocity(t, **_no_Omega(kw))
        # trim the endpoints, where np.gradient drops to first order
        worst = max(worst, np.max(np.abs(v_fd - v_closed)[5:-5]) / np.max(np.abs(v_closed)))
    ok = worst < 1e-6
    return ok, f"max relative |d(z_src)/dt - closed-form v_z| over e in [0,0.4] = {worst:.2e} (tol 1e-6)"


def check_retarded_time_is_the_whole_first_order_doppler():
    """d t_em / d t_obs must equal 1/(1+beta_los) evaluated AT THE EMISSION
    TIME. This is what licenses the cases to apply Doppler by evaluating the
    waveform at `emission_time` and multiplying by nothing else: if the
    identity did not hold, a separate frequency-shift factor would be
    missing, and Case B's 90-cycle phase modulation would be wrong.

    Evaluated at t_em, not t_obs: those differ by up to a_out/c ~ 900 s,
    over which beta itself changes by ~2.7e-4 -- numerically the same size
    as beta^2, so comparing against beta(t_obs) instead looks exactly like a
    spurious second-order disagreement and invites 'fixing' it with a gamma
    factor that does not belong here (doppler.doppler_factor's docstring).
    The message below reports both numbers so that trap stays visible."""
    kw = _orbit_kwargs()
    t_obs = np.linspace(0.0, kw["period_out"], 200001)
    t_em, n_iter = dp.emission_time(t_obs, **kw)
    d_num = np.gradient(t_em, t_obs)[5:-5]
    at_em = np.max(np.abs(d_num - dp.classical_doppler_factor(t_em, **_no_Omega(kw))[5:-5]))
    at_obs = np.max(np.abs(d_num - dp.classical_doppler_factor(t_obs, **_no_Omega(kw))[5:-5]))
    residual = np.max(np.abs(t_em + dp.roemer_delay(t_em, **kw) - t_obs))
    ok = at_em < 1e-8 and residual < 1e-9
    return ok, (f"max |dt_em/dt_obs - 1/(1+beta(t_em))| = {at_em:.2e} (tol 1e-8); "
                f"the same against beta(t_obs) would be {at_obs:.2e}. "
                f"Retarded-time residual |t_em + z/c - t_obs| = {residual:.2e} s "
                f"after {n_iter} iterations.")


def check_roemer_dwarfs_the_image_delay():
    """The finding this module exists for, recorded as a number rather than
    asserted in prose: at Case B's frequency the outer orbit's light-travel
    modulation is hundreds of times the image time delay Delta_T that the
    case is built to show. Also fixes the Roemer amplitude against its own
    closed form 2 a_src sin(i)/c, so the headline number in report/ has a
    second, independent derivation behind it."""
    kw = _orbit_kwargs()
    t = np.linspace(0.0, kw["period_out"], 200001)
    roemer_ptp = float(np.ptp(dp.roemer_delay(t, **kw)))
    a_src = dp.source_semimajor_axis(kw["a_out_m"], kw["M_lens_msun"], kw["M_binary_msun"])
    closed_form = 2.0 * a_src * np.sin(kw["i_out"]) / units.C_SI
    y_t, lensed, _ = geo.impact_parameter_of_time(
        t, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT,
        system.LITTLE_OMEGA_OUT, system.T_PERI_OUT, system.M_LENS_MSUN,
        system.D_L_PC)
    y_min = float(np.min(y_t[lensed]))
    dT = float(wo.time_delay_difference(y_min)) * 4.0 * units.msun_to_seconds(system.M_LENS_MSUN)
    ok = abs(roemer_ptp - closed_form) / closed_form < 1e-7 and roemer_ptp / dT > 100.0
    return ok, (f"Roemer p-p = {roemer_ptp:.1f} s (closed form 2 a_src sin i/c = {closed_form:.1f} s) "
                f"vs image delay Delta_T(y_min={y_min:.3f}) = {dT:.3f} s -> {roemer_ptp/dT:.0f}x; "
                f"{system.F_B_HZ * roemer_ptp:.1f} GW cycles p-p at f_B={system.F_B_HZ} Hz")


def check_envelope_unaffected_by_roemer():
    """The claim every unchanged Case B amplitude result rests on: the Roemer
    substitution is pure TIMING, so it must leave the inspiral strain
    ENVELOPE alone (the envelope evolves on the 240-day time-to-merger scale;
    the delay only reshuffles it by ~900 s). If this failed, Case B's |F|^2
    pulses, rms amplification and every diffraction figure would need
    regenerating rather than being demonstrably untouched."""
    kw = _orbit_kwargs()
    t = np.linspace(0.0, system.T_OBS_B_S, 20001)
    t_em, _ = dp.emission_time(t, **kw)
    t_c = chirp.time_to_merger(system.F_B_HZ, system.MCHIRP_MSUN)
    d_eff = system.D_L_PC / 1.0e6
    A_no = chirp.restricted_pn_amplitude_td(
        chirp.freq_of_time(t, t_c, system.MCHIRP_MSUN), system.MCHIRP_MSUN, d_eff)
    A_yes = chirp.restricted_pn_amplitude_td(
        chirp.freq_of_time(t_em, t_c, system.MCHIRP_MSUN), system.MCHIRP_MSUN, d_eff)
    rel = float(np.max(np.abs(A_yes - A_no) / A_no))
    ok = rel < 1e-3
    return ok, f"max fractional change in the strain envelope from the Roemer delay = {rel:.2e} (tol 1e-3)"


def check_case_A_roemer_is_negligible():
    """Case A's entire justification for treating the source as kinematically
    static, computed rather than asserted. Over the 15.2 s chirp a constant
    delay and a constant slope are exactly degenerate with t_c and with the
    chirp mass, so only the CURVATURE (line-of-sight acceleration) is
    observable. Evaluated in the window Case A actually uses -- ending at the
    orbital phase of closest approach, t/P=0.25, where run.py places the
    merger and where the line-of-sight velocity happens to cross zero.
    Also reports the constant orbital redshift, which is not a waveform
    feature at all but a bias on any inferred chirp mass."""
    kw = _orbit_kwargs()
    tau = chirp.time_to_merger(system.F_A_START_HZ, system.MCHIRP_MSUN)
    t = np.linspace(0.25 * kw["period_out"] - tau, 0.25 * kw["period_out"], 4001)
    residual, slope, _ = dp.roemer_residual(t, **kw)
    phase_rad = 2.0 * np.pi * system.F_ISCO_HZ * float(np.ptp(residual))
    z_orb = dp.orbital_redshift_factor(system.A_OUT_M, system.M_LENS_MSUN) - 1.0
    ok = phase_rad < 0.05
    return ok, (f"over the {tau:.1f} s chirp: degenerate slope beta_los={slope:.2e} "
                f"(absorbed into M_chirp), non-degenerate residual {np.ptp(residual):.2e} s "
                f"-> {phase_rad:.4f} rad at f_isco={system.F_ISCO_HZ:.1f} Hz (require < 0.05 rad). "
                f"Constant orbital redshift {z_orb:.3e}: an M_chirp bias, not a waveform feature.")


def check_orbital_redshift_is_the_schwarzschild_factor():
    """`orbital_redshift_factor` returns the EXACT Schwarzschild dtau/dt of a
    circular geodesic, 1/sqrt(1 - 3GM/(r c^2)). Rebuilt here from its two
    pieces instead of from that formula: the static metric term 2GM/(r c^2),
    and the transverse-Doppler term (r Omega/c)^2, which for the exact
    Schwarzschild coordinate angular velocity Omega^2 = GM/r^3 is a further
    GM/(r c^2). Their sum is the 3 in the factor, and the leading term of the
    result must then be 3GM/(2 r c^2) -- the split the module docstring states
    and that cases/ report as a chirp-mass bias. Worth its own check because
    the redshift is constant: it produces no waveform feature, so no other
    check in this suite would notice a wrong coefficient here."""
    r = system.A_OUT_M
    r_g = units.msun_to_meters(system.M_LENS_MSUN)      # GM_L/c^2, in metres
    Omega = np.sqrt(r_g * units.C_SI ** 2 / r ** 3)     # Omega^2 = GM/r^3
    dtau_dt = np.sqrt(1.0 - 2.0 * r_g / r - (r * Omega / units.C_SI) ** 2)
    exact = dp.orbital_redshift_factor(r, system.M_LENS_MSUN)
    rebuilt_err = abs(exact - 1.0 / dtau_dt) * dtau_dt

    lead = 1.5 * r_g / r        # grav GM/(r c^2) + transverse GM/(2 r c^2)
    lead_err = abs((exact - 1.0) - lead) / lead

    ok = rebuilt_err < 1e-13 and lead_err < 1e-3
    return ok, (f"1+z_orb = {exact:.12f}; rebuilt from 2GM/rc^2 and (r Omega/c)^2 "
                f"it agrees to {rebuilt_err:.1e}; its leading term "
                f"3GM/(2 r c^2) = {lead:.6e} reproduces z_orb = {exact - 1.0:.6e} "
                f"to {lead_err:.1e} relative, as 2GM/rc^2 = {2 * r_g / r:.3e}")


CHECKS = [
    ("los_velocity_matches_the_orbit_code", check_los_velocity_matches_the_orbit_code),
    ("retarded_time_is_the_whole_first_order_doppler", check_retarded_time_is_the_whole_first_order_doppler),
    ("roemer_dwarfs_the_image_delay", check_roemer_dwarfs_the_image_delay),
    ("envelope_unaffected_by_roemer", check_envelope_unaffected_by_roemer),
    ("case_A_roemer_is_negligible", check_case_A_roemer_is_negligible),
    ("orbital_redshift_is_the_schwarzschild_factor",
     check_orbital_redshift_is_the_schwarzschild_factor),
]


def main():
    results = {}
    all_ok = True
    for name, fn in CHECKS:
        try:
            ok, msg = fn()
        except Exception as exc:  # noqa: BLE001 - report, don't hide
            ok, msg = False, f"EXCEPTION: {exc!r}"
        all_ok &= ok
        results[name] = {"passed": bool(ok), "message": msg}
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")
    out = Path(__file__).parent / "CHECKS_doppler.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
