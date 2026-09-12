"""Line-of-sight kinematics of the OUTER orbit: Roemer delay, first-order
Doppler, and the constant orbital redshift.

The source is not at rest: it orbits the lens at v/c ~ 1.6e-2 with a
light-crossing time a_out/c ~ 900 s, so the same outer orbit that produces
the lensing also imprints its own kinematics on the observed waveform. At
Case B's f_B = 0.05 Hz the Roemer modulation is 1810 s peak to peak, i.e.
90.5 GW cycles of phase, against the 0.017 cycles the lens itself writes:
this is the larger of the two things the orbit does, by orders of
magnitude, and it is pure timing, so it leaves every amplitude result
alone. It is kinematics of the source rather than optics of the lens, so
the static-lens formalism has nothing to say about it; D'Orazio & Loeb
(2020) model the two together in their Appendix A. Added 2026-09-07 after
being missing entirely -- see wiki/log.md.

Three physically distinct pieces, kept separate because they are observable
in completely different ways:

  1. `roemer_delay`      the classical light-travel-time term z_src(t)/c.
                         This is the whole first-order effect: to O(beta) the
                         observed waveform is simply the emitted one
                         evaluated at the retarded emission time
                         (`emission_time` below). It is pure TIMING -- the
                         strain ENVELOPE is untouched, because the envelope
                         evolves on the inspiral timescale (240 d here) and
                         the delay only reshuffles it by ~900 s (relative
                         change ~6e-5, checked in tests/test_doppler.py).
                         That is why adding it changes none of this repo's
                         amplitude results.
  2. `classical_doppler_factor`
                         1/(1 + beta_los), the FIRST-ORDER frequency shift.
                         This one is not applied separately anywhere either:
                         differentiating the retarded-time relation gives
                         d t_em/d t_obs = 1/(1 + beta_los) identically, so
                         evaluating the emitted waveform at `emission_time`
                         already produces exactly this shift. That identity
                         is checked to machine precision in
                         tests/test_doppler.py, and it is what makes the
                         retarded-time substitution the *complete*
                         first-order treatment rather than one term of it.
                         (It is specifically NOT the special-relativistic
                         D = [gamma(1+beta)]^-1; see the note under 3 for
                         why the retarded-time construction does not
                         reproduce that one and must not be expected to.)
  3. `orbital_redshift_factor`
                         everything beyond first order, which for a CIRCULAR
                         orbit is exactly constant: the gravitational
                         redshift in the lens's potential plus the transverse
                         (second-order, 1/gamma) Doppler shift. For a
                         circular geodesic these combine into the exact
                         Schwarzschild result
                         dtau/dt = sqrt(1 - 3 G M_L/(r c^2)), i.e.
                         1+z_orb = (1 - 3GM_L/(r c^2))^{-1/2}.
                         So the complete frequency relation splits cleanly,
                             f_obs = f_em,proper / [ (1+beta_los)(1+z_orb) ],
                         with the first bracket time-varying (and carried by
                         `emission_time`) and the second constant. Constant
                         means EXACTLY degenerate with a rescaling of the
                         source's chirp mass and of t_c, so it is not
                         separately observable in a single event -- it is a
                         bias on the inferred M_chirp (4.1e-4 here), not a
                         waveform feature, and is reported rather than
                         applied. This split is also why the 1/gamma of the
                         special-relativistic Doppler factor must NOT be
                         applied on top of `emission_time`: it is already
                         inside 1+z_orb, and applying both double-counts it.

Degeneracy, and why Case A and Case B use this module differently:
over a SHORT observation (Case A: 15 s out of a 4-day orbit) a constant
offset and a constant slope in the delay are themselves degenerate -- a
constant delay is reabsorbed into t_c, a constant slope into the observed
chirp mass. The only non-degenerate, genuinely observable part there is the
CURVATURE, i.e. the line-of-sight acceleration; `roemer_residual` returns
exactly that (the delay with its best-fit constant + linear part removed
over the observation window). Over Case B's 6 full outer periods, by
contrast, the delay is periodic with zero mean and zero trend, nothing is
degenerate, and the full `roemer_delay` applies.
"""
import numpy as np

from . import units, geometry


def source_semimajor_axis(a_out_m, M_lens_msun, M_binary_msun):
    """The SOURCE's own orbit about the triple's barycentre, from the
    relative orbit: a_src = a_out * M_L/(M_L+M_b).

    `geometry.relative_separation` returns the RELATIVE separation (source
    minus lens); the Roemer delay is set by the source's motion about the
    barycentre, which is that relative vector scaled by the lens's mass
    fraction. Here M_L/(M_L+M_b)=0.99930, so this is a 7e-4 correction --
    kept explicit anyway rather than silently approximated by 1, since
    nothing else in this module needs it to be small.
    """
    return a_out_m * M_lens_msun / (M_lens_msun + M_binary_msun)


def source_los_offset(t, a_out_m, e_out, period_out, i_out, Omega_out,
                      omega_out, t_peri, M_lens_msun, M_binary_msun):
    """z_src(t): the source's line-of-sight displacement from the barycentre,
    in metres, positive AWAY from the observer -- the same axis sense as
    `geometry.orbital_plane_to_sky` (whose docstring fixes it, and whose
    z_los>0 means "source behind the lens" in `impact_parameter_of_time`).

    `a_out_m` is in METRES here, not the parsecs `geometry` takes: this is a
    light-travel time, so metres is the only unit that does not invite a
    silent 3e16 error. Callers pass `system.A_OUT_M` directly.
    """
    r, nu = geometry.relative_separation(t, a_out_m, e_out, period_out, t_peri)
    _, _, z_rel = geometry.orbital_plane_to_sky(r, nu, omega_out, i_out, Omega_out)
    return z_rel * M_lens_msun / (M_lens_msun + M_binary_msun)


def roemer_delay(t, a_out_m, e_out, period_out, i_out, Omega_out, omega_out,
                 t_peri, M_lens_msun, M_binary_msun):
    """Light-travel-time delay z_src(t)/c, in seconds, relative to the
    barycentre. Positive = source farther from the observer = signal arrives
    later."""
    return source_los_offset(t, a_out_m, e_out, period_out, i_out, Omega_out,
                             omega_out, t_peri, M_lens_msun,
                             M_binary_msun) / units.C_SI


def los_velocity(t, a_out_m, e_out, period_out, i_out, omega_out, t_peri,
                 M_lens_msun, M_binary_msun):
    """dz_src/dt in m/s, positive = receding, from the standard closed-form
    radial-velocity equation rather than by differentiating
    `source_los_offset` numerically:

        v_z(t) = K [ cos(nu+omega) + e cos(omega) ],
        K = 2 pi a_src sin(i) / (P sqrt(1-e^2))

    (the spectroscopic-binary radial-velocity curve; e.g. Murray & Dermott,
    "Solar System Dynamics", Sec. 2.5-2.8). Note `Omega` does not appear --
    correct, and a useful consistency check on
    `geometry.orbital_plane_to_sky`, whose z_los = r sin(nu+omega) sin(i) is
    likewise Omega-independent. Cross-checked against a finite-difference
    derivative of `source_los_offset` in tests/test_doppler.py, which is what
    actually ties this closed form to the orbit code the rest of the repo
    uses.
    """
    a_src = source_semimajor_axis(a_out_m, M_lens_msun, M_binary_msun)
    K = 2.0 * np.pi * a_src * np.sin(i_out) / (period_out * np.sqrt(1.0 - e_out ** 2))
    _, nu = geometry.relative_separation(t, a_out_m, e_out, period_out, t_peri)
    return K * (np.cos(nu + omega_out) + e_out * np.cos(omega_out))


def classical_doppler_factor(t, a_out_m, e_out, period_out, i_out, omega_out,
                             t_peri, M_lens_msun, M_binary_msun):
    """1/(1 + beta_los): the first-order (light-travel-time / Roemer)
    frequency shift, beta_los = v_z/c positive for a receding source.

    This is exactly d t_em / d t_obs for the retarded-time relation solved by
    `emission_time` -- differentiate t_obs = t_em + z_src(t_em)/c and the
    result is immediate -- so a caller that evaluates its waveform at
    `emission_time` has already applied this and must not apply it again.
    The identity is checked to machine precision in tests/test_doppler.py.

    Note what this is NOT: the special-relativistic Doppler factor
    D = [gamma(1+beta)]^-1. The retarded-time construction is carried out in
    the barycentric frame and so yields the classical factor; the extra
    1/gamma is the source clock's own time dilation, which for this circular
    orbit is constant and is already inside `orbital_redshift_factor`
    (together with the gravitational redshift, which the special-relativistic
    expression does not contain at all). Multiplying by D on top of
    `emission_time` would double-count the 1/gamma and still miss the
    gravitational piece. A `doppler_factor` returning D used to live here for
    "reporting"; nothing ever called it, so it is gone -- the warning it
    carried is what mattered and it is this paragraph.
    """
    beta = los_velocity(t, a_out_m, e_out, period_out, i_out, omega_out,
                        t_peri, M_lens_msun, M_binary_msun) / units.C_SI
    return 1.0 / (1.0 + beta)


def emission_time(t_obs, a_out_m, e_out, period_out, i_out, Omega_out,
                  omega_out, t_peri, M_lens_msun, M_binary_msun,
                  tol=1e-9, max_iter=50):
    """Solve t_em from t_obs = t_em + z_src(t_em)/c by fixed-point iteration.

    The observed waveform is the emitted one evaluated here:
    h_obs(t_obs) = h_em(t_em(t_obs)). That single substitution IS the whole
    first-order Doppler/Roemer effect -- both the arrival-time modulation and
    the frequency shift f_obs = f_em/(1+beta_los) fall out of it, with no
    separate factor applied anywhere, since differentiating the relation
    above gives d t_em/d t_obs = 1/(1+beta_los) =
    `classical_doppler_factor` identically (checked in
    tests/test_doppler.py). What it does NOT contain is the second-order and
    gravitational piece, which is constant here and lives in
    `orbital_redshift_factor`; see the module docstring for the split.

    The map contracts at rate |dz/dt|/c = beta ~ 1.6e-2 per iteration, so
    tol=1e-9 s is reached in ~5 iterations; `max_iter` is a guard, not a
    working limit, and the iteration count is returned so a caller can assert
    it converged rather than assume it.

    Returns (t_em, n_iter).
    """
    t_obs = np.atleast_1d(np.asarray(t_obs, dtype=float))
    t_em = t_obs.copy()
    for n in range(1, max_iter + 1):
        delay = roemer_delay(t_em, a_out_m, e_out, period_out, i_out,
                             Omega_out, omega_out, t_peri, M_lens_msun,
                             M_binary_msun)
        t_new = t_obs - delay
        if np.max(np.abs(t_new - t_em)) < tol:
            return t_new, n
        t_em = t_new
    raise RuntimeError(f"emission_time did not converge in {max_iter} iterations")


def roemer_residual(t, a_out_m, e_out, period_out, i_out, Omega_out,
                    omega_out, t_peri, M_lens_msun, M_binary_msun):
    """The Roemer delay with its best-fit constant and linear parts removed
    over the sampled window -- the only part observable at all in a SHORT
    observation (Case A).

    A constant delay is exactly reabsorbed into the coalescence time t_c; a
    constant slope (a constant line-of-sight velocity) is exactly reabsorbed
    into the observed chirp mass, being an unmeasurable constant redshift.
    What survives is the curvature, i.e. the line-of-sight ACCELERATION, and
    that is genuinely non-degenerate: it makes the binary appear to chirp at
    a slightly wrong rate. Returns (residual_seconds, slope_s_per_s,
    offset_seconds) so a caller can report the degenerate pieces too rather
    than silently discarding them.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    d = roemer_delay(t, a_out_m, e_out, period_out, i_out, Omega_out,
                     omega_out, t_peri, M_lens_msun, M_binary_msun)
    slope, offset = np.polyfit(t, d, 1)
    return d - (slope * t + offset), float(slope), float(offset)


def orbital_redshift_factor(a_out_m, M_lens_msun):
    """1+z_orb = (1 - 3 G M_L / (r c^2))^{-1/2} for a circular orbit of radius
    r=a_out about the lens: the EXACT Schwarzschild combination of the
    gravitational redshift in the lens's potential and the transverse
    (second-order) Doppler shift, for a clock on a circular geodesic as read
    by a distant static observer (dtau/dt = sqrt(1-3GM/(r c^2))).

    Being constant for a circular orbit, this is exactly degenerate with a
    rescaling of the source's chirp mass -- it biases an inferred M_chirp by
    this factor and produces no waveform feature at all, which is why nothing
    in cases/ applies it. Reported, not applied. (Its two pieces separately:
    gravitational GM/(r c^2), transverse-Doppler GM/(2 r c^2), since
    v^2 = GM/r on a circular orbit -- their sum 3GM/(2 r c^2) is the leading
    term of the expression above, checked in tests/test_doppler.py.)
    """
    r_g = units.msun_to_meters(M_lens_msun)  # GM_L/c^2, metres
    return 1.0 / np.sqrt(1.0 - 3.0 * r_g / a_out_m)
