"""Leading-order (Newtonian/quadrupole, restricted) inspiral waveform for the
inner binary -- the GW source in this project.

Deliberately NOT full PN (see wiki/conventions.md: the lensing physics is the
point, not waveform accuracy). Every formula here is derived from one input
that is textbook-safe and stated explicitly: the quadrupole-formula frequency
evolution

    df/dt = (96/5) * pi^{8/3} * (G*Mchirp/c^3)^{5/3} * f^{11/3}      (*)

(Peters 1964, Eq. (5.6)'s da/dt rewritten in terms of the GW frequency via
Kepler's third law; also standard textbook material, e.g. Maggiore 2008
Vol. 1 Ch. 4, cited by chapter rather than by equation number, which this
project has no copy of to check).
Time-to-coalescence and f(t) below are obtained by directly
integrating (*), not copied from a table -- reproduced in
tests/test_chirp.py by numerically integrating (*) with
scipy.integrate.solve_ivp and checking it agrees with the closed form.
"""
import numpy as np
from scipy.integrate import solve_ivp

from . import units


def chirp_mass(m1, m2):
    return (m1 * m2) ** (3.0 / 5.0) / (m1 + m2) ** (1.0 / 5.0)


def isco_frequency(m_total_msun):
    """GW frequency at the Schwarzschild ISCO (r=6GM/c^2), leading order:
    f_isco = c^3 / (6^{3/2} * pi * G * M_total) (GW frequency = 2x orbital)."""
    M_sec = units.msun_to_seconds(m_total_msun)
    return 1.0 / (6.0 ** 1.5 * np.pi * M_sec)


def time_to_merger(f_hz, mchirp_msun):
    """tau(f) = t_c - t, from integrating df/dt = k f^{11/3}:
    tau(f) = (5/256) * (pi f)^{-8/3} * (G*Mchirp/c^3)^{-5/3}."""
    Mc_sec = units.msun_to_seconds(mchirp_msun)
    return (5.0 / 256.0) * (np.pi * f_hz) ** (-8.0 / 3.0) * Mc_sec ** (-5.0 / 3.0)


def freq_of_time(t, t_c, mchirp_msun):
    """f(t) = (1/pi) * (5/256)^{3/8} * (t_c-t)^{-3/8} * Mc_sec^{-5/8} for t<t_c,
    the closed-form solution of (*) -- algebraically re-derived from
    `time_to_merger` by inverting tau(f) for f (see module docstring); cross
    -checked against direct numerical ODE integration in tests/test_chirp.py.
    """
    Mc_sec = units.msun_to_seconds(mchirp_msun)
    tau = np.clip(t_c - np.asarray(t, dtype=float), 1e-12, None)
    # f = (1/pi) * (5/256)^(3/8) * tau^(-3/8) * Mc_sec^(-5/8) -- solved for f
    # from time_to_merger(f)=tau by direct algebraic inversion (see the
    # module docstring derivation). Written as separate factors on purpose:
    # folded into one bracket, the exponents on (5/256) and Mc_sec are easy
    # to get wrong, and tests/test_chirp.py exists partly to catch that.
    return (1.0 / np.pi) * (5.0 / 256.0) ** (3.0 / 8.0) * tau ** (-3.0 / 8.0) * Mc_sec ** (-5.0 / 8.0)


def _k_of_mchirp(mchirp_msun):
    Mc_sec = units.msun_to_seconds(mchirp_msun)
    return (96.0 / 5.0) * np.pi ** (8.0 / 3.0) * Mc_sec ** (5.0 / 3.0)


def freq_of_time_numeric(t_array, f0, mchirp_msun):
    """Same f(t), but by direct numerical integration of df/dt=(*) with
    solve_ivp -- an independent check on `freq_of_time`, not a replacement
    for it (used only in tests/, not in the case scripts)."""
    k = _k_of_mchirp(mchirp_msun)

    def rhs(t, f):
        return [k * f[0] ** (11.0 / 3.0)]

    sol = solve_ivp(rhs, (t_array[0], t_array[-1]), [f0], t_eval=t_array,
                     rtol=1e-10, atol=1e-14, method="RK45")
    return sol.y[0]


def phase_of_time(f_func, t0, t1, n=200_000):
    """Phi(t) = 2*pi*int_{t0}^t f(t') dt', by cumulative trapezoid
    over a fine grid (numerical -- avoids trusting a closed-form phase
    formula that this project does not need elsewhere)."""
    tt = np.linspace(t0, t1, n)
    ff = f_func(tt)
    phase = 2.0 * np.pi * np.concatenate(([0.0], np.cumsum(
        0.5 * (ff[1:] + ff[:-1]) * np.diff(tt))))
    return tt, phase


def restricted_pn_amplitude_td(f_hz, mchirp_msun, d_eff_mpc):
    """Restricted (leading-order, quadrupole) TIME-DOMAIN strain envelope:

        h(t) = (4/D_eff) * (G*Mchirp/c^2)^(5/3) * (pi*f_gw(t)/c)^(2/3)

    (e.g. Maggiore, "Gravitational Waves" Vol. 1, restricted-PN (2,2) mode;
    O(1) inclination/polarization prefactors folded into the leading 4, not
    separately calibrated -- logged as a choice in wiki/log.md, along with
    everywhere else this project does not claim an absolutely calibrated
    strain). It GROWS with f, as a real inspiral's amplitude must. Do not
    substitute `taylorf2.spa_amplitude` here: that is the frequency-domain
    f^-7/6 scaling, which decreases with f, and using it as a time-domain
    envelope makes the plotted chirp shrink towards merger.
    """
    Mc_sec = units.msun_to_seconds(mchirp_msun)
    R_c_m = Mc_sec * units.C_SI  # G*Mchirp/c^2, metres
    D_eff_m = d_eff_mpc * units.MPC_SI
    return (4.0 / D_eff_m) * R_c_m ** (5.0 / 3.0) * (np.pi * np.asarray(f_hz, dtype=float) / units.C_SI) ** (2.0 / 3.0)
