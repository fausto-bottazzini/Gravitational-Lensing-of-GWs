"""Leading-order (Newtonian/quadrupole, restricted) inspiral waveform for the
inner binary -- the GW source in this project.

Deliberately NOT full PN (see wiki/conventions.md: the lensing physics is the
point, not waveform accuracy). Every formula here is derived from one input
that is textbook-safe and stated explicitly: the quadrupole-formula frequency
evolution

    df/dt = (96/5) * pi^{8/3} * (G*Mchirp/c^3)^{5/3} * f^{11/3}      (*)

(e.g. Peters 1964; Maggiore, "Gravitational Waves" Vol.1 Eq. 4.??-class
result). Time-to-coalescence and f(t) below are obtained by directly
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


def _k_of_mchirp(mchirp_msun):
    Mc_sec = units.msun_to_seconds(mchirp_msun)
    return (96.0 / 5.0) * np.pi ** (8.0 / 3.0) * Mc_sec ** (5.0 / 3.0)


def time_to_merger(f_hz, mchirp_msun):
    """tau(f) = t_c - t, from integrating df/dt = k f^{11/3}:
    tau(f) = (5/256) * (pi f)^{-8/3} * (G*Mchirp/c^3)^{-5/3}."""
    Mc_sec = units.msun_to_seconds(mchirp_msun)
    return (5.0 / 256.0) * (np.pi * f_hz) ** (-8.0 / 3.0) * Mc_sec ** (-5.0 / 3.0)


def freq_of_time(t, t_c, mchirp_msun):
    """f(t) = (1/pi) * [(5/256) * (Mc_sec)^{-5/3} * (t_c-t)]^{-3/8} for t<t_c,
    the closed-form solution of (*) -- algebraically re-derived from
    `time_to_merger` by inverting tau(f) for f (see module docstring); cross
    -checked against direct numerical ODE integration in tests/test_chirp.py.
    """
    Mc_sec = units.msun_to_seconds(mchirp_msun)
    tau = np.clip(t_c - np.asarray(t, dtype=float), 1e-12, None)
    # f = (1/pi) * (5/256)^(3/8) * tau^(-3/8) * Mc_sec^(-5/8) -- solved for f
    # from time_to_merger(f)=tau by direct algebraic inversion (see the
    # module docstring derivation). An earlier version got the exponents on
    # (5/256) and Mc_sec backwards by folding them into one bracket; caught
    # because f(t=0) did not reproduce the input f0 in tests/test_chirp.py.
    return (1.0 / np.pi) * (5.0 / 256.0) ** (3.0 / 8.0) * tau ** (-3.0 / 8.0) * Mc_sec ** (-5.0 / 8.0)


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


def phase_of_time(t, f_func, t0, t1, n=200_000):
    """Phi(t) = 2*pi*int_{t0}^t f(t') dt', by cumulative Simpson/trapezoid
    over a fine grid (numerical -- avoids trusting a closed-form phase
    formula that this project does not need elsewhere)."""
    tt = np.linspace(t0, t1, n)
    ff = f_func(tt)
    phase = 2.0 * np.pi * np.concatenate(([0.0], np.cumsum(
        0.5 * (ff[1:] + ff[:-1]) * np.diff(tt))))
    return tt, phase


def restricted_pn_amplitude(f_hz, mchirp_msun, d_eff_mpc):
    """Restricted (leading-order, quadrupole) PN amplitude prefactor, the
    standard textbook stationary-phase-approximation scaling
    |h(f)| ~ A0 f^{-7/6}, A0 = sqrt(5/24)/pi^{2/3} * (G Mchirp/c^3)^{5/6} *
    c / (D_eff), so that the plotted time- and frequency-domain amplitudes
    have the right *relative* scaling; the absolute normalization is not the
    point of this project (no calibration against a detector's absolute
    strain units is claimed) and is logged as a choice in wiki/log.md.
    """
    Mc_sec = units.msun_to_seconds(mchirp_msun)
    D_sec = d_eff_mpc * units.MPC_SI / units.C_SI  # Mpc -> light-seconds
    A0 = np.sqrt(5.0 / 24.0) / np.pi ** (2.0 / 3.0) * Mc_sec ** (5.0 / 6.0) / D_sec
    return A0 * np.asarray(f_hz, dtype=float) ** (-7.0 / 6.0)
