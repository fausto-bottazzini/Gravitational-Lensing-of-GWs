"""TaylorF2: the standard frequency-domain, stationary-phase-approximation
(SPA) post-Newtonian inspiral waveform -- the actual functional form used
for the inspiral part of LIGO/Virgo low-mass compact-binary template banks
(e.g. the PyCBC/LALSimulation "TaylorF2" approximant). Added on request
(the leading-order-only, hand-built time-domain chirp in `chirp.py` is not
what a real search would template against) -- see wiki/log.md.

Truncated at 2PN (v^4) here, not the full 3.5PN bank-standard order:
the 2PN coefficients below are well-established and independently
cross-checked (see `check_leading_order_group_delay` in
tests/test_taylorf2.py); the 2.5PN-3.5PN coefficients involve enough extra
rational/log terms that reproducing them from memory without a verified
reference implementation to check against was judged NOT worth the risk of
a silent numerical error, given this project's emphasis on not asserting a
number that has not been checked. Logged as a scope choice, not an
oversight -- see wiki/log.md and theory/theory.tex.

Reference for the phase coefficients: Buonanno, Iyer, Ochsner, Pan &
Sathyaprakash (2009), PRD 80, 084043, Eq. (3.10)-(3.13) (and consistent
with e.g. Cutler & Flanagan 1994; Blanchet, Living Rev. Relativity).
"""
import numpy as np

from . import units


def _v_of_f(f_hz, mtot_msun):
    M_sec = units.msun_to_seconds(mtot_msun)
    return (np.pi * M_sec * np.asarray(f_hz, dtype=float)) ** (1.0 / 3.0)


def spa_phase(f_hz, m1_msun, m2_msun, t_c, phi_c=0.0):
    """Psi(f), the TaylorF2 SPA phase to 2PN (v^4), in THIS project's
    Fourier convention (wiki/conventions.md: h(f)=int h(t) e^{-2pi i f t}
    dt, h(t)=int h(f) e^{+2pi i f t} df -- same sign as numpy's ifft):

        Psi(f) = -[ 2*pi*f*t_c - phi_c - pi/4
                    + (3/(128*eta*v^5)) * (1 + phi2*v^2 + phi3*v^3 + phi4*v^4) ]

    v = (pi*M*f)^(1/3), M = m1+m2 (geometrized seconds), eta = m1*m2/M^2.

    The overall minus sign is NOT the sign printed in most PN-literature
    statements of this formula (which assume the opposite FT convention,
    h(f)=int h(t) e^{+2pi i f t} dt). Using the literature sign directly
    with numpy's ifft/irfft put the reconstructed chirp's peak at
    t~17-245 s depending on zero-padding (it should sit at t_c=15.2 s,
    independent of padding) -- caught exactly because it MOVED with the
    padding length, the signature of a sign error aliasing power to the
    wrong end of the (periodic) IFFT window, not of a real physical
    effect. Fixed by flipping the sign and confirming the reconstructed
    peak lands at t_c to <0.2% and stays there as padding is increased;
    see wiki/log.md and tests/test_taylorf2.py::
    check_ifft_reconstructs_chirp_at_tc.
    """
    mtot = m1_msun + m2_msun
    eta = m1_msun * m2_msun / mtot ** 2
    v = _v_of_f(f_hz, mtot)

    phi2 = 3715.0 / 756.0 + 55.0 * eta / 9.0                              # 1PN
    phi3 = -16.0 * np.pi                                                   # 1.5PN (tail)
    phi4 = 15293365.0 / 508032.0 + 27145.0 * eta / 504.0 + 3085.0 * eta ** 2 / 72.0  # 2PN

    pn_series = 1.0 + phi2 * v ** 2 + phi3 * v ** 3 + phi4 * v ** 4
    psi_literature_sign = (2.0 * np.pi * f_hz * t_c - phi_c - np.pi / 4.0
                            + (3.0 / (128.0 * eta * v ** 5)) * pn_series)
    return -psi_literature_sign


def spa_amplitude(f_hz, m1_msun, m2_msun, d_eff_mpc):
    """|h(f)| ~ f^(-7/6) SPA amplitude -- the same formula as
    chirp.restricted_pn_amplitude_fd, re-derived here in terms of (m1,m2)
    directly rather than the chirp mass, for a self-contained module; the
    two are cross-checked for consistency in tests/test_taylorf2.py."""
    mchirp = (m1_msun * m2_msun) ** 0.6 / (m1_msun + m2_msun) ** 0.2
    Mc_sec = units.msun_to_seconds(mchirp)
    D_sec = d_eff_mpc * units.MPC_SI / units.C_SI
    A0 = np.sqrt(5.0 / 24.0) / np.pi ** (2.0 / 3.0) * Mc_sec ** (5.0 / 6.0) / D_sec
    return A0 * np.asarray(f_hz, dtype=float) ** (-7.0 / 6.0)


def htilde(f_hz, m1_msun, m2_msun, t_c, d_eff_mpc, phi_c=0.0):
    """The restricted TaylorF2 frequency-domain waveform,
    h(f) = A(f) * exp(i*Psi(f)), f>0 only (apply Hermitian symmetry for
    the negative-frequency half when building a real time series)."""
    amp = spa_amplitude(f_hz, m1_msun, m2_msun, d_eff_mpc)
    psi = spa_phase(f_hz, m1_msun, m2_msun, t_c, phi_c)
    return amp * np.exp(1j * psi)
