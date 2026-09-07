"""Wave-optics gravitational lensing of a scalar field, point-mass lens.

Notation follows wiki/conventions.md (dimensionless lens-plane position x,
dimensionless source position/impact parameter y, dimensionless frequency
w = 8 pi G M_L f / c^3). Full derivation: theory/theory.tex.

Everything here is evaluated at LEAST TWO independent ways and cross-checked
in tests/test_waveoptics.py — see wiki/todo.md. Three evaluators of F(w,y):

  1. `F_bruteforce_2d`   direct 2D numerical quadrature of the Kirchhoff
                         diffraction integral, regularized by a convergence
                         factor. Slow, first-principles, no closed-form input
                         at all -- the ground truth.
  2. `F_radial_1d`       the same integral after doing the azimuthal part
                         analytically (Bessel J0 identity), leaving a 1D
                         radial integral. Faster, still first-principles.
  3. `F_point_lens`      the closed-form confluent-hypergeometric expression
                         (Deguchi & Watson 1986; Takahashi & Nakamura 2003),
                         evaluated at arbitrary precision with mpmath. Fast;
                         this is what the case scripts actually call.

Geometric-optics (ray) quantities -- image positions, magnifications, time
delay -- are derived independently below and used both as w -> infinity
checks on F and, via the Paczynski (1986) total-magnification formula, as an
internal cross-check on the magnification formula itself.
"""
import numpy as np
from scipy import integrate
import mpmath as mp


# ---------------------------------------------------------------------
# Geometric optics (point-mass / Schwarzschild lens, y along the x-axis
# by symmetry; x here is a SIGNED scalar along that axis unless noted)
# ---------------------------------------------------------------------

def image_positions(y):
    """Signed 1D positions of the two images, y = x - 1/x solved for x.

    x_+ > 1 > 0 > x_- > -1 (x_+ x_- = -1 always; |x_-| < 1 < x_+ for y > 0).
    """
    y = np.asarray(y, dtype=float)
    disc = np.sqrt(y**2 + 4.0)
    x_plus = 0.5 * (y + disc)
    x_minus = 0.5 * (y - disc)
    return x_plus, x_minus


def magnification(x):
    """Point-lens magnification of a single image at signed position x.

    mu(x) = [1-(1/x)^4]^-1 = x^4/(x^4-1), from the Jacobian of the radial
    lens map y(x) = x - 1/x: det(dy/dx) = (y(x)/x) * y'(x) for an
    axisymmetric map, giving (1-1/x^2)(1+1/x^2) = 1-1/x^4.
    """
    x = np.asarray(x, dtype=float)
    return x**4 / (x**4 - 1.0)


def total_magnification_paczynski(y):
    """A(y) = (y^2+2) / (y*sqrt(y^2+4)) -- Paczynski (1986) microlensing
    total (unsigned) magnification. Independent closed form used ONLY to
    cross-check `magnification(image_positions(y))` below."""
    y = np.asarray(y, dtype=float)
    return (y**2 + 2.0) / (y * np.sqrt(y**2 + 4.0))


def fermat_potential(x, y):
    """phi(x,y) = (1/2)(x-y)^2 - ln|x|, x,y signed scalars (source on-axis)."""
    return 0.5 * (x - y)**2 - np.log(np.abs(x))


def time_delay_difference(y):
    """Delta T(y) = phi(x_-,y) - phi(x_+,y) > 0 (x_- is the delayed image)."""
    x_plus, x_minus = image_positions(y)
    dT = fermat_potential(x_minus, y) - fermat_potential(x_plus, y)
    return dT


def F_geometric_optics(w, y):
    """Leading (infinite-w) geometric-optics amplification factor.

    F -> exp(i*w*phi_+) * [sqrt(mu_+) - i*sqrt(|mu_-|)*exp(i*w*DeltaT)], the
    standard two-image stationary-phase (Morse-index) result: minimum image
    x_+ carries phase w*phi_+(y) (its own, UN-subtracted Fermat potential --
    `F_point_lens`/`F_radial_1d`/`F_bruteforce_2d` below do not subtract a
    reference phase either, so this term must be kept for the comparison to
    close), saddle image x_- carries an extra Morse phase -pi/2 (factor -i)
    relative to it, in the e^{-i2*pi*f*t} convention fixed in
    wiki/conventions.md. The exp(i*w*phi_+) piece was missing in an earlier
    version of this function; restored after `F_point_lens` disagreed with
    it in phase (not amplitude) at large w -- see wiki/log.md, and
    tests/test_waveoptics.py::test_geometric_optics_limit for the numbers
    (relative error ~3e-4 at w=10, ~3e-5 at w=1000, shrinking as w grows).
    """
    w = np.asarray(w, dtype=float)
    x_plus, x_minus = image_positions(y)
    mu_plus = magnification(x_plus)
    mu_minus = magnification(x_minus)
    dT = time_delay_difference(y)
    phi_plus = fermat_potential(x_plus, y)
    envelope = np.sqrt(mu_plus) - 1j * np.sqrt(np.abs(mu_minus)) * np.exp(1j * w * dT)
    return np.exp(1j * w * phi_plus) * envelope


# ---------------------------------------------------------------------
# Wave optics: the diffraction integral, three ways
# ---------------------------------------------------------------------

def _phi2d(x1, x2, y):
    """2D Fermat potential, source at (y,0): (1/2)|x-y|^2 - ln|x|."""
    r2 = x1**2 + x2**2
    return 0.5 * ((x1 - y)**2 + x2**2) - 0.5 * np.log(r2)


def F_bruteforce_2d(w, y, x_max=40.0, n=1200, eta=0.0):
    """Direct 2D quadrature of F(w,y) = (w/2pi i) int d^2x exp[i w phi(x,y)].

    Regularized with an optional Gaussian convergence factor exp(-eta*x^2)
    (eta=0 -> unregularized; increase eta slightly if the raw sum does not
    look converged for a given (w, x_max, n) -- checked in tests/, not just
    asserted here). O(n^2) points on a square grid; slow, spot-check only.
    """
    xs = np.linspace(-x_max, x_max, n)
    dx = xs[1] - xs[0]
    x1, x2 = np.meshgrid(xs, xs, indexing="ij")
    r2 = x1**2 + x2**2
    mask = r2 > 1e-12  # exclude the point-lens singularity itself (measure zero)
    phase = w * _phi2d(x1, x2, y)
    integrand = np.where(mask, np.exp(1j * phase - eta * r2), 0.0)
    integral = np.sum(integrand) * dx * dx
    return (w / (2j * np.pi)) * integral


def F_radial_1d(w, y, eta=0.05, safety=25.0, min_breakpoints=60, max_breakpoints=3000):
    """1D radial reduction of the same integral via
    int_0^2pi dtheta exp(i a cos theta) = 2 pi J0(a) (Abramowitz & Stegun
    9.1.21), applied to the cross term in |x-y|^2 = x^2+y^2-2xy*cos(theta):

        F(w,y) = -i w exp(i w y^2/2) int_0^inf dx x J0(w x y)
                 exp[i w (x^2/2 - ln x)]

    The integral is only conditionally convergent; regularized with a small
    positive `eta` added to the exponent's imaginary part (x^2 -> x^2*(1-i
    eta)), i.e. w -> w(1-i*eta), the standard i-epsilon prescription
    (Ulmer & Goodman 1995): the integrand acquires a factor
    exp(-eta*w*x^2/2), which is negligible once eta*w*x^2/2 > `safety`.
    That fixes an x_max = sqrt(2*safety/(eta*w)) automatically, and the
    breakpoint count is set from the number of oscillations
    (~w*x_max^2/(4pi)) expected below it, so this stays accurate without a
    hand-tuned cutoff. Convergence as eta -> 0 (at fixed accuracy) is
    checked in tests/test_waveoptics.py against `F_point_lens`.
    """
    from scipy.special import j0

    x_max = np.sqrt(2.0 * safety / (eta * w))
    n_osc = w * x_max**2 / (4.0 * np.pi)
    n_periods = int(np.clip(4 * n_osc, min_breakpoints, max_breakpoints))

    def integrand(x):
        if x <= 0:
            return 0.0 + 0.0j
        # exp(i * w*(1+i*eta) * x^2/2) = exp(i w x^2/2) * exp(-eta*w*x^2/2):
        # the (1+i*eta) factor -- NOT (1-i*eta) -- is what damps at large x;
        # the sign was wrong in an earlier version and blew the integral up
        # (see wiki/log.md).
        decay = np.exp(-eta * w * x**2 / 2.0)
        oscillation = np.exp(1j * w * (0.5 * x**2 - np.log(x)))
        return x * j0(w * x * y) * oscillation * decay

    breakpoints = np.linspace(0, x_max, n_periods)
    real_part = 0.0
    imag_part = 0.0
    for a, b in zip(breakpoints[:-1], breakpoints[1:]):
        re, _ = integrate.quad(lambda x: integrand(x).real, a, b, limit=200)
        im, _ = integrate.quad(lambda x: integrand(x).imag, a, b, limit=200)
        real_part += re
        imag_part += im
    integral = real_part + 1j * imag_part
    return -1j * w * np.exp(1j * w * y**2 / 2.0) * integral


def F_point_lens(w, y, dps=30):
    """Closed-form point-lens amplification factor (Deguchi & Watson 1986;
    Takahashi & Nakamura 2003):

        F(w,y) = exp[ pi w/4 + i (w/2) ln(w/2) ] * Gamma(1 - i w/2)
                 * 1F1(i w/2, 1; i w y^2 / 2)

    Evaluated at `dps` decimal digits of precision with mpmath (1F1 with
    complex parameters is not in scipy.special). This is the fast evaluator
    the case scripts use; cross-checked against `F_radial_1d` and
    `F_bruteforce_2d` in tests/test_waveoptics.py, and against the w->0 and
    w->infinity limits derived independently above.
    """
    mp.mp.dps = dps
    w = mp.mpf(w)
    y = mp.mpf(y)
    i = mp.mpc(0, 1)
    prefac = mp.e ** (mp.pi * w / 4 + i * (w / 2) * mp.log(w / 2))
    gam = mp.gamma(1 - i * w / 2)
    hyp = mp.hyp1f1(i * w / 2, 1, i * w * y**2 / 2)
    return complex(prefac * gam * hyp)


def F_point_lens_array(w_array, y):
    """Vectorized convenience wrapper around F_point_lens for an array of w
    at fixed y (the common case: sweeping frequency at one impact parameter,
    or one frequency at many impact-parameter samples along an orbit)."""
    return np.array([F_point_lens(w, y) for w in np.atleast_1d(w_array)])


def w_of_frequency(f_hz, M_lens_msun):
    """w = 8 pi G M_L f / c^3, restoring physical units.

    Using geometrized time unit GM/c^3 (seconds) from units.py so this file
    has no hidden constants of its own.
    """
    from . import units
    M_sec = units.msun_to_seconds(M_lens_msun)
    return 8.0 * np.pi * M_sec * f_hz
