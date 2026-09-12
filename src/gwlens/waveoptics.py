"""Wave-optics gravitational lensing of a scalar field, point-mass lens.

Notation: dimensionless lens-plane position x, source position / impact
parameter y, frequency w = 8 pi G M_L f / c^3. Derivations in
theory/theory.pdf chapters 2-3; conventions in wiki/conventions.md.

Two conventions are built into every F(w,y) this module returns, and both
change the answer:

  * Fourier sign. The literal Fresnel-Kirchhoff integral,
    F_raw = (w/2 pi i) int d^2x exp(i w phi(x,y)), is written in the lensing
    literature's convention, which is the opposite of this project's
    (h(f) = int h(t) e^{-2 pi i f t} dt, matching numpy.fft). Every function
    below computes F_raw and returns its CONJUGATE. Without that, the weak
    image arrives BEFORE the strong one. |F| is conjugation-invariant, so
    comparing evaluators against each other cannot detect a mistake here;
    only a time-domain reconstruction can, which is what
    tests/test_waveoptics.py::check_causality_of_lensed_pulse does.

  * Reference phase. psi(x) = ln|x| solves the deflection Poisson equation
    only up to an additive constant, so F's ABSOLUTE phase is not an
    observable. Every F below is referenced to the strong (minimum) image's
    own arrival, via a factor exp(i w phi_+(y)): that image then carries no
    phase, and only the physical delay Delta_T(y) between images survives.

Both of these were wrong once and were fixed on 2026-09-07; wiki/log.md has
the account, including why every evaluator agreeing with every other one was
not evidence that they were right.

Three evaluators of F(w,y), in decreasing order of cost:

  1. F_bruteforce_2d   direct 2D quadrature of the diffraction integral,
                       regularized by a convergence factor. No closed form
                       anywhere in it.
  2. F_radial_1d       the same integral with the azimuthal part done
                       analytically (Bessel J0 identity). Still from first
                       principles, much faster.
  3. F_point_lens      the closed-form confluent-hypergeometric expression
                       (Deguchi & Watson 1986; Takahashi & Nakamura 2003),
                       evaluated with mpmath. This is what the case scripts
                       use, through F_hybrid.

The geometric-optics quantities below -- image positions, magnifications,
time delay -- are derived independently of all three, and are used both as
the w -> infinity limit of F and, through the Paczynski (1986) total
magnification, as a cross-check on the magnification formula itself.
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

    F -> sqrt(mu_+) + i*sqrt(|mu_-|)*exp(-i*w*DeltaT), the standard two-image
    stationary-phase (Morse-index) result, REFERENCED TO THE STRONG (minimum,
    x_+) IMAGE'S OWN ARRIVAL: its own Fermat potential phi_+(y) has been
    subtracted (see the module docstring, "Reference phase"),
    so it carries zero phase here by construction and the saddle image x_-
    carries only the physical, convention-independent relative delay
    Delta_T(y) plus its Morse phase (factor i, in the e^{-i2*pi*f*t}
    convention fixed in wiki/conventions.md).
    """
    w = np.asarray(w, dtype=float)
    x_plus, x_minus = image_positions(y)
    mu_plus = magnification(x_plus)
    mu_minus = magnification(x_minus)
    dT = time_delay_difference(y)
    F_raw = np.sqrt(mu_plus) - 1j * np.sqrt(np.abs(mu_minus)) * np.exp(1j * w * dT)
    return np.conjugate(F_raw)  # sign-convention + reference-phase fixes, see module docstring


# ---------------------------------------------------------------------
# Wave optics: the diffraction integral, three ways
# ---------------------------------------------------------------------

def _phi2d(x1, x2, y):
    """2D Fermat potential, source at (y,0): (1/2)|x-y|^2 - ln|x|."""
    r2 = x1**2 + x2**2
    return 0.5 * ((x1 - y)**2 + x2**2) - 0.5 * np.log(r2)


def F_bruteforce_2d(w, y, x_max=25.0, n=1200, eta=0.005):
    """Direct 2D quadrature of F(w,y) = (w/2pi i) int d^2x exp[i w phi(x,y)],
    referenced to the strong image (see module docstring).

    Regularized with a Gaussian convergence factor exp(-eta*x^2). `eta` is
    NOT optional in practice and its default is not zero, because the
    unregularized integral does not converge under a hard cutoff: the tail
    of int d^2x exp(i w x^2/2) beyond radius R contributes
    (2 pi/i w)[exp(i w R^2/2) - 1], which the (w/2 pi i) prefactor turns into
    an O(1) oscillation in R that never decays. Measured, at w=y=1: the
    relative error is 4.6e-2 at (x_max=25, n=1200, eta=0) and still 4.4e-2 at
    n=2000 -- it does not shrink with resolution, only wanders with x_max
    (1.0e-2 at x_max=40). With eta=0.005 it is 1.1e-2 and behaves as the
    O(eta) regularization error it should be, the same way `F_radial_1d`'s
    does. The defaults here are exactly the arguments
    tests/test_waveoptics.py::check_bruteforce_2d_agrees passes, so the
    default call is the one that has actually been validated -- and it is
    validated ONLY near w~1. `x_max` and `n` fix a grid while the
    integrand oscillates faster with w, so with these same defaults the
    relative error is 2.3e-2 at (w=1, y=1.589), 1.2e-1 at (w=6.19,
    y=1.589) and 7.5e-1 at (w=10, y=1). Raise `n` and lower `eta`
    together for larger w, or use `F_radial_1d`, which is far cheaper.

    O(n^2) points on a square grid; slow, spot-check only.
    """
    xs = np.linspace(-x_max, x_max, n)
    dx = xs[1] - xs[0]
    x1, x2 = np.meshgrid(xs, xs, indexing="ij")
    r2 = x1**2 + x2**2
    mask = r2 > 1e-12  # exclude the point-lens singularity itself (measure zero)
    phase = w * _phi2d(x1, x2, y)
    integrand = np.where(mask, np.exp(1j * phase - eta * r2), 0.0)
    integral = np.sum(integrand) * dx * dx
    F_raw = (w / (2j * np.pi)) * integral
    x_plus, _ = image_positions(y)
    phi_plus = fermat_potential(x_plus, y)
    # sign-convention + reference-phase fixes, see module docstring:
    return np.conjugate(F_raw) * np.exp(1j * w * phi_plus)


def F_radial_1d(w, y, eta=0.01, safety=25.0, min_breakpoints=60, max_breakpoints=3000):
    """1D radial reduction of the same integral via
    int_0^2pi dtheta exp(i a cos theta) = 2 pi J0(a) (Abramowitz & Stegun
    9.1.21), applied to the cross term in |x-y|^2 = x^2+y^2-2xy*cos(theta):

        F(w,y) = -i w exp(i w y^2/2) int_0^inf dx x J0(w x y)
                 exp[i w (x^2/2 - ln x)]

    The integral is only conditionally convergent; regularized with a small
    positive `eta` added to the exponent's imaginary part (x^2 -> x^2*(1+i
    eta)), i.e. w -> w(1+i*eta), the standard i-epsilon prescription
    (Ulmer & Goodman 1995): the integrand acquires a factor
    exp(-eta*w*x^2/2), which is negligible once eta*w*x^2/2 > `safety`.
    That fixes an x_max = sqrt(2*safety/(eta*w)) automatically, and the
    breakpoint count is set from the number of oscillations
    (~w*x_max^2/(4pi)) expected below it, so this stays accurate without a
    hand-tuned cutoff. Convergence as eta -> 0 (at fixed accuracy) is
    checked in tests/test_waveoptics.py against `F_point_lens`.

    The residual regularization error is O(eta) at FIXED (w,y), but its
    coefficient is not 1. The regulator multiplies the integrand by
    exp(-eta*w*x^2/2), so the small parameter is eta*w*x_+^2/2, and the
    measured relative error tracks 1-exp(-eta*w*x_+^2/2) across the plane:

        w=1     y=1.000  eta=0.01  ->  1.3e-2
        w=25    y=1.589  eta=0.01  ->  4.9e-1   (just below F_hybrid's cut)
        w=62    y=1.589  eta=0.01  ->  6.5e-1

    So `eta` has to be chosen for the w in hand: the default eta=0.01 is
    ~1.3% only near w~1, and is worthless above w~10. At fixed (w,y) the
    error does fall linearly with eta -- at w=25, y=1.589 it runs 6.1e-2
    then 1.6e-2 for eta = 1e-3 then 2.5e-4 -- which is what makes the
    eta->0 extrapolation in tests/test_waveoptics.py an independent check
    on `F_point_lens` rather than a coincidence at one point.

    `max_breakpoints` caps a request of 4*n_osc = 2*safety/(pi*eta)
    subdivisions, which does not depend on w, so the cap binds for
    eta < 2*safety/(pi*max_breakpoints) = 5.3e-3 with the defaults --
    including at the eta values quoted just above. Measured, raising the
    cap to 4e5 at eta=1e-4 changes the answer by nothing at w=6.19: the
    quadrature is converged well before the cap, so the cap costs time
    rather than accuracy. It does mean the subdivision count is not a
    convergence indicator, and that eta, not the cap, is the knob.
    """
    from scipy.special import j0

    x_max = np.sqrt(2.0 * safety / (eta * w))
    n_osc = w * x_max**2 / (4.0 * np.pi)
    n_periods = int(np.clip(4 * n_osc, min_breakpoints, max_breakpoints))

    def integrand(x):
        if x <= 0:
            return 0.0 + 0.0j
        # exp(i * w*(1+i*eta) * x^2/2) = exp(i w x^2/2) * exp(-eta*w*x^2/2):
        # it is (1+i*eta), NOT (1-i*eta), that damps at large x. The other
        # sign makes the integral diverge.
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
    F_raw = -1j * w * np.exp(1j * w * y**2 / 2.0) * integral
    x_plus, _ = image_positions(y)
    phi_plus = fermat_potential(x_plus, y)
    # sign-convention + reference-phase fixes, see module docstring:
    return np.conjugate(F_raw) * np.exp(1j * w * phi_plus)


def F_point_lens(w, y, dps=30):
    """Closed-form point-lens amplification factor (Deguchi & Watson 1986;
    Takahashi & Nakamura 2003):

        F(w,y) = exp[ pi w/4 + i (w/2) ln(w/2) ] * Gamma(1 - i w/2)
                 * 1F1(i w/2, 1; i w y^2 / 2)

    Evaluated at `dps` decimal digits of precision with mpmath (1F1 with
    complex parameters is not in scipy.special). This is the fast evaluator
    the case scripts use; cross-checked against `F_radial_1d` and
    `F_bruteforce_2d` in tests/test_waveoptics.py, and against the w->0 and
    w->infinity limits derived independently above. The closed-form
    expression above is the literal published one, carrying no reference
    phase; it is conjugated and then multiplied by exp(i*w*phi_+(y)) to
    apply the two conventions described in the module docstring.
    """
    mp.mp.dps = dps
    w = mp.mpf(w)
    y = mp.mpf(y)
    i = mp.mpc(0, 1)
    prefac = mp.e ** (mp.pi * w / 4 + i * (w / 2) * mp.log(w / 2))
    gam = mp.gamma(1 - i * w / 2)
    hyp = mp.hyp1f1(i * w / 2, 1, i * w * y**2 / 2)
    F_raw = complex(prefac * gam * hyp)
    y_f = float(y)
    x_plus, _ = image_positions(y_f)
    phi_plus = fermat_potential(x_plus, y_f)
    # sign-convention + reference-phase fixes, see module docstring:
    return F_raw.conjugate() * complex(np.exp(1j * float(w) * phi_plus))


def F_hybrid(w, y, w_geo_threshold=30.0):
    """F(w,y), fast: `F_point_lens` (exact, mpmath) below `w_geo_threshold`,
    `F_geometric_optics` (closed-form trig/algebra, no special functions)
    above it.

    Why this is needed, not just convenient: mpmath's confluent
    hypergeometric series for `F_point_lens` converges slowly once its
    argument `i*w*y^2/2` gets large (many terms before the series settles),
    which happens exactly in the regime -- large w -- where
    `F_geometric_optics` is independently validated to agree with it to
    <0.1% (`tests/test_waveoptics.py::check_geometric_optics_limit`, errors
    3e-4 at w=10 falling to 3e-5 at w=1000). `w_geo_threshold=30` sits
    comfortably inside that validated agreement (3.1% at w=30 itself, per
    `check_hybrid_matches_at_threshold`, well under its 4% tolerance and
    shrinking fast on either side of that point). The cost is not marginal:
    one Case A frequency sweep through `F_point_lens` takes ~150 s, and
    through `F_hybrid` no mpmath call is needed at all.
    """
    if w < w_geo_threshold:
        return F_point_lens(w, y)
    return complex(F_geometric_optics(w, y))


def w_of_frequency(f_hz, M_lens_msun):
    """w = 8 pi G M_L f / c^3, restoring physical units.

    Using geometrized time unit GM/c^3 (seconds) from units.py so this file
    has no hidden constants of its own.
    """
    from . import units
    M_sec = units.msun_to_seconds(M_lens_msun)
    return 8.0 * np.pi * M_sec * f_hz
