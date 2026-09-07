"""Wave-optics gravitational lensing of a scalar field, point-mass lens.

Notation follows wiki/conventions.md (dimensionless lens-plane position x,
dimensionless source position/impact parameter y, dimensionless frequency
w = 8 pi G M_L f / c^3). Full derivation: theory/theory.tex.

*** SIGN-CONVENTION FIX (2026-09-07, see wiki/log.md) ***
Every F(w,y) below is the COMPLEX CONJUGATE of the literal Fresnel integral
F_raw(w,y) = (w/2*pi*i) * int d^2x exp(i*w*phi(x,y)) (theory.tex Eq. 2.5) --
each function computes F_raw internally (unchanged from the original
derivation) and returns its conjugate. This was found to be necessary by an
independent review's causality test, confirmed here independently: with
F_raw applied directly (i.e. WITHOUT this conjugation) as
h_lensed(f)=F_raw(f,y)*h_unlensed(f) and inverse-transformed with
numpy.fft.irfft (which implements THIS project's own declared convention,
wiki/conventions.md: h(t)=int h(f) e^{+2pi i f t} df), a short test pulse's
weak (saddle-point) image comes out BEFORE the strong (minimum) image --
acausal. Conjugating fixes it: re-running the same test, the weak image
lands AFTER the strong one, at the correct delay Delta_t, with the correct
amplitude ratio sqrt(|mu_-/mu_+|). The root cause: wiki/conventions.md
asserted "a later-arriving image gets a +i*w*Delta_t phase" in the
h(f)=int h(t) e^{-2*pi*i*f*t} dt convention -- that assertion is simply
wrong (the standard Fourier delay theorem gives e^{-i*2*pi*f*Delta_t} for a
delay of Delta_t>0 in THAT convention, not e^{+i*2*pi*f*Delta_t}), and
every formula below (all independently re-derived or cited under that one
wrong premise) inherited the same overall sign error, consistently, which
is exactly why they kept passing their OWN cross-checks against each other
without it ever showing up: |F| is conjugation-invariant, so every
magnitude-only check (Paczynski magnification, the w->0 and w->infinity
limits, the three-evaluator F(w,y) agreement, every Case B pulse-height and
Case A fringe-amplitude plot) already gave the same answer either way.
Only phase-sensitive, time-domain reconstruction (Case A's lensed
waveform/peak amplification, and the new causality regression test below)
could catch this, and none of the original tests did.

*** REFERENCE-PHASE NORMALIZATION FIX (2026-09-07, see wiki/log.md) ***
Every F(w,y) below now also carries an extra factor exp(i*w*phi_+(y)),
phi_+(y) = fermat_potential(x_plus(y), y) the MINIMUM (strong) image's own
Fermat potential. Without it (the state of this file until this fix), F's
absolute phase was the raw, un-subtracted phi(x,y) of theory.tex Eq. 2.5 --
mathematically self-consistent (all four evaluators agreed with each other,
since none of them subtracted anything, so every existing magnitude- or
relative-phase-only test kept passing), but physically not what "F(w,y)"
conventionally means in this literature (Takahashi & Nakamura 2003, their
eq. 17-ish reference-phase subtraction): since ln|x| solves the deflection
Poisson equation only up to an arbitrary additive constant, the ABSOLUTE
phase w*phi_+(y) an un-normalized F carries is not a physical observable --
only the phase DIFFERENCE between images (i.e. Delta_T(y), which this
constant cancels out of) is. An independent review found this by asking
why Case A's lensed merger peak lands 0.60s "before" the unlensed one: that
number turned out to be EXACTLY 4*G*M_L*phi_+(y_A)/c^3 (agreement <2e-4s,
under half a sample) -- an arbitrary convention artifact, not the
interference effect an earlier version of this project's own analysis
(RESULTS.md, report.tex, wiki/log.md) mistakenly concluded it was. With
this fix, F is referenced to the strong image's OWN arrival: the
geometric-optics limit becomes exactly F -> sqrt(mu_+) (zero phase) as
w -> infinity along the strong-image-only direction, the lensed merger
peaks at EXACTLY the unlensed merger time (not 0.6s before it), and the
second-image echo lands at EXACTLY t_peak_unlensed + Delta_T(y) (not offset
by the same arbitrary 0.6s). |F(w,y)| is unaffected (a unit-magnitude
phase factor), so every Case B result, every diffraction-pattern figure,
and rms_strain_amplification (Parseval) are unchanged; only Case A's
absolute waveform TIMING changes.


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

    F -> sqrt(mu_+) + i*sqrt(|mu_-|)*exp(-i*w*DeltaT), the standard two-image
    stationary-phase (Morse-index) result, REFERENCED TO THE STRONG (minimum,
    x_+) IMAGE'S OWN ARRIVAL: its own Fermat potential phi_+(y) has been
    subtracted (see module docstring, "REFERENCE-PHASE NORMALIZATION FIX"),
    so it carries zero phase here by construction and the saddle image x_-
    carries only the physical, convention-independent relative delay
    Delta_T(y) plus its Morse phase (factor i, in the e^{-i2*pi*f*t}
    convention fixed in wiki/conventions.md). An earlier version of this
    function kept the un-subtracted exp(i*w*phi_+) prefactor instead
    (matching `F_point_lens` before ITS normalization fix, so the two still
    agreed -- see wiki/log.md for why that agreement was not, by itself,
    evidence of correctness).
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


def F_bruteforce_2d(w, y, x_max=40.0, n=1200, eta=0.0):
    """Direct 2D quadrature of F(w,y) = (w/2pi i) int d^2x exp[i w phi(x,y)],
    referenced to the strong image (see module docstring).

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
    F_raw = (w / (2j * np.pi)) * integral
    x_plus, _ = image_positions(y)
    phi_plus = fermat_potential(x_plus, y)
    # sign-convention + reference-phase fixes, see module docstring:
    return np.conjugate(F_raw) * np.exp(1j * w * phi_plus)


def F_radial_1d(w, y, eta=0.05, safety=25.0, min_breakpoints=60, max_breakpoints=3000):
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
    expression above is the literal published one (no reference phase);
    after conjugating (sign-convention fix) an extra factor
    exp(i*w*phi_+(y)) is applied to reference it to the strong image's own
    arrival (reference-phase fix) -- see module docstring, both dated
    2026-09-07.
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
    shrinking fast on either side of that point). First noticed as a stall:
    evaluating `F_point_lens` at the actual (w,y) used in Case A took ~150 s
    for one frequency sweep; `F_hybrid` needs no mpmath calls at all there
    (see wiki/log.md).
    """
    if w < w_geo_threshold:
        return F_point_lens(w, y)
    return complex(F_geometric_optics(w, y))


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
