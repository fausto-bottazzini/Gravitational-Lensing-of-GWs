"""Lensing geometry: physical distances/masses <-> the dimensionless (w, y)
used everywhere in waveoptics.py.

Local (Galactic-scale) triple, no cosmological redshift (wiki/conventions.md):
D_L, D_S ordinary Euclidean distances, D_LS = D_S - D_L, z_L = 0.
"""
import numpy as np

from . import units


def einstein_radius_angle(M_lens_msun, D_L_pc, D_S_pc):
    """theta_E = sqrt(4 G M_L/c^2 * D_LS/(D_L D_S)) (point lens, z_L=0)."""
    D_LS_pc = D_S_pc - D_L_pc
    if D_LS_pc <= 0:
        raise ValueError("lens must be closer than the source (D_L < D_S)")
    R_L_m = units.msun_to_meters(M_lens_msun)  # GM/c^2, metres
    D_L_m = D_L_pc * units.PC_SI
    D_S_m = D_S_pc * units.PC_SI
    D_LS_m = D_LS_pc * units.PC_SI
    return np.sqrt(4.0 * R_L_m * D_LS_m / (D_L_m * D_S_m))


def y_of_angular_offset(theta_rad, theta_E_rad):
    return theta_rad / theta_E_rad


# ---------------------------------------------------------------------
# Kepler orbit: relative separation vector of the binary CM around the lens
# ---------------------------------------------------------------------

def solve_kepler_equation(M_anom, e, tol=1e-13, max_iter=100):
    """Eccentric anomaly E from mean anomaly M via Newton-Raphson on
    M = E - e sin E (standard; e.g. Murray & Dermott, "Solar System
    Dynamics" Sec. 2.4). Vectorized over M_anom."""
    M_anom = np.atleast_1d(np.asarray(M_anom, dtype=float))
    E = M_anom.copy() if e < 0.8 else np.full_like(M_anom, np.pi)
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M_anom
        fp = 1.0 - e * np.cos(E)
        dE = -f / fp
        E += dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def true_anomaly(E, e):
    return 2.0 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2.0),
                             np.sqrt(1 - e) * np.cos(E / 2.0))


def relative_separation(t, a, e, period, t_peri=0.0):
    """Relative separation r(t) and true anomaly nu(t) of a Keplerian orbit,
    a=semi-major axis, e=eccentricity, period=orbital period, t_peri=time of
    periapsis passage. Returns (r, nu), same units as `a`."""
    n = 2.0 * np.pi / period
    M_anom = n * (np.asarray(t, dtype=float) - t_peri)
    E = solve_kepler_equation(M_anom, e)
    nu = true_anomaly(E, e)
    r = a * (1.0 - e * np.cos(E))
    return r, nu


def orbital_plane_to_sky(r, nu, omega, i, Omega):
    """Rotate a position (r,nu) in the orbital plane into the observer
    frame: z = line of sight (toward the observer), (x,y) = sky plane.
    Standard R_z(Omega) R_x(i) R_z(omega) rotation (e.g. Murray & Dermott
    Sec. 2.8, adapted: i=0 is face-on/orbital plane == sky plane, i=pi/2 is
    edge-on with the line of sight in the orbital plane).

    Returns (x_sky, y_sky, z_los), same shape as r/nu.
    """
    theta = nu + omega  # argument of latitude
    ct, st = np.cos(theta), np.sin(theta)
    cO, sO = np.cos(Omega), np.sin(Omega)
    ci, si = np.cos(i), np.sin(i)

    x_sky = r * (cO * ct - sO * st * ci)
    y_sky = r * (sO * ct + cO * st * ci)
    z_los = r * st * si
    return x_sky, y_sky, z_los


def impact_parameter_of_time(t, a_out, e_out, period_out, i_out, Omega_out,
                              omega_out, t_peri, M_lens_msun, D_L_pc):
    """y(t): dimensionless source-plane impact parameter from the outer
    Keplerian orbit of the source (binary CM) around the lens (M_lens,
    effectively fixed at the system barycenter since M_lens >> M_binary).

    IMPORTANT and easy to get wrong: unlike textbook cosmological lensing,
    here the lens and source are in the SAME physical system, so D_LS is NOT
    a large, ~fixed distance -- it IS the (time-dependent) line-of-sight
    separation between them, `z_los(t)` from `orbital_plane_to_sky`, at the
    AU-to-sub-pc scale of the outer orbit itself, while D_L (source-to-
    observer) is set by D_L (kpc-to-pc scale, effectively also D_S to
    superb approximation: |z_los|/D_L ~ 1e-8 here, dropped everywhere
    EXCEPT in D_LS itself, where it is the whole story):

        theta_E(t)^2 = (4 G M_lens/c^2) * D_LS(t) / (D_L * D_S(t))
                     ~ (4 G M_lens/c^2) * z_los(t) / D_L^2   (D_S ~ D_L)

    Two direct physical consequences, both real, not approximation
    artefacts:
      * z_los(t) < 0 (source instantaneously nearer Earth than the lens,
        i.e. in FRONT of it) -> D_LS<0 -> no lensing geometry exists for
        that half of the orbit. Returned as `lensed=False`, y=inf, F=1
        should be used by the caller.
      * z_los(t) -> 0+ (source crossing the plane through the lens
        perpendicular to the line of sight) -> theta_E -> 0 -> y -> infinity
        for any nonzero transverse offset -> effectively unlensed there too,
        with NO special-casing needed: the formula does this on its own.
    This is exactly the origin of the "repeated lensing" pulses once per
    outer orbit found in hierarchical-triple GW lensing (only the far-side
    half-orbit is lensed at all) -- see D'Orazio & Loeb (2020) and
    theory/theory.tex Sec. 5.

    Returns (y, lensed_mask, theta_E_rad) -- y is np.inf where lensed_mask
    is False.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    r, nu = relative_separation(t, a_out, e_out, period_out, t_peri)
    x_sky, y_sky, z_los = orbital_plane_to_sky(r, nu, omega_out, i_out, Omega_out)
    rho_pc = np.hypot(x_sky, y_sky)          # projected transverse offset, pc
    D_LS_pc = z_los                          # exact in this approximation
    lensed_mask = D_LS_pc > 0.0

    theta_E_rad = np.zeros_like(t)
    theta_E_rad[lensed_mask] = np.sqrt(
        4.0 * units.msun_to_meters(M_lens_msun) * (D_LS_pc[lensed_mask] * units.PC_SI)
        / (D_L_pc * units.PC_SI) ** 2
    )
    theta_rad = rho_pc / (D_L_pc)  # D_S ~ D_L; rho_pc/D_L_pc is already the (dimensionless) small angle in rad since both in pc

    y = np.full_like(t, np.inf)
    y[lensed_mask] = theta_rad[lensed_mask] / theta_E_rad[lensed_mask]
    return y, lensed_mask, theta_E_rad
