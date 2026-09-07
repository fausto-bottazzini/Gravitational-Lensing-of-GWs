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
                              omega_out, t_peri, M_lens_msun, D_L_pc, D_S_pc):
    """y(t): dimensionless source-plane impact parameter from the outer
    Keplerian orbit, projected onto the sky and scaled by the angular
    Einstein radius. `a_out` in the same physical length unit as the
    distances (pc here, for consistency with D_L/D_S)."""
    r, nu = relative_separation(t, a_out, e_out, period_out, t_peri)
    x_sky, y_sky, _ = orbital_plane_to_sky(r, nu, omega_out, i_out, Omega_out)
    rho_pc = np.hypot(x_sky, y_sky)  # projected transverse offset, pc
    theta_rad = rho_pc / D_S_pc      # small-angle: physical offset / D_S
    theta_E_rad = einstein_radius_angle(M_lens_msun, D_L_pc, D_S_pc)
    return theta_rad / theta_E_rad, theta_E_rad
