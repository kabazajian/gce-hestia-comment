"""Line-of-sight annihilation map, J(l,b) = int rho^2 ds, for a generalized NFW.

Vendored 2026-09-17 from the `gcereduce.dm_templates` module of a separate
project so that this repository stands alone; only the pieces the Comment's
gNFW^2 comparison template needs are kept (spherical or ellipsoidal gNFW; no
cored profiles, no ROI J-factors). Numerically identical to the original:
same node placement, same quadrature.

PROFILE  (Abazajian et al. 2020, PRD 102, 023023, Eq. 2), normalised so that
rho(R_sun) = rho_sun exactly for any (gamma, r_s):

    rho(r) = rho_sun (R_sun/r)^gamma [(r_s + R_sun)/(r_s + r)]^(3-gamma)

Ellipsoidal variant: r -> r' with r'^2 = x^2 + (y/q_b)^2 + (z/q_c)^2, z along
the disk normal. Defaults are spherical.

QUADRATURE. A uniform grid in s cannot resolve a cusp: the integrand of
int rho^2 ds peaks over a width ~p about the tangent point, p = R_sun sin(psi)
the impact parameter, and p falls below any fixed step for psi ~< 1 deg. The
substitution s = s0 + p sinh(u) makes r = p cosh(u) exactly for a spherical
profile, so the integrand is smooth in u and a uniform u grid resolves it at
every angle. Converged to 1e-6 between n_s = 400 and 3200 for gamma = 1.2.

Units: kpc, GeV cm^-3; j_map returns GeV^2 cm^-5 sr^-1 (only the shape matters
for the Comment, which normalises every template).
"""
from __future__ import annotations

import numpy as np

KPC_CM = 3.0856775814913673e21

#: The parameter set used for the Comment's gNFW^2 template: the centres of
#: the priors marginalised in Abazajian et al. 2020 (rho_sun = 0.28 GeV cm^-3,
#: r_s = 26 kpc), gamma = 1.2, R_sun = 8.18 kpc, spherical.
DEFAULTS = dict(gamma=1.2, rs=26.0, rho_sun=0.28, R_sun=8.18, q_b=1.0, q_c=1.0)


def rho_gnfw(r, gamma=1.2, rs=26.0, rho_sun=0.28, R_sun=8.18):
    r = np.maximum(np.asarray(r, dtype=float), 1e-8)
    return rho_sun * (R_sun / r) ** gamma * ((rs + R_sun) / (rs + r)) ** (3.0 - gamma)


def _los_nodes(grid_l, grid_b, R_sun, s_max, n_s, q_b, q_c):
    """Sample points along each sight line, concentrated at closest approach.
    Returns (r, jacobian, t) broadcasting to (nb, nl, n_s)."""
    L = np.deg2rad(np.asarray(grid_l, dtype=float))[None, :]
    B = np.deg2rad(np.asarray(grid_b, dtype=float))[:, None]
    cos_psi = np.cos(B) * np.cos(L)
    s0 = R_sun * cos_psi
    p_imp = R_sun * np.sqrt(np.maximum(1.0 - cos_psi ** 2, 1e-12))
    u_lo = np.arcsinh((0.0 - s0) / p_imp)
    u_hi = np.arcsinh((s_max - s0) / p_imp)
    t = np.linspace(0.0, 1.0, n_s)
    u = u_lo[..., None] + (u_hi - u_lo)[..., None] * t
    S = s0[..., None] + p_imp[..., None] * np.sinh(u)
    jac = p_imp[..., None] * np.cosh(u) * (u_hi - u_lo)[..., None]
    Lc, Bc = L[..., None], B[..., None]
    x = R_sun - S * np.cos(Bc) * np.cos(Lc)
    y = -S * np.cos(Bc) * np.sin(Lc)
    z = S * np.sin(Bc)
    r = np.sqrt(x ** 2 + (y / q_b) ** 2 + (z / q_c) ** 2)
    return r, jac, t


def j_map(grid_l, grid_b, s_max=100.0, n_s=800, **kw):
    """int rho^2 ds on the (grid_b x grid_l) mesh, GeV^2 cm^-5 sr^-1.
    Keyword overrides: gamma, rs, rho_sun, R_sun, q_b, q_c. Unknown keys raise
    rather than being silently ignored."""
    p = dict(DEFAULTS)
    unknown = set(kw) - set(p)
    if unknown:
        raise TypeError(f"unknown parameter(s) {sorted(unknown)}; expected {sorted(p)}")
    p.update(kw)
    r, jac, t = _los_nodes(grid_l, grid_b, p["R_sun"], s_max, n_s, p["q_b"], p["q_c"])
    rho = rho_gnfw(r, p["gamma"], p["rs"], p["rho_sun"], p["R_sun"])
    return np.trapezoid(rho ** 2 * jac, t, axis=-1) * KPC_CM
