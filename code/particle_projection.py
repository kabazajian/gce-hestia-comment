"""Annihilation (J-factor) maps from simulation particles — done correctly.

Written 2026-07-30 ahead of possible HESTIA data access, so that the corrected
calculation is ready and *validated* before the particles arrive rather than
after. Supports `validation/muru_comment_prep.md` item 4.3: do the corrected
calculation ourselves BEFORE submitting anything.

THE ESTIMATOR, AND WHY IT DIFFERS FROM A MASS SUM
--------------------------------------------------
The annihilation signal along a sight line is ``J(Ω) = ∫ ρ²(s) ds``. For a beam
of solid angle ΔΩ, using ``dV = s² ds dΩ``,

    J(Ω)·ΔΩ = ∫∫ ρ² ds dΩ = ∫ (ρ / s²) · ρ dV

and a particle realization replaces ``∫ f ρ dV → Σ_i m_i f_i``, giving

    **J(Ω)·ΔΩ  =  Σ_{i ∈ beam}  m_i · ρ_i / s_i²**

where ``ρ_i`` is the *local density at particle i* (needs a density estimator —
here kNN) and ``s_i`` is its distance from the observer.

Contrast that with the statistic used by Muru et al. (arXiv:2508.06314, PRL),
whose published code sums bare masses in an angular aperture and squares the
result afterwards:

    their map  =  ( Σ_{i ∈ beam} m_i )²

The correct estimator carries **two factors theirs does not**: the extra ``ρ_i``
(this is the difference between ∫ρ² ds and a projection of ρ) and the ``1/s_i²``
(the beam's volume grows as s², and the physical J-factor has no such
weighting). Both are visible in a single line of code below, which is the
cleanest way to state the objection.

`project_muru` reproduces their statistic deliberately, so the two can be run on
identical particles and plotted side by side — which is what a Comment figure
needs. It is NOT an endorsement of that statistic.
"""

from __future__ import annotations

import numpy as np

__all__ = ["knn_density", "project_annihilation", "project_muru", "nfw_density",
           "analytic_j_profile"]


def nfw_density(r, rs=20.0, rho_s=1.0, gamma=1.0):
    """Generalized NFW, ρ ∝ (r/rs)^-γ (1 + r/rs)^(γ-3)."""
    x = np.maximum(np.asarray(r, float) / rs, 1e-12)
    return rho_s * x ** (-gamma) * (1.0 + x) ** (gamma - 3.0)


def knn_density(pos: np.ndarray, mass: np.ndarray, k: int = 32,
                block: int = 2_000_000) -> np.ndarray:
    """Local density at each particle from its k nearest neighbours.

    ρ_i = (Σ of the k neighbour masses) / (4/3 π r_k³). Simple and adequate:
    the J-factor is an integral, so per-particle density noise averages down
    along a sight line. Bias matters more than scatter, which is why the test
    suite checks this against an analytic profile rather than trusting it.

    Queries run in blocks of `block` particles: an all-at-once query on a
    ~10M-particle HESTIA halo allocates ~10 GB for (d, idx); chunking caps
    peak memory at ~block*(k+1)*16 bytes with identical results (the tree is
    built once over all particles — only the queries are batched).
    """
    from scipy.spatial import cKDTree

    pos = np.asarray(pos, float)
    mass = np.asarray(mass, float)
    tree = cKDTree(pos)
    rho = np.empty(len(pos))
    for lo in range(0, len(pos), block):
        hi = min(lo + block, len(pos))
        d, idx = tree.query(pos[lo:hi], k=k + 1, workers=-1)  # 1st nb = self
        r_k = d[:, -1]
        enclosed = mass[idx[:, 1:]].sum(axis=1)
        vol = (4.0 / 3.0) * np.pi * np.maximum(r_k, 1e-12) ** 3
        rho[lo:hi] = enclosed / vol
    return rho


def _angles_from_observer(pos, observer):
    """(ℓ, b) in degrees and distance s, with the GC at the origin.

    Observer sits at `observer`; the GC direction is −observer normalized, and
    ℓ increases toward +y of the frame in which the disc plane is z = 0.
    """
    rel = np.asarray(pos, float) - np.asarray(observer, float)
    s = np.linalg.norm(rel, axis=1)
    ez = np.array([0.0, 0.0, 1.0])
    ex = -np.asarray(observer, float)
    ex = ex / np.linalg.norm(ex)               # toward the GC
    ey = np.cross(ez, ex)
    ey /= np.linalg.norm(ey)
    x, y, z = rel @ ex, rel @ ey, rel @ ez
    b = np.degrees(np.arcsin(np.clip(z / np.maximum(s, 1e-12), -1, 1)))
    l = np.degrees(np.arctan2(y, x))
    return l, b, s


def _beam_sum(l, b, s, weight, grid_l, grid_b, aperture_deg, metric="flat_lcosb"):
    """Sum `weight` over particles within `aperture_deg` of each grid direction.

    metric="flat_lcosb": flat (ℓcos b, b) — a true small-angle aperture,
    adequate for the inner ROI. metric="flat_lb": flat (ℓ, b) with NO cos b —
    NOT a constant angular aperture (it widens as 1/cos b in true angle), but
    it is exactly what Muru et al.'s published code does (KDTree built on raw
    (l, b) degrees, `BulgeDensityProjection.calc_KDTree`), so the anchor
    comparison must use it. Measured on G3.2: the two metrics differ by up to
    ~6% in summed mass at |b| = 20°.
    """
    from scipy.spatial import cKDTree

    L, B = np.meshgrid(grid_l, grid_b)
    if metric == "flat_lb":
        pts = np.column_stack([l, b])
        q = np.column_stack([L.ravel(), B.ravel()])
    elif metric == "flat_lcosb":
        pts = np.column_stack([l * np.cos(np.deg2rad(b)), b])
        q = np.column_stack([(L * np.cos(np.deg2rad(B))).ravel(), B.ravel()])
    else:
        raise ValueError(f"unknown metric {metric!r}")
    tree = cKDTree(pts)
    out = np.zeros(q.shape[0])
    for j, nb in enumerate(tree.query_ball_point(q, aperture_deg)):
        if nb:
            out[j] = weight[nb].sum()
    return out.reshape(L.shape)


def project_annihilation(pos, mass, dens, observer, grid_l, grid_b,
                         aperture_deg=1.0, max_dist=None):
    """CORRECT annihilation map: J(Ω) = Σ m_i ρ_i / s_i² per beam, / ΔΩ.

    Returns J in [mass²/length⁵/sr] for whatever units `mass`, `pos` are in.
    """
    l, b, s = _angles_from_observer(pos, observer)
    w = np.asarray(mass, float) * np.asarray(dens, float) / np.maximum(s, 1e-12) ** 2
    if max_dist is not None:
        w = np.where(s <= max_dist, w, 0.0)
    dOmega = 2.0 * np.pi * (1.0 - np.cos(np.deg2rad(aperture_deg)))
    return _beam_sum(l, b, s, w, grid_l, grid_b, aperture_deg) / dOmega


def project_muru(pos, mass, observer, grid_l, grid_b,
                 aperture_deg=3.0, max_dist=None, exponent=2,
                 metric="flat_lcosb"):
    """Reproduce the Muru et al. statistic: (Σ masses in the aperture)^exponent.

    Provided ONLY so the two can be compared on identical particles. Note it
    lacks both the ρ_i weighting and the 1/s_i² of the correct estimator.
    For an exact match to their published code pass metric="flat_lb" (their
    KDTree ignores the cos b factor — see `_beam_sum`).
    """
    l, b, s = _angles_from_observer(pos, observer)
    w = np.asarray(mass, float).copy()
    if max_dist is not None:
        w = np.where(s <= max_dist, w, 0.0)
    return _beam_sum(l, b, s, w, grid_l, grid_b, aperture_deg, metric) ** exponent


def analytic_j_profile(theta_deg, observer_dist=8.0, s_max=50.0, n=20000, **kw):
    """J(θ) = ∫ρ² ds by direct quadrature, for a SPHERICAL profile. Test truth."""
    theta = np.atleast_1d(np.deg2rad(theta_deg))
    s = np.linspace(1e-4, s_max, n)
    # law of cosines: r² = R0² + s² − 2 R0 s cos θ
    r = np.sqrt(observer_dist ** 2 + s[None, :] ** 2
                - 2 * observer_dist * s[None, :] * np.cos(theta)[:, None])
    return np.trapezoid(nfw_density(r, **kw) ** 2, s, axis=1)
