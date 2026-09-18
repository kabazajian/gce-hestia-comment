"""Validate the corrected J-factor estimator against an analytic profile.

The point of these tests is that the estimator is trustworthy BEFORE any real
simulation data arrives — `validation/muru_comment_prep.md` item 4.3. An
estimator validated after the fact, on data whose answer we do not independently
know, would be worth much less.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

pytest.importorskip("scipy")

from particle_projection import (  # noqa: E402
    analytic_j_profile,
    knn_density,
    nfw_density,
    project_annihilation,
    project_muru,
)

RS, RHO_S, R0 = 20.0, 1.0, 8.0


def sample_nfw(n, rmax=60.0, rs=RS, seed=0):
    """Draw n equal-mass particles from a spherical NFW out to rmax."""
    rng = np.random.default_rng(seed)
    r = np.linspace(1e-3, rmax, 40000)
    m = np.cumsum(nfw_density(r, rs=rs, rho_s=RHO_S) * r ** 2)
    m /= m[-1]
    u = rng.random(n)
    rr = np.interp(u, m, r)
    cos_t = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    st = np.sqrt(1 - cos_t ** 2)
    pos = np.column_stack([rr * st * np.cos(phi), rr * st * np.sin(phi), rr * cos_t])
    # total mass so that the mean density matches the analytic normalization
    shell = np.trapezoid(4 * np.pi * r ** 2 * nfw_density(r, rs=rs, rho_s=RHO_S), r)
    return pos, np.full(n, shell / n)


def test_knn_density_recovers_the_analytic_profile():
    """Bias, not scatter, is what would corrupt an integral — so check bias."""
    pos, mass = sample_nfw(200_000)
    rho = knn_density(pos, mass, k=32)
    r = np.linalg.norm(pos, axis=1)
    sel = (r > 1.0) & (r < 30.0)
    ratio = rho[sel] / nfw_density(r[sel], rs=RS, rho_s=RHO_S)
    assert 0.85 < np.median(ratio) < 1.15, f"kNN density biased: median {np.median(ratio):.3f}"


def test_annihilation_estimator_matches_analytic_j():
    """🚨 THE validation: Σ m_i ρ_i / s_i² per beam must reproduce ∫ρ² ds."""
    pos, mass = sample_nfw(400_000)
    rho = knn_density(pos, mass, k=32)

    grid_l = np.array([2.0, 5.0, 10.0])          # off-centre, where MC noise is manageable
    grid_b = np.array([0.0])
    got = project_annihilation(pos, mass, rho, np.array([R0, 0.0, 0.0]),
                               grid_l, grid_b, aperture_deg=2.0, max_dist=50.0)[0]
    want = analytic_j_profile(grid_l, observer_dist=R0, s_max=50.0, rs=RS, rho_s=RHO_S)

    ratio = got / want
    assert np.all((ratio > 0.6) & (ratio < 1.6)), (
        f"J estimator off: got {got}, analytic {want}, ratio {ratio}")


def test_their_statistic_is_a_different_functional():
    """The two estimators must not be proportional — that is the whole point.

    If they were, the objection would be about normalization rather than
    physics. They are not: the correct J falls off far faster with angle,
    because ρ² weights the inner halo.
    """
    pos, mass = sample_nfw(300_000)
    rho = knn_density(pos, mass, k=32)
    obs = np.array([R0, 0.0, 0.0])
    grid_l = np.array([1.0, 4.0, 16.0])
    grid_b = np.array([0.0])

    j = project_annihilation(pos, mass, rho, obs, grid_l, grid_b,
                             aperture_deg=2.0, max_dist=15.0)[0]
    m = project_muru(pos, mass, obs, grid_l, grid_b,
                     aperture_deg=2.0, max_dist=15.0, exponent=2)[0]

    j_fall = j[0] / j[-1]
    m_fall = m[0] / m[-1]
    assert j_fall > m_fall * 1.5, (
        f"correct J should fall off much faster with angle: "
        f"J {j_fall:.1f}x vs their statistic {m_fall:.1f}x over 1->16 deg")


def test_flat_lb_metric_matches_their_aperture_not_ours():
    """metric="flat_lb" is Muru et al.'s KDTree convention (raw (l,b) degrees,
    no cos b): agrees with our flat_lcosb metric on the plane, diverges off it.

    The anchor (scripts/muru_redo/anchor/) pins the flat_lb branch against
    their actual Julia pipeline to 5e-7; this test only keeps the two branches
    from being silently conflated in refactors.
    """
    pos, mass = sample_nfw(150_000)
    obs = np.array([R0, 0.0, 0.0])
    gl = np.array([2.0, 5.0])
    b_plane, b_high = np.array([0.0]), np.array([18.0])

    for gb, tol, must_differ in ((b_plane, 2e-3, False), (b_high, None, True)):
        a = project_muru(pos, mass, obs, gl, gb, aperture_deg=3.0,
                         exponent=1, metric="flat_lcosb")
        m = project_muru(pos, mass, obs, gl, gb, aperture_deg=3.0,
                         exponent=1, metric="flat_lb")
        rel = np.abs(a - m) / np.maximum(a, 1e-30)
        if must_differ:
            assert rel.max() > 5e-3, (
                f"metrics indistinguishable at |b|=18: {rel.max():.2e}")
        else:
            assert rel.max() < tol, (
                f"metrics should agree on the plane: {rel.max():.2e}")


def test_muru_statistic_squaring_is_monotonic_relabeling():
    """Their exponent cannot change isophote SHAPE — the Comment's core point."""
    pos, mass = sample_nfw(150_000)
    obs = np.array([R0, 0.0, 0.0])
    gl = np.linspace(-10, 10, 41)
    gb = np.linspace(-10, 10, 41)
    a = project_muru(pos, mass, obs, gl, gb, aperture_deg=3.0, exponent=1)
    b = project_muru(pos, mass, obs, gl, gb, aperture_deg=3.0, exponent=2)

    # the level sets of b at c² are exactly those of a at c
    lvl = 0.5
    ca = a / a.max() >= np.sqrt(lvl)
    cb = b / b.max() >= lvl
    assert np.array_equal(ca, cb), "squaring must be a pure relabeling of contours"
