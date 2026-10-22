from __future__ import annotations

import numpy as np

from characterization.forward.dda import _ray_box_intersection, traverse
from characterization.forward.system_matrix import (
    build_F,
    expected_counts,
    flux_integral_above,
    opacity,
)
from characterization.forward.synthetic_scenarios import add_cavity, homogeneous_slab
from characterization.forward.voxel_grid import VoxelGrid
from characterization.geometry import DetectorGeometry


def test_dda_chord_length_conservation_random_rays() -> None:
    grid = VoxelGrid(
        origin=(-5.0, -5.0, -5.0), voxel_size=(1.0, 1.0, 1.0), n_voxels=(10, 10, 10)
    )
    rng = np.random.default_rng(10)
    for _ in range(50):
        origin = rng.uniform(-20.0, 20.0, 3)
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        hit = _ray_box_intersection(
            origin, direction, np.asarray(grid.origin), grid.max_corner
        )
        indices, lengths = traverse(origin, direction, grid)
        if hit is None or hit[1] <= max(hit[0], 0.0):
            assert len(indices) == 0
            continue
        expected_cm = (hit[1] - max(hit[0], 0.0)) * 0.1
        assert abs(float(np.sum(lengths)) - expected_cm) < 1e-9


def test_translation_invariance_of_F() -> None:
    geom = DetectorGeometry(n_channels=2, D=20.0)
    grid = VoxelGrid(
        origin=(-1000.0, -1000.0, 0.0),
        voxel_size=(500.0, 500.0, 500.0),
        n_voxels=(4, 4, 4),
    )
    shift = (2000.0, -1500.0, 300.0)
    shifted_grid = VoxelGrid(
        origin=tuple(np.asarray(grid.origin) + np.asarray(shift)),
        voxel_size=grid.voxel_size,
        n_voxels=grid.n_voxels,
    )
    F = build_F(geom, grid)
    F_shift = build_F(geom, shifted_grid, detector_offset_mm=shift)
    diff = (F - F_shift).data
    assert np.allclose(diff, 0.0, atol=1e-12)


def test_homogeneous_slab_opacity_matches_analytic() -> None:
    geom = DetectorGeometry(n_channels=1, D=100.0)
    grid = VoxelGrid(
        origin=(-5000.0, -5000.0, 0.0),
        voxel_size=(1000.0, 1000.0, 1000.0),
        n_voxels=(10, 10, 10),
    )
    F = build_F(geom, grid)
    rho = homogeneous_slab(grid, 2.7).ravel()
    gamma = opacity(F, rho)
    assert abs(gamma[0] - 2.7 * 1000.0) / (2.7 * 1000.0) < 1e-10

    geom2 = DetectorGeometry(n_channels=2, D=100.0)
    F2 = build_F(geom2, grid)
    gamma2 = opacity(F2, homogeneous_slab(grid, 1.0).ravel())
    trajectories = np.asarray([row.sum() for row in F2])
    nonzero = trajectories > 0
    assert np.max(np.abs(gamma2[nonzero] - trajectories[nonzero])) < 1e-6


def test_cavity_contrast_central_pair() -> None:
    geom = DetectorGeometry(n_channels=1, D=100.0)
    grid = VoxelGrid(
        origin=(-5000.0, -5000.0, 0.0),
        voxel_size=(500.0, 500.0, 500.0),
        n_voxels=(20, 20, 40),
    )
    F = build_F(geom, grid)
    background = homogeneous_slab(grid, 2.7)
    cavity = add_cavity(background, grid, center=(0.0, 0.0, 10_000.0), radius=1000.0)
    gamma_bg = opacity(F, background.ravel())[0]
    gamma_cavity = opacity(F, cavity.ravel())[0]
    expected_delta = 2.7 * 200.0
    assert abs((gamma_bg - gamma_cavity) - expected_delta) / expected_delta < 0.3
    _, n_bg = expected_counts(geom, grid, F, background.ravel(), 3600.0)
    _, n_cavity = expected_counts(geom, grid, F, cavity.ravel(), 3600.0)
    assert n_cavity[0] > n_bg[0]


def test_flux_model_sanity_vertical_rate() -> None:
    rate_cm2 = flux_integral_above(
        np.asarray([0.0]), np.asarray([0.0]), 1.0, model="reyna"
    )[0]
    rate_m2 = rate_cm2 * 10_000.0
    assert abs(rate_m2 - 70.0) / 70.0 < 1e-12
