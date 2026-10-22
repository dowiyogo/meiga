from __future__ import annotations

import time

import numpy as np
import pytest

from characterization.forward.dda import (
    HAS_NUMBA,
    traverse,
    traverse_many,
    traverse_njit,
    traverse_python,
)
from characterization.forward.system_matrix import build_F
from characterization.forward.voxel_grid import VoxelGrid
from characterization.geometry import DetectorGeometry


def _random_regression_rays(n_rays: int = 1000) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(12345)
    origins = rng.uniform([-250.0, -250.0, -50.0], [250.0, 250.0, 450.0], (n_rays, 3))
    directions = rng.normal(size=(n_rays, 3))
    return origins, directions


def test_traverse_njit_matches_python_reference_random_rays() -> None:
    grid = VoxelGrid(
        origin=(-123.0, -80.0, 10.0),
        voxel_size=(17.0, 23.0, 31.0),
        n_voxels=(9, 7, 11),
    )
    origins, directions = _random_regression_rays()
    for origin, direction in zip(origins, directions):
        indices_py, lengths_py = traverse_python(origin, direction, grid)
        if HAS_NUMBA:
            indices_nb, lengths_nb = traverse_njit(origin, direction, grid)
        else:
            indices_nb, lengths_nb = traverse(origin, direction, grid)
        assert np.array_equal(indices_nb, indices_py)
        assert np.array_equal(lengths_nb, lengths_py)


def test_traverse_many_matches_per_ray_calls() -> None:
    grid = VoxelGrid(
        origin=(-100.0, -100.0, 0.0),
        voxel_size=(25.0, 25.0, 50.0),
        n_voxels=(8, 8, 8),
    )
    origins, directions = _random_regression_rays(64)
    batched = traverse_many(origins, directions, grid)
    expected = [
        traverse(origin, direction, grid)
        for origin, direction in zip(origins, directions)
    ]
    assert len(batched) == len(expected)
    for (indices_b, lengths_b), (indices_e, lengths_e) in zip(batched, expected):
        assert np.array_equal(indices_b, indices_e)
        assert np.array_equal(lengths_b, lengths_e)


def test_build_F_batched_matches_legacy_reference() -> None:
    geom = DetectorGeometry(n_channels=3, D=40.0)
    grid = VoxelGrid(
        origin=(-1000.0, -1000.0, 0.0),
        voxel_size=(500.0, 500.0, 500.0),
        n_voxels=(4, 4, 4),
    )
    f_batched = build_F(geom, grid, use_batched=True)
    f_legacy = build_F(geom, grid, use_batched=False)
    diff = f_batched - f_legacy
    assert diff.nnz == 0 or np.max(np.abs(diff.data)) < 1.0e-12


@pytest.mark.skipif(not HAS_NUMBA, reason="numba is unavailable")
def test_build_F_numba_standard_benchmark_under_20s() -> None:
    geom = DetectorGeometry()
    grid = VoxelGrid(
        origin=(-50_000.0, -50_000.0, 0.0),
        voxel_size=(1_000.0, 1_000.0, 1_000.0),
        n_voxels=(100, 100, 50),
    )
    warm_grid = VoxelGrid(
        origin=(-1_000.0, -1_000.0, 0.0),
        voxel_size=(1_000.0, 1_000.0, 1_000.0),
        n_voxels=(2, 2, 2),
    )
    build_F(DetectorGeometry(n_channels=1), warm_grid)
    start = time.perf_counter()
    matrix = build_F(geom, grid)
    elapsed = time.perf_counter() - start
    assert matrix.shape == (65_536, 500_000)
    assert elapsed < 20.0
