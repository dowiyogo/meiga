"""Sparse forward matrix and flux model for muon transmission tomography.

Minimal example:

```python
import numpy as np
from characterization.geometry import DetectorGeometry
from characterization.forward.system_matrix import build_F, opacity
from characterization.forward.voxel_grid import VoxelGrid

geom = DetectorGeometry(n_channels=2)
grid = VoxelGrid((-1000, -1000, 0), (1000, 1000, 1000), (2, 2, 2))
F = build_F(geom, grid)
gamma = opacity(F, np.ones(grid.n_total) * 2.7)
```
"""

from __future__ import annotations

import logging
import time
from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray
from scipy import sparse

from characterization.forward.dda import traverse_python, traverse_many
from characterization.forward.frames import detector_bottom_z
from characterization.forward.voxel_grid import VoxelGrid
from characterization.geometry import DetectorGeometry

LOGGER = logging.getLogger(__name__)

CSDA_G_PER_CM2_PER_GEV = 500.0
DEFAULT_P_MIN_GEV = 1.0


def _fast_pair_kinematics(
    geom: DetectorGeometry,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    centers = geom.pixel_centers()
    n_pairs = geom.n_pixels * geom.n_pixels
    origins = np.zeros((n_pairs, 3), dtype=float)
    directions = np.zeros((n_pairs, 3), dtype=float)
    theta = np.zeros(n_pairs, dtype=float)
    phi = np.zeros(n_pairs, dtype=float)
    omega = np.zeros(n_pairs, dtype=float)
    z_bot = detector_bottom_z(geom.D)
    row = 0
    area_mm2 = geom.pixel_area
    for top in range(geom.n_pixels):
        top_xy = centers[top]
        for bot in range(geom.n_pixels):
            bot_xy = centers[bot]
            origin = np.asarray([bot_xy[0], bot_xy[1], z_bot], dtype=float)
            top_point = np.asarray([top_xy[0], top_xy[1], 0.0], dtype=float)
            vector = top_point - origin
            distance = float(np.linalg.norm(vector))
            direction = vector / distance
            origins[row] = origin
            directions[row] = direction
            theta[row] = np.arccos(np.clip(direction[2], -1.0, 1.0))
            phi[row] = np.arctan2(direction[1], direction[0])
            omega[row] = area_mm2 * geom.D * geom.D / (distance**4)
            row += 1
    return origins, directions, theta, phi, omega


def _oversample_offsets(
    geom: DetectorGeometry, oversample_n: int
) -> NDArray[np.float64]:
    if oversample_n == 1:
        return np.zeros((1, 4), dtype=float)
    indices = np.arange(oversample_n, dtype=float) + 0.5
    fractions = np.column_stack(
        [
            (indices * 0.7548776662466927) % 1.0,
            (indices * 0.5698402909980532) % 1.0,
            (indices * 0.4385794979158758) % 1.0,
            (indices * 0.3275021914881652) % 1.0,
        ]
    )
    return (fractions - 0.5) * geom.bar_width


def _pair_rays(
    geom: DetectorGeometry, oversample_n: int
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.int64]]:
    centers = geom.pixel_centers()
    n_pairs = geom.n_pixels * geom.n_pixels
    n_rays = n_pairs * oversample_n
    origins = np.zeros((n_rays, 3), dtype=float)
    directions = np.zeros((n_rays, 3), dtype=float)
    ray_rows = np.zeros(n_rays, dtype=np.int64)
    z_bot = detector_bottom_z(geom.D)
    offsets = _oversample_offsets(geom, oversample_n)
    cursor = 0
    for top in range(geom.n_pixels):
        top_xy = centers[top]
        for bot in range(geom.n_pixels):
            bot_xy = centers[bot]
            row = top * geom.n_pixels + bot
            for offset in offsets:
                origin = np.asarray(
                    [bot_xy[0] + offset[0], bot_xy[1] + offset[1], z_bot],
                    dtype=float,
                )
                top_point = np.asarray(
                    [top_xy[0] + offset[2], top_xy[1] + offset[3], 0.0],
                    dtype=float,
                )
                vector = top_point - origin
                origins[cursor] = origin
                directions[cursor] = vector / float(np.linalg.norm(vector))
                ray_rows[cursor] = row
                cursor += 1
    return origins, directions, ray_rows


def build_F(
    geom: DetectorGeometry,
    grid: VoxelGrid,
    oversample_n: int = 1,
    detector_offset_mm: tuple[float, float, float] | None = None,
    use_batched: bool = True,
) -> sparse.csr_matrix:
    """Build sparse system matrix ``F`` with lengths in cm.

    Args:
        geom: Detector geometry in mm.
        grid: VoxelGrid in world frame.
        oversample_n: Number of pencil beams per pixel pair. ``1`` uses the
            center ray. Larger values deterministically sample the bottom/top
            pixel footprints with a low-discrepancy sequence and average
            lengths.
        detector_offset_mm: Optional rigid detector translation in mm. This is
            useful for frame-invariance tests; the default keeps the detector
            top module at ``z = 0``.
        use_batched: When true, use the batched DDA path. When false, use the
            pure-Python reference loop for regression checks.

    Returns:
        CSR matrix of shape ``(n_pairs, n_voxels)``. Entries are intersection
        lengths ``l_ij,k`` in cm.

    Edge cases:
        The dense matrix is never materialized. Rays missing the grid contribute
        empty rows.
    """

    if oversample_n < 1:
        raise ValueError("oversample_n must be >= 1")
    start = time.perf_counter()
    origins, directions, ray_rows = _pair_rays(geom, oversample_n)
    offset = np.asarray(detector_offset_mm or (0.0, 0.0, 0.0), dtype=float)
    if use_batched:
        rows, cols, data = traverse_many(
            origins + offset, directions, grid, return_coo=True
        )
        rows = ray_rows[rows]
        if oversample_n > 1:
            data = data / float(oversample_n)
    else:
        row_parts = []
        col_parts = []
        data_parts = []
        for ray, (origin, direction) in enumerate(zip(origins + offset, directions)):
            indices, lengths = traverse_python(origin, direction, grid)
            if len(indices) == 0:
                continue
            row_parts.append(np.full(len(indices), ray_rows[ray], dtype=np.int64))
            col_parts.append(indices)
            data_parts.append(lengths / float(oversample_n))
        if row_parts:
            rows = np.concatenate(row_parts)
            cols = np.concatenate(col_parts)
            data = np.concatenate(data_parts)
        else:
            rows = np.asarray([], dtype=np.int64)
            cols = np.asarray([], dtype=np.int64)
            data = np.asarray([], dtype=float)
    matrix = sparse.csr_matrix(
        (data, (rows, cols)),
        shape=(geom.n_pixels * geom.n_pixels, grid.n_total),
    )
    LOGGER.info(
        "Built F shape=%s nnz=%s in %.3f s",
        matrix.shape,
        matrix.nnz,
        time.perf_counter() - start,
    )
    return matrix


def opacity(
    F: sparse.csr_matrix, rho_vector: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Return opacity ``gamma = F @ rho`` in g/cm^2.

    Args:
        F: Sparse matrix with entries in cm.
        rho_vector: Density vector in g/cm^3.

    Returns:
        Opacity per pixel pair in g/cm^2.
    """

    rho = np.asarray(rho_vector, dtype=float)
    if rho.ndim != 1 or rho.shape[0] != F.shape[1]:
        raise ValueError("rho_vector must be 1D with length F.shape[1]")
    return np.asarray(F @ rho, dtype=float)


def _power_law_integral(
    p_min: NDArray[np.float64], gamma: float, p0: float, integral_at_1: float
) -> NDArray[np.float64]:
    norm = integral_at_1 * (gamma - 1.0) * (1.0 + p0) ** (gamma - 1.0)
    return norm / (gamma - 1.0) * (p_min + p0) ** (1.0 - gamma)


def flux_integral_above(
    theta: NDArray[np.float64],
    phi: NDArray[np.float64],
    p_min_gev: NDArray[np.float64] | float = DEFAULT_P_MIN_GEV,
    model: str | Callable[..., NDArray[np.float64]] = "reyna",
) -> NDArray[np.float64]:
    """Integrated sea-level muon flux above momentum threshold.

    Units are muons / (cm^2 s sr). The default Reyna-like model is normalized so
    the vertical flux above 1 GeV/c is 70 / (m^2 s sr) = 0.007 / (cm^2 s sr).
    """

    theta_arr = np.asarray(theta, dtype=float)
    p_min = np.asarray(p_min_gev, dtype=float)
    angular = np.clip(np.cos(theta_arr), 0.0, 1.0) ** 2
    if callable(model):
        return np.asarray(model(theta_arr, phi, p_min), dtype=float)
    if model == "reyna":
        return angular * _power_law_integral(
            p_min, gamma=2.7, p0=0.6, integral_at_1=0.007
        )
    if model == "gaisser":
        return angular * _power_law_integral(
            p_min, gamma=2.65, p0=1.0, integral_at_1=0.007
        )
    raise ValueError("model must be 'reyna', 'gaisser', or a callable")


def flux_model(
    theta: NDArray[np.float64] | float,
    phi: NDArray[np.float64] | float,
    momentum: NDArray[np.float64] | float | None = None,
    model: str | Callable[..., NDArray[np.float64]] = "reyna",
) -> NDArray[np.float64]:
    """Differential cosmic muon flux at sea level.

    Args:
        theta: Zenith angle in rad.
        phi: Azimuth angle in rad; currently azimuth-symmetric.
        momentum: Momentum in GeV/c. If ``None``, returns the integral above
            1 GeV/c.
        model: ``"reyna"``, ``"gaisser"``, or a callable.

    Returns:
        If momentum is supplied, differential flux in
        muons/(cm^2 s sr GeV). Otherwise integrated flux in muons/(cm^2 s sr).

    Edge cases:
        These are smooth parameterizations for Phase C-1. They are not a full
        atmospheric transport calculation.
    """

    theta_arr = np.asarray(theta, dtype=float)
    phi_arr = np.asarray(phi, dtype=float)
    if momentum is None:
        return flux_integral_above(theta_arr, phi_arr, DEFAULT_P_MIN_GEV, model=model)
    p = np.asarray(momentum, dtype=float)
    angular = np.clip(np.cos(theta_arr), 0.0, 1.0) ** 2
    if callable(model):
        return np.asarray(model(theta_arr, phi_arr, p), dtype=float)
    if model == "reyna":
        gamma, p0 = 2.7, 0.6
    elif model == "gaisser":
        gamma, p0 = 2.65, 1.0
    else:
        raise ValueError("model must be 'reyna', 'gaisser', or a callable")
    norm = 0.007 * (gamma - 1.0) * (1.0 + p0) ** (gamma - 1.0)
    return angular * norm * (p + p0) ** (-gamma)


def expected_counts(
    geom: DetectorGeometry,
    grid: VoxelGrid,
    F: sparse.csr_matrix,
    rho_vector: NDArray[np.float64],
    exposure_s: float,
    flux_kwargs: dict | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return no-absorber and absorber expected counts per pixel pair.

    Args:
        geom: Detector geometry in mm.
        grid: Voxel grid; included to keep frame/unit metadata explicit.
        F: Sparse system matrix in cm.
        rho_vector: Density vector in g/cm^3.
        exposure_s: Exposure in seconds.
        flux_kwargs: Optional arguments. Recognized keys include ``model`` and
            ``p_min_gev``.

    Returns:
        ``(N0, N)`` arrays. ``N0`` uses the no-absorber threshold. ``N`` uses a
        CSDA cutoff ``p_cut = max(p_min, gamma / 500)`` GeV/c, where
        ``500 g/cm^2/GeV`` approximates minimum-ionizing muon range.

    Edge cases:
        This Phase C-1 simplification ignores stochastic losses, multiple
        scattering, and detector response. It should be validated against Geant4
        before quantitative inversion.
    """

    del grid
    if exposure_s <= 0:
        raise ValueError("exposure_s must be positive")
    kwargs = flux_kwargs or {}
    model = kwargs.get("model", "reyna")
    p_min = float(kwargs.get("p_min_gev", DEFAULT_P_MIN_GEV))
    gamma_opacity = opacity(F, rho_vector)
    _, _, theta, phi, omega = _fast_pair_kinematics(geom)
    area_cm2 = geom.pixel_area * 0.01
    base_flux = flux_integral_above(theta, phi, p_min, model=model)
    p_cut = np.maximum(p_min, gamma_opacity / CSDA_G_PER_CM2_PER_GEV)
    attenuated_flux = flux_integral_above(theta, phi, p_cut, model=model)
    n0 = exposure_s * area_cm2 * omega * base_flux
    n = exposure_s * area_cm2 * omega * attenuated_flux
    return n0, n


def weight_matrix_gamma(
    counts: NDArray[np.float64], floor: float = 1.0
) -> NDArray[np.float64]:
    """Return diagonal weights for opacity-space inverse problems.

    Args:
        counts: Expected or observed counts per ray.
        floor: Minimum variance stabilizer in counts.

    Returns:
        One-dimensional diagonal entries for ``W_gamma``. This utility is built
        in C-1 and consumed by the later inverse phase.
    """

    return 1.0 / np.sqrt(np.maximum(np.asarray(counts, dtype=float), floor))
