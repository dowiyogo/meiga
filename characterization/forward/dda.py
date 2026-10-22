"""Deterministic Amanatides-Woo DDA traversal for voxel grids.

Minimal example:

```python
import numpy as np
from characterization.forward.dda import traverse
from characterization.forward.voxel_grid import VoxelGrid

grid = VoxelGrid((0, 0, 0), (10, 10, 10), (2, 2, 2))
indices, lengths_cm = traverse(np.array([-5, 5, 5]), np.array([1, 0, 0]), grid)
```
"""

from __future__ import annotations

import logging

import numpy as np
from numpy.typing import ArrayLike, NDArray

from characterization.forward.frames import CM_PER_MM
from characterization.forward.voxel_grid import VoxelGrid

LOGGER = logging.getLogger(__name__)

try:
    from numba import njit

    HAS_NUMBA = True
except ImportError:  # pragma: no cover - environment dependent
    HAS_NUMBA = False
    LOGGER.warning("numba is unavailable; DDA will use the pure-Python path")

    def njit(*args, **kwargs):  # type: ignore[no-redef]
        def wrapper(func):
            return func

        return wrapper


def _ray_box_intersection(
    origin: NDArray[np.float64],
    direction: NDArray[np.float64],
    bounds_min: NDArray[np.float64],
    bounds_max: NDArray[np.float64],
) -> tuple[float, float] | None:
    t_min = -np.inf
    t_max = np.inf
    for axis in range(3):
        if abs(direction[axis]) < 1e-15:
            if origin[axis] < bounds_min[axis] or origin[axis] > bounds_max[axis]:
                return None
            continue
        inv = 1.0 / direction[axis]
        t0 = (bounds_min[axis] - origin[axis]) * inv
        t1 = (bounds_max[axis] - origin[axis]) * inv
        if t0 > t1:
            t0, t1 = t1, t0
        t_min = max(t_min, t0)
        t_max = min(t_max, t1)
        if t_max < t_min:
            return None
    return float(t_min), float(t_max)


def traverse_python(
    origin: ArrayLike,
    direction: ArrayLike,
    grid: VoxelGrid,
    t_max: float | None = None,
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Pure-Python reference traversal through ``grid``.

    Args:
        origin: Ray origin in mm.
        direction: Ray direction. It is normalized internally.
        grid: VoxelGrid in the same world frame.
        t_max: Optional maximum ray parameter in mm. If omitted, traversal ends
            at grid exit.

    Returns:
        ``(voxel_indices, segment_lengths)`` where voxel indices are flat
        C-order indices and lengths are in cm.

    Edge cases:
        Rays missing the grid return empty arrays. Rays starting inside the grid
        work. Boundary grazing is made deterministic by nudging the entry point
        by a tiny distance along the ray direction.
    """

    ray_origin = np.asarray(origin, dtype=float)
    ray_direction = np.asarray(direction, dtype=float)
    norm = float(np.linalg.norm(ray_direction))
    if norm <= 0:
        raise ValueError("direction must be non-zero")
    ray_direction = ray_direction / norm

    bounds_min = np.asarray(grid.origin, dtype=float)
    bounds_max = grid.max_corner
    hit = _ray_box_intersection(ray_origin, ray_direction, bounds_min, bounds_max)
    if hit is None:
        return np.asarray([], dtype=np.int64), np.asarray([], dtype=float)
    enter, exit_ = hit
    enter = max(enter, 0.0)
    if t_max is not None:
        exit_ = min(exit_, float(t_max))
    if exit_ <= enter:
        return np.asarray([], dtype=np.int64), np.asarray([], dtype=float)

    eps = 1e-9
    point = ray_origin + (enter + eps) * ray_direction
    current = grid.position_to_index(point)
    if current is None:
        point = ray_origin + (enter + 1e-6) * ray_direction
        current = grid.position_to_index(point)
    if current is None:
        return np.asarray([], dtype=np.int64), np.asarray([], dtype=float)

    voxel_size = np.asarray(grid.voxel_size, dtype=float)
    grid_origin = np.asarray(grid.origin, dtype=float)
    n_voxels = np.asarray(grid.n_voxels, dtype=int)
    step = np.sign(ray_direction).astype(int)
    step[step == 0] = 0
    t_delta = np.full(3, np.inf, dtype=float)
    t_next = np.full(3, np.inf, dtype=float)
    current_arr = np.asarray(current, dtype=int)
    for axis in range(3):
        if ray_direction[axis] > 0:
            boundary = grid_origin[axis] + (current_arr[axis] + 1) * voxel_size[axis]
            t_next[axis] = (boundary - ray_origin[axis]) / ray_direction[axis]
            t_delta[axis] = voxel_size[axis] / ray_direction[axis]
        elif ray_direction[axis] < 0:
            boundary = grid_origin[axis] + current_arr[axis] * voxel_size[axis]
            t_next[axis] = (boundary - ray_origin[axis]) / ray_direction[axis]
            t_delta[axis] = -voxel_size[axis] / ray_direction[axis]

    indices: list[int] = []
    lengths: list[float] = []
    t_current = enter
    while (
        np.all(current_arr >= 0)
        and np.all(current_arr < n_voxels)
        and t_current < exit_
    ):
        axis = int(np.argmin(t_next))
        t_segment_end = min(float(t_next[axis]), exit_)
        segment = max(0.0, t_segment_end - t_current)
        if segment > 0:
            indices.append(grid.flat_index(tuple(int(x) for x in current_arr)))
            lengths.append(segment * CM_PER_MM)
        t_current = t_segment_end
        if t_current >= exit_:
            break
        tied = np.where(np.isclose(t_next, t_next[axis], rtol=0.0, atol=1e-10))[0]
        for tied_axis in tied:
            current_arr[tied_axis] += step[tied_axis]
            t_next[tied_axis] += t_delta[tied_axis]

    return np.asarray(indices, dtype=np.int64), np.asarray(lengths, dtype=float)


@njit(cache=True)
def _traverse_one_njit(
    ray_origin,
    ray_direction,
    grid_origin,
    voxel_size,
    n_voxels,
    t_max_value,
    use_t_max,
    out_indices,
    out_lengths,
):
    norm = 0.0
    for axis in range(3):
        norm += ray_direction[axis] * ray_direction[axis]
    norm = np.sqrt(norm)
    if norm <= 0.0:
        return 0
    direction = np.empty(3, dtype=np.float64)
    for axis in range(3):
        direction[axis] = ray_direction[axis] / norm

    bounds_max = np.empty(3, dtype=np.float64)
    for axis in range(3):
        bounds_max[axis] = grid_origin[axis] + voxel_size[axis] * n_voxels[axis]

    t_enter = -1.0e300
    t_exit = 1.0e300
    for axis in range(3):
        if abs(direction[axis]) < 1.0e-15:
            if (
                ray_origin[axis] < grid_origin[axis]
                or ray_origin[axis] > bounds_max[axis]
            ):
                return 0
        else:
            inv = 1.0 / direction[axis]
            t0 = (grid_origin[axis] - ray_origin[axis]) * inv
            t1 = (bounds_max[axis] - ray_origin[axis]) * inv
            if t0 > t1:
                tmp = t0
                t0 = t1
                t1 = tmp
            if t0 > t_enter:
                t_enter = t0
            if t1 < t_exit:
                t_exit = t1
            if t_exit < t_enter:
                return 0

    if t_enter < 0.0:
        t_enter = 0.0
    if use_t_max and t_max_value < t_exit:
        t_exit = t_max_value
    if t_exit <= t_enter:
        return 0

    point = np.empty(3, dtype=np.float64)
    for axis in range(3):
        point[axis] = ray_origin[axis] + (t_enter + 1.0e-9) * direction[axis]

    current = np.empty(3, dtype=np.int64)
    for axis in range(3):
        rel = (point[axis] - grid_origin[axis]) / voxel_size[axis]
        current[axis] = int(np.floor(rel))
        if current[axis] < 0 or current[axis] >= n_voxels[axis]:
            point[axis] = ray_origin[axis] + (t_enter + 1.0e-6) * direction[axis]
            rel = (point[axis] - grid_origin[axis]) / voxel_size[axis]
            current[axis] = int(np.floor(rel))
            if current[axis] < 0 or current[axis] >= n_voxels[axis]:
                return 0

    step = np.zeros(3, dtype=np.int64)
    t_delta = np.empty(3, dtype=np.float64)
    t_next = np.empty(3, dtype=np.float64)
    for axis in range(3):
        t_delta[axis] = 1.0e300
        t_next[axis] = 1.0e300
        if direction[axis] > 0.0:
            step[axis] = 1
            boundary = grid_origin[axis] + (current[axis] + 1) * voxel_size[axis]
            t_next[axis] = (boundary - ray_origin[axis]) / direction[axis]
            t_delta[axis] = voxel_size[axis] / direction[axis]
        elif direction[axis] < 0.0:
            step[axis] = -1
            boundary = grid_origin[axis] + current[axis] * voxel_size[axis]
            t_next[axis] = (boundary - ray_origin[axis]) / direction[axis]
            t_delta[axis] = -voxel_size[axis] / direction[axis]

    count = 0
    t_current = t_enter
    ny = n_voxels[1]
    nz = n_voxels[2]
    while (
        current[0] >= 0
        and current[0] < n_voxels[0]
        and current[1] >= 0
        and current[1] < n_voxels[1]
        and current[2] >= 0
        and current[2] < n_voxels[2]
        and t_current < t_exit
    ):
        axis_min = 0
        if t_next[1] < t_next[axis_min]:
            axis_min = 1
        if t_next[2] < t_next[axis_min]:
            axis_min = 2
        t_segment_end = t_next[axis_min]
        if t_segment_end > t_exit:
            t_segment_end = t_exit
        segment = t_segment_end - t_current
        if segment > 0.0:
            out_indices[count] = (current[0] * ny + current[1]) * nz + current[2]
            out_lengths[count] = segment * 0.1
            count += 1
        t_current = t_segment_end
        if t_current >= t_exit:
            break
        target = t_next[axis_min]
        for axis in range(3):
            if abs(t_next[axis] - target) <= 1.0e-10:
                current[axis] += step[axis]
                t_next[axis] += t_delta[axis]
    return count


@njit(cache=True)
def _traverse_many_coo_njit(
    origins, directions, grid_origin, voxel_size, n_voxels, t_max_value, use_t_max
):
    max_segments = n_voxels[0] + n_voxels[1] + n_voxels[2] + 8
    n_rays = origins.shape[0]
    rows = np.empty(n_rays * max_segments, dtype=np.int64)
    cols = np.empty(n_rays * max_segments, dtype=np.int64)
    data = np.empty(n_rays * max_segments, dtype=np.float64)
    tmp_indices = np.empty(max_segments, dtype=np.int64)
    tmp_lengths = np.empty(max_segments, dtype=np.float64)
    cursor = 0
    for row in range(n_rays):
        count = _traverse_one_njit(
            origins[row],
            directions[row],
            grid_origin,
            voxel_size,
            n_voxels,
            t_max_value,
            use_t_max,
            tmp_indices,
            tmp_lengths,
        )
        for j in range(count):
            rows[cursor] = row
            cols[cursor] = tmp_indices[j]
            data[cursor] = tmp_lengths[j]
            cursor += 1
    return rows[:cursor], cols[:cursor], data[:cursor]


def traverse_njit(
    origin: ArrayLike,
    direction: ArrayLike,
    grid: VoxelGrid,
    t_max: float | None = None,
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Numba-accelerated single-ray traversal, equivalent to ``traverse_python``."""

    max_segments = sum(grid.n_voxels) + 8
    indices = np.empty(max_segments, dtype=np.int64)
    lengths = np.empty(max_segments, dtype=float)
    count = _traverse_one_njit(
        np.asarray(origin, dtype=float),
        np.asarray(direction, dtype=float),
        np.asarray(grid.origin, dtype=float),
        np.asarray(grid.voxel_size, dtype=float),
        np.asarray(grid.n_voxels, dtype=np.int64),
        0.0 if t_max is None else float(t_max),
        t_max is not None,
        indices,
        lengths,
    )
    return indices[:count].copy(), lengths[:count].copy()


def traverse(
    origin: ArrayLike,
    direction: ArrayLike,
    grid: VoxelGrid,
    t_max: float | None = None,
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Traverse a ray, dispatching to numba when available."""

    if HAS_NUMBA:
        return traverse_njit(origin, direction, grid, t_max=t_max)
    return traverse_python(origin, direction, grid, t_max=t_max)


def traverse_many(
    origins: NDArray[np.float64],
    directions: NDArray[np.float64],
    grid: VoxelGrid,
    t_max: float | None = None,
    return_coo: bool = False,
) -> (
    list[tuple[NDArray[np.int64], NDArray[np.float64]]]
    | tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.float64]]
):
    """Traverse many rays through a grid.

    Args:
        origins: Array of shape ``(n_rays, 3)`` in mm.
        directions: Array of shape ``(n_rays, 3)``.
        grid: VoxelGrid.
        t_max: Optional maximum path length in mm.
        return_coo: When true, return batched ``(rows, cols, lengths_cm)`` COO
            triplets instead of a per-ray list.

    Returns:
        By default, a list of ``traverse`` outputs. With ``return_coo=True``,
        returns ``(rows, cols, lengths_cm)`` suitable for sparse matrix
        assembly. The traversal itself is deterministic; the list or row order
        follows the input order.
    """

    origins_arr = np.asarray(origins, dtype=float)
    directions_arr = np.asarray(directions, dtype=float)
    if HAS_NUMBA and return_coo:
        return _traverse_many_coo_njit(
            origins_arr,
            directions_arr,
            np.asarray(grid.origin, dtype=float),
            np.asarray(grid.voxel_size, dtype=float),
            np.asarray(grid.n_voxels, dtype=np.int64),
            0.0 if t_max is None else float(t_max),
            t_max is not None,
        )
    results = [
        traverse(origin, direction, grid, t_max=t_max)
        for origin, direction in zip(origins_arr, directions_arr)
    ]
    if return_coo:
        row_parts = []
        col_parts = []
        data_parts = []
        for row, (indices, lengths) in enumerate(results):
            if len(indices) == 0:
                continue
            row_parts.append(np.full(len(indices), row, dtype=np.int64))
            col_parts.append(indices)
            data_parts.append(lengths)
        if not row_parts:
            return (
                np.asarray([], dtype=np.int64),
                np.asarray([], dtype=np.int64),
                np.asarray([], dtype=float),
            )
        return (
            np.concatenate(row_parts),
            np.concatenate(col_parts),
            np.concatenate(data_parts),
        )
    return results
