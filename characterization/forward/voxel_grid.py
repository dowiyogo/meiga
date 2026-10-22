"""Voxel grid definitions for the forward muography model.

Minimal example:

```python
from characterization.forward.voxel_grid import VoxelGrid

grid = VoxelGrid(origin=(-5000, -5000, 0), voxel_size=(1000, 1000, 1000), n_voxels=(10, 10, 5))
center = grid.index_to_center((0, 0, 0))
```
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class VoxelGrid:
    """Axis-aligned voxel grid in the detector world frame.

    Attributes:
        origin: Minimum ``(x, y, z)`` corner in mm.
        voxel_size: Voxel side lengths ``(dx, dy, dz)`` in mm.
        n_voxels: Number of voxels along each axis.

    Edge cases:
        ``z`` must be allowed to start at zero for geology above the top module.
        Non-positive sizes or counts raise ``ValueError``.
    """

    origin: tuple[float, float, float]
    voxel_size: tuple[float, float, float]
    n_voxels: tuple[int, int, int]

    def __post_init__(self) -> None:
        origin = tuple(float(x) for x in self.origin)
        voxel_size = tuple(float(x) for x in self.voxel_size)
        n_voxels = tuple(int(x) for x in self.n_voxels)
        if any(size <= 0 for size in voxel_size):
            raise ValueError("voxel_size entries must be positive")
        if any(count <= 0 for count in n_voxels):
            raise ValueError("n_voxels entries must be positive")
        object.__setattr__(self, "origin", origin)
        object.__setattr__(self, "voxel_size", voxel_size)
        object.__setattr__(self, "n_voxels", n_voxels)

    @property
    def shape(self) -> tuple[int, int, int]:
        """Return ``(nx, ny, nz)``."""

        return self.n_voxels

    @property
    def n_total(self) -> int:
        """Return total number of voxels."""

        return int(np.prod(self.n_voxels))

    @property
    def max_corner(self) -> NDArray[np.float64]:
        """Return maximum ``(x, y, z)`` grid corner in mm."""

        return np.asarray(self.origin) + np.asarray(self.voxel_size) * np.asarray(
            self.n_voxels
        )

    def flat_index(self, index: tuple[int, int, int]) -> int:
        """Return flat index for ``(ix, iy, iz)`` using C-order layout."""

        ix, iy, iz = index
        nx, ny, nz = self.n_voxels
        if not (0 <= ix < nx and 0 <= iy < ny and 0 <= iz < nz):
            raise IndexError(index)
        return int((ix * ny + iy) * nz + iz)

    def unravel_index(self, flat_index: int) -> tuple[int, int, int]:
        """Return ``(ix, iy, iz)`` for a flat C-order voxel index."""

        nx, ny, nz = self.n_voxels
        if flat_index < 0 or flat_index >= self.n_total:
            raise IndexError(flat_index)
        ix, rem = divmod(int(flat_index), ny * nz)
        iy, iz = divmod(rem, nz)
        return ix, iy, iz

    def index_to_center(self, index: tuple[int, int, int] | int) -> NDArray[np.float64]:
        """Return voxel center in mm."""

        if isinstance(index, int):
            index = self.unravel_index(index)
        return np.asarray(self.origin) + (
            np.asarray(index, dtype=float) + 0.5
        ) * np.asarray(self.voxel_size)

    def bounds(
        self, index: tuple[int, int, int] | int
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return ``(min_corner, max_corner)`` for one voxel in mm."""

        if isinstance(index, int):
            index = self.unravel_index(index)
        lower = np.asarray(self.origin) + np.asarray(index, dtype=float) * np.asarray(
            self.voxel_size
        )
        return lower, lower + np.asarray(self.voxel_size)

    def position_to_index(self, point: ArrayLike) -> tuple[int, int, int] | None:
        """Return voxel index containing ``point`` in mm, or ``None`` if outside."""

        relative = (
            np.asarray(point, dtype=float) - np.asarray(self.origin)
        ) / np.asarray(self.voxel_size)
        index = np.floor(relative).astype(int)
        if np.any(index < 0) or np.any(index >= np.asarray(self.n_voxels)):
            return None
        return tuple(int(x) for x in index)


def from_topography(
    dem_array: NDArray[np.float64] | None = None,
    origin: tuple[float, float, float] = (-50_000.0, -50_000.0, 0.0),
    voxel_size: tuple[float, float, float] = (1_000.0, 1_000.0, 1_000.0),
    surface_density: float = 0.0,
    basement_density: float = 2.7,
) -> tuple[VoxelGrid, NDArray[np.float64]]:
    """Build a grid and default density from optional DEM topography.

    Args:
        dem_array: Optional surface height map in mm above detector top. If
            omitted, a flat 50 m slab is generated.
        origin: Grid minimum corner in mm.
        voxel_size: Voxel dimensions in mm.
        surface_density: Density above topography in g/cm^3.
        basement_density: Density below topography in g/cm^3.

    Returns:
        ``(grid, rho_default)`` where ``rho_default`` has shape
        ``grid.n_voxels`` and units g/cm^3.

    Edge cases:
        This is a coarse Phase C-1 helper. DEM cells are assumed to align with
        x-y voxels; non-aligned DEMs should be resampled before calling.
    """

    if dem_array is None:
        n_voxels = (100, 100, 50)
        grid = VoxelGrid(origin=origin, voxel_size=voxel_size, n_voxels=n_voxels)
        return grid, np.full(n_voxels, basement_density, dtype=float)

    dem = np.asarray(dem_array, dtype=float)
    if dem.ndim != 2:
        raise ValueError("dem_array must be 2D")
    max_height = float(np.nanmax(dem))
    nz = max(1, int(np.ceil(max_height / voxel_size[2])))
    grid = VoxelGrid(
        origin=origin,
        voxel_size=voxel_size,
        n_voxels=(dem.shape[0], dem.shape[1], nz),
    )
    rho = np.full(grid.shape, surface_density, dtype=float)
    z_centers = origin[2] + (np.arange(nz) + 0.5) * voxel_size[2]
    for iz, z_center in enumerate(z_centers):
        rho[:, :, iz] = np.where(z_center <= dem, basement_density, surface_density)
    return grid, rho
