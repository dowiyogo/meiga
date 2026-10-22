"""Mining-relevant density scenarios for the compact hodoscope.

Minimal example:

```python
from characterization.forward.mining_scenarios import scenario_borehole_cavity

grid, rho, label = scenario_borehole_cavity()
```

Density references:
    Common host-rock densities around 2.6-2.8 g/cm^3 are typical of crystalline
    silicate rocks. Magnetite-rich or massive sulfide mineralization can reach
    roughly 4-5 g/cm^3. These values are engineering inputs for synthetic
    detectability studies, not site-specific claims.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from characterization.forward.synthetic_scenarios import add_cavity, add_vein
from characterization.forward.voxel_grid import VoxelGrid


def _grid_for_box(
    half_width_m: float, height_m: float, resolution_m: float
) -> VoxelGrid:
    n_xy = int(np.ceil(2.0 * half_width_m / resolution_m))
    n_z = int(np.ceil(height_m / resolution_m))
    return VoxelGrid(
        origin=(-half_width_m * 1000.0, -half_width_m * 1000.0, 0.0),
        voxel_size=(
            resolution_m * 1000.0,
            resolution_m * 1000.0,
            resolution_m * 1000.0,
        ),
        n_voxels=(n_xy, n_xy, n_z),
    )


def scenario_borehole_cavity(
    grid_resolution_m: float = 0.25,
    cavity_radius_m: float = 2.0,
    cavity_distance_m: float = 8.0,
    overburden_m: float = 15.0,
    rho_rock: float = 2.7,
) -> tuple[VoxelGrid, NDArray[np.float64], str]:
    """Return a close-range borehole void scenario.

    Args:
        grid_resolution_m: Cubic voxel pitch in m.
        cavity_radius_m: Spherical void radius in m.
        cavity_distance_m: Cavity-center height above detector in m.
        overburden_m: Rock thickness above the top of the cavity in m.
        rho_rock: Host-rock density in g/cm^3.

    Returns:
        ``(grid, rho_field, label)``. The cavity has density 0 g/cm^3.
    """

    height = cavity_distance_m + cavity_radius_m + overburden_m
    half_width = max(4.0, cavity_radius_m + 2.5)
    grid = _grid_for_box(half_width, height, grid_resolution_m)
    rho = np.full(grid.shape, rho_rock, dtype=float)
    rho = add_cavity(
        rho,
        grid,
        center=(0.0, 0.0, cavity_distance_m * 1000.0),
        radius=cavity_radius_m * 1000.0,
        rho_inside=0.0,
    )
    return grid, rho, "borehole_cavity"


def scenario_gallery_vein(
    grid_resolution_m: float = 0.25,
    vein_thickness_m: float = 0.5,
    vein_dip_deg: float = 60.0,
    vein_distance_m: float = 10.0,
    rho_vein: float = 4.5,
    rho_host: float = 2.7,
    overburden_m: float = 10.0,
) -> tuple[VoxelGrid, NDArray[np.float64], str]:
    """Return a dense dipping vein in host rock.

    The vein density is representative of magnetite/sulfide-rich ore, while the
    host density is typical hard rock.
    """

    height = vein_distance_m + overburden_m
    grid = _grid_for_box(5.0, height, grid_resolution_m)
    rho = np.full(grid.shape, rho_host, dtype=float)
    dip = np.deg2rad(vein_dip_deg)
    length = 10.0
    start = (-0.5 * length * np.cos(dip) * 1000.0, -3000.0, vein_distance_m * 1000.0)
    end = (
        0.5 * length * np.cos(dip) * 1000.0,
        3000.0,
        (vein_distance_m + length * np.sin(dip)) * 1000.0,
    )
    rho = add_vein(
        rho,
        grid,
        start=start,
        end=end,
        radius=0.5 * vein_thickness_m * 1000.0,
        rho_inside=rho_vein,
    )
    return grid, rho, "gallery_vein"


def scenario_tunnel_pillar(
    grid_resolution_m: float = 0.5,
    pillar_radius_m: float = 3.0,
    pillar_height_m: float = 6.0,
    pillar_distance_m: float = 12.0,
    overburden_m: float = 8.0,
) -> tuple[VoxelGrid, NDArray[np.float64], str]:
    """Return a dense rock pillar embedded in an air-filled gallery background."""

    height = pillar_distance_m + pillar_height_m + overburden_m
    grid = _grid_for_box(6.0, height, grid_resolution_m)
    rho = np.zeros(grid.shape, dtype=float)
    z0 = pillar_distance_m * 1000.0
    z1 = (pillar_distance_m + pillar_height_m) * 1000.0
    xs = grid.origin[0] + (np.arange(grid.n_voxels[0]) + 0.5) * grid.voxel_size[0]
    ys = grid.origin[1] + (np.arange(grid.n_voxels[1]) + 0.5) * grid.voxel_size[1]
    zs = grid.origin[2] + (np.arange(grid.n_voxels[2]) + 0.5) * grid.voxel_size[2]
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    mask = (xx**2 + yy**2 <= (pillar_radius_m * 1000.0) ** 2) & (zz >= z0) & (zz <= z1)
    rho[mask] = 2.7
    rho[zz > z1] = 2.7
    return grid, rho, "tunnel_pillar"
