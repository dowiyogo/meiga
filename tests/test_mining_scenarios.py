from __future__ import annotations

import numpy as np

from characterization.forward.detectability import exposure_curve
from characterization.forward.mining_scenarios import (
    scenario_borehole_cavity,
    scenario_gallery_vein,
    scenario_tunnel_pillar,
)
from characterization.forward.voxel_grid import VoxelGrid
from characterization.geometry import DetectorGeometry


def test_mining_scenarios_return_valid_density_fields() -> None:
    scenarios = [
        scenario_borehole_cavity(grid_resolution_m=1.0),
        scenario_gallery_vein(grid_resolution_m=1.0),
        scenario_tunnel_pillar(grid_resolution_m=1.0),
    ]
    for grid, rho, label in scenarios:
        assert isinstance(grid, VoxelGrid)
        assert rho.shape == grid.shape
        assert rho.dtype.kind == "f"
        assert label in {"borehole_cavity", "gallery_vein", "tunnel_pillar"}

    _, cavity, _ = scenarios[0]
    assert np.isclose(cavity.min(), 0.0)
    assert np.isclose(cavity.max(), 2.7)

    _, vein, _ = scenarios[1]
    assert np.isclose(vein.max(), 4.5)
    assert np.isclose(vein.min(), 2.7)

    _, pillar, _ = scenarios[2]
    assert np.isclose(pillar.min(), 0.0)
    assert np.isclose(pillar.max(), 2.7)


def test_borehole_cavity_design_point_reaches_five_sigma() -> None:
    curve = exposure_curve(
        scenario_borehole_cavity,
        [30.0 * 86400.0],
        geom=DetectorGeometry(),
        grid_resolution_m=1.0,
        cavity_radius_m=2.0,
        cavity_distance_m=8.0,
        overburden_m=5.0,
    )
    z_30d = float(curve["Z"].iloc[0])
    assert z_30d >= 5.0, (
        "borehole cavity design point should be detectable at 30 days; "
        f"got {z_30d:.3g} sigma"
    )
