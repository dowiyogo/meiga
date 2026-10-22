"""Detectability utilities for mining scenarios.

Minimal example:

```python
from characterization.forward.detectability import exposure_curve
from characterization.forward.mining_scenarios import scenario_borehole_cavity

curve = exposure_curve(scenario_borehole_cavity, [3600, 86400])
```
"""

from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from characterization.forward.system_matrix import build_F, expected_counts
from characterization.geometry import DetectorGeometry


def angular_signature(
    F,
    rho_field_with_anomaly: NDArray[np.float64],
    rho_field_baseline: NDArray[np.float64],
    exposure_s: float,
    geom: DetectorGeometry | None = None,
    grid=None,
) -> dict[str, NDArray[np.float64]]:
    """Return per-pixel-pair anomaly counts and Poisson significance.

    Args:
        F: Sparse forward matrix in cm.
        rho_field_with_anomaly: Density field including anomaly, g/cm^3.
        rho_field_baseline: Baseline density field, g/cm^3.
        exposure_s: Exposure in seconds.
        geom: Detector geometry. Defaults to nominal 16-channel geometry.
        grid: Voxel grid; required by ``expected_counts`` but not used directly.

    Returns:
        Dict with ``Delta_N_i``, ``sigma_Poisson`` and ``Z_i`` arrays.
    """

    geometry = geom or DetectorGeometry()
    _, n_anomaly = expected_counts(
        geometry,
        grid,
        F,
        rho_field_with_anomaly.ravel(),
        exposure_s,
        {"model": "reyna"},
    )
    _, n_baseline = expected_counts(
        geometry, grid, F, rho_field_baseline.ravel(), exposure_s, {"model": "reyna"}
    )
    delta = n_anomaly - n_baseline
    sigma = np.sqrt(np.maximum(n_baseline, 1.0e-12))
    return {"Delta_N_i": delta, "sigma_Poisson": sigma, "Z_i": delta / sigma}


def aggregated_significance(
    Z_map: NDArray[np.float64], aggregation: str = "quadrature"
) -> float:
    """Return global anomaly significance.

    ``quadrature`` computes ``sqrt(sum Z_i^2)`` and is an optimistic upper bound
    that ignores correlated systematic errors. Phase B covariance enters in the
    next stage.
    """

    z = np.asarray(Z_map, dtype=float)
    if aggregation == "quadrature":
        return float(np.sqrt(np.sum(z * z)))
    if aggregation == "max":
        return float(np.max(np.abs(z)))
    raise ValueError("aggregation must be 'quadrature' or 'max'")


def _baseline_for(label: str, rho: NDArray[np.float64]) -> NDArray[np.float64]:
    if label == "tunnel_pillar":
        return np.zeros_like(rho)
    if label == "gallery_vein":
        base = np.array(rho, copy=True)
        base[base > 2.7] = 2.7
        return base
    base = np.array(rho, copy=True)
    base[base < 2.7] = 2.7
    return base


def exposure_curve(
    scenario_func: Callable[..., tuple],
    exposure_grid_s: Iterable[float],
    geom: DetectorGeometry | None = None,
    **scenario_kwargs,
) -> pd.DataFrame:
    """Return aggregated significance versus exposure time."""

    geometry = geom or DetectorGeometry()
    grid, rho, label = scenario_func(**scenario_kwargs)
    baseline = _baseline_for(label, rho)
    F = build_F(geometry, grid)
    rows = []
    for exposure_s in exposure_grid_s:
        signature = angular_signature(
            F, rho, baseline, float(exposure_s), geometry, grid
        )
        rows.append(
            {
                "scenario": label,
                "exposure_s": float(exposure_s),
                "exposure_days": float(exposure_s) / 86400.0,
                "Z": aggregated_significance(signature["Z_i"]),
            }
        )
    return pd.DataFrame(rows)


def distance_curve(
    scenario_func: Callable[..., tuple],
    distance_grid_m: Iterable[float],
    exposure_s: float = 30.0 * 86400.0,
    geom: DetectorGeometry | None = None,
    distance_kwarg: str = "cavity_distance_m",
    **scenario_kwargs,
) -> pd.DataFrame:
    """Return aggregated significance versus anomaly distance."""

    rows = []
    for distance_m in distance_grid_m:
        kwargs = dict(scenario_kwargs)
        kwargs[distance_kwarg] = float(distance_m)
        curve = exposure_curve(scenario_func, [exposure_s], geom=geom, **kwargs)
        rows.append(
            {
                "distance_m": float(distance_m),
                "exposure_s": exposure_s,
                "Z": float(curve["Z"].iloc[0]),
            }
        )
    return pd.DataFrame(rows)
