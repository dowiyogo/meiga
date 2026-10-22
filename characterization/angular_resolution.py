"""Angular-resolution estimators for the hodoscope.

Minimal reproducible example:

```python
from characterization.angular_resolution import sigma_theta_geom
from characterization.geometry import DetectorGeometry

sigma = sigma_theta_geom(DetectorGeometry(D=98.4))
```
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from characterization.geometry import DetectorGeometry

LOGGER = logging.getLogger(__name__)


def sigma_theta_geom(geom: DetectorGeometry) -> float:
    """Return the analytical pixel-limited zenith resolution in rad.

    Args:
        geom: Detector geometry. ``D`` and ``pixel_pitch`` are in mm.

    Returns:
        ``sigma_theta ~= d / (sqrt(6) * D)`` in rad, where ``d`` is the
        effective pixel pitch in mm.

    Edge cases:
        Raises ``ValueError`` if ``D`` is not positive.
    """

    if geom.D <= 0:
        raise ValueError("DetectorGeometry.D must be positive")
    return float(geom.pixel_pitch / (np.sqrt(6.0) * geom.D))


def sigma_theta_mc(
    events: pd.DataFrame, geom: DetectorGeometry, truth_window_rad: float = 0.01
) -> float:
    """Return empirical ``theta_rec`` spread for a narrow truth-angle window.

    Args:
        events: DataFrame from ``load_events``. ``theta_rec`` is in rad. If
            ``events.attrs['theta_true']`` exists, it is used as truth.
        geom: Detector geometry in mm. Present for API symmetry and future
            tolerance studies; it is not modified.
        truth_window_rad: Half-width around the median true theta in rad.

    Returns:
        Standard deviation of ``theta_rec`` in rad. Returns ``nan`` when no
        truth-level information is available.

    Edge cases:
        Current Meiga output may not include truth direction for every event.
        In that case this function logs the limitation and returns ``nan``.
    """

    del geom
    theta_true = events.attrs.get("theta_true")
    if theta_true is None:
        LOGGER.warning(
            "Truth-level theta is unavailable; sigma_theta_mc is a documented limitation"
        )
        return float("nan")
    theta_true = np.asarray(theta_true, dtype=float)
    if len(theta_true) != len(events) or len(events) == 0:
        return float("nan")
    center = float(np.median(theta_true))
    mask = np.abs(theta_true - center) <= truth_window_rad
    if np.count_nonzero(mask) < 2:
        return float("nan")
    return float(np.std(events.loc[mask, "theta_rec"].to_numpy(dtype=float), ddof=1))
