"""Frame conventions for the forward muography model.

Minimal example:

```python
from characterization.forward.frames import detector_bottom_z

z_bottom = detector_bottom_z(98.4)
```

World frame:
    Coordinates are in mm. ``z`` increases upward. The detector top module is
    centered at ``z = 0``. The bottom module is centered at ``z = -D``. Geology
    occupies ``z > 0``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

MM_PER_CM = 10.0
CM_PER_MM = 0.1
MM_PER_M = 1000.0

TOP_MODULE_Z_MM = 0.0


def detector_bottom_z(D: float) -> float:
    """Return bottom-module z coordinate in mm.

    Args:
        D: Inter-plane distance in mm.

    Returns:
        ``-D`` in mm.
    """

    return -float(D)


def mm_to_cm(length_mm: ArrayLike) -> NDArray[np.float64]:
    """Convert lengths from mm to cm."""

    return np.asarray(length_mm, dtype=float) * CM_PER_MM


def cm_to_mm(length_cm: ArrayLike) -> NDArray[np.float64]:
    """Convert lengths from cm to mm."""

    return np.asarray(length_cm, dtype=float) * MM_PER_CM


def world_to_geology(point_mm: ArrayLike) -> NDArray[np.float64]:
    """Return geology-frame coordinates from world coordinates.

    The Phase C-1 geology frame is currently identical to the world frame. This
    helper exists to prevent silent convention drift in later phases.
    """

    return np.asarray(point_mm, dtype=float)


def geology_to_world(point_mm: ArrayLike) -> NDArray[np.float64]:
    """Return world-frame coordinates from geology coordinates in mm."""

    return np.asarray(point_mm, dtype=float)
