"""Detector pixel-pair trajectories for forward tomography.

Minimal example:

```python
from characterization.geometry import DetectorGeometry
from characterization.forward.trajectories import pair_trajectories

traj = pair_trajectories(DetectorGeometry())
```
"""

from __future__ import annotations

import numpy as np

from characterization.geometry import DetectorGeometry
from characterization.forward.frames import detector_bottom_z
from characterization.solid_angle_general import Rectangle3D, omega_general


TRAJECTORY_DTYPE = np.dtype(
    [
        ("i_top", "i4"),
        ("j_top", "i4"),
        ("i_bot", "i4"),
        ("j_bot", "i4"),
        ("origin", "f8", (3,)),
        ("direction", "f8", (3,)),
        ("theta", "f8"),
        ("phi", "f8"),
        ("omega", "f8"),
        ("acceptance_T", "f8"),
    ]
)


def fast_pair_kinematics(
    geom: DetectorGeometry,
    n_quad: int = 6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return public per-pixel-pair origins, directions, angles, and omega.

    The row order is ``top_pixel * geom.n_pixels + bottom_pixel`` and matches
    ``pair_trajectories``. Lengths are in mm, angles are in rad, and omega is
    in sr.
    """

    traj = pair_trajectories(geom, n_quad=n_quad)
    return (
        np.asarray(traj["origin"], dtype=float),
        np.asarray(traj["direction"], dtype=float),
        np.asarray(traj["theta"], dtype=float),
        np.asarray(traj["phi"], dtype=float),
        np.asarray(traj["omega"], dtype=float),
    )


def _fast_pair_kinematics(
    geom: DetectorGeometry,
    n_quad: int = 6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Deprecated alias for ``fast_pair_kinematics``."""

    import warnings

    warnings.warn(
        "_fast_pair_kinematics is deprecated; use fast_pair_kinematics",
        DeprecationWarning,
        stacklevel=2,
    )
    return fast_pair_kinematics(geom, n_quad=n_quad)


def pair_trajectories(
    geom: DetectorGeometry,
    oversample_n: int = 1,
    n_quad: int = 6,
) -> np.ndarray:
    """Return representative trajectories for all top-bottom pixel pairs.

    Args:
        geom: Detector geometry in mm.
        oversample_n: Number of pencil beams per pixel pair. Phase C-1 returns
            center rays for ``oversample_n = 1``; oversampling is implemented in
            ``build_F`` where ray sampling is needed.
        n_quad: Gauss-Legendre order for the general-orientation solid angle.

    Returns:
        Structured array with one row per pixel pair. ``origin`` is the center
        of the bottom pixel in mm and ``direction`` points from bottom to top
        and then upward into geology. ``theta`` and ``phi`` are rad.

    Edge cases:
        ``acceptance_T`` is set to 1.0 in this deterministic geometric forward
        model. Phase A angular bin acceptances can be multiplied externally if
        desired.
    """

    if oversample_n < 1:
        raise ValueError("oversample_n must be >= 1")
    centers = geom.pixel_centers()
    n_pairs = geom.n_pixels * geom.n_pixels
    result = np.zeros(n_pairs, dtype=TRAJECTORY_DTYPE)
    row = 0
    z_bot = detector_bottom_z(geom.D)
    for top in range(geom.n_pixels):
        top_xy = centers[top]
        for bot in range(geom.n_pixels):
            bot_xy = centers[bot]
            origin = np.asarray([bot_xy[0], bot_xy[1], z_bot], dtype=float)
            top_point = np.asarray([top_xy[0], top_xy[1], 0.0], dtype=float)
            direction = top_point - origin
            direction = direction / np.linalg.norm(direction)
            theta = np.arccos(np.clip(direction[2], -1.0, 1.0))
            phi = np.arctan2(direction[1], direction[0])
            half = 0.5 * geom.bar_width
            rect_top = Rectangle3D(
                center=top_point,
                normal=(0.0, 0.0, -1.0),
                u_axis=(1.0, 0.0, 0.0),
                v_axis=(0.0, 1.0, 0.0),
                half_extent_u=half,
                half_extent_v=half,
            )
            rect_bot = Rectangle3D(
                center=origin,
                normal=(0.0, 0.0, 1.0),
                u_axis=(1.0, 0.0, 0.0),
                v_axis=(0.0, 1.0, 0.0),
                half_extent_u=half,
                half_extent_v=half,
            )
            result[row] = (
                top % geom.n_channels,
                top // geom.n_channels,
                bot % geom.n_channels,
                bot // geom.n_channels,
                origin,
                direction,
                theta,
                phi,
                omega_general(
                    rect_top, rect_bot, method="gauss_legendre", n_quad=n_quad
                ),
                1.0,
            )
            row += 1
    return result
