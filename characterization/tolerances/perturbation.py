"""Perturbed detector geometry for Phase-B tolerance studies.

Minimal reproducible example:

```python
from characterization.geometry import DetectorGeometry
from characterization.tolerances.perturbation import apply_perturbation

geom_p = apply_perturbation(DetectorGeometry(), {"dD": 0.5, "tilt_x": 0.005})
top_xyz, bot_xyz = geom_p.pixel_centers_perturbed()
```
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping

import numpy as np
from numpy.typing import NDArray
from scipy.spatial.transform import Rotation

from characterization.geometry import DetectorGeometry


@dataclass(frozen=True)
class PerturbedGeometry(DetectorGeometry):
    """Immutable geometry with first-order tolerance fields.

    Attributes:
        dD: Change in top-bottom module separation in mm. Positive values move
            the bottom module farther from the top module.
        tilt_x: Bottom-plane rotation around the common-frame x axis in rad.
        tilt_y: Bottom-plane rotation around the common-frame y axis in rad.
        dx: Bottom-module lateral shift along +x in mm.
        dy: Bottom-module lateral shift along +y in mm.

    Edge cases:
        The top module is never perturbed. Tilts are applied to bottom-module
        local coordinates around the bottom-module center, then the center is
        translated to ``(dx, dy, D + dD)``.
    """

    dD: float = 0.0
    tilt_x: float = 0.0
    tilt_y: float = 0.0
    dx: float = 0.0
    dy: float = 0.0

    @property
    def effective_D(self) -> float:
        """Return the perturbed module-center separation in mm."""

        return self.D + self.dD

    def pixel_centers_perturbed(
        self,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return perturbed top and bottom pixel centers in a common frame.

        Returns:
            ``(top_xyz, bottom_xyz)`` arrays, each of shape ``(n_pixels, 3)``.
            Coordinates are in mm.

        Edge cases:
            Raises ``ValueError`` if the effective separation is not positive.
            Rotations use full 3D rotation matrices; no small-angle
            approximation is introduced in the geometry itself.
        """

        return pixel_centers_perturbed(self)


def rotation_matrix(tilt_x: float, tilt_y: float) -> NDArray[np.float64]:
    """Return the bottom-plane rotation matrix for tilts in rad.

    Args:
        tilt_x: Right-hand-rule rotation around common-frame x in rad.
        tilt_y: Right-hand-rule rotation around common-frame y in rad.

    Returns:
        ``(3, 3)`` matrix using ``scipy.spatial.transform.Rotation``. The
        composition is y-after-x, matching the original Phase-B convention.

    Edge cases:
        No small-angle approximation is used.
    """

    return (
        Rotation.from_rotvec([0.0, tilt_y, 0.0])
        * Rotation.from_rotvec([tilt_x, 0.0, 0.0])
    ).as_matrix()


def apply_perturbation(
    geom: DetectorGeometry | PerturbedGeometry,
    perturbation: Mapping[str, float] | PerturbedGeometry,
) -> PerturbedGeometry:
    """Return a new immutable perturbed geometry.

    Args:
        geom: Nominal or already perturbed geometry. Lengths are mm.
        perturbation: Mapping with any of ``dD``, ``tilt_x``, ``tilt_y``,
            ``dx`` or ``dy``; angles are rad.

    Returns:
        ``PerturbedGeometry`` with the same nominal detector dimensions and the
        requested tolerance fields.

    Sign convention:
        ``dD > 0`` increases module separation. ``dx > 0`` and ``dy > 0`` shift
        the bottom module in +x/+y. Positive tilts follow the right-hand rule
        around the common-frame x or y axis.

    Edge cases:
        Unknown mapping keys are ignored to keep CLI-generated dictionaries
        convenient. Existing perturbation values are replaced, not added.
    """

    if isinstance(perturbation, PerturbedGeometry):
        values = {
            "dD": perturbation.dD,
            "tilt_x": perturbation.tilt_x,
            "tilt_y": perturbation.tilt_y,
            "dx": perturbation.dx,
            "dy": perturbation.dy,
        }
    else:
        values = {key: float(value) for key, value in perturbation.items()}
    if isinstance(geom, PerturbedGeometry):
        base = replace(geom, dD=0.0, tilt_x=0.0, tilt_y=0.0, dx=0.0, dy=0.0)
    else:
        base = geom
    return PerturbedGeometry(
        D=base.D,
        bar_length=base.bar_length,
        bar_width=base.bar_width,
        bar_thickness=base.bar_thickness,
        n_channels=base.n_channels,
        pixel_pitch=base.pixel_pitch,
        sub_plane_offset=base.sub_plane_offset,
        module_size=base.module_size,
        dD=values.get("dD", 0.0),
        tilt_x=values.get("tilt_x", 0.0),
        tilt_y=values.get("tilt_y", 0.0),
        dx=values.get("dx", 0.0),
        dy=values.get("dy", 0.0),
    )


def pixel_centers_perturbed(
    geom_p: PerturbedGeometry,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return top and bottom pixel centers after perturbation.

    Args:
        geom_p: Perturbed geometry. Lengths are mm, angles are rad.

    Returns:
        Tuple ``(top_xyz, bottom_xyz)`` in mm.

    Edge cases:
        The top plane lies at ``z = 0``. The bottom plane center lies at
        ``z = D + dD``. A non-positive effective distance raises ``ValueError``.
    """

    if geom_p.effective_D <= 0:
        raise ValueError("Effective module separation D + dD must be positive")
    xy = geom_p.pixel_centers()
    zeros = np.zeros((geom_p.n_pixels, 1), dtype=float)
    top = np.column_stack([xy, zeros])
    bottom_local = np.column_stack([xy, zeros])
    rotation = rotation_matrix(geom_p.tilt_x, geom_p.tilt_y)
    bottom = bottom_local @ rotation.T
    bottom += np.asarray([geom_p.dx, geom_p.dy, geom_p.effective_D], dtype=float)
    return top, bottom
