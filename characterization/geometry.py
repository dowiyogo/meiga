"""Detector geometry helpers for the NA64-mini-style hodoscope.

Minimal reproducible example:

```python
from characterization.geometry import DetectorGeometry

geom = DetectorGeometry(D=98.4, bar_length=75.0, bar_width=3.0)
top_xy = geom.pixel_center(108)
bot_xy = geom.pixel_center(101)
```
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class DetectorGeometry:
    """Geometry constants for the two-module hodoscope.

    Attributes:
        D: Center-to-center distance between top and bottom modules in mm.
        bar_length: Active scintillator bar length in mm.
        bar_width: Active scintillator bar width in mm.
        bar_thickness: Active scintillator bar thickness in mm.
        n_channels: Number of reconstructed X or Y channels per module.
        pixel_pitch: Effective pitch between sorted pixel centers in mm.
        sub_plane_offset: Offset between interleaved sub-planes in mm.
        module_size: Active module envelope as ``(x_size, y_size)`` in mm.
    """

    D: float = 98.4
    bar_length: float = 75.0
    bar_width: float = 3.0
    bar_thickness: float = 1.0
    n_channels: int = 16
    pixel_pitch: float = 2.0
    sub_plane_offset: float = 2.0
    module_size: tuple[float, float] = (75.0, 75.0)

    @property
    def n_pixels(self) -> int:
        """Return the number of X-Y pixels in one module."""

        return self.n_channels * self.n_channels

    @property
    def pixel_area(self) -> float:
        """Return the effective pixel area in mm^2."""

        return self.bar_width * self.bar_width

    def channel_centers(self) -> NDArray[np.float64]:
        """Return reconstructed channel centers in mm.

        The ordering follows the simulator channel IDs, including the
        interleaved A/B sub-plane convention. Edge cases: if ``n_channels`` is
        not divisible by two, the method still returns evenly split channels
        using integer division, but the NA64-mini default is 16.
        """

        if self.n_channels == 1:
            return np.asarray([0.0], dtype=float)
        half = self.n_channels // 2
        channels = np.arange(self.n_channels, dtype=float)
        subplanes = np.floor_divide(channels.astype(int), half)
        local = np.mod(channels.astype(int), half)
        nominal_pitch = 2.0 * self.pixel_pitch
        origin = 0.5 * ((half - 1) * nominal_pitch + self.sub_plane_offset)
        return local * nominal_pitch + subplanes * self.sub_plane_offset - origin

    def pixel_centers(self) -> NDArray[np.float64]:
        """Return all pixel centers as an array of shape ``(n_pixels, 2)``.

        Units are mm. Pixel IDs follow ``pixel_id = x + n_channels * y``.
        """

        centers = self.channel_centers()
        x_index, y_index = np.meshgrid(
            np.arange(self.n_channels), np.arange(self.n_channels), indexing="xy"
        )
        xy = np.column_stack(
            [centers[x_index.ravel()], centers[y_index.ravel()]]
        ).astype(float)
        return xy

    def pixel_center(self, pixel_id: int) -> tuple[float, float]:
        """Return one pixel center ``(x, y)`` in mm.

        Args:
            pixel_id: Integer pixel ID in ``[0, n_pixels)``.

        Raises:
            ValueError: If ``pixel_id`` is outside the module.
        """

        if pixel_id < 0 or pixel_id >= self.n_pixels:
            raise ValueError(f"pixel_id={pixel_id} outside [0, {self.n_pixels})")
        centers = self.channel_centers()
        x_index = pixel_id % self.n_channels
        y_index = pixel_id // self.n_channels
        return float(centers[x_index]), float(centers[y_index])

    def rectangle_for_pixel(self, pixel_id: int) -> tuple[float, float, float, float]:
        """Return ``(xmin, xmax, ymin, ymax)`` for a pixel rectangle in mm.

        The effective finite rectangle is the overlap of two bars, modeled as a
        square of side ``bar_width`` centered at the reconstructed pixel center.
        """

        x, y = self.pixel_center(pixel_id)
        half = 0.5 * self.bar_width
        return x - half, x + half, y - half, y + half
