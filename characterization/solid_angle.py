"""Solid-angle estimators for pairs of finite hodoscope pixels.

Minimal reproducible example:

```python
from characterization.geometry import DetectorGeometry
from characterization.solid_angle import omega_analytical

geom = DetectorGeometry(D=1.0, bar_width=1.0, pixel_pitch=1.0, n_channels=1)
omega_sr = omega_analytical(0, 0, geom)
```

``omega_analytical`` returns the top-pixel-area-averaged mutual solid angle in
sr. Equivalently, ``pixel_area * omega`` is the usual two-aperture geometrical
factor in ``mm^2 sr``.
"""

from __future__ import annotations

import logging
import math

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.integrate import dblquad

from characterization.geometry import DetectorGeometry

LOGGER = logging.getLogger(__name__)


def _overlap_length(
    delta: float, top_min: float, top_max: float, bot_min: float, bot_max: float
) -> float:
    shifted_min = bot_min - delta
    shifted_max = bot_max - delta
    return max(0.0, min(top_max, shifted_max) - max(top_min, shifted_min))


def omega_analytical(pixel_top: int, pixel_bot: int, geom: DetectorGeometry) -> float:
    """Compute the finite-rectangle mutual solid angle in sr.

    Args:
        pixel_top: Top-module pixel ID.
        pixel_bot: Bottom-module pixel ID.
        geom: Detector geometry. Lengths are mm and ``D`` is the inter-module
            distance in mm.

    Returns:
        Area-averaged solid angle in sr. The dimensional kernel is
        ``D^2 / r^4``; after integrating over both rectangles and dividing by
        the top-pixel area, the result has units sr.

    Edge cases:
        Raises ``ValueError`` for non-positive ``D`` or invalid pixel IDs. The
        implementation evaluates the Khadjavi/Gardner-Verghese convolution form
        over rectangle-edge difference limits with adaptive quadrature.
    """

    if geom.D <= 0:
        raise ValueError("DetectorGeometry.D must be positive")
    top = geom.rectangle_for_pixel(pixel_top)
    bot = geom.rectangle_for_pixel(pixel_bot)
    top_area = (top[1] - top[0]) * (top[3] - top[2])
    if top_area <= 0:
        raise ValueError("Top pixel area must be positive")

    ux_min = bot[0] - top[1]
    ux_max = bot[1] - top[0]
    uy_min = bot[2] - top[3]
    uy_max = bot[3] - top[2]

    def integrand(v_delta: float, u_delta: float) -> float:
        wx = _overlap_length(u_delta, top[0], top[1], bot[0], bot[1])
        wy = _overlap_length(v_delta, top[2], top[3], bot[2], bot[3])
        r2 = geom.D * geom.D + u_delta * u_delta + v_delta * v_delta
        return wx * wy * geom.D * geom.D / (r2 * r2)

    integral, _ = dblquad(
        integrand,
        ux_min,
        ux_max,
        lambda _u: uy_min,
        lambda _u: uy_max,
        epsabs=1e-10,
        epsrel=1e-8,
    )
    return float(integral / top_area)


def _pixel_ids(
    events: pd.DataFrame, geom: DetectorGeometry
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    top = events["pixel_top_x"].to_numpy(dtype=np.int64) + geom.n_channels * events[
        "pixel_top_y"
    ].to_numpy(dtype=np.int64)
    bot = events["pixel_bot_x"].to_numpy(dtype=np.int64) + geom.n_channels * events[
        "pixel_bot_y"
    ].to_numpy(dtype=np.int64)
    return top, bot


def _cos2_pdf(theta: NDArray[np.float64]) -> NDArray[np.float64]:
    return (3.0 / (2.0 * np.pi)) * np.cos(theta) ** 2


def omega_mc(events: pd.DataFrame, geom: DetectorGeometry) -> NDArray[np.float64]:
    """Estimate solid angle per pixel pair from accepted events.

    Args:
        events: Four-plane coincidence events. Angles are rad. For synthetic
            validation, set ``events.attrs['sampling_pdf']`` to
            ``'cosine_hemisphere'`` and ``events.attrs['n_generated_per_top_pixel']``
            to the number of generated rays per top pixel.
        geom: Detector geometry in mm.

    Returns:
        ``(n_pix_top, n_pix_bot)`` array in sr.

    Edge cases:
        For real EcoMug-like samples without per-top-pixel incident counters,
        the denominator is approximated as total incident events divided evenly
        over top pixels. This is suitable for Phase-A diagnostics, not a final
        flux unfolding.
    """

    omega = np.zeros((geom.n_pixels, geom.n_pixels), dtype=float)
    if events.empty:
        return omega

    top, bot = _pixel_ids(events, geom)
    valid = (top >= 0) & (top < geom.n_pixels) & (bot >= 0) & (bot < geom.n_pixels)
    counts = np.zeros_like(omega)
    np.add.at(counts, (top[valid], bot[valid]), 1.0)

    generated = float(
        events.attrs.get(
            "n_generated_per_top_pixel",
            events.attrs.get("n_incident_top", len(events)) / geom.n_pixels,
        )
    )
    if generated <= 0:
        return omega

    if events.attrs.get("sampling_pdf") == "uniform_hemisphere":
        return (2.0 * np.pi) * counts / generated
    if events.attrs.get("sampling_pdf") == "cosine_hemisphere":
        return np.pi * counts / generated

    centers = geom.pixel_centers()
    for top_id in np.unique(top[valid]):
        pair_mask = valid & (top == top_id)
        bot_ids = np.unique(bot[pair_mask])
        top_xy = centers[top_id]
        for bot_id in bot_ids:
            bot_xy = centers[bot_id]
            theta = math.atan2(float(np.linalg.norm(bot_xy - top_xy)), geom.D)
            pdf = float(_cos2_pdf(np.asarray([theta]))[0])
            if pdf > 0:
                omega[top_id, bot_id] = counts[top_id, bot_id] / (generated * pdf)
    return omega


def synthetic_pair_events(
    pixel_top: int,
    pixel_bot: int,
    geom: DetectorGeometry,
    n_rays: int,
    seed: int = 12345,
) -> pd.DataFrame:
    """Generate accepted events for one pixel pair using cosine-hemisphere rays.

    Args:
        pixel_top: Top pixel ID.
        pixel_bot: Bottom pixel ID.
        geom: Detector geometry in mm.
        n_rays: Number of generated rays from the top pixel.
        seed: Random seed.

    Returns:
        DataFrame compatible with ``omega_mc``. Units are mm and rad.

    Edge cases:
        This is a geometry-only generator; energy and timing columns are filled
        with deterministic placeholders.
    """

    rng = np.random.default_rng(seed)
    top_rect = geom.rectangle_for_pixel(pixel_top)
    bot_rect = geom.rectangle_for_pixel(pixel_bot)
    x0 = rng.uniform(top_rect[0], top_rect[1], n_rays)
    y0 = rng.uniform(top_rect[2], top_rect[3], n_rays)
    cos_theta = np.sqrt(rng.uniform(0.0, 1.0, n_rays))
    sin_theta = np.sqrt(1.0 - cos_theta * cos_theta)
    phi = rng.uniform(-np.pi, np.pi, n_rays)
    dx = geom.D * sin_theta * np.cos(phi) / np.maximum(cos_theta, 1e-15)
    dy = geom.D * sin_theta * np.sin(phi) / np.maximum(cos_theta, 1e-15)
    xb = x0 + dx
    yb = y0 + dy
    hit = (
        (xb >= bot_rect[0])
        & (xb <= bot_rect[1])
        & (yb >= bot_rect[2])
        & (yb <= bot_rect[3])
    )
    top_x = pixel_top % geom.n_channels
    top_y = pixel_top // geom.n_channels
    bot_x = pixel_bot % geom.n_channels
    bot_y = pixel_bot // geom.n_channels
    rows = {
        "event_id": np.arange(int(np.count_nonzero(hit)), dtype=np.int64),
        "pixel_top_x": np.full(np.count_nonzero(hit), top_x, dtype=np.int64),
        "pixel_top_y": np.full(np.count_nonzero(hit), top_y, dtype=np.int64),
        "pixel_bot_x": np.full(np.count_nonzero(hit), bot_x, dtype=np.int64),
        "pixel_bot_y": np.full(np.count_nonzero(hit), bot_y, dtype=np.int64),
        "theta_rec": np.full(np.count_nonzero(hit), 0.0, dtype=float),
        "phi_rec": np.full(np.count_nonzero(hit), 0.0, dtype=float),
        "energy_top": np.ones(np.count_nonzero(hit), dtype=float),
        "energy_bot": np.ones(np.count_nonzero(hit), dtype=float),
        "t_top": np.zeros(np.count_nonzero(hit), dtype=float),
        "t_bot": np.zeros(np.count_nonzero(hit), dtype=float),
    }
    frame = pd.DataFrame(rows)
    frame.attrs["sampling_pdf"] = "cosine_hemisphere"
    frame.attrs["n_generated_per_top_pixel"] = n_rays
    frame.attrs["n_incident_top"] = n_rays
    return frame


def compare_omega(geom: DetectorGeometry, n_samples: int) -> pd.DataFrame:
    """Compare analytical and MC solid angle for representative pixel pairs.

    Args:
        geom: Detector geometry in mm.
        n_samples: Number of MC rays per representative pair.

    Returns:
        DataFrame with top/bottom pixel IDs, analytical ``omega`` in sr, MC
        ``omega`` in sr, and relative absolute error.

    Edge cases:
        To keep Phase-A runs lightweight, the table covers a representative set
        of central, edge, and diagonal pairs rather than all 65536 pairs.
    """

    representatives = sorted(
        {
            geom.n_pixels // 2,
            geom.n_channels * (geom.n_channels // 2) + geom.n_channels // 2,
            0,
            geom.n_channels - 1,
            geom.n_pixels - geom.n_channels,
            geom.n_pixels - 1,
        }
    )
    rows: list[dict[str, float | int]] = []
    for i, top_id in enumerate(representatives):
        for bot_id in representatives:
            analytical = omega_analytical(top_id, bot_id, geom)
            synthetic = synthetic_pair_events(
                top_id, bot_id, geom, n_samples, seed=1000 + i * 100 + bot_id
            )
            mc_value = omega_mc(synthetic, geom)[top_id, bot_id]
            rel = abs(analytical - mc_value) / analytical if analytical > 0 else np.nan
            rows.append(
                {
                    "pixel_top": top_id,
                    "pixel_bot": bot_id,
                    "omega_analytical_sr": analytical,
                    "omega_mc_sr": mc_value,
                    "relative_error": rel,
                }
            )
    table = pd.DataFrame(rows)
    table.attrs["mean_relative_error"] = float(table["relative_error"].mean())
    table.attrs["max_relative_error"] = float(table["relative_error"].max())
    return table
