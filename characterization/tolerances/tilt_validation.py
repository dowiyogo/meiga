"""General-kernel validation for tilted hodoscope pixels.

Minimal reproducible example:

```python
import numpy as np
from characterization.geometry import DetectorGeometry
from characterization.tolerances.tilt_validation import compare_kernels

geom = DetectorGeometry()
table = compare_kernels(geom, np.linspace(-20, 20, 9), {"center": (136, 136)})
```
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from characterization.geometry import DetectorGeometry
from characterization.solid_angle import omega_analytical
from characterization.solid_angle_general import Rectangle3D, omega_general
from characterization.tolerances.finite_diff import omega_perturbed_pair
from characterization.tolerances.perturbation import (
    PerturbedGeometry,
    apply_perturbation,
    rotation_matrix,
)


def build_rectangles(
    geom_perturbed: PerturbedGeometry,
    i_top: int,
    j_top: int,
    i_bot: int,
    j_bot: int,
) -> tuple[Rectangle3D, Rectangle3D]:
    """Construct top and bottom ``Rectangle3D`` pixels from perturbed geometry.

    Args:
        geom_perturbed: Perturbed geometry. Lengths are mm and tilts are rad.
        i_top: Top pixel x index.
        j_top: Top pixel y index.
        i_bot: Bottom pixel x index.
        j_bot: Bottom pixel y index.

    Returns:
        ``(rect_top, rect_bot)`` with normals facing each other. Units are mm.

    Edge cases:
        Pixel indices outside the module raise ``ValueError`` through the
        nominal geometry pixel lookup. The same ``rotation_matrix`` utility as
        ``perturbation.py`` is used to avoid sign-convention drift.
    """

    top_id = i_top + geom_perturbed.n_channels * j_top
    bot_id = i_bot + geom_perturbed.n_channels * j_bot
    top_x, top_y = geom_perturbed.pixel_center(top_id)
    bot_x, bot_y = geom_perturbed.pixel_center(bot_id)
    rotation = rotation_matrix(geom_perturbed.tilt_x, geom_perturbed.tilt_y)
    bottom_center = rotation @ np.asarray([bot_x, bot_y, 0.0], dtype=float)
    bottom_center += np.asarray(
        [geom_perturbed.dx, geom_perturbed.dy, geom_perturbed.effective_D],
        dtype=float,
    )
    half = 0.5 * geom_perturbed.bar_width
    rect_top = Rectangle3D(
        center=np.asarray([top_x, top_y, 0.0], dtype=float),
        normal=np.asarray([0.0, 0.0, 1.0], dtype=float),
        u_axis=np.asarray([1.0, 0.0, 0.0], dtype=float),
        v_axis=np.asarray([0.0, 1.0, 0.0], dtype=float),
        half_extent_u=half,
        half_extent_v=half,
    )
    rect_bot = Rectangle3D(
        center=bottom_center,
        normal=rotation @ np.asarray([0.0, 0.0, -1.0], dtype=float),
        u_axis=rotation @ np.asarray([1.0, 0.0, 0.0], dtype=float),
        v_axis=rotation @ np.asarray([0.0, 1.0, 0.0], dtype=float),
        half_extent_u=half,
        half_extent_v=half,
    )
    return rect_top, rect_bot


def _ids_to_indices(
    geom: DetectorGeometry, pixel_pair: tuple[int, int]
) -> tuple[int, int, int, int]:
    top, bot = pixel_pair
    return (
        top % geom.n_channels,
        top // geom.n_channels,
        bot % geom.n_channels,
        bot // geom.n_channels,
    )


def omega_general_for_pair(
    geom_perturbed: PerturbedGeometry,
    pixel_pair: tuple[int, int],
    n_quad: int = 12,
) -> float:
    """Evaluate the general kernel for one pixel pair.

    Args:
        geom_perturbed: Perturbed geometry in mm/rad.
        pixel_pair: ``(top_pixel_id, bottom_pixel_id)``.
        n_quad: Gauss-Legendre order.

    Returns:
        Effective solid angle in sr.

    Edge cases:
        Propagates clear geometry errors from ``Rectangle3D`` and
        ``omega_general``.
    """

    rect_top, rect_bot = build_rectangles(
        geom_perturbed, *_ids_to_indices(geom_perturbed, pixel_pair)
    )
    return omega_general(rect_top, rect_bot, method="gauss_legendre", n_quad=n_quad)


def compare_kernels(
    geom: DetectorGeometry,
    alpha_grid: np.ndarray,
    pixel_pairs: Mapping[str, tuple[int, int]],
    axis: str = "x",
    n_quad: int = 12,
    n_samples: int = 200_000,
) -> pd.DataFrame:
    """Compare parallel, Gauss-Legendre general, and MC general kernels.

    Args:
        geom: Nominal geometry in mm.
        alpha_grid: Tilt values in mrad.
        pixel_pairs: Mapping from label to ``(top_pixel_id, bottom_pixel_id)``.
        axis: ``"x"`` or ``"y"`` tilt axis.
        n_quad: Gauss-Legendre order.
        n_samples: Monte Carlo samples for the independent general-kernel check.

    Returns:
        DataFrame with columns ``pixel_pair_label``, ``alpha_mrad``,
        ``omega_parallel``, ``omega_gl``, ``omega_mc``,
        ``rel_diff_gl_vs_parallel``, and ``rel_diff_mc_vs_gl``.

    Edge cases:
        The Monte Carlo column is intentionally stochastic and uses deterministic
        seeds. It is for bias scale, not for high-precision derivatives.
    """

    if axis not in {"x", "y"}:
        raise ValueError("axis must be 'x' or 'y'")
    rows: list[dict[str, float | str]] = []
    for label, pair in pixel_pairs.items():
        for index, alpha_mrad in enumerate(alpha_grid):
            perturbation = {
                "tilt_x" if axis == "x" else "tilt_y": float(alpha_mrad) * 1e-3
            }
            geom_p = apply_perturbation(geom, perturbation)
            rect_top, rect_bot = build_rectangles(geom_p, *_ids_to_indices(geom, pair))
            omega_gl = omega_general(
                rect_top, rect_bot, method="gauss_legendre", n_quad=n_quad
            )
            omega_mc = omega_general(
                rect_top,
                rect_bot,
                method="monte_carlo",
                n_samples=n_samples,
                seed=17_000 + 100 * index + pair[0] + pair[1],
            )
            omega_parallel = omega_perturbed_pair(pair[0], pair[1], geom_p)
            rows.append(
                {
                    "pixel_pair_label": label,
                    "alpha_mrad": float(alpha_mrad),
                    "omega_parallel": omega_parallel,
                    "omega_gl": omega_gl,
                    "omega_mc": omega_mc,
                    "rel_diff_gl_vs_parallel": (omega_gl - omega_parallel) / omega_gl,
                    "rel_diff_mc_vs_gl": (omega_mc - omega_gl) / omega_gl,
                }
            )
    return pd.DataFrame(rows)


def dOmega_dalpha_general(
    geom: DetectorGeometry,
    pixel_pair: tuple[int, int],
    axis: str,
    h: float = 0.1,
    n_quad: int = 12,
) -> float:
    """Return central finite-difference tilt derivative using ``omega_general``.

    Args:
        geom: Nominal geometry in mm.
        pixel_pair: ``(top_pixel_id, bottom_pixel_id)``.
        axis: ``"x"`` or ``"y"``.
        h: Step size in mrad.
        n_quad: Gauss-Legendre order.

    Returns:
        ``dOmega/dalpha`` in sr/mrad.

    Edge cases:
        Raises ``ValueError`` for non-positive ``h`` or invalid axis.
    """

    if h <= 0:
        raise ValueError("h must be positive")
    if axis not in {"x", "y"}:
        raise ValueError("axis must be 'x' or 'y'")
    field = "tilt_x" if axis == "x" else "tilt_y"
    plus = apply_perturbation(geom, {field: h * 1e-3})
    minus = apply_perturbation(geom, {field: -h * 1e-3})
    omega_plus = omega_general_for_pair(plus, pixel_pair, n_quad=n_quad)
    omega_minus = omega_general_for_pair(minus, pixel_pair, n_quad=n_quad)
    return float((omega_plus - omega_minus) / (2.0 * h))


def alpha_threshold_for_bias(
    table: pd.DataFrame, threshold: float = 0.01
) -> pd.DataFrame:
    """Find where the parallel-kernel bias exceeds a threshold.

    Args:
        table: Output from ``compare_kernels``.
        threshold: Absolute relative bias threshold.

    Returns:
        DataFrame with one row per pixel-pair label and threshold in mrad. If no
        grid point exceeds the threshold, ``alpha_mrad_at_threshold`` is ``nan``.
    """

    rows = []
    for label, group in table.groupby("pixel_pair_label"):
        ordered = group.assign(abs_alpha=group["alpha_mrad"].abs()).sort_values(
            "abs_alpha"
        )
        exceeded = ordered[ordered["rel_diff_gl_vs_parallel"].abs() > threshold]
        rows.append(
            {
                "pixel_pair_label": label,
                "alpha_mrad_at_1pct_bias": (
                    float(exceeded["abs_alpha"].iloc[0])
                    if not exceeded.empty
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)
