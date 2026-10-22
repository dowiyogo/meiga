"""Finite-difference tolerance derivatives for Phase B.

Minimal reproducible example:

```python
from characterization.geometry import DetectorGeometry
from characterization.tolerances.finite_diff import fd_dOmega

geom = DetectorGeometry(n_channels=2)
d_omega_dD = fd_dOmega(geom, "D", 1e-3)
```
"""

from __future__ import annotations

import hashlib
import logging
import pickle
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.integrate import dblquad

from characterization.acceptance import pixel_acceptance
from characterization.geometry import DetectorGeometry
from characterization.solid_angle import _overlap_length
from characterization.tolerances.perturbation import (
    PerturbedGeometry,
    apply_perturbation,
    pixel_centers_perturbed,
)

LOGGER = logging.getLogger(__name__)

PixelPair = tuple[int, int]


def _perturb_dict(field: str, value: float) -> dict[str, float]:
    if field == "D":
        return {"dD": value}
    if field in {"tilt_x", "tilt_y", "dx", "dy"}:
        return {field: value}
    raise ValueError("perturbation_field must be one of D, tilt_x, tilt_y, dx, dy")


def _omega_from_rectangles(
    top_rect: tuple[float, float, float, float],
    bot_rect: tuple[float, float, float, float],
    distance: float,
) -> float:
    top_area = (top_rect[1] - top_rect[0]) * (top_rect[3] - top_rect[2])
    ux_min = bot_rect[0] - top_rect[1]
    ux_max = bot_rect[1] - top_rect[0]
    uy_min = bot_rect[2] - top_rect[3]
    uy_max = bot_rect[3] - top_rect[2]

    def integrand(v_delta: float, u_delta: float) -> float:
        wx = _overlap_length(
            u_delta, top_rect[0], top_rect[1], bot_rect[0], bot_rect[1]
        )
        wy = _overlap_length(
            v_delta, top_rect[2], top_rect[3], bot_rect[2], bot_rect[3]
        )
        r2 = distance * distance + u_delta * u_delta + v_delta * v_delta
        return wx * wy * distance * distance / (r2 * r2)

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


def omega_perturbed_pair(
    pixel_top: int, pixel_bot: int, geom_p: PerturbedGeometry
) -> float:
    """Evaluate Phase-A ``omega_analytical`` for a perturbed pixel pair.

    Args:
        pixel_top: Top pixel ID.
        pixel_bot: Bottom pixel ID.
        geom_p: Perturbed geometry with lengths in mm and angles in rad.

    Returns:
        Solid angle in sr.

    Edge cases:
        Tilts are mapped to the exact 3D perturbed pixel-center separation, then
        evaluated with the Phase-A parallel-rectangle kernel. This keeps the
        finite-difference path compatible with Phase A while documenting that
        finite rectangle foreshortening from tilt is not yet modeled.
    """

    top_xyz, bot_xyz = pixel_centers_perturbed(geom_p)
    top_center = top_xyz[pixel_top]
    bot_center = bot_xyz[pixel_bot]
    delta = bot_center - top_center
    half = 0.5 * geom_p.bar_width
    top_rect = (-half, half, -half, half)
    bot_rect = (
        float(delta[0] - half),
        float(delta[0] + half),
        float(delta[1] - half),
        float(delta[1] + half),
    )
    return _omega_from_rectangles(top_rect, bot_rect, abs(float(delta[2])))


def _all_pairs(geom: DetectorGeometry) -> list[PixelPair]:
    return [(top, bot) for top in range(geom.n_pixels) for bot in range(geom.n_pixels)]


def fd_dOmega(
    geom: DetectorGeometry,
    perturbation_field: str,
    magnitude: float,
    scheme: str = "central",
    pixel_pairs: Iterable[PixelPair] | None = None,
) -> NDArray[np.float64]:
    """Return finite-difference ``dOmega/dxi`` per pixel pair.

    Args:
        geom: Nominal detector geometry in mm.
        perturbation_field: One of ``"D"``, ``"tilt_x"``, ``"tilt_y"``,
            ``"dx"``, ``"dy"``. Tilt magnitudes are rad.
        magnitude: Step size in the physical unit of the selected field.
        scheme: ``"central"`` or ``"forward"``. Default is central.
        pixel_pairs: Optional iterable of ``(top, bottom)`` pairs to evaluate.

    Returns:
        ``(n_pixels, n_pixels)`` derivative array. Unrequested pairs are zero
        when ``pixel_pairs`` is supplied. Units are sr/mm or sr/rad.

    Edge cases:
        The perturbed geometry uses exact 3D rotations, but the solid-angle
        kernel remains the Phase-A parallel-rectangle approximation.
    """

    if magnitude <= 0:
        raise ValueError("magnitude must be positive")
    if scheme not in {"forward", "central"}:
        raise ValueError("scheme must be 'forward' or 'central'")
    pairs = list(pixel_pairs) if pixel_pairs is not None else _all_pairs(geom)
    derivative = np.zeros((geom.n_pixels, geom.n_pixels), dtype=float)

    for top, bot in pairs:
        if scheme == "central":
            plus = apply_perturbation(
                geom, _perturb_dict(perturbation_field, magnitude)
            )
            minus = apply_perturbation(
                geom, _perturb_dict(perturbation_field, -magnitude)
            )
            derivative[top, bot] = (
                omega_perturbed_pair(top, bot, plus)
                - omega_perturbed_pair(top, bot, minus)
            ) / (2.0 * magnitude)
        else:
            plus = apply_perturbation(
                geom, _perturb_dict(perturbation_field, magnitude)
            )
            nominal = apply_perturbation(geom, {})
            derivative[top, bot] = (
                omega_perturbed_pair(top, bot, plus)
                - omega_perturbed_pair(top, bot, nominal)
            ) / magnitude
    return derivative


def _geometry_hash(geom: PerturbedGeometry) -> str:
    payload = pickle.dumps(geom)
    return hashlib.sha256(payload).hexdigest()[:16]


def _events_cached(
    events_func: Callable[[PerturbedGeometry], pd.DataFrame],
    geom_p: PerturbedGeometry,
    cache_dir: Path,
) -> pd.DataFrame:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"events_{_geometry_hash(geom_p)}.pkl"
    if path.exists():
        return pd.read_pickle(path)
    events = events_func(geom_p)
    events.to_pickle(path)
    return events


def fd_dT(
    events_func: Callable[[PerturbedGeometry], pd.DataFrame],
    geom: DetectorGeometry,
    perturbation_field: str,
    magnitude: float,
    theta_bins: NDArray[np.float64] | None = None,
    phi_bins: NDArray[np.float64] | None = None,
    cache_dir: str | Path = "characterization/results/tolerances/cache",
) -> NDArray[np.float64]:
    """Return finite-difference derivative of ``T(theta, phi)``.

    Args:
        events_func: Callable receiving a ``PerturbedGeometry`` and returning a
            Monte Carlo event table for that geometry.
        geom: Nominal geometry in mm.
        perturbation_field: One of ``"D"``, ``"tilt_x"``, ``"tilt_y"``,
            ``"dx"``, ``"dy"``.
        magnitude: Central-difference step in mm or rad.
        theta_bins: Optional theta bin edges in rad.
        phi_bins: Optional phi bin edges in rad.
        cache_dir: Directory for small pickle caches. Do not commit large MC
            caches.

    Returns:
        Derivative of Phase-A ``T`` in inverse physical units.

    Edge cases:
        Cache keys are content hashes of the complete perturbed geometry. Any
        exception raised by ``events_func`` is propagated.
    """

    theta = theta_bins if theta_bins is not None else np.linspace(0.0, 0.5 * np.pi, 9)
    phi = phi_bins if phi_bins is not None else np.linspace(-np.pi, np.pi, 17)
    plus = apply_perturbation(geom, _perturb_dict(perturbation_field, magnitude))
    minus = apply_perturbation(geom, _perturb_dict(perturbation_field, -magnitude))
    cache_path = Path(cache_dir)
    plus_events = _events_cached(events_func, plus, cache_path)
    minus_events = _events_cached(events_func, minus, cache_path)
    return (
        pixel_acceptance(plus_events, plus, theta, phi)
        - pixel_acceptance(minus_events, minus, theta, phi)
    ) / (2.0 * magnitude)
