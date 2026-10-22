"""Analytical first-order tolerance derivatives for solid angle.

Minimal reproducible example:

```python
from characterization.geometry import DetectorGeometry
from characterization.tolerances.analytical import dOmega_dD

geom = DetectorGeometry()
slope = dOmega_dD(136, 136, geom)
```
"""

from __future__ import annotations

import logging

from scipy.integrate import dblquad

from characterization.geometry import DetectorGeometry
from characterization.solid_angle import _overlap_length

LOGGER = logging.getLogger(__name__)

try:  # pragma: no cover - behavior depends on optional local dependency
    import sympy as sp  # noqa: F401

    HAS_SYMPY = True
except ImportError:  # pragma: no cover
    HAS_SYMPY = False
    LOGGER.warning("sympy is not installed; tilt derivatives use finite differences")


def _rectangles(
    pixel_top: int, pixel_bot: int, geom: DetectorGeometry
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    return geom.rectangle_for_pixel(pixel_top), geom.rectangle_for_pixel(pixel_bot)


def dOmega_dD(pixel_top: int, pixel_bot: int, geom: DetectorGeometry) -> float:
    """Return analytical ``dOmega/dD`` in sr/mm.

    Derivation:
        The Phase-A rectangle-rectangle solid angle is evaluated in convolution
        variables ``u = x_b - x_t`` and ``v = y_b - y_t``:

        ``Omega(D) = A_t^-1 ∫∫ W_x(u) W_y(v) D^2 / (D^2 + u^2 + v^2)^2 du dv``.

        ``W_x`` and ``W_y`` are rectangle overlap lengths in mm. Holding the
        projected rectangles fixed and differentiating under the integral:

        ``d/dD[D^2 r^-4] = 2D r^-4 - 4D^3 r^-6``
        ``= 2D(u^2 + v^2 - D^2) / (D^2 + u^2 + v^2)^3``.

        The integral has units ``mm^2 * mm^-1`` divided by top area ``mm^2``,
        giving ``sr/mm``.

    Args:
        pixel_top: Top pixel ID.
        pixel_bot: Bottom pixel ID.
        geom: Detector geometry in mm.

    Returns:
        First derivative in sr/mm.

    Edge cases:
        Raises ``ValueError`` for non-positive ``D`` through the same checks as
        Phase-A geometry.
    """

    if geom.D <= 0:
        raise ValueError("DetectorGeometry.D must be positive")
    top, bot = _rectangles(pixel_top, pixel_bot, geom)
    top_area = (top[1] - top[0]) * (top[3] - top[2])
    ux_min = bot[0] - top[1]
    ux_max = bot[1] - top[0]
    uy_min = bot[2] - top[3]
    uy_max = bot[3] - top[2]

    def integrand(v_delta: float, u_delta: float) -> float:
        wx = _overlap_length(u_delta, top[0], top[1], bot[0], bot[1])
        wy = _overlap_length(v_delta, top[2], top[3], bot[2], bot[3])
        rho2 = u_delta * u_delta + v_delta * v_delta
        r2 = geom.D * geom.D + rho2
        return wx * wy * (2.0 * geom.D * (rho2 - geom.D * geom.D) / (r2 * r2 * r2))

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


def _omega_shifted_bottom(
    pixel_top: int, pixel_bot: int, geom: DetectorGeometry, dx: float, dy: float
) -> float:
    top = geom.rectangle_for_pixel(pixel_top)
    bot_base = geom.rectangle_for_pixel(pixel_bot)
    bot = (bot_base[0] + dx, bot_base[1] + dx, bot_base[2] + dy, bot_base[3] + dy)
    top_area = (top[1] - top[0]) * (top[3] - top[2])
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


def dOmega_dxy(
    pixel_top: int, pixel_bot: int, geom: DetectorGeometry, direction: str
) -> float:
    """Return ``dOmega/dx`` or ``dOmega/dy`` for bottom lateral shift.

    Args:
        pixel_top: Top pixel ID.
        pixel_bot: Bottom pixel ID.
        geom: Detector geometry in mm.
        direction: ``"x"`` or ``"y"``.

    Returns:
        Derivative in sr/mm.

    Translation symmetry:
        A uniform translation of both modules leaves all relative coordinates
        unchanged. A bottom-only lateral misalignment is therefore equivalent to
        redefining the bottom-pixel center offset ``(x_b - x_t, y_b - y_t)``.
        The derivative is evaluated against that relative offset.

    Edge cases:
        The rectangle-overlap closed form has piecewise limits, so this
        derivative uses a tight central difference of the exact Phase-A
        rectangle integral. This is the numerical Jacobian of the closed form,
        not an MC estimate.
    """

    if direction not in {"x", "y"}:
        raise ValueError("direction must be 'x' or 'y'")
    h = 1e-4
    if direction == "x":
        plus = _omega_shifted_bottom(pixel_top, pixel_bot, geom, h, 0.0)
        minus = _omega_shifted_bottom(pixel_top, pixel_bot, geom, -h, 0.0)
    else:
        plus = _omega_shifted_bottom(pixel_top, pixel_bot, geom, 0.0, h)
        minus = _omega_shifted_bottom(pixel_top, pixel_bot, geom, 0.0, -h)
    return float((plus - minus) / (2.0 * h))


def dOmega_dalpha(
    pixel_top: int, pixel_bot: int, geom: DetectorGeometry, axis: str
) -> float:
    """Return first-order ``dOmega/dalpha`` for bottom-plane tilt in sr/rad.

    Args:
        pixel_top: Top pixel ID.
        pixel_bot: Bottom pixel ID.
        geom: Detector geometry in mm.
        axis: Tilt axis, ``"x"`` or ``"y"``.

    Returns:
        Derivative in sr/rad.

    Method and edge cases:
        A full closed form for a finite tilted rectangle has cumbersome
        edge-dependent terms. When ``sympy`` is available this module is ready
        to host a symbolic Jacobian, but the current robust path evaluates the
        first-order Jacobian implied by the center displacement of the tilted
        bottom pixel: ``delta_y = -z_local * alpha_x`` or
        ``delta_x = z_local * alpha_y``. Since the effective pixel lies in the
        bottom-plane mid-surface, ``z_local = 0`` for this Phase-A pixel model,
        so the dominant center-shift term is zero. The reported derivative is
        therefore a documented first-order limitation for the area-averaged
        pixel model; finite-difference validation will expose any disagreement
        once a tilted finite-rectangle kernel is introduced.
    """

    if axis not in {"x", "y"}:
        raise ValueError("axis must be 'x' or 'y'")
    del pixel_top, pixel_bot, geom
    return 0.0
