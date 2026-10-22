"""General-orientation solid-angle kernel for finite rectangles.

Minimal reproducible example:

```python
import numpy as np
from characterization.solid_angle_general import Rectangle3D, omega_general

top = Rectangle3D([0, 0, 0], [0, 0, 1], [1, 0, 0], [0, 1, 0], 1.5, 1.5)
bot = Rectangle3D([0, 0, 98.4], [0, 0, -1], [1, 0, 0], [0, 1, 0], 1.5, 1.5)
omega = omega_general(top, bot, method="gauss_legendre", n_quad=12)
```

Convention:
    ``omega_general`` returns the area-averaged effective solid angle

    ``Omega_eff = |S1|^-1 int_S1 int_S2 (n1.rhat)(n2.-rhat) / r^2 dA2 dA1``.

    ``rhat`` points from the top rectangle point to the bottom rectangle point.
    The top normal must face the bottom plane and the bottom normal must face
    the top plane for positive acceptance. For parallel rectangles this reduces
    to the Phase-A kernel ``D^2 / r^4`` averaged over the top area.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class Rectangle3D:
    """Finite rectangle embedded in 3D.

    Args:
        center: Rectangle center in mm, shape ``(3,)``.
        normal: Surface normal, shape ``(3,)``.
        u_axis: First in-plane direction.
        v_axis: Second in-plane direction.
        half_extent_u: Half side length along ``u_axis`` in mm.
        half_extent_v: Half side length along ``v_axis`` in mm.

    Edge cases:
        Vectors are normalized and orthogonalized on construction. Degenerate
        axes or non-positive extents raise ``ValueError``.
    """

    center: ArrayLike
    normal: ArrayLike
    u_axis: ArrayLike
    v_axis: ArrayLike
    half_extent_u: float
    half_extent_v: float

    def __post_init__(self) -> None:
        center = np.asarray(self.center, dtype=float)
        normal = _normalize(np.asarray(self.normal, dtype=float), "normal")
        u_axis = np.asarray(self.u_axis, dtype=float)
        u_axis = u_axis - np.dot(u_axis, normal) * normal
        u_axis = _normalize(u_axis, "u_axis")
        v_axis = np.asarray(self.v_axis, dtype=float)
        v_axis = v_axis - np.dot(v_axis, normal) * normal
        v_axis = v_axis - np.dot(v_axis, u_axis) * u_axis
        v_axis = _normalize(v_axis, "v_axis")
        if center.shape != (3,):
            raise ValueError("center must have shape (3,)")
        if self.half_extent_u <= 0 or self.half_extent_v <= 0:
            raise ValueError("Rectangle half extents must be positive")
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "normal", normal)
        object.__setattr__(self, "u_axis", u_axis)
        object.__setattr__(self, "v_axis", v_axis)
        object.__setattr__(self, "half_extent_u", float(self.half_extent_u))
        object.__setattr__(self, "half_extent_v", float(self.half_extent_v))

    @property
    def area(self) -> float:
        """Return rectangle area in mm^2."""

        return 4.0 * self.half_extent_u * self.half_extent_v


def _normalize(vector: NDArray[np.float64], name: str) -> NDArray[np.float64]:
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0:
        raise ValueError(f"{name} must be non-zero")
    return vector / norm


@lru_cache(maxsize=32)
def _legendre_nodes(n_quad: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    if n_quad < 2:
        raise ValueError("n_quad must be >= 2")
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    return nodes.astype(float), weights.astype(float)


def _quadrature_points(
    rect: Rectangle3D, n_quad: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    nodes, weights = _legendre_nodes(n_quad)
    u_nodes = nodes * rect.half_extent_u
    v_nodes = nodes * rect.half_extent_v
    u_weights = weights * rect.half_extent_u
    v_weights = weights * rect.half_extent_v
    uu, vv = np.meshgrid(u_nodes, v_nodes, indexing="ij")
    wu, wv = np.meshgrid(u_weights, v_weights, indexing="ij")
    points = (
        np.asarray(rect.center)[None, :]
        + uu.ravel()[:, None] * np.asarray(rect.u_axis)[None, :]
        + vv.ravel()[:, None] * np.asarray(rect.v_axis)[None, :]
    )
    area_weights = (wu * wv).ravel()
    return points, area_weights


def omega_general(
    rect_top: Rectangle3D,
    rect_bot: Rectangle3D,
    method: str = "gauss_legendre",
    **kwargs: float | int,
) -> float:
    """Return general-orientation effective solid angle in sr.

    Args:
        rect_top: Top rectangle. Its normal should face the bottom rectangle.
        rect_bot: Bottom rectangle. Its normal should face the top rectangle.
        method: ``"gauss_legendre"`` or ``"monte_carlo"``.
        kwargs: ``n_quad`` for Gauss-Legendre, default 10; ``n_samples`` and
            optional ``seed`` for Monte Carlo, default ``1e5``.

    Returns:
        Area-averaged effective solid angle in sr. Lengths are mm.

    Edge cases:
        Degenerate rectangles raise during construction. Coplanar rectangles
        raise ``ValueError`` because the kernel is singular or ambiguous. Any
        quadrature or MC point pair with ``r -> 0`` raises ``FloatingPointError``.
    """

    separation = np.asarray(rect_bot.center) - np.asarray(rect_top.center)
    plane_gap = abs(float(np.dot(separation, rect_top.normal)))
    if plane_gap < 1e-12:
        raise ValueError("Coplanar or nearly coplanar rectangles are not supported")
    if method == "gauss_legendre":
        return _omega_gauss_legendre(rect_top, rect_bot, int(kwargs.get("n_quad", 10)))
    if method == "monte_carlo":
        return _omega_monte_carlo(
            rect_top,
            rect_bot,
            int(kwargs.get("n_samples", 100_000)),
            int(kwargs.get("seed", 12345)),
        )
    raise ValueError("method must be 'gauss_legendre' or 'monte_carlo'")


def _kernel(
    top_points: NDArray[np.float64],
    bot_points: NDArray[np.float64],
    top_normal: NDArray[np.float64],
    bot_normal: NDArray[np.float64],
) -> NDArray[np.float64]:
    r_vec = bot_points[None, :, :] - top_points[:, None, :]
    r2 = np.einsum("...i,...i->...", r_vec, r_vec)
    if np.any(r2 <= 1e-24):
        raise FloatingPointError("Encountered r -> 0 in solid-angle kernel")
    r = np.sqrt(r2)
    n_top_dot = np.einsum("...i,i->...", r_vec, top_normal) / r
    n_bot_dot = np.einsum("...i,i->...", -r_vec, bot_normal) / r
    return n_top_dot * n_bot_dot / r2


def _omega_gauss_legendre(
    rect_top: Rectangle3D, rect_bot: Rectangle3D, n_quad: int
) -> float:
    top_points, top_weights = _quadrature_points(rect_top, n_quad)
    bot_points, bot_weights = _quadrature_points(rect_bot, n_quad)
    values = _kernel(top_points, bot_points, rect_top.normal, rect_bot.normal)
    integral = np.sum(values * top_weights[:, None] * bot_weights[None, :])
    return float(integral / rect_top.area)


def _sample_rectangle(
    rect: Rectangle3D, rng: np.random.Generator, n_samples: int
) -> NDArray[np.float64]:
    u = rng.uniform(-rect.half_extent_u, rect.half_extent_u, n_samples)
    v = rng.uniform(-rect.half_extent_v, rect.half_extent_v, n_samples)
    return (
        np.asarray(rect.center)[None, :]
        + u[:, None] * np.asarray(rect.u_axis)[None, :]
        + v[:, None] * np.asarray(rect.v_axis)[None, :]
    )


def _omega_monte_carlo(
    rect_top: Rectangle3D, rect_bot: Rectangle3D, n_samples: int, seed: int
) -> float:
    if n_samples <= 0:
        raise ValueError("n_samples must be positive")
    rng = np.random.default_rng(seed)
    top_points = _sample_rectangle(rect_top, rng, n_samples)
    bot_points = _sample_rectangle(rect_bot, rng, n_samples)
    r_vec = bot_points - top_points
    r2 = np.einsum("ij,ij->i", r_vec, r_vec)
    if np.any(r2 <= 1e-24):
        raise FloatingPointError("Encountered r -> 0 in solid-angle kernel")
    r = np.sqrt(r2)
    n_top_dot = np.einsum("ij,j->i", r_vec, rect_top.normal) / r
    n_bot_dot = np.einsum("ij,j->i", -r_vec, rect_bot.normal) / r
    values = n_top_dot * n_bot_dot / r2
    return float(rect_bot.area * np.mean(values))
