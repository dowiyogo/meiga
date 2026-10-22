"""First-order tolerance error budget for hodoscope solid angle.

Minimal reproducible example:

```python
from characterization.geometry import DetectorGeometry
from characterization.tolerances.error_budget import propagate_to_acceptance

geom = DetectorGeometry(n_channels=2)
variance = propagate_to_acceptance(geom, {"D": 0.5, "dx": 0.5, "dy": 0.5})
```
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from characterization.geometry import DetectorGeometry
from characterization.tolerances.analytical import dOmega_dD, dOmega_dalpha, dOmega_dxy
from characterization.tolerances.tilt_validation import dOmega_dalpha_general


DEFAULT_FIELDS = ("D", "tilt_x", "tilt_y", "dx", "dy")


def _jacobian_field(
    field: str,
    top: int,
    bot: int,
    geom: DetectorGeometry,
    use_general_kernel_for_tilt: bool,
) -> float:
    if field == "D":
        return dOmega_dD(top, bot, geom)
    if field == "tilt_x":
        if use_general_kernel_for_tilt:
            return 1000.0 * dOmega_dalpha_general(geom, (top, bot), "x")
        return dOmega_dalpha(top, bot, geom, "x")
    if field == "tilt_y":
        if use_general_kernel_for_tilt:
            return 1000.0 * dOmega_dalpha_general(geom, (top, bot), "y")
        return dOmega_dalpha(top, bot, geom, "y")
    if field == "dx":
        return dOmega_dxy(top, bot, geom, "x")
    if field == "dy":
        return dOmega_dxy(top, bot, geom, "y")
    raise ValueError(f"Unsupported tolerance field: {field}")


def jacobian_matrix(
    geom: DetectorGeometry,
    fields: Sequence[str] = DEFAULT_FIELDS,
    use_general_kernel_for_tilt: bool = True,
) -> dict[str, NDArray[np.float64]]:
    """Return ``dOmega/dxi`` matrices for requested tolerance fields.

    Args:
        geom: Detector geometry in mm.
        fields: Tolerance fields. Units are mm for ``D``, ``dx``, ``dy`` and
            rad for tilts.
        use_general_kernel_for_tilt: If true, tilt derivatives use the
            general-orientation Gauss-Legendre kernel. If false, they use the
            legacy Phase-B parallel-kernel chain.

    Returns:
        Mapping from field name to ``(n_pixels, n_pixels)`` derivative arrays.

    Edge cases:
        This exact scalar path is intended for small grids and report maps.
        Large full 16x16 pair matrices can be slow because derivatives call the
        Phase-A rectangle integral.
    """

    result = {
        field: np.zeros((geom.n_pixels, geom.n_pixels), dtype=float) for field in fields
    }
    for top in range(geom.n_pixels):
        for bot in range(geom.n_pixels):
            for field in fields:
                result[field][top, bot] = _jacobian_field(
                    field, top, bot, geom, use_general_kernel_for_tilt
                )
    return result


def propagate_to_acceptance(
    geom: DetectorGeometry,
    sigma_dict: Mapping[str, float],
    covariance: NDArray[np.float64] | None = None,
    fields: Sequence[str] | None = None,
    use_general_kernel_for_tilt: bool = True,
) -> NDArray[np.float64]:
    """Propagate tolerance uncertainty into ``Var(Omega)`` per pixel pair.

    Args:
        geom: Detector geometry in mm.
        sigma_dict: Mapping from tolerance field to one-sigma uncertainty. Units
            are mm for ``D``, ``dx``, ``dy`` and rad for tilts.
        covariance: Optional covariance matrix ``Sigma`` in the same field
            order. If omitted, uncorrelated tolerances are assumed and
            ``Var(Omega) = sum_i (dOmega/dxi)^2 sigma_i^2``.
        fields: Optional explicit field order. Required to interpret
            ``covariance``. Defaults to keys in ``sigma_dict``.
        use_general_kernel_for_tilt: Default true. Tilt derivatives are more
            expensive but physically correct because they use the general
            non-parallel rectangle kernel. Set false for backward-compatible
            Phase-B behavior or very fast exploratory checks.

    Returns:
        ``(n_pixels, n_pixels)`` variance of ``Omega`` in sr^2.

    Edge cases:
        If ``covariance`` is supplied, ``sigma_dict`` is only used to infer the
        default field order unless ``fields`` is provided. The correlated case
        evaluates ``Var(Omega) = J Sigma J^T`` independently for each pixel pair.
    """

    selected = tuple(fields if fields is not None else sigma_dict.keys())
    jac = jacobian_matrix(
        geom, selected, use_general_kernel_for_tilt=use_general_kernel_for_tilt
    )
    if covariance is None:
        variance = np.zeros((geom.n_pixels, geom.n_pixels), dtype=float)
        for field in selected:
            sigma = float(sigma_dict.get(field, 0.0))
            variance += jac[field] ** 2 * sigma * sigma
        return variance

    sigma = np.asarray(covariance, dtype=float)
    if sigma.shape != (len(selected), len(selected)):
        raise ValueError("covariance shape must match selected fields")
    stacked = np.stack([jac[field] for field in selected], axis=-1)
    return np.einsum("...i,ij,...j->...", stacked, sigma, stacked)


def dominant_term(
    geom: DetectorGeometry,
    sigma_dict: Mapping[str, float],
    fields: Sequence[str] | None = None,
    use_general_kernel_for_tilt: bool = True,
) -> NDArray[np.str_]:
    """Return the dominant tolerance source per pixel pair.

    Args:
        geom: Detector geometry in mm.
        sigma_dict: One-sigma uncertainties in mm or rad.
        fields: Optional explicit field list.
        use_general_kernel_for_tilt: Use the general-orientation kernel for
            tilt derivatives when true. Set false for legacy Phase-B behavior.

    Returns:
        ``(n_pixels, n_pixels)`` array of field names.

    Edge cases:
        Ties are resolved by field order. Fields with zero sigma are allowed and
        contribute zero variance.
    """

    selected = tuple(fields if fields is not None else sigma_dict.keys())
    jac = jacobian_matrix(
        geom, selected, use_general_kernel_for_tilt=use_general_kernel_for_tilt
    )
    contributions = np.stack(
        [(jac[field] * float(sigma_dict.get(field, 0.0))) ** 2 for field in selected],
        axis=-1,
    )
    indices = np.argmax(contributions, axis=-1)
    labels = np.asarray(selected, dtype=str)
    return labels[indices]
