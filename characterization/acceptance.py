"""Angular acceptance and geometric efficiency estimators.

Minimal reproducible example:

```python
import numpy as np
from characterization.acceptance import geometric_efficiency, pixel_acceptance
from characterization.io import load_events, load_geometry

geom = load_geometry("output.json")
events = load_events("output.json")
theta_bins = np.linspace(0.0, np.pi / 2.0, 10)
phi_bins = np.linspace(-np.pi, np.pi, 17)
T = pixel_acceptance(events, geom, theta_bins, phi_bins)
epsilon = geometric_efficiency(events, geom, theta_bins)
```

The normalization assumes EcoMug's default downward cosmic-muon angular law
``p(theta, phi) dtheta dphi = (3 / 2pi) cos(theta)^2 sin(theta) dtheta dphi``.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray

from characterization.geometry import DetectorGeometry

LOGGER = logging.getLogger(__name__)


def cos2_theta_bin_prob(theta_bins: ArrayLike) -> NDArray[np.float64]:
    """Return normalized ``cos^2(theta)`` probabilities per theta bin.

    Args:
        theta_bins: Bin edges in rad. Values outside ``[0, pi/2]`` are clipped
            to the physical downward hemisphere.

    Returns:
        Array of probabilities, one per theta bin. Dimensionless.

    Edge cases:
        Non-overlapping bins outside the downward hemisphere receive zero
        probability.
    """

    theta = np.asarray(theta_bins, dtype=float)
    low = np.clip(theta[:-1], 0.0, 0.5 * np.pi)
    high = np.clip(theta[1:], 0.0, 0.5 * np.pi)
    prob = np.cos(low) ** 3 - np.cos(high) ** 3
    return np.where(high > low, prob, 0.0)


def cos2_bin_prob(theta_bins: ArrayLike, phi_bins: ArrayLike) -> NDArray[np.float64]:
    """Return normalized ``cos^2(theta)`` probabilities per theta-phi bin.

    Args:
        theta_bins: Zenith-angle bin edges in rad.
        phi_bins: Azimuth bin edges in rad.

    Returns:
        Array of shape ``(n_theta, n_phi)`` with dimensionless probabilities.

    Edge cases:
        Phi bins are interpreted by their width and are not wrapped; use
        ``[-pi, pi]`` for full azimuth.
    """

    theta_prob = cos2_theta_bin_prob(theta_bins)
    phi = np.asarray(phi_bins, dtype=float)
    phi_prob = np.diff(phi) / (2.0 * np.pi)
    return theta_prob[:, None] * phi_prob[None, :]


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


def pixel_acceptance(
    events: pd.DataFrame,
    geom: DetectorGeometry,
    theta_bins: ArrayLike,
    phi_bins: ArrayLike,
) -> NDArray[np.float64]:
    """Estimate per-pixel-pair angular transmission ``T(theta, phi)``.

    Args:
        events: DataFrame from ``load_events`` with one four-plane coincidence
            per row. Angles are in rad.
        geom: Detector geometry with lengths in mm.
        theta_bins: Zenith-angle bin edges in rad.
        phi_bins: Azimuth bin edges in rad.

    Returns:
        Array with shape ``(n_pix_top, n_pix_bot, n_theta, n_phi)``. Values are
        dimensionless ratios ``N_pair_bin / (N_incident_top * P_cos2_bin)``.

    Edge cases:
        If ``events.attrs['n_incident_top']`` is absent, the number of loaded
        four-plane events is used as a fallback and a warning is logged. Empty
        angular bins with zero incident probability are returned as zero.
    """

    n_incident = float(events.attrs.get("n_incident_top", len(events)))
    if "n_incident_top" not in events.attrs:
        LOGGER.warning("events.attrs['n_incident_top'] missing; using len(events)")
    n_theta = len(theta_bins) - 1
    n_phi = len(phi_bins) - 1
    counts = np.zeros((geom.n_pixels, geom.n_pixels, n_theta, n_phi), dtype=float)
    if events.empty or n_incident <= 0:
        return counts

    top, bot = _pixel_ids(events, geom)
    theta_idx = np.digitize(events["theta_rec"].to_numpy(dtype=float), theta_bins) - 1
    phi_idx = np.digitize(events["phi_rec"].to_numpy(dtype=float), phi_bins) - 1
    valid = (
        (top >= 0)
        & (top < geom.n_pixels)
        & (bot >= 0)
        & (bot < geom.n_pixels)
        & (theta_idx >= 0)
        & (theta_idx < n_theta)
        & (phi_idx >= 0)
        & (phi_idx < n_phi)
    )
    np.add.at(counts, (top[valid], bot[valid], theta_idx[valid], phi_idx[valid]), 1.0)

    bin_prob = cos2_bin_prob(theta_bins, phi_bins)
    denom = n_incident * bin_prob
    with np.errstate(divide="ignore", invalid="ignore"):
        acceptance = counts / denom[None, None, :, :]
    acceptance[~np.isfinite(acceptance)] = 0.0
    return acceptance


def geometric_efficiency(
    events: pd.DataFrame, geom: DetectorGeometry, theta_bins: ArrayLike
) -> pd.DataFrame:
    """Estimate ``epsilon(theta) = N_4plane(theta) / N_incident_top(theta)``.

    Args:
        events: DataFrame from ``load_events`` with reconstructed ``theta_rec``
            in rad.
        geom: Detector geometry in mm. Present for API symmetry and future
            tolerance studies; it is not modified.
        theta_bins: Zenith-angle bin edges in rad.

    Returns:
        DataFrame with ``theta_low``, ``theta_high``, ``theta_center``,
        ``n_4plane``, ``n_incident_top`` and ``epsilon``. Angles are rad.

    Edge cases:
        Empty input returns zero counts. Missing incident-count metadata falls
        back to ``len(events)`` and logs a warning.
    """

    del geom
    n_incident_total = float(events.attrs.get("n_incident_top", len(events)))
    if "n_incident_top" not in events.attrs:
        LOGGER.warning("events.attrs['n_incident_top'] missing; using len(events)")
    theta = np.asarray(theta_bins, dtype=float)
    counts, _ = np.histogram(events["theta_rec"].to_numpy(dtype=float), bins=theta)
    p_theta = cos2_theta_bin_prob(theta)
    n_incident = n_incident_total * p_theta
    with np.errstate(divide="ignore", invalid="ignore"):
        epsilon = counts / n_incident
    epsilon[~np.isfinite(epsilon)] = 0.0
    return pd.DataFrame(
        {
            "theta_low": theta[:-1],
            "theta_high": theta[1:],
            "theta_center": 0.5 * (theta[:-1] + theta[1:]),
            "n_4plane": counts.astype(int),
            "n_incident_top": n_incident,
            "epsilon": epsilon,
        }
    )
