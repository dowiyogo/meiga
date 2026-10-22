"""Scan inter-plane distance D for Phase-A hodoscope characterization.

Minimal reproducible example:

```bash
python -m characterization.scan_D --d-values 40 80 120 --synthetic-events 2000
```

By default the script uses a lightweight synthetic geometry-only sample. It can
also load precomputed Meiga outputs with ``--input-template 'runs/D_{D}/output.json'``
or invoke a simulator command with ``--simulator-cmd``.
"""

from __future__ import annotations

import argparse
import logging
import shlex
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from characterization.acceptance import geometric_efficiency, pixel_acceptance
from characterization.angular_resolution import sigma_theta_geom, sigma_theta_mc
from characterization.geometry import DetectorGeometry
from characterization.io import load_events, load_geometry
from characterization.solid_angle import compare_omega
from characterization.solid_angle_general import Rectangle3D, omega_general

LOGGER = logging.getLogger(__name__)


def _nearest_channel(
    values: np.ndarray, centers: np.ndarray, half_width: float
) -> np.ndarray:
    distances = np.abs(values[:, None] - centers[None, :])
    index = np.argmin(distances, axis=1)
    valid = distances[np.arange(len(values)), index] <= half_width
    return np.where(valid, index, -1)


def synthetic_events_for_D(
    geom: DetectorGeometry, n_events: int, seed: int = 42
) -> pd.DataFrame:
    """Generate a small EcoMug-like ``cos^2(theta)`` synthetic event sample.

    Args:
        geom: Detector geometry in mm.
        n_events: Number of incident top-module rays to generate.
        seed: Random seed.

    Returns:
        DataFrame with accepted four-plane events. Angles are rad, energies are
        MeV, times are ns.

    Edge cases:
        This generator is geometry-only and does not model energy loss,
        scattering, optical response, or electronics.
    """

    rng = np.random.default_rng(seed)
    centers = geom.channel_centers()
    top_x_index = rng.integers(0, geom.n_channels, n_events)
    top_y_index = rng.integers(0, geom.n_channels, n_events)
    top_x = centers[top_x_index] + rng.uniform(-0.5, 0.5, n_events) * geom.bar_width
    top_y = centers[top_y_index] + rng.uniform(-0.5, 0.5, n_events) * geom.bar_width

    cos_theta = (1.0 - rng.uniform(0.0, 1.0, n_events)) ** (1.0 / 3.0)
    sin_theta = np.sqrt(1.0 - cos_theta * cos_theta)
    phi = rng.uniform(-np.pi, np.pi, n_events)
    dx = geom.D * sin_theta * np.cos(phi) / np.maximum(cos_theta, 1e-12)
    dy = geom.D * sin_theta * np.sin(phi) / np.maximum(cos_theta, 1e-12)
    bot_x = top_x + dx
    bot_y = top_y + dy
    bot_x_index = _nearest_channel(bot_x, centers, 0.5 * geom.bar_width)
    bot_y_index = _nearest_channel(bot_y, centers, 0.5 * geom.bar_width)
    accepted = (bot_x_index >= 0) & (bot_y_index >= 0)

    rows = pd.DataFrame(
        {
            "event_id": np.arange(np.count_nonzero(accepted), dtype=np.int64),
            "pixel_top_x": top_x_index[accepted].astype(np.int64),
            "pixel_top_y": top_y_index[accepted].astype(np.int64),
            "pixel_bot_x": bot_x_index[accepted].astype(np.int64),
            "pixel_bot_y": bot_y_index[accepted].astype(np.int64),
            "theta_rec": np.arctan2(
                np.hypot(
                    centers[bot_x_index[accepted]] - centers[top_x_index[accepted]],
                    centers[bot_y_index[accepted]] - centers[top_y_index[accepted]],
                ),
                geom.D,
            ),
            "phi_rec": np.arctan2(
                centers[bot_y_index[accepted]] - centers[top_y_index[accepted]],
                centers[bot_x_index[accepted]] - centers[top_x_index[accepted]],
            ),
            "energy_top": np.ones(np.count_nonzero(accepted), dtype=float),
            "energy_bot": np.ones(np.count_nonzero(accepted), dtype=float),
            "t_top": np.zeros(np.count_nonzero(accepted), dtype=float),
            "t_bot": geom.D
            / 299.792458
            * np.ones(np.count_nonzero(accepted), dtype=float),
        }
    )
    rows.attrs["n_incident_top"] = n_events
    rows.attrs["n_total_events"] = n_events
    rows.attrs["n_4plane"] = len(rows)
    rows.attrs["geometry"] = geom
    rows.attrs["theta_true"] = np.arccos(cos_theta[accepted])
    rows.attrs["phi_true"] = phi[accepted]
    return rows


def _load_or_generate_events(
    D: float,
    args: argparse.Namespace,
    base_geom: DetectorGeometry,
) -> tuple[pd.DataFrame, DetectorGeometry]:
    geom = DetectorGeometry(
        D=D,
        bar_length=base_geom.bar_length,
        bar_width=base_geom.bar_width,
        bar_thickness=base_geom.bar_thickness,
        n_channels=base_geom.n_channels,
        pixel_pitch=base_geom.pixel_pitch,
        sub_plane_offset=base_geom.sub_plane_offset,
        module_size=base_geom.module_size,
    )
    if args.simulator_cmd:
        command = args.simulator_cmd.format(D=D)
        LOGGER.info("Running simulator command for D=%s: %s", D, command)
        subprocess.run(shlex.split(command), check=True)
        if not args.input_template:
            raise ValueError(
                "--simulator-cmd requires --input-template for output loading"
            )
    if args.input_template:
        path = Path(args.input_template.format(D=D))
        loaded_geom = load_geometry(path)
        geom = DetectorGeometry(
            D=D,
            bar_length=loaded_geom.bar_length,
            bar_width=loaded_geom.bar_width,
            bar_thickness=loaded_geom.bar_thickness,
            n_channels=loaded_geom.n_channels,
            pixel_pitch=loaded_geom.pixel_pitch,
            sub_plane_offset=loaded_geom.sub_plane_offset,
            module_size=loaded_geom.module_size,
        )
        events = load_events(path)
        events.attrs["geometry"] = geom
        return events, geom
    return (
        synthetic_events_for_D(
            geom, args.synthetic_events, seed=args.seed + int(D * 10)
        ),
        geom,
    )


def _omega_pair_distribution(geom: DetectorGeometry) -> np.ndarray:
    """Return omega for the full pixel-pair grid using offset multiplicities.

    The uncertainty bars in the D-scan figures are the standard deviation of
    this distribution across pixel pairs. They therefore represent geometrical
    pair-to-pair scatter, not Monte Carlo noise.
    """

    centers = geom.pixel_centers()
    half = 0.5 * geom.bar_width
    top_rect = Rectangle3D(
        center=(0.0, 0.0, 0.0),
        normal=(0.0, 0.0, -1.0),
        u_axis=(1.0, 0.0, 0.0),
        v_axis=(0.0, 1.0, 0.0),
        half_extent_u=half,
        half_extent_v=half,
    )
    offsets: dict[tuple[float, float], int] = {}
    for top_xy in centers:
        for bot_xy in centers:
            key = tuple(float(v) for v in np.round(bot_xy - top_xy, 12))
            offsets[key] = offsets.get(key, 0) + 1

    values: list[float] = []
    for (dx, dy), count in offsets.items():
        bot_rect = Rectangle3D(
            center=(dx, dy, -geom.D),
            normal=(0.0, 0.0, 1.0),
            u_axis=(1.0, 0.0, 0.0),
            v_axis=(0.0, 1.0, 0.0),
            half_extent_u=half,
            half_extent_v=half,
        )
        omega = omega_general(top_rect, bot_rect, method="gauss_legendre", n_quad=6)
        values.extend([omega] * count)
    return np.asarray(values, dtype=float)


def _fit_power_law(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    valid = (x > 0.0) & (y > 0.0)
    if np.count_nonzero(valid) < 2:
        return float("nan"), float("nan")
    slope, intercept = np.polyfit(np.log(x[valid]), np.log(y[valid]), 1)
    return float(np.exp(intercept)), float(slope)


def _plot_observables(table: pd.DataFrame, outdir: Path) -> None:
    observables = [
        (
            "mean_Omega",
            "std_Omega",
            "Mean per-pair solid angle Ω [sr]",
            True,
        ),
        (
            "mean_T",
            "std_T",
            "Mean per-pair acceptance T [dimensionless]",
            True,
        ),
        ("sigma_theta_rec", None, "sigma(theta_rec) [rad]", True),
        ("n_4plane_per_event", None, "Four-plane events / incident event", False),
    ]
    for column, std_column, ylabel, loglog in observables:
        fig, ax = plt.subplots(figsize=(6.8, 4.6))
        x = table["D"].to_numpy(dtype=float)
        y = table[column].to_numpy(dtype=float)
        yerr = (
            table[std_column].to_numpy(dtype=float)
            if std_column and std_column in table
            else None
        )
        if yerr is not None:
            ax.errorbar(
                x,
                y,
                yerr=yerr,
                marker="o",
                capsize=3,
                linewidth=1.4,
                label="mean ± 1 std across pairs",
            )
            amplitude, exponent = _fit_power_law(x, y)
            if np.isfinite(amplitude):
                x_fit = np.geomspace(float(np.min(x)), float(np.max(x)), 200)
                ax.plot(
                    x_fit,
                    amplitude * x_fit**exponent,
                    linestyle="--",
                    color="black",
                    label=rf"fit $A D^n$, $n={exponent:.2f}$",
                )
        else:
            ax.plot(x, y, marker="o")
        ax.set_xlabel("Inter-plane separation D [mm]")
        ax.set_ylabel(ylabel)
        if column in {"mean_Omega", "mean_T"}:
            ax.set_title(
                "NA64-mini telescope, 16×16 pixel grid, far-field scaling"
            )
        ax.grid(True, alpha=0.3)
        if loglog:
            ax.set_xscale("log")
            ax.set_yscale("log")
        if yerr is not None:
            ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"{column}_vs_D.png", dpi=220)
        fig.savefig(outdir / f"{column}_vs_D.pdf")
        plt.close(fig)


def write_report(
    outdir: Path,
    omega_table: pd.DataFrame,
    epsilon_table: pd.DataFrame,
    scan_table: pd.DataFrame,
) -> None:
    """Write the Phase-A Markdown report.

    Args:
        outdir: Output directory.
        omega_table: Analytical-vs-MC solid-angle comparison in sr.
        epsilon_table: Geometric-efficiency curve versus theta in rad.
        scan_table: D-scan observables with lengths in mm.

    Returns:
        ``None``. Writes ``REPORT_phaseA.md``.

    Edge cases:
        Existing reports are overwritten intentionally by the characterization
        pipeline.
    """

    def markdown_table(frame: pd.DataFrame) -> str:
        columns = list(frame.columns)
        lines = [
            "| " + " | ".join(columns) + " |",
            "| " + " | ".join(["---"] * len(columns)) + " |",
        ]
        for _, row in frame.iterrows():
            values = []
            for column in columns:
                value = row[column]
                if isinstance(value, float):
                    values.append(f"{value:.6g}")
                else:
                    values.append(str(value))
            lines.append("| " + " | ".join(values) + " |")
        return "\n".join(lines)

    report = outdir / "REPORT_phaseA.md"
    omega_preview = markdown_table(omega_table.head(12))
    epsilon_preview = markdown_table(epsilon_table)
    sigma_preview = markdown_table(
        scan_table[["D", "sigma_theta_geom", "sigma_theta_rec"]]
    )
    report.write_text(
        "\n".join(
            [
                "# Phase A Characterization Report",
                "",
                "## Analytical vs MC Solid Angle",
                "",
                omega_preview,
                "",
                f"- mean relative error: {omega_table['relative_error'].mean():.6g}",
                f"- max relative error: {omega_table['relative_error'].max():.6g}",
                "",
                "## Geometric Efficiency epsilon(theta)",
                "",
                epsilon_preview,
                "",
                "## Angular Resolution vs D",
                "",
                sigma_preview,
                "",
                "## Open Issues and Known Limitations",
                "",
                "- MC solid-angle estimates use lightweight synthetic sampling in the default pipeline.",
                "- Real Meiga output does not yet expose generated top-module incidence per pixel.",
                "- Passive material, optical attenuation, SiPM response, and electronics are outside Phase A.",
                "- Tolerance studies for Delta D, tilt, and lateral misalignment are intentionally deferred.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    """Run the D scan and write CSV, figures, and report.

    Args:
        argv: Optional CLI argument list.

    Returns:
        Process exit code. Units are mm, rad, ns, and MeV as documented by the
        individual modules.

    Edge cases:
        Without precomputed outputs or simulator command, a synthetic dataset is
        used so the pipeline remains reproducible on machines without Geant4.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--d-values", nargs="+", type=float)
    parser.add_argument(
        "--input-template", help="Precomputed output path template using {D}"
    )
    parser.add_argument(
        "--simulator-cmd", help="Command template to run before loading {D}"
    )
    parser.add_argument("--synthetic-events", type=int, default=5000)
    parser.add_argument("--omega-samples", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", type=Path, default=Path("characterization/results"))
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    base_geom = DetectorGeometry()
    theta_bins = np.linspace(0.0, 0.5 * np.pi, 9)
    phi_bins = np.linspace(-np.pi, np.pi, 17)
    rows: list[dict[str, float]] = []
    first_epsilon: pd.DataFrame | None = None
    first_omega: pd.DataFrame | None = None

    d_values = (
        np.asarray(args.d_values, dtype=float)
        if args.d_values
        else np.geomspace(30.0, 300.0, 20)
    )

    for D in d_values:
        events, geom = _load_or_generate_events(D, args, base_geom)
        epsilon = geometric_efficiency(events, geom, theta_bins)
        acceptance = pixel_acceptance(events, geom, theta_bins, phi_bins)
        omega_table = compare_omega(geom, args.omega_samples)
        omega_pairs = _omega_pair_distribution(geom)
        t_pairs = omega_pairs / (2.0 * np.pi)
        if first_epsilon is None:
            first_epsilon = epsilon
            first_omega = omega_table
        positive_acceptance = acceptance[acceptance > 0]
        rows.append(
            {
                "D": D,
                "mean_Omega": float(np.mean(omega_pairs)),
                "std_Omega": float(np.std(omega_pairs, ddof=1)),
                "mean_T": float(np.mean(t_pairs)),
                "std_T": float(np.std(t_pairs, ddof=1)),
                "mean_T_synthetic_acceptance": (
                    float(positive_acceptance.mean()) if positive_acceptance.size else 0.0
                ),
                "sigma_theta_rec": (
                    float(events["theta_rec"].std(ddof=1)) if len(events) > 1 else 0.0
                ),
                "sigma_theta_geom": sigma_theta_geom(geom),
                "sigma_theta_mc": sigma_theta_mc(events, geom),
                "n_4plane_per_event": len(events)
                / float(events.attrs.get("n_incident_top", len(events) or 1)),
            }
        )

    table = pd.DataFrame(rows)
    table.to_csv(args.outdir / "scan_D_summary.csv", index=False)
    _plot_observables(table, args.outdir)
    assert first_epsilon is not None and first_omega is not None
    first_omega.to_csv(args.outdir / "omega_comparison_representative.csv", index=False)
    first_epsilon.to_csv(args.outdir / "epsilon_theta_representative.csv", index=False)
    write_report(args.outdir, first_omega, first_epsilon, table)
    LOGGER.info("Wrote characterization results to %s", args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
