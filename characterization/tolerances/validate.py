"""Validate Phase-B analytical tolerance derivatives against finite differences.

Minimal reproducible example:

```bash
python -m characterization.tolerances.validate
```

The script writes figures and ``REPORT_phaseB.md`` under
``characterization/results/tolerances``.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from characterization.geometry import DetectorGeometry
from characterization.solid_angle import omega_analytical
from characterization.tolerances.analytical import dOmega_dD, dOmega_dalpha, dOmega_dxy
from characterization.tolerances.error_budget import dominant_term
from characterization.tolerances.finite_diff import fd_dOmega, omega_perturbed_pair
from characterization.tolerances.perturbation import apply_perturbation
from characterization.tolerances.tilt_validation import (
    alpha_threshold_for_bias,
    compare_kernels,
    dOmega_dalpha_general,
)

LOGGER = logging.getLogger(__name__)


def _representative_pairs(geom: DetectorGeometry) -> dict[str, tuple[int, int]]:
    center = (geom.n_channels // 2) + geom.n_channels * (geom.n_channels // 2)
    edge = geom.n_channels // 2
    corner = geom.n_pixels - 1
    return {
        "center": (center, center),
        "edge": (edge, edge),
        "corner_to_corner": (0, corner),
    }


def _pair_omega(
    field: str, value: float, top: int, bot: int, geom: DetectorGeometry
) -> float:
    if field == "D":
        return omega_perturbed_pair(top, bot, apply_perturbation(geom, {"dD": value}))
    return omega_perturbed_pair(top, bot, apply_perturbation(geom, {field: value}))


def _analytical_derivative(
    field: str, top: int, bot: int, geom: DetectorGeometry
) -> float:
    if field == "D":
        return dOmega_dD(top, bot, geom)
    if field == "tilt_x":
        return dOmega_dalpha(top, bot, geom, "x")
    if field == "tilt_y":
        return dOmega_dalpha(top, bot, geom, "y")
    if field == "dx":
        return dOmega_dxy(top, bot, geom, "x")
    if field == "dy":
        return dOmega_dxy(top, bot, geom, "y")
    raise ValueError(field)


def _plot_sweep(
    field: str,
    top: int,
    bot: int,
    geom: DetectorGeometry,
    values: np.ndarray,
    outdir: Path,
    label: str,
) -> pd.DataFrame:
    nominal = omega_analytical(top, bot, geom)
    analytical = _analytical_derivative(field, top, bot, geom)
    omega = np.asarray([_pair_omega(field, value, top, bot, geom) for value in values])
    fd = np.gradient(omega, values)
    line = nominal + analytical * values
    with np.errstate(divide="ignore", invalid="ignore"):
        discrepancy = (analytical - fd) / fd
    discrepancy[~np.isfinite(discrepancy)] = np.nan

    fig, axes = plt.subplots(3, 1, figsize=(6.4, 8.0), sharex=True)
    axes[0].plot(values, omega, marker="o", label="perturbed omega")
    axes[0].plot(values, line, linestyle="--", label="analytical slope")
    axes[0].set_ylabel("Omega [sr]")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(values, fd, marker="o", label="finite difference")
    axes[1].axhline(analytical, linestyle="--", color="black", label="analytical")
    axes[1].set_ylabel("dOmega/dxi")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[2].plot(values, discrepancy, marker="o")
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_xlabel(
        f"{field} perturbation [{'mm' if field in {'D', 'dx', 'dy'} else 'rad'}]"
    )
    axes[2].set_ylabel("(analytical - FD) / FD")
    axes[2].grid(True, alpha=0.3)
    fig.suptitle(f"{label}: {field}")
    fig.tight_layout()
    stem = f"{label}_{field}"
    fig.savefig(outdir / f"{stem}.png", dpi=160)
    fig.savefig(outdir / f"{stem}.pdf")
    plt.close(fig)
    return pd.DataFrame(
        {
            "pair": label,
            "field": field,
            "value": values,
            "omega_sr": omega,
            "fd_derivative": fd,
            "analytical_derivative": analytical,
            "relative_discrepancy": discrepancy,
        }
    )


def _diagonal_derivative_maps(geom: DetectorGeometry) -> dict[str, np.ndarray]:
    fields = ["D", "tilt_x", "tilt_y", "dx", "dy"]
    maps = {
        field: np.zeros((geom.n_channels, geom.n_channels), dtype=float)
        for field in fields
    }
    for pixel in range(geom.n_pixels):
        x = pixel % geom.n_channels
        y = pixel // geom.n_channels
        maps["D"][y, x] = dOmega_dD(pixel, pixel, geom)
        maps["tilt_x"][y, x] = 1000.0 * dOmega_dalpha_general(
            geom, (pixel, pixel), "x", h=0.1, n_quad=12
        )
        maps["tilt_y"][y, x] = 1000.0 * dOmega_dalpha_general(
            geom, (pixel, pixel), "y", h=0.1, n_quad=12
        )
        maps["dx"][y, x] = dOmega_dxy(pixel, pixel, geom, "x")
        maps["dy"][y, x] = dOmega_dxy(pixel, pixel, geom, "y")
    return maps


def _plot_maps(maps: dict[str, np.ndarray], outdir: Path) -> None:
    labels = {
        "D": r"$\partial\Omega/\partial D$ [sr/mm]",
        "tilt_x": r"$\partial\Omega/\partial\alpha_x$ [sr/mrad]",
        "tilt_y": r"$\partial\Omega/\partial\alpha_y$ [sr/mrad]",
        "dx": r"$\partial\Omega/\partial x$ [sr/mm]",
        "dy": r"$\partial\Omega/\partial y$ [sr/mm]",
    }
    titles = {
        "D": r"$\partial\Omega/\partial D$: nearly uniform; D affects all pairs equally",
        "tilt_x": r"$\partial\Omega/\partial\alpha_x$: stripes; y is the lever arm",
        "tilt_y": r"$\partial\Omega/\partial\alpha_y$: stripes; x is the lever arm",
        "dx": r"$\partial\Omega/\partial x$: near zero; centered pairs cancel",
        "dy": r"$\partial\Omega/\partial y$: near zero; centered pairs cancel",
    }
    for name, matrix in maps.items():
        vmax = float(np.nanmax(np.abs(matrix))) or 1.0
        fig, ax = plt.subplots(figsize=(5.8, 4.8))
        image = ax.imshow(
            matrix, origin="lower", cmap="coolwarm", vmin=-vmax, vmax=vmax
        )
        ax.set_title(titles[name], fontsize=10)
        ax.set_xlabel("pixel x index")
        ax.set_ylabel("pixel y index")
        fig.colorbar(image, ax=ax, label=labels[name])
        fig.tight_layout()
        fig.savefig(outdir / f"map_dOmega_d{name}.png", dpi=220)
        fig.savefig(outdir / f"map_dOmega_d{name}.pdf")
        plt.close(fig)

    order = ["D", "tilt_x", "tilt_y", "dx", "dy"]
    fig, axes = plt.subplots(1, 5, figsize=(20.0, 4.2), constrained_layout=True)
    for ax, name in zip(axes, order):
        matrix = maps[name]
        vmax = float(np.nanmax(np.abs(matrix))) or 1.0
        image = ax.imshow(
            matrix, origin="lower", cmap="coolwarm", vmin=-vmax, vmax=vmax
        )
        ax.set_title(titles[name], fontsize=8)
        ax.set_xlabel("pixel x index")
        ax.set_ylabel("pixel y index")
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label=labels[name])
    fig.suptitle(
        "Sensitivity of per-pair solid angle to geometric tolerances at D = 150 mm"
    )
    fig.savefig(outdir / "derivative_maps_composite.png", dpi=220)
    fig.savefig(outdir / "derivative_maps_composite.pdf")
    plt.close(fig)


def _dominant_diagonal_map(
    geom: DetectorGeometry, maps: dict[str, np.ndarray], outdir: Path
) -> tuple[np.ndarray, pd.DataFrame]:
    sigmas = {"D": 0.5, "alpha": 0.005, "xy": 0.5}
    contributions = {
        "D": np.abs(maps["D"]) * sigmas["D"],
        "alpha": np.sqrt(maps["tilt_x"] ** 2 + maps["tilt_y"] ** 2) * sigmas["alpha"],
        "xy": np.sqrt(maps["dx"] ** 2 + maps["dy"] ** 2) * sigmas["xy"],
    }
    labels = np.asarray(list(contributions.keys()))
    stack = np.stack([contributions[field] ** 2 for field in labels], axis=-1)
    dominant = labels[np.argmax(stack, axis=-1)]
    colors = {"D": 0, "alpha": 1, "xy": 2}
    coded = np.vectorize(colors.get)(dominant)
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    image = ax.imshow(coded, origin="lower", cmap="tab10", vmin=0, vmax=2)
    ax.set_title("Dominant tolerance source")
    ax.set_xlabel("pixel x")
    ax.set_ylabel("pixel y")
    cbar = fig.colorbar(image, ax=ax, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(["D", "alpha", "xy"])
    fig.tight_layout()
    fig.savefig(outdir / "dominant_tolerance_map.png", dpi=160)
    fig.savefig(outdir / "dominant_tolerance_map.pdf")
    plt.close(fig)

    center = (geom.n_channels // 2, geom.n_channels // 2)
    corner = (0, 0)
    rows = []
    for field, contribution in contributions.items():
        rows.append(
            {
                "field": field,
                "sigma_input": sigmas[field],
                "sigmaOmega_center_sr": contribution[center[1], center[0]],
                "sigmaOmega_corner_sr": contribution[corner[1], corner[0]],
                "sigmaOmega_global_mean_sr": float(np.mean(contribution)),
            }
        )
    return dominant, pd.DataFrame(rows)


def _legacy_tilt_budget(geom: DetectorGeometry) -> pd.DataFrame:
    maps_x = np.zeros((geom.n_channels, geom.n_channels), dtype=float)
    maps_y = np.zeros_like(maps_x)
    for pixel in range(geom.n_pixels):
        x = pixel % geom.n_channels
        y = pixel // geom.n_channels
        maps_x[y, x] = fd_dOmega(geom, "tilt_x", 1e-4, pixel_pairs=[(pixel, pixel)])[
            pixel, pixel
        ]
        maps_y[y, x] = fd_dOmega(geom, "tilt_y", 1e-4, pixel_pairs=[(pixel, pixel)])[
            pixel, pixel
        ]
    sigma = 0.005
    contribution = np.sqrt(maps_x**2 + maps_y**2) * sigma
    center = (geom.n_channels // 2, geom.n_channels // 2)
    corner = (0, 0)
    return pd.DataFrame(
        [
            {
                "field": "alpha_legacy_parallel",
                "sigma_input": sigma,
                "sigmaOmega_center_sr": contribution[center[1], center[0]],
                "sigmaOmega_corner_sr": contribution[corner[1], corner[0]],
                "sigmaOmega_global_mean_sr": float(np.mean(contribution)),
            }
        ]
    )


def _plot_tilt_kernel_comparison(table: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    for label, group in table.groupby("pixel_pair_label"):
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        ax.plot(
            group["alpha_mrad"].to_numpy(),
            group["omega_parallel"].to_numpy(),
            marker="o",
            label="parallel kernel",
        )
        ax.plot(
            group["alpha_mrad"].to_numpy(),
            group["omega_gl"].to_numpy(),
            marker="s",
            label="general GL kernel",
        )
        ax.set_xlabel("tilt alpha_x [mrad]")
        ax.set_ylabel("Omega [sr]")
        ax.set_title(f"Tilt kernel comparison: {label}")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"tilt_kernel_{label}.png", dpi=220)
        fig.savefig(outdir / f"tilt_kernel_{label}.pdf")
        plt.close(fig)

    _plot_tilt_kernel_validation(table, outdir)
    return table.loc[np.isclose(table["alpha_mrad"], 0.0)].copy()


def _fit_bias_power(abs_alpha: np.ndarray, abs_delta: np.ndarray) -> tuple[float, float]:
    valid = (abs_alpha > 0.0) & (abs_delta > 0.0)
    if np.count_nonzero(valid) < 2:
        return float("nan"), float("nan")
    exponent, intercept = np.polyfit(np.log(abs_alpha[valid]), np.log(abs_delta[valid]), 1)
    return float(np.exp(intercept)), float(exponent)


def _plot_tilt_kernel_validation(table: pd.DataFrame, outdir: Path) -> None:
    """Plot kernel tilt curves and bias growth.

    The bottom-panel uncertainty scale is not a Monte Carlo error bar. It is
    the absolute model difference between the general-orientation kernel and
    the parallel approximation, evaluated for the same tilted geometry.
    """

    styles = {
        "center": "-",
        "edge": "--",
        "corner_to_corner": "-.",
    }
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(7.2, 7.0), sharex=False, constrained_layout=True
    )
    for label, group in table.groupby("pixel_pair_label"):
        group = group.sort_values("alpha_mrad")
        style = styles.get(label, "-")
        ax_top.plot(
            group["alpha_mrad"].to_numpy(dtype=float),
            group["omega_gl"].to_numpy(dtype=float),
            linestyle=style,
            marker="o",
            label=label.replace("_", " "),
        )
        grouped = (
            group.assign(
                abs_alpha=lambda frame: np.abs(frame["alpha_mrad"]),
                abs_delta=lambda frame: np.abs(
                    frame["omega_gl"] - frame["omega_parallel"]
                ),
            )
            .groupby("abs_alpha", as_index=False)["abs_delta"]
            .mean()
        )
        nonzero = grouped[grouped["abs_alpha"] > 0.0]
        ax_bot.loglog(
            nonzero["abs_alpha"].to_numpy(dtype=float),
            nonzero["abs_delta"].to_numpy(dtype=float),
            linestyle=style,
            marker="o",
            label=label.replace("_", " "),
        )

    center = table[table["pixel_pair_label"] == "center"].assign(
        abs_alpha=lambda frame: np.abs(frame["alpha_mrad"]),
        abs_delta=lambda frame: np.abs(frame["omega_gl"] - frame["omega_parallel"]),
    )
    fit_amp, fit_exp = _fit_bias_power(
        center["abs_alpha"].to_numpy(dtype=float),
        center["abs_delta"].to_numpy(dtype=float),
    )
    if np.isfinite(fit_amp):
        x_fit = np.geomspace(2.5, 20.0, 100)
        ax_bot.loglog(
            x_fit,
            fit_amp * x_fit**fit_exp,
            color="black",
            linestyle=":",
            label=rf"center fit $|\Delta\Omega|=c\alpha^n$, $n={fit_exp:.2f}$",
        )

    corner = table[table["pixel_pair_label"] == "corner_to_corner"].copy()
    for alpha_mrad in (5.0, 20.0):
        rows = corner[np.isclose(np.abs(corner["alpha_mrad"]), alpha_mrad)]
        if rows.empty:
            continue
        row = rows.iloc[-1]
        abs_delta = abs(float(row["omega_gl"] - row["omega_parallel"]))
        rel_pct = 100.0 * abs(float(row["rel_diff_gl_vs_parallel"]))
        ax_bot.annotate(
            f"{alpha_mrad:.0f} mrad: {rel_pct:.2f}%",
            xy=(alpha_mrad, abs_delta),
            xytext=(1.08 * alpha_mrad, 1.4 * abs_delta),
            arrowprops={"arrowstyle": "->", "linewidth": 0.8},
            fontsize=8,
        )
    rows_20 = corner[np.isclose(np.abs(corner["alpha_mrad"]), 20.0)]
    if not rows_20.empty:
        row_20 = rows_20.iloc[-1]
        quadratic_ref = abs(float(row_20["omega_gl"] - row_20["omega_parallel"])) / 16.0
        ax_bot.scatter([5.0], [quadratic_ref], marker="x", color="black", zorder=5)
        ax_bot.annotate(
            "0.04% quadratic reference",
            xy=(5.0, quadratic_ref),
            xytext=(6.0, 0.65 * quadratic_ref),
            fontsize=8,
            arrowprops={"arrowstyle": "->", "linewidth": 0.8},
        )

    ax_top.set_ylabel(r"$\Omega$ [sr]")
    ax_top.set_title(
        "General-orientation kernel: parallel-limit recovery and quadratic bias growth"
    )
    ax_top.grid(True, alpha=0.3)
    ax_top.legend()
    ax_bot.set_xlabel(r"Tilt magnitude $|\alpha_x|$ [mrad]")
    ax_bot.set_ylabel(r"$|\Omega_{\rm GL}-\Omega_{\rm parallel}|$ [sr]")
    ax_bot.grid(True, alpha=0.3, which="both")
    ax_bot.legend(fontsize=8)
    fig.savefig(outdir / "tilt_kernel_validation.png", dpi=220)
    fig.savefig(outdir / "tilt_kernel_validation.pdf")
    plt.close(fig)


def _max_tilt_bias(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, group in table.groupby("pixel_pair_label"):
        index = group["rel_diff_gl_vs_parallel"].abs().idxmax()
        row = group.loc[index]
        rows.append(
            {
                "pixel_pair_label": label,
                "max_abs_rel_bias": abs(float(row["rel_diff_gl_vs_parallel"])),
                "alpha_mrad": float(row["alpha_mrad"]),
            }
        )
    return pd.DataFrame(rows)


def _markdown_table(frame: pd.DataFrame) -> str:
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


def _agreement_table(geom: DetectorGeometry) -> pd.DataFrame:
    rows = []
    for label, (top, bot) in _representative_pairs(geom).items():
        for field, step in [
            ("D", 1e-3),
            ("dx", 1e-3),
            ("dy", 1e-3),
            ("tilt_x", 1e-4),
            ("tilt_y", 1e-4),
        ]:
            if field == "tilt_x":
                analytical = 1000.0 * dOmega_dalpha_general(
                    geom, (top, bot), "x", h=0.1, n_quad=12
                )
                fd = 1000.0 * dOmega_dalpha_general(
                    geom, (top, bot), "x", h=0.2, n_quad=12
                )
            elif field == "tilt_y":
                analytical = 1000.0 * dOmega_dalpha_general(
                    geom, (top, bot), "y", h=0.1, n_quad=12
                )
                fd = 1000.0 * dOmega_dalpha_general(
                    geom, (top, bot), "y", h=0.2, n_quad=12
                )
            else:
                analytical = _analytical_derivative(field, top, bot, geom)
                fd = fd_dOmega(geom, field, step, pixel_pairs=[(top, bot)])[top, bot]
            rel = abs(analytical - fd) / abs(fd) if fd != 0 else np.nan
            rows.append(
                {
                    "pair": label,
                    "field": field,
                    "analytical": analytical,
                    "finite_difference": fd,
                    "relative_discrepancy": rel,
                }
            )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    """Run tolerance validation and rebuild ``REPORT_phaseB.md``.

    Args:
        argv: Optional CLI arguments.

    Returns:
        Process exit code.

    Edge cases:
        The default report uses the Phase-A parallel-rectangle kernel for
        finite differences and flags tilt foreshortening as an open limitation.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outdir", type=Path, default=Path("characterization/results/tolerances")
    )
    parser.add_argument("--n-channels", type=int, default=16)
    parser.add_argument("--D-mm", type=float, default=150.0)
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)
    geom = DetectorGeometry(D=args.D_mm, n_channels=args.n_channels)

    sweeps = []
    for label, (top, bot) in _representative_pairs(geom).items():
        sweeps.append(
            _plot_sweep(
                "D", top, bot, geom, np.linspace(-3.0, 3.0, 13), args.outdir, label
            )
        )
        sweeps.append(
            _plot_sweep(
                "dx", top, bot, geom, np.linspace(-1.0, 1.0, 11), args.outdir, label
            )
        )
        sweeps.append(
            _plot_sweep(
                "dy", top, bot, geom, np.linspace(-1.0, 1.0, 11), args.outdir, label
            )
        )
        sweeps.append(
            _plot_sweep(
                "tilt_x",
                top,
                bot,
                geom,
                np.linspace(-0.01, 0.01, 11),
                args.outdir,
                label,
            )
        )
        sweeps.append(
            _plot_sweep(
                "tilt_y",
                top,
                bot,
                geom,
                np.linspace(-0.01, 0.01, 11),
                args.outdir,
                label,
            )
        )
    pd.concat(sweeps, ignore_index=True).to_csv(
        args.outdir / "tolerance_sweeps.csv", index=False
    )

    agreement = _agreement_table(geom)
    agreement.to_csv(args.outdir / "derivative_agreement.csv", index=False)
    alpha_grid = np.linspace(-20.0, 20.0, 17)
    kernel_table = compare_kernels(
        geom,
        alpha_grid,
        _representative_pairs(geom),
        axis="x",
        n_quad=12,
        n_samples=200_000,
    )
    kernel_table.to_csv(args.outdir / "tilt_kernel_comparison.csv", index=False)
    alpha_zero = _plot_tilt_kernel_comparison(kernel_table, args.outdir)
    alpha_thresholds = alpha_threshold_for_bias(kernel_table)
    max_tilt_bias = _max_tilt_bias(kernel_table)
    alpha_thresholds.to_csv(
        args.outdir / "tilt_parallel_bias_thresholds.csv", index=False
    )
    max_tilt_bias.to_csv(args.outdir / "tilt_parallel_max_bias.csv", index=False)
    maps = _diagonal_derivative_maps(geom)
    _plot_maps(maps, args.outdir)
    _, budget = _dominant_diagonal_map(geom, maps, args.outdir)
    legacy_alpha = _legacy_tilt_budget(geom)
    budget.to_csv(args.outdir / "error_budget_default_tolerances.csv", index=False)
    legacy_alpha.to_csv(args.outdir / "legacy_tilt_error_budget.csv", index=False)
    alpha_budget = budget[budget["field"] == "alpha"].iloc[0]
    legacy_row = legacy_alpha.iloc[0]
    alpha_delta_center = (
        alpha_budget["sigmaOmega_center_sr"] - legacy_row["sigmaOmega_center_sr"]
    )
    alpha_delta_global = (
        alpha_budget["sigmaOmega_global_mean_sr"]
        - legacy_row["sigmaOmega_global_mean_sr"]
    )
    alpha_percent_center = (
        100.0 * alpha_delta_center / legacy_row["sigmaOmega_center_sr"]
        if legacy_row["sigmaOmega_center_sr"] != 0
        else np.nan
    )
    alpha_percent_global = (
        100.0 * alpha_delta_global / legacy_row["sigmaOmega_global_mean_sr"]
        if legacy_row["sigmaOmega_global_mean_sr"] != 0
        else np.nan
    )

    report = args.outdir.parent / "REPORT_phaseB.md"
    report.write_text(
        "\n".join(
            [
                "# Phase B Tolerance Propagation Report",
                "",
                "## Analytical vs Finite-Difference Agreement",
                "",
                _markdown_table(agreement),
                "",
                "Figures are written under `characterization/results/tolerances/`.",
                "",
                "## Tilt Kernel Validation",
                "",
                "At `alpha = 0`, the general Gauss-Legendre kernel is required to match the parallel kernel.",
                "",
                _markdown_table(
                    alpha_zero[
                        [
                            "pixel_pair_label",
                            "alpha_mrad",
                            "omega_parallel",
                            "omega_gl",
                            "rel_diff_gl_vs_parallel",
                        ]
                    ]
                ),
                "",
                "Tilt comparison curves:",
                "",
                "- `tilt_kernel_center.png/pdf`",
                "- `tilt_kernel_edge.png/pdf`",
                "- `tilt_kernel_corner_to_corner.png/pdf`",
                "",
                "Parallel-kernel 1% bias thresholds:",
                "",
                _markdown_table(alpha_thresholds),
                "",
                "Maximum absolute parallel-kernel bias over `[-20, 20] mrad`:",
                "",
                _markdown_table(max_tilt_bias),
                "",
                "## Per-Pixel Derivative Maps",
                "",
                "- `map_dOmega_dD.png/pdf`",
                "- `map_dOmega_dtilt_x.png/pdf`",
                "- `map_dOmega_dtilt_y.png/pdf`",
                "- `map_dOmega_ddx.png/pdf`",
                "- `map_dOmega_ddy.png/pdf`",
                "",
                "The maps use diagonal top-bottom pixel pairs as the nominal pixel grid.",
                "",
                "## Error Budget with Engineering Inputs",
                "",
                "Inputs, not fit results: `sigma_D = 0.5 mm`, `sigma_alpha = 5 mrad`, `sigma_xy = 0.5 mm`.",
                "",
                _markdown_table(budget),
                "",
                "Legacy parallel-kernel tilt contribution:",
                "",
                _markdown_table(legacy_alpha),
                "",
                f"Corrected center alpha contribution differs from the previous value by {alpha_delta_center:.6g} sr ({alpha_percent_center:.3g}%).",
                f"Corrected global-mean alpha contribution differs from the previous value by {alpha_delta_global:.6g} sr ({alpha_percent_global:.3g}%).",
                "",
                "## Dominant Tolerance Map",
                "",
                "- `dominant_tolerance_map.png/pdf` identifies which tolerance dominates each diagonal pixel.",
                "",
                "## Open Issues and Known Limitations",
                "",
                "- Tilt derivatives in the error budget now use the general non-parallel rectangle kernel.",
                "- The parallel kernel remains adequate only where `tilt_parallel_bias_thresholds.csv` shows less than 1% bias over the engineering alpha range.",
                "- Full per-pixel-pair 65536-entry exact maps are intentionally avoided in the default report to keep Phase B reproducible on a laptop.",
                "- Pixel deformation `d -> d + delta_d` remains a Phase-B stretch TODO.",
                "",
                "## Suggested Phase C",
                "",
                "Propagate the tolerance covariance through the acceptance matrix `F` and into the reconstructed density estimator `rho_hat`.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    LOGGER.info("Wrote Phase-B tolerance report to %s", report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
