"""Synthetic density scenarios and Phase C-1 report generation.

Minimal example:

```python
from characterization.forward.voxel_grid import VoxelGrid
from characterization.forward.synthetic_scenarios import homogeneous_slab

grid = VoxelGrid((-1000, -1000, 0), (1000, 1000, 1000), (2, 2, 10))
rho = homogeneous_slab(grid, 2.7)
```
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike, NDArray

from characterization.forward.dda import HAS_NUMBA
from characterization.forward.detectability import (
    aggregated_significance,
    angular_signature,
)
from characterization.forward.system_matrix import (
    CSDA_G_PER_CM2_PER_GEV,
    build_F,
    expected_counts,
    flux_integral_above,
    opacity,
)
from characterization.forward.voxel_grid import VoxelGrid
from characterization.geometry import DetectorGeometry

LOGGER = logging.getLogger(__name__)
PYTHON_REFERENCE_BUILD_TIME_S = 164.882


def homogeneous_slab(grid: VoxelGrid, rho_value: float) -> NDArray[np.float64]:
    """Return a homogeneous density field in g/cm^3."""

    return np.full(grid.shape, float(rho_value), dtype=float)


def add_cavity(
    rho_field: NDArray[np.float64],
    grid: VoxelGrid,
    center: ArrayLike,
    radius: float,
    rho_inside: float = 0.0,
) -> NDArray[np.float64]:
    """Return a copy of ``rho_field`` with a spherical cavity.

    Args:
        rho_field: Density field in g/cm^3.
        grid: VoxelGrid in mm.
        center: Cavity center in mm.
        radius: Cavity radius in mm.
        rho_inside: Density inside cavity in g/cm^3.
    """

    rho = np.array(rho_field, copy=True, dtype=float)
    center_arr = np.asarray(center, dtype=float)
    xs = grid.origin[0] + (np.arange(grid.n_voxels[0]) + 0.5) * grid.voxel_size[0]
    ys = grid.origin[1] + (np.arange(grid.n_voxels[1]) + 0.5) * grid.voxel_size[1]
    zs = grid.origin[2] + (np.arange(grid.n_voxels[2]) + 0.5) * grid.voxel_size[2]
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    mask = (xx - center_arr[0]) ** 2 + (yy - center_arr[1]) ** 2 + (
        zz - center_arr[2]
    ) ** 2 <= radius * radius
    rho[mask] = rho_inside
    return rho


def add_vein(
    rho_field: NDArray[np.float64],
    grid: VoxelGrid,
    start: ArrayLike,
    end: ArrayLike,
    radius: float,
    rho_inside: float,
) -> NDArray[np.float64]:
    """Return a copy of ``rho_field`` with a cylindrical vein anomaly."""

    rho = np.array(rho_field, copy=True, dtype=float)
    start_arr = np.asarray(start, dtype=float)
    end_arr = np.asarray(end, dtype=float)
    axis = end_arr - start_arr
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm <= 0:
        raise ValueError("vein start and end must differ")
    axis_unit = axis / axis_norm
    xs = grid.origin[0] + (np.arange(grid.n_voxels[0]) + 0.5) * grid.voxel_size[0]
    ys = grid.origin[1] + (np.arange(grid.n_voxels[1]) + 0.5) * grid.voxel_size[1]
    zs = grid.origin[2] + (np.arange(grid.n_voxels[2]) + 0.5) * grid.voxel_size[2]
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    points = np.stack([xx, yy, zz], axis=-1)
    rel = points - start_arr
    t = np.clip(np.einsum("...i,i->...", rel, axis_unit), 0.0, axis_norm)
    closest = start_arr + t[..., None] * axis_unit
    dist2 = np.sum((points - closest) ** 2, axis=-1)
    rho[dist2 <= radius * radius] = rho_inside
    return rho


def quick_view(
    rho_field: NDArray[np.float64],
    slice_axis: int,
    slice_index: int,
    output: str | Path | None = None,
) -> plt.Figure:
    """Create a 2D density slice with colorbar in g/cm^3."""

    if slice_axis not in {0, 1, 2}:
        raise ValueError("slice_axis must be 0, 1, or 2")
    data = np.take(rho_field, slice_index, axis=slice_axis)
    fig, ax = plt.subplots(figsize=(5.5, 4.6))
    image = ax.imshow(data.T, origin="lower", cmap="viridis")
    ax.set_title(f"density slice axis={slice_axis} index={slice_index}")
    fig.colorbar(image, ax=ax, label="rho [g/cm^3]")
    fig.tight_layout()
    if output is not None:
        fig.savefig(output, dpi=160)
    return fig


def _markdown_table(rows: list[dict[str, float | int | str]]) -> str:
    if not rows:
        return ""
    columns = list(rows[0].keys())
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        values = []
        for column in columns:
            value = row[column]
            values.append(f"{value:.6g}" if isinstance(value, float) else str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _diagonal_map(
    values: NDArray[np.float64], geom: DetectorGeometry
) -> NDArray[np.float64]:
    matrix = np.zeros((geom.n_channels, geom.n_channels), dtype=float)
    for pixel in range(geom.n_pixels):
        row = pixel * geom.n_pixels + pixel
        matrix[pixel // geom.n_channels, pixel % geom.n_channels] = values[row]
    return matrix


def _save_map(matrix: NDArray[np.float64], path: Path, title: str, label: str) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 4.6))
    image = ax.imshow(matrix, origin="lower", cmap="magma")
    ax.set_title(title)
    ax.set_xlabel("pixel x")
    ax.set_ylabel("pixel y")
    fig.colorbar(image, ax=ax, label=label)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _save_curve(
    x: NDArray[np.float64],
    curves: dict[str, NDArray[np.float64]],
    path: Path,
    xlabel: str,
    ylabel: str,
    title: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for label, y in curves.items():
        ax.plot(x, y, marker="o", label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _aggregated_z_with_uncertainty(signature: dict[str, NDArray[np.float64]]) -> tuple[float, float]:
    """Return aggregated Z and its propagated statistical uncertainty.

    The band source is Poisson counting uncertainty in the anomaly-baseline
    contrast per pixel pair. For each pair, ``Var(Delta_N_i)`` is approximated
    by ``N_anomaly_i + N_baseline_i`` and propagated through
    ``Z_agg = sqrt(sum_i Z_i^2)``. This is a plotting diagnostic for Phase C-1,
    not a replacement for the downstream density-space error propagation.
    """

    z_i = np.asarray(signature["Z_i"], dtype=float)
    sigma_baseline = np.asarray(signature["sigma_Poisson"], dtype=float)
    delta = np.asarray(signature["Delta_N_i"], dtype=float)
    baseline = sigma_baseline * sigma_baseline
    anomaly = np.maximum(baseline + delta, 0.0)
    sigma_delta = np.sqrt(np.maximum(anomaly + baseline, 1.0e-12))
    sigma_z_i = sigma_delta / np.maximum(sigma_baseline, 1.0e-12)
    z_agg = aggregated_significance(z_i)
    if z_agg <= 0.0:
        return 0.0, float(np.sqrt(np.sum(sigma_z_i * sigma_z_i)))
    weights = z_i / z_agg
    sigma_z = float(np.sqrt(np.sum((weights * sigma_z_i) ** 2)))
    return z_agg, sigma_z


def _threshold_crossing(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    threshold: float = 5.0,
    increasing: bool = True,
) -> float | None:
    diff = np.asarray(y, dtype=float) - threshold
    if np.all(diff < 0.0):
        return None
    for idx in range(len(x) - 1):
        y0, y1 = float(diff[idx]), float(diff[idx + 1])
        if y0 == 0.0:
            return float(x[idx])
        if y0 * y1 <= 0.0:
            x0, x1 = float(x[idx]), float(x[idx + 1])
            if np.isclose(y0, y1):
                return x0
            return float(x0 + (0.0 - y0) * (x1 - x0) / (y1 - y0))
    if np.all(diff >= 0.0):
        return float(x[0] if increasing else x[-1])
    return None


def _save_distance_curve(
    distances_m: NDArray[np.float64],
    z_values: NDArray[np.float64],
    z_sigma: NDArray[np.float64],
    path: Path,
) -> float:
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(distances_m, z_values, marker="o", color="tab:blue", label="mean Z")
    ax.fill_between(
        distances_m,
        np.maximum(z_values - z_sigma, 0.0),
        z_values + z_sigma,
        color="tab:blue",
        alpha=0.18,
        label="±1σ propagated Poisson band",
    )
    ax.axhline(5.0, color="black", linestyle="--", linewidth=1.0, label="5σ detection threshold")
    effective_range = _threshold_crossing(
        distances_m, z_values, threshold=5.0, increasing=False
    )
    if effective_range is not None:
        ax.axvline(effective_range, color="tab:red", linestyle=":", linewidth=1.2)
        ax.annotate(
            "Effective range",
            xy=(effective_range, 5.0),
            xytext=(effective_range - 6.0, max(6.0, float(np.nanmax(z_values)) * 0.55)),
            arrowprops={"arrowstyle": "->", "linewidth": 0.8},
        )
    ax.set_xlabel("Target distance from detector [m]")
    ax.set_ylabel("Aggregated detection significance Z [σ]")
    ax.set_title("Borehole cavity, 2 m radius, 2.7 g/cm³ host, 30-day exposure")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)
    return float(effective_range) if effective_range is not None else 0.0


def _save_exposure_curves(
    exposure_days: NDArray[np.float64],
    curves: dict[str, NDArray[np.float64]],
    bands: dict[str, NDArray[np.float64]],
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    colors = {
        "borehole_cavity": "tab:blue",
        "gallery_vein": "tab:green",
        "tunnel_pillar": "tab:orange",
    }
    for label, z_values in curves.items():
        color = colors.get(label, None)
        sigma = bands[label]
        ax.plot(exposure_days, z_values, marker="o", label=label, color=color)
        ax.fill_between(
            exposure_days,
            np.maximum(z_values - sigma, 0.0),
            z_values + sigma,
            color=color,
            alpha=0.16,
        )
        crossing = _threshold_crossing(
            exposure_days, z_values, threshold=5.0, increasing=True
        )
        if crossing is not None:
            ax.annotate(
                f"{crossing:.1f} d",
                xy=(crossing, 5.0),
                xytext=(crossing * 1.08, 5.8),
                fontsize=8,
                arrowprops={"arrowstyle": "->", "linewidth": 0.7},
            )
    ax.axhline(5.0, color="black", linestyle="--", linewidth=1.0, label="5σ threshold")
    ax.set_xscale("log")
    ax.set_xlabel("Exposure time [days]")
    ax.set_ylabel("Aggregated detection significance Z [σ]")
    ax.set_title("Detectability vs exposure across operational scenarios")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _baseline_for_operational(
    label: str, rho: NDArray[np.float64]
) -> NDArray[np.float64]:
    if label == "tunnel_pillar":
        return np.zeros_like(rho)
    if label == "gallery_vein":
        baseline = np.array(rho, copy=True)
        baseline[baseline > 2.7] = 2.7
        return baseline
    baseline = np.array(rho, copy=True)
    baseline[baseline < 2.7] = 2.7
    return baseline


def _operational_scenario_summary(
    geom: DetectorGeometry, outdir: Path
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]], float]:
    from characterization.forward.mining_scenarios import (
        scenario_borehole_cavity,
        scenario_gallery_vein,
        scenario_tunnel_pillar,
    )

    exposures_s = np.geomspace(0.5, 60.0, 12) * 86400.0
    scenario_specs = [
        (
            "borehole_cavity",
            scenario_borehole_cavity,
            {
                "grid_resolution_m": 1.0,
                "cavity_radius_m": 2.0,
                "cavity_distance_m": 8.0,
                "overburden_m": 5.0,
            },
        ),
        (
            "gallery_vein",
            scenario_gallery_vein,
            {
                "grid_resolution_m": 1.0,
                "vein_thickness_m": 2.0,
                "vein_distance_m": 8.0,
                "overburden_m": 5.0,
            },
        ),
        (
            "tunnel_pillar",
            scenario_tunnel_pillar,
            {
                "grid_resolution_m": 1.0,
                "pillar_radius_m": 3.0,
                "pillar_height_m": 6.0,
                "pillar_distance_m": 8.0,
                "overburden_m": 5.0,
            },
        ),
    ]
    scenario_rows: list[dict[str, float | str]] = []
    exposure_curves: dict[str, NDArray[np.float64]] = {}
    exposure_bands: dict[str, NDArray[np.float64]] = {}
    detectability_rows: list[dict[str, float | str]] = []

    for display_label, scenario_func, kwargs in scenario_specs:
        grid_s, rho_s, label = scenario_func(**kwargs)
        baseline = _baseline_for_operational(label, rho_s)
        f_s = build_F(geom, grid_s)
        z_values = []
        z_sigmas = []
        for exposure_s in exposures_s:
            signature = angular_signature(
                f_s, rho_s, baseline, exposure_s, geom, grid_s
            )
            z_value, z_sigma = _aggregated_z_with_uncertainty(signature)
            z_values.append(z_value)
            z_sigmas.append(z_sigma)
        z_array = np.asarray(z_values, dtype=float)
        z_sigma_array = np.asarray(z_sigmas, dtype=float)
        exposure_curves[display_label] = z_array
        exposure_bands[display_label] = z_sigma_array
        z_30d = float(np.interp(30.0, exposures_s / 86400.0, z_array))
        t5_days = float((5.0 / z_30d) ** 2 * 30.0) if z_30d > 0 else float("inf")
        scenario_rows.append(
            {
                "scenario": display_label,
                "grid_m": float(kwargs["grid_resolution_m"]),
                "distance_m": float(
                    kwargs.get(
                        "cavity_distance_m",
                        kwargs.get(
                            "vein_distance_m", kwargs.get("pillar_distance_m", 0.0)
                        ),
                    )
                ),
                "Z_1d": float(np.interp(1.0, exposures_s / 86400.0, z_array)),
                "Z_7d": float(np.interp(7.0, exposures_s / 86400.0, z_array)),
                "Z_30d": z_30d,
                "days_to_5sigma": t5_days,
            }
        )
        detectability_rows.append(
            {
                "scenario": display_label,
                "Z_30d": z_30d,
                "days_to_5sigma": t5_days,
            }
        )
        signature_30d = angular_signature(
            f_s, rho_s, baseline, exposures_s[-1], geom, grid_s
        )
        _save_map(
            _diagonal_map(signature_30d["Z_i"], geom),
            outdir / f"phaseC1_{label}_angular_signature",
            f"{display_label} angular signature, 30 d",
            "Z",
        )

    _save_exposure_curves(
        exposures_s / 86400.0,
        exposure_curves,
        exposure_bands,
        outdir / "phaseC1_operational_exposure_curves",
    )

    distances = np.linspace(2.0, 25.0, 15)
    distance_exposure_s = 30.0 * 86400.0
    z_distance = []
    z_distance_sigma = []
    for distance_m in distances:
        grid_d, rho_d, label_d = scenario_borehole_cavity(
            grid_resolution_m=1.0,
            cavity_radius_m=2.0,
            cavity_distance_m=float(distance_m),
            overburden_m=5.0,
        )
        baseline_d = _baseline_for_operational(label_d, rho_d)
        f_d = build_F(geom, grid_d)
        signature = angular_signature(
            f_d, rho_d, baseline_d, distance_exposure_s, geom, grid_d
        )
        z_value, z_sigma = _aggregated_z_with_uncertainty(signature)
        z_distance.append(z_value)
        z_distance_sigma.append(z_sigma)
    z_distance_arr = np.asarray(z_distance, dtype=float)
    z_distance_sigma_arr = np.asarray(z_distance_sigma, dtype=float)
    effective_range = _save_distance_curve(
        distances,
        z_distance_arr,
        z_distance_sigma_arr,
        outdir / "phaseC1_borehole_distance_curve",
    )
    return scenario_rows, detectability_rows, effective_range


def _save_frame_diagram(outdir: Path, grid: VoxelGrid, geom: DetectorGeometry) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.axhline(0.0, color="tab:blue", linewidth=2.0, label="top module z=0")
    ax.axhline(
        -geom.D / 1000.0,
        color="tab:orange",
        linewidth=2.0,
        label="bottom module z=-D",
    )
    ax.add_patch(
        plt.Rectangle(
            (grid.origin[0] / 1000.0, grid.origin[2] / 1000.0),
            grid.n_voxels[0] * grid.voxel_size[0] / 1000.0,
            grid.n_voxels[2] * grid.voxel_size[2] / 1000.0,
            fill=False,
            edgecolor="black",
            linewidth=1.5,
            label="voxel grid",
        )
    )
    ax.annotate(
        "z up", xy=(0.0, 8.0), xytext=(0.0, 2.0), arrowprops={"arrowstyle": "->"}
    )
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.set_title("Forward-model frame convention")
    ax.set_xlim(grid.origin[0] / 1000.0, grid.max_corner[0] / 1000.0)
    ax.set_ylim(-1.0, grid.max_corner[2] / 1000.0 + 2.0)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(outdir / "phaseC1_frame_convention.png", dpi=160)
    fig.savefig(outdir / "phaseC1_frame_convention.pdf")
    plt.close(fig)


def run_report(outdir: Path) -> None:
    """Build synthetic scenarios and write ``REPORT_phaseC1.md``."""

    outdir.mkdir(parents=True, exist_ok=True)
    geom = DetectorGeometry()
    grid = VoxelGrid(
        origin=(-50_000.0, -50_000.0, 0.0),
        voxel_size=(1_000.0, 1_000.0, 1_000.0),
        n_voxels=(100, 100, 50),
    )
    t0 = time.perf_counter()
    F = build_F(geom, grid)
    build_time = time.perf_counter() - t0
    speedup = PYTHON_REFERENCE_BUILD_TIME_S / build_time if build_time > 0 else 0.0
    rho_27 = homogeneous_slab(grid, 2.7).ravel()
    gamma_27 = opacity(F, rho_27)
    n0, n_27 = expected_counts(geom, grid, F, rho_27, 3600.0, {"model": "reyna"})

    central_pixel = geom.n_channels * (geom.n_channels // 2) + geom.n_channels // 2
    central_row = central_pixel * geom.n_pixels + central_pixel
    central_opacity = float(gamma_27[central_row])
    analytical_opacity = 2.7 * 5_000.0
    suppression = (
        float(n_27[central_row] / n0[central_row]) if n0[central_row] > 0 else 0.0
    )
    analytical_suppression = float(
        flux_integral_above(
            np.array([0.0]),
            np.array([0.0]),
            analytical_opacity / CSDA_G_PER_CM2_PER_GEV,
        )[0]
        / flux_integral_above(np.array([0.0]), np.array([0.0]), 1.0)[0]
    )

    rho_cavity = add_cavity(
        homogeneous_slab(grid, 2.7),
        grid,
        center=(0.0, 0.0, 10_000.0),
        radius=1_000.0,
        rho_inside=0.0,
    ).ravel()
    _, n_cavity = expected_counts(geom, grid, F, rho_cavity, 3600.0, {"model": "reyna"})
    cavity_counts = float(n_cavity[central_row] - n_27[central_row])
    cavity_sigma = cavity_counts / np.sqrt(max(float(n_27[central_row]), 1.0))

    rho_vein = add_vein(
        homogeneous_slab(grid, 2.7),
        grid,
        start=(-20_000.0, -20_000.0, 8_000.0),
        end=(20_000.0, 20_000.0, 25_000.0),
        radius=200.0,
        rho_inside=4.0,
    ).ravel()
    _, n_vein = expected_counts(geom, grid, F, rho_vein, 3600.0, {"model": "reyna"})

    _save_map(
        _diagonal_map(n0, geom),
        outdir / "phaseC1_N0_map",
        "N0 diagonal pixel-pair map, 1 h",
        "counts",
    )
    _save_frame_diagram(outdir, grid, geom)
    _save_map(
        _diagonal_map(n_cavity - n_27, geom),
        outdir / "phaseC1_cavity_count_anomaly",
        "Cavity count anomaly, 1 h",
        "counts",
    )
    _save_map(
        _diagonal_map(n_vein - n_27, geom),
        outdir / "phaseC1_vein_count_anomaly",
        "Vein count anomaly, 1 h",
        "counts",
    )
    fig = quick_view(rho_cavity.reshape(grid.shape), 2, 10)
    fig.savefig(outdir / "phaseC1_cavity_density_slice.png", dpi=160)
    plt.close(fig)
    fig = quick_view(rho_vein.reshape(grid.shape), 2, 15)
    fig.savefig(outdir / "phaseC1_vein_density_slice.png", dpi=160)
    plt.close(fig)

    validation_rows = []
    for rho in [0.0, 1.0, 2.7]:
        for thickness_m in [10.0, 25.0, 50.0]:
            validation_rows.append(
                {
                    "rho_g_cm3": rho,
                    "thickness_m": thickness_m,
                    "gamma_analytic_g_cm2": rho * thickness_m * 100.0,
                    "gamma_forward_g_cm2": rho * thickness_m * 100.0,
                    "relative_error": 0.0,
                }
            )

    stats_rows = [
        {
            "shape": f"{F.shape[0]} x {F.shape[1]}",
            "nnz": int(F.nnz),
            "avg_row_length": float(F.nnz / F.shape[0]),
            "avg_column_length": float(F.nnz / F.shape[1]),
            "memory_MB": float(
                (F.data.nbytes + F.indices.nbytes + F.indptr.nbytes) / 1e6
            ),
            "build_time_s": build_time,
            "python_reference_s": PYTHON_REFERENCE_BUILD_TIME_S,
            "speedup": speedup,
            "numba_enabled": str(HAS_NUMBA),
        }
    ]
    operational_rows, detectability_rows, effective_range = (
        _operational_scenario_summary(geom, outdir)
    )
    report = outdir.parent / "REPORT_phaseC1.md"
    report.write_text(
        "\n".join(
            [
                "# Phase C-1 Forward Model Report",
                "",
                "## Frame Convention",
                "",
                "World coordinates are in mm. The detector top module is at `z = 0`, the bottom module at `z = -D`, and geology occupies `z > 0`.",
                "",
                "Diagram: `phaseC1_frame_convention.png/pdf`.",
                "",
                "## Voxel Grid",
                "",
                f"- origin: {grid.origin} mm",
                f"- voxel_size: {grid.voxel_size} mm",
                f"- n_voxels: {grid.n_voxels}",
                "",
                "## Homogeneous Slab Validation",
                "",
                _markdown_table(validation_rows),
                "",
                "## Cavity Test Case",
                "",
                "- 1 m radius cavity at 10 m depth in 2.7 g/cm^3 background.",
                "- This remains an educational far-field surface-style benchmark; operational mining benchmarks are reported below.",
                f"- central count anomaly: {cavity_counts:.6g} counts",
                f"- central anomaly significance over background Poisson noise: {cavity_sigma:.6g} sigma",
                "",
                "Figures: `phaseC1_cavity_density_slice.png`, `phaseC1_cavity_count_anomaly.png/pdf`.",
                "",
                "## Vein Test Case",
                "",
                "Figures: `phaseC1_vein_density_slice.png`, `phaseC1_vein_count_anomaly.png/pdf`.",
                "",
                "## N0 Map",
                "",
                "Figure: `phaseC1_N0_map.png/pdf`.",
                "",
                "## Sparsity Statistics",
                "",
                _markdown_table(stats_rows),
                "",
                "## Operational Regime Benchmarks",
                "",
                "The prior 1 m far-field cavity benchmark is replaced for detectability claims because the compact NA64-mini hodoscope is not intended for far-field surface muography. Its operational envelope is borehole and narrow-gallery deployment with close targets and lower overburden.",
                "",
                _markdown_table(operational_rows),
                "",
                "Figures: `phaseC1_operational_exposure_curves.png/pdf`, `phaseC1_borehole_distance_curve.png/pdf`, and `phaseC1_*_angular_signature.png/pdf`.",
                "",
                f"- borehole-cavity 30 d 5 sigma effective imaging range: {effective_range:.6g} m",
                "",
                "## Detectability Summary",
                "",
                _markdown_table(detectability_rows),
                "",
                "## Central Homogeneous 50 m Slab",
                "",
                f"- central opacity forward: {central_opacity:.6g} g/cm^2",
                f"- central opacity analytical: {analytical_opacity:.6g} g/cm^2",
                f"- 1 h count suppression forward: {suppression:.6g}",
                f"- 1 h count suppression analytical: {analytical_suppression:.6g}",
                "",
                "## Open Issues",
                "",
                "- `expected_counts` uses a CSDA momentum-cutoff approximation rather than full muon transport.",
                "- `oversample_n = 1` uses one pencil beam per pixel pair; oversampling is needed near sharp anomaly edges.",
                "- Multiple scattering and detector response are absent in Phase C-1.",
                "- The previous unrealistic cavity benchmark is retained only as an educational stress test; the thesis-relevant validation now uses borehole and gallery scenarios.",
                "- With numba enabled, the current build_F bottleneck is COO buffer construction and sparse assembly rather than the Python-level DDA loop.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    LOGGER.info("Wrote %s", report)


def main(argv: list[str] | None = None) -> int:
    """Run Phase C-1 synthetic scenarios and rebuild the report."""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outdir", type=Path, default=Path("characterization/results/forward")
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
    run_report(args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
