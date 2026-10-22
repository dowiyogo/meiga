"""Aggregate detector figures of merit without changing upstream physics.

This module reads existing Phase A/B/C outputs and optional serialized forward
or reconstruction artifacts. Missing inputs are warnings, not fatal errors.

Minimal example:

```bash
python -m characterization.figures_of_merit.run_fom
```
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from characterization.geometry import DetectorGeometry

LOGGER = logging.getLogger(__name__)
TOLERANCE_FIELDS = ("D", "alpha_x", "alpha_y", "dx", "dy")


@dataclass
class Metric:
    """One long-form summary value with explicit provenance and maturity."""

    section: str
    scope: str
    metric: str
    value: Any
    unit: str
    status: str
    source: str
    note: str = ""


class FoMCollector:
    """Collect metrics and non-fatal warnings while reading upstream outputs."""

    def __init__(self, root: Path | None = None) -> None:
        self.root = root
        self.metrics: list[Metric] = []
        self.warnings: list[str] = []

    def warn(self, message: str) -> None:
        LOGGER.warning(message)
        self.warnings.append(message)

    def add(
        self,
        section: str,
        scope: str,
        metric: str,
        value: Any,
        unit: str,
        status: str,
        source: str | Path,
        note: str = "",
    ) -> None:
        source_path = Path(source) if source else None
        if source_path is not None and self.root is not None:
            try:
                source = source_path.resolve().relative_to(self.root.resolve())
            except ValueError:
                source = source_path
        self.metrics.append(
            Metric(
                section=section,
                scope=scope,
                metric=metric,
                value=_json_value(value),
                unit=unit,
                status=status,
                source=str(source),
                note=note,
            )
        )


def _json_value(value: Any) -> Any:
    """Return JSON-safe scalar values, representing non-finite values as null."""

    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _search_roots(root: Path) -> list[Path]:
    candidates = [root / "data" / "results", root / "characterization" / "results"]
    candidates.extend(sorted((root / "characterization").glob("**/outputs")))
    unique: list[Path] = []
    for candidate in candidates:
        if candidate.exists() and candidate not in unique:
            unique.append(candidate)
    return unique


def _find_first(roots: Iterable[Path], names: Iterable[str]) -> Path | None:
    for directory in roots:
        for name in names:
            matches = sorted(directory.rglob(name))
            if matches:
                return matches[0]
    return None


def _load_yaml_config(root: Path, collector: FoMCollector) -> dict[str, Any]:
    path = root / "configs" / "thesis_scenario.yaml"
    if not path.exists():
        collector.warn(
            "Optional configs/thesis_scenario.yaml is absent; using DetectorGeometry "
            "defaults only for documented detector constants."
        )
        return {}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        collector.warn(f"Could not read {path}: {exc}")
        return {}
    return data if isinstance(data, dict) else {}


def _find_nested(config: dict[str, Any], keys: Iterable[str]) -> Any | None:
    wanted = {key.lower() for key in keys}
    stack: list[Any] = [config]
    while stack:
        item = stack.pop()
        if not isinstance(item, dict):
            continue
        for key, value in item.items():
            if str(key).lower() in wanted:
                return value
            stack.append(value)
    return None


def _geometry_from_config(config: dict[str, Any]) -> DetectorGeometry:
    defaults = DetectorGeometry()

    def number(keys: Iterable[str], default: float) -> float:
        value = _find_nested(config, keys)
        try:
            return float(value) if value is not None else default
        except (TypeError, ValueError):
            return default

    return DetectorGeometry(
        D=number(("D", "inter_plane_distance_mm", "detector_separation_mm"), defaults.D),
        bar_length=number(("bar_length", "bar_length_mm"), defaults.bar_length),
        bar_width=number(("bar_width", "bar_width_mm"), defaults.bar_width),
        bar_thickness=number(("bar_thickness", "bar_thickness_mm"), defaults.bar_thickness),
        n_channels=int(number(("n_channels", "channels_per_axis"), defaults.n_channels)),
        pixel_pitch=number(("pixel_pitch", "pixel_pitch_mm"), defaults.pixel_pitch),
        sub_plane_offset=number(
            ("sub_plane_offset", "sub_plane_offset_mm"), defaults.sub_plane_offset
        ),
        module_size=defaults.module_size,
    )


def summarize_fractional_uncertainty(
    omega: np.ndarray, sigma_omega: np.ndarray
) -> dict[str, float]:
    """Summarize ``sigma_omega / omega`` where both inputs are finite and positive."""

    omega_arr = np.asarray(omega, dtype=float)
    sigma_arr = np.asarray(sigma_omega, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        fractional = sigma_arr / omega_arr
    valid = np.isfinite(fractional) & (omega_arr > 0.0) & (sigma_arr >= 0.0)
    values = fractional[valid]
    if values.size == 0:
        return {name: float("nan") for name in ("median", "p90", "p95", "max")}
    return {
        "median": float(np.percentile(values, 50.0)),
        "p90": float(np.percentile(values, 90.0)),
        "p95": float(np.percentile(values, 95.0)),
        "max": float(np.max(values)),
    }


def robustness_index(omega: np.ndarray, sigma_omega: np.ndarray) -> float:
    """Return median ``omega / sigma_omega`` over finite strictly-positive pairs.

    Pairs with zero omega or zero sigma are omitted because their ratio is not a
    finite comparison value. This keeps the index stable when sparse or
    placeholder arrays contain zeros.
    """

    omega_arr = np.asarray(omega, dtype=float)
    sigma_arr = np.asarray(sigma_omega, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = omega_arr / sigma_arr
    valid = (
        np.isfinite(ratio)
        & np.isfinite(omega_arr)
        & np.isfinite(sigma_arr)
        & (omega_arr > 0.0)
        & (sigma_arr > 0.0)
    )
    return float(np.median(ratio[valid])) if np.any(valid) else float("nan")


def compute_sparsity_stats(matrix: sparse.spmatrix) -> dict[str, Any]:
    """Return sparse forward-operator metrics from an existing matrix."""

    csr = sparse.csr_matrix(matrix)
    n_rows, n_cols = csr.shape
    total = n_rows * n_cols
    row_nnz = np.diff(csr.indptr)
    col_nnz = np.diff(csr.tocsc().indptr)
    return {
        "shape": [int(n_rows), int(n_cols)],
        "nnz": int(csr.nnz),
        "sparsity": float(1.0 - csr.nnz / total) if total else float("nan"),
        "row_nnz": row_nnz,
        "column_nnz": col_nnz,
        "row_nnz_mean": float(np.mean(row_nnz)) if row_nnz.size else float("nan"),
        "row_nnz_median": float(np.median(row_nnz)) if row_nnz.size else float("nan"),
        "row_nnz_p95": float(np.percentile(row_nnz, 95.0)) if row_nnz.size else float("nan"),
        "column_nnz_mean": float(np.mean(col_nnz)) if col_nnz.size else float("nan"),
        "column_nnz_median": float(np.median(col_nnz)) if col_nnz.size else float("nan"),
        "column_nnz_p95": (
            float(np.percentile(col_nnz, 95.0)) if col_nnz.size else float("nan")
        ),
        "voxel_coverage_fraction": (
            float(np.count_nonzero(col_nnz) / n_cols) if n_cols else float("nan")
        ),
        "track_coverage_fraction": (
            float(np.count_nonzero(row_nnz) / n_rows) if n_rows else float("nan")
        ),
    }


def _read_csv(path: Path | None, collector: FoMCollector, label: str) -> pd.DataFrame:
    if path is None:
        collector.warn(f"Missing {label}; related figures of merit were skipped.")
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.ParserError) as exc:
        collector.warn(f"Could not read {path}: {exc}")
        return pd.DataFrame()


def _markdown_tables(path: Path) -> dict[str, pd.DataFrame]:
    """Parse simple Markdown tables and map them to their preceding heading."""

    lines = path.read_text(encoding="utf-8").splitlines()
    heading = ""
    tables: dict[str, pd.DataFrame] = {}
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if line.startswith("#"):
            heading = line.lstrip("#").strip()
        if (
            line.startswith("|")
            and index + 1 < len(lines)
            and re.match(r"^\|\s*:?-+", lines[index + 1].strip())
        ):
            columns = [cell.strip() for cell in line.strip("|").split("|")]
            rows: list[list[str]] = []
            index += 2
            while index < len(lines) and lines[index].strip().startswith("|"):
                rows.append(
                    [cell.strip() for cell in lines[index].strip().strip("|").split("|")]
                )
                index += 1
            tables[heading] = pd.DataFrame(rows, columns=columns)
            continue
        index += 1
    return tables


def _numeric(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    copy = frame.copy()
    for column in columns:
        if column in copy:
            copy[column] = pd.to_numeric(copy[column], errors="coerce")
    return copy


def _save_figure(fig: plt.Figure, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def _add_geometric_fom(
    roots: list[Path],
    geom: DetectorGeometry,
    collector: FoMCollector,
    figures: Path,
) -> pd.DataFrame:
    path = _find_first(roots, ("scan_D_summary.csv",))
    frame = _read_csv(path, collector, "Phase A scan_D_summary.csv")
    if frame.empty or path is None:
        return pd.DataFrame()
    frame = _numeric(
        frame,
        ("D", "mean_Omega", "mean_T", "sigma_theta_geom", "sigma_theta_rec"),
    )
    pairs = geom.n_pixels * geom.n_pixels
    rows: list[dict[str, Any]] = []
    for _, item in frame.iterrows():
        D = float(item["D"])
        mean_omega = float(item["mean_Omega"])
        mean_t = float(item["mean_T"]) if "mean_T" in frame else float("nan")
        scope = f"D={D:g} mm"
        row = {
            "D_mm": D,
            "mean_omega_sr": mean_omega,
            "median_omega_sr": None,
            "sum_omega_sr": mean_omega * pairs,
            "mean_T_dimensionless": _json_value(mean_t),
            "median_T_dimensionless": None,
            "effective_acceptance_sum_omega_sr": mean_omega * pairs,
            "angular_pitch_mrad": geom.pixel_pitch / D * 1000.0,
            "sigma_theta_proxy_mrad": (
                float(item["sigma_theta_geom"]) * 1000.0
                if "sigma_theta_geom" in frame
                else geom.pixel_pitch / (np.sqrt(6.0) * D) * 1000.0
            ),
        }
        rows.append(row)
        for metric, value, unit in (
            ("mean_omega", row["mean_omega_sr"], "sr"),
            ("median_omega", None, "sr"),
            ("sum_omega", row["sum_omega_sr"], "sr"),
            ("mean_T", row["mean_T_dimensionless"], "dimensionless"),
            ("median_T", None, "dimensionless"),
            ("effective_acceptance_sum_omega", row["effective_acceptance_sum_omega_sr"], "sr"),
            ("angular_pitch", row["angular_pitch_mrad"], "mrad"),
            ("sigma_theta_proxy", row["sigma_theta_proxy_mrad"], "mrad"),
        ):
            note = ""
            if metric == "effective_acceptance_sum_omega":
                note = "sum(omega); no detector area factor is included"
            if value is None:
                note = "not stored by current Phase A outputs"
            collector.add("geometric", scope, metric, value, unit, "SIMULATION", path, note)
    table = pd.DataFrame(rows)
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.8), sharex=True)
    axes[0].plot(
        table["D_mm"].to_numpy(dtype=float),
        table["mean_omega_sr"].to_numpy(dtype=float),
        marker="o",
        color="tab:blue",
    )
    axes[0].set_ylabel("Mean per-pair Omega [sr]")
    axes[0].set_yscale("log")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(
        table["D_mm"].to_numpy(dtype=float),
        table["angular_pitch_mrad"].to_numpy(dtype=float),
        marker="o",
        color="tab:orange",
    )
    axes[1].set_xlabel("Inter-plane separation D [mm]")
    axes[1].set_ylabel("Angular pitch d / D [mrad]")
    axes[1].grid(True, alpha=0.3)
    axes[2].text(
        0.5,
        0.5,
        "robustness_index unavailable:\nper-pair sigma_omega is not serialized",
        ha="center",
        va="center",
        transform=axes[2].transAxes,
    )
    axes[2].set_xlabel("Inter-plane separation D [mm]")
    axes[2].set_ylabel("robustness_index")
    axes[2].grid(True, alpha=0.3)
    _save_figure(fig, figures / "fom_vs_D")
    return table


def _add_robustness_fom(
    roots: list[Path], collector: FoMCollector, figures: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pair_path = _find_first(
        roots,
        (
            "per_pair_tolerance_budget.csv",
            "tolerance_budget_per_pair.csv",
            "fractional_omega_uncertainty.csv",
        ),
    )
    pair_table = _read_csv(pair_path, collector, "per-pair tolerance budget CSV")
    summary_rows: list[dict[str, Any]] = []
    contributions = pd.DataFrame()
    if not pair_table.empty and {"omega_sr", "sigma_omega_sr"} <= set(pair_table.columns):
        omega = pair_table["omega_sr"].to_numpy(dtype=float)
        sigma = pair_table["sigma_omega_sr"].to_numpy(dtype=float)
        stats = summarize_fractional_uncertainty(omega, sigma)
        stats["robustness_index"] = robustness_index(omega, sigma)
        for metric, value in stats.items():
            collector.add(
                "robustness",
                "all pixel pairs",
                metric if metric == "robustness_index" else f"fractional_omega_uncertainty_{metric}",
                value,
                "dimensionless",
                "SIMULATION",
                pair_path or "",
            )
        summary_rows.append(stats)
        split_rows = []
        for field in TOLERANCE_FIELDS:
            candidates = (
                f"sigma_omega_{field}_sr",
                f"sigmaOmega_{field}_sr",
                f"contribution_{field}_sr",
            )
            column = next((name for name in candidates if name in pair_table), None)
            if column is None:
                continue
            values = pd.to_numeric(pair_table[column], errors="coerce").to_numpy(dtype=float)
            split_rows.append(
                {
                    "field": field,
                    "variance_sum": float(np.nansum(values**2)),
                    "source_resolution": "split",
                }
            )
        if len(split_rows) == len(TOLERANCE_FIELDS):
            contributions = pd.DataFrame(split_rows)
            total = float(contributions["variance_sum"].sum())
            contributions["relative_variance_contribution"] = (
                contributions["variance_sum"] / total if total > 0 else np.nan
            )
            for _, row in contributions.iterrows():
                collector.add(
                    "robustness",
                    "tolerance budget",
                    f"relative_variance_contribution_{row['field']}",
                    row["relative_variance_contribution"],
                    "fraction",
                    "SIMULATION",
                    pair_path or "",
                )
    else:
        collector.warn(
            "Per-pair omega and sigma_omega are unavailable; fractional uncertainty "
            "percentiles and robustness_index cannot be reported."
        )

    budget_path = _find_first(roots, ("error_budget_default_tolerances.csv",))
    budget = (
        _read_csv(budget_path, collector, "Phase B tolerance error budget CSV")
        if contributions.empty
        else pd.DataFrame()
    )
    if not budget.empty and budget_path is not None:
        budget = _numeric(budget, ("sigmaOmega_global_mean_sr",))
        variance = budget["sigmaOmega_global_mean_sr"].to_numpy(dtype=float) ** 2
        total = float(np.nansum(variance))
        contributions = budget[["field"]].copy()
        contributions["relative_variance_contribution"] = variance / total if total > 0 else np.nan
        contributions["source_resolution"] = "grouped"
        for _, row in contributions.iterrows():
            collector.add(
                "robustness",
                "tolerance budget",
                f"relative_variance_contribution_{row['field']}",
                row["relative_variance_contribution"],
                "fraction",
                "SIMULATION",
                budget_path,
                "Current Phase B CSV groups alpha_x/alpha_y as alpha and dx/dy as xy.",
            )
        collector.warn(
            "Current Phase B output stores grouped alpha and xy contributions; "
            "alpha_x, alpha_y, dx and dy cannot be split without a richer upstream table."
        )
    if not contributions.empty:
        fig, ax = plt.subplots(figsize=(6.8, 4.4))
        ax.bar(
            contributions["field"].astype(str).to_numpy(),
            contributions["relative_variance_contribution"].to_numpy(dtype=float),
            color="tab:purple",
        )
        ax.set_xlabel("Tolerance source (stored grouping)")
        ax.set_ylabel("Relative contribution to variance [fraction]")
        ax.grid(True, alpha=0.3, axis="y")
        _save_figure(fig, figures / "tolerance_budget_bar")
    return pd.DataFrame(summary_rows), contributions


def _phase_c_report(roots: list[Path], collector: FoMCollector) -> tuple[Path | None, dict[str, pd.DataFrame]]:
    path = _find_first(roots, ("REPORT_phaseC1.md", "REPORT_phaseC*.md"))
    if path is None:
        collector.warn("Missing Phase C report; sensitivity and report-level forward FoMs were skipped.")
        return None, {}
    try:
        return path, _markdown_tables(path)
    except OSError as exc:
        collector.warn(f"Could not read {path}: {exc}")
        return None, {}


def _add_sensitivity_fom(
    report_path: Path | None,
    tables: dict[str, pd.DataFrame],
    collector: FoMCollector,
    figures: Path,
) -> pd.DataFrame:
    scenarios = tables.get("Operational Regime Benchmarks", pd.DataFrame())
    if scenarios.empty or report_path is None:
        collector.warn("Phase C operational scenario table is unavailable.")
        return pd.DataFrame()
    scenarios = _numeric(
        scenarios,
        ("distance_m", "Z_1d", "Z_7d", "Z_30d", "days_to_5sigma"),
    )
    rows: list[dict[str, Any]] = []
    for _, scenario in scenarios.iterrows():
        for days, column in ((1.0, "Z_1d"), (7.0, "Z_7d"), (30.0, "Z_30d")):
            zagg = float(scenario[column])
            row = {
                "scenario": str(scenario["scenario"]),
                "exposure_days": days,
                "target_distance_m": float(scenario["distance_m"]),
                "Zagg_sigma": zagg,
                "reaches_5sigma": bool(zagg >= 5.0),
                "minimum_time_to_5sigma_days": float(scenario["days_to_5sigma"]),
            }
            rows.append(row)
            for metric, value, unit in (
                ("Zagg", zagg, "sigma"),
                ("target_distance", row["target_distance_m"], "m"),
                ("reaches_5sigma", row["reaches_5sigma"], "boolean"),
                ("minimum_time_to_5sigma", row["minimum_time_to_5sigma_days"], "days"),
            ):
                collector.add(
                    "sensitivity",
                    f"{row['scenario']}, exposure={days:g} days",
                    metric,
                    value,
                    unit,
                    "ILLUSTRATIVE",
                    report_path,
                    "Synthetic Phase C-1 CSDA detectability benchmark.",
                )
    table = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for label, group in table.groupby("scenario"):
        ax.plot(
            group["exposure_days"].to_numpy(dtype=float),
            group["Zagg_sigma"].to_numpy(dtype=float),
            marker="o",
            label=label,
        )
    ax.axhline(5.0, color="black", linestyle="--", linewidth=1.0, label="5 sigma threshold")
    ax.set_xscale("log")
    ax.set_xlabel("Exposure [days]")
    ax.set_ylabel("Aggregated significance Zagg [sigma]")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend()
    _save_figure(fig, figures / "sensitivity_range_exposure")
    collector.warn(
        "Phase C does not serialize the distance sweep values. "
        "sensitivity_range_exposure therefore plots stored Zagg versus exposure; "
        "it does not infer an exposure-versus-range curve."
    )

    report_text = report_path.read_text(encoding="utf-8")
    match = re.search(r"30 d 5 sigma effective imaging range:\s*([0-9.eE+-]+)\s*m", report_text)
    if match:
        collector.add(
            "sensitivity",
            "borehole_cavity, exposure=30 days",
            "maximum_range_at_5sigma",
            float(match.group(1)),
            "m",
            "ILLUSTRATIVE",
            report_path,
            "Interpolated Phase C-1 synthetic benchmark.",
        )
    else:
        collector.warn("Phase C report does not store a maximum 5-sigma range.")
    return table


def _load_sparse_matrix(roots: list[Path], collector: FoMCollector) -> tuple[Path | None, sparse.csr_matrix | None]:
    path = _find_first(roots, ("F.npz", "forward_operator.npz", "system_matrix.npz"))
    if path is None:
        return None, None
    try:
        return path, sparse.load_npz(path).tocsr()
    except (OSError, ValueError) as exc:
        collector.warn(f"Could not load sparse forward operator {path}: {exc}")
        return None, None


def _add_forward_fom(
    roots: list[Path],
    report_path: Path | None,
    tables: dict[str, pd.DataFrame],
    collector: FoMCollector,
    figures: Path,
) -> dict[str, Any]:
    matrix_path, matrix = _load_sparse_matrix(roots, collector)
    stats: dict[str, Any] = {}
    source: Path | str = matrix_path or report_path or ""
    if matrix is not None:
        stats = compute_sparsity_stats(matrix)
    else:
        report_stats = tables.get("Sparsity Statistics", pd.DataFrame())
        if report_stats.empty:
            collector.warn("No serialized F matrix or Phase C sparsity table was found.")
            return stats
        row = report_stats.iloc[0]
        match = re.match(r"\s*(\d+)\s*x\s*(\d+)\s*", str(row["shape"]))
        if not match:
            collector.warn("Could not parse F shape from Phase C report.")
            return stats
        n_rows, n_cols = int(match.group(1)), int(match.group(2))
        nnz = int(float(row["nnz"]))
        stats = {
            "shape": [n_rows, n_cols],
            "nnz": nnz,
            "sparsity": 1.0 - nnz / (n_rows * n_cols),
            "row_nnz_mean": float(row["avg_row_length"]),
            "column_nnz_mean": float(row["avg_column_length"]),
            "build_time_s": float(row["build_time_s"]),
            "speedup": float(row["speedup"]),
        }
        collector.warn(
            "F is not serialized; p95 nnz and voxel/track coverage fractions cannot "
            "be recovered from the Phase C report alone."
        )
    units = {
        "shape": "rows x voxels",
        "nnz": "elements",
        "sparsity": "fraction",
        "row_nnz_mean": "elements/row",
        "row_nnz_median": "elements/row",
        "row_nnz_p95": "elements/row",
        "column_nnz_mean": "elements/column",
        "column_nnz_median": "elements/column",
        "column_nnz_p95": "elements/column",
        "voxel_coverage_fraction": "fraction",
        "track_coverage_fraction": "fraction",
        "build_time_s": "s",
        "speedup": "dimensionless",
    }
    for metric, unit in units.items():
        collector.add(
            "forward_operator",
            "F",
            metric,
            stats.get(metric),
            unit,
            "ILLUSTRATIVE",
            source,
            "Phase C-1 synthetic forward operator.",
        )
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    if matrix is not None:
        row_nnz = stats["row_nnz"]
        col_nnz = stats["column_nnz"]
        ax.hist(row_nnz, bins=30, alpha=0.65, label="nnz per row")
        ax.hist(col_nnz, bins=30, alpha=0.65, label="nnz per column")
        ax.set_xlabel("Non-zero elements [count]")
        ax.set_ylabel("Rows or columns [count]")
    else:
        ax.bar(
            ["mean nnz / row", "mean nnz / column"],
            [stats["row_nnz_mean"], stats["column_nnz_mean"]],
            color=["tab:blue", "tab:orange"],
        )
        ax.set_ylabel("Non-zero elements [count]")
    ax.grid(True, alpha=0.3, axis="y")
    ax.legend() if matrix is not None else None
    _save_figure(fig, figures / "forward_operator_sparsity")
    return stats


def _load_reconstruction_npz(roots: list[Path], collector: FoMCollector) -> tuple[Path | None, dict[str, np.ndarray]]:
    path = _find_first(roots, ("reconstruction*.npz", "rho_reco*.npz"))
    if path is None:
        collector.warn("No reconstruction artifact was found; reconstruction FoMs remain future work.")
        return None, {}
    try:
        with np.load(path) as archive:
            return path, {key: archive[key] for key in archive.files}
    except (OSError, ValueError) as exc:
        collector.warn(f"Could not read reconstruction artifact {path}: {exc}")
        return None, {}


def _add_reconstruction_fom(
    roots: list[Path], collector: FoMCollector, figures: Path
) -> dict[str, Any]:
    path, arrays = _load_reconstruction_npz(roots, collector)
    if path is None:
        return {}
    rho_true = arrays.get("rho_true")
    rho_reco = arrays.get("rho_reco")
    if rho_true is None or rho_reco is None or rho_true.shape != rho_reco.shape:
        collector.warn("Reconstruction artifact must contain shape-compatible rho_true and rho_reco arrays.")
        return {}
    true = np.asarray(rho_true, dtype=float)
    reco = np.asarray(rho_reco, dtype=float)
    valid = np.isfinite(true) & np.isfinite(reco)
    metrics: dict[str, Any] = {
        "RMSE": float(np.sqrt(np.mean((reco[valid] - true[valid]) ** 2))) if np.any(valid) else float("nan"),
        "MAE": float(np.mean(np.abs(reco[valid] - true[valid]))) if np.any(valid) else float("nan"),
    }
    coverage_mask = arrays.get("coverage_mask")
    metrics["coverage_fraction"] = (
        float(np.count_nonzero(coverage_mask) / coverage_mask.size)
        if coverage_mask is not None and coverage_mask.size
        else float(np.count_nonzero(valid) / valid.size)
    )
    anomaly_mask = arrays.get("anomaly_mask")
    if anomaly_mask is not None:
        anomaly = np.asarray(anomaly_mask, dtype=bool) & valid
        host_mask = arrays.get("host_mask")
        host = (np.asarray(host_mask, dtype=bool) if host_mask is not None else ~np.asarray(anomaly_mask, dtype=bool)) & valid
        if np.any(anomaly) and np.any(host):
            true_contrast = float(np.mean(true[anomaly]) - np.mean(true[host]))
            reco_contrast = float(np.mean(reco[anomaly]) - np.mean(reco[host]))
            metrics["CRF"] = reco_contrast / true_contrast if true_contrast != 0 else float("nan")
            host_sigma = float(np.std(reco[host], ddof=1)) if np.count_nonzero(host) > 1 else float("nan")
            metrics["CNR"] = abs(reco_contrast) / host_sigma if host_sigma > 0 else float("nan")
    for metric, value in metrics.items():
        unit = (
            "fraction"
            if metric in {"CRF", "coverage_fraction"}
            else "dimensionless"
            if metric == "CNR"
            else "g/cm^3"
        )
        collector.add(
            "reconstruction",
            "rho reconstruction",
            metric,
            value,
            unit,
            "ILLUSTRATIVE",
            path,
            "Prototype reconstruction metric; not a closed posterior-uncertainty result.",
        )
    fig, ax = plt.subplots(figsize=(5.8, 5.0))
    ax.scatter(true[valid], reco[valid], s=8, alpha=0.45)
    if np.any(valid):
        low = float(min(np.min(true[valid]), np.min(reco[valid])))
        high = float(max(np.max(true[valid]), np.max(reco[valid])))
        ax.plot([low, high], [low, high], linestyle="--", color="black")
    ax.set_xlabel("rho_true [g/cm^3]")
    ax.set_ylabel("rho_reco [g/cm^3]")
    ax.grid(True, alpha=0.3)
    _save_figure(fig, figures / "reconstruction_fom")
    return metrics


def _markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "_No available rows._"
    columns = list(frame.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        values: list[str] = []
        for column in columns:
            value = row[column]
            if value is None or (isinstance(value, float) and not math.isfinite(value)):
                values.append("not available")
            elif isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _write_outputs(
    outdir: Path,
    collector: FoMCollector,
    tables: dict[str, Any],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    records = [asdict(metric) for metric in collector.metrics]
    pd.DataFrame(records).to_csv(outdir / "fom_summary.csv", index=False)
    payload = {
        "metadata": {
            "purpose": "Figures of merit aggregated from existing characterization outputs",
            "posterior_uncertainty_status": "illustrative/prototype; not a closed result",
        },
        "warnings": collector.warnings,
        "metrics": records,
        "tables": {name: _records(frame) for name, frame in tables.items()},
    }
    (outdir / "fom_summary.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    simulation = pd.DataFrame(
        [record for record in records if record["status"] == "SIMULATION"]
    )
    illustrative = pd.DataFrame(
        [record for record in records if record["status"] == "ILLUSTRATIVE"]
    )
    warning_lines = [f"- {warning}" for warning in collector.warnings] or ["- None."]
    report = "\n".join(
        [
            "# Figures of Merit Report",
            "",
            "This report aggregates existing outputs. It does not alter the detector physics, "
            "reinterpret experimental data as cosmic-muon data, or claim closed posterior uncertainty.",
            "",
            "## [SIMULATION]",
            "",
            "Geometry and tolerance values derived from the current Phase A/B simulation outputs.",
            "",
            _markdown_table(simulation),
            "",
            "The `effective_acceptance_sum_omega` metric is `sum(omega)` in sr. It is not a physical "
            "area-weighted acceptance because no detector area factor is included.",
            "",
            "## [ILLUSTRATIVE]",
            "",
            "Phase C-1 sensitivity and forward-operator values are synthetic CSDA benchmarks. "
            "They are useful for comparisons, not experimental cosmic-muon measurements.",
            "",
            _markdown_table(illustrative),
            "",
            "## [FUTURE WORK]",
            "",
            "- Serialize per-pixel-pair `omega` and `sigma_omega` to report fractional uncertainty percentiles and `robustness_index`.",
            "- Serialize split tolerance contributions for `D`, `alpha_x`, `alpha_y`, `dx`, and `dy`.",
            "- Serialize `F` to report p95 nnz and voxel/track coverage fractions directly.",
            "- Serialize the Phase C distance sweep before reporting an exposure-versus-range curve.",
            "- Add consistent `rho_true`, `rho_reco`, and anomaly masks before reporting reconstruction FoMs.",
            "- Posterior uncertainty remains illustrative/prototype until shared-grid posterior covariance is implemented and validated.",
            "",
            "## Warnings",
            "",
            *warning_lines,
        ]
    )
    (outdir / "fom_report.md").write_text(report + "\n", encoding="utf-8")


def _records(value: Any) -> Any:
    if isinstance(value, pd.DataFrame):
        return [
            {key: _json_value(item) for key, item in row.items()}
            for row in value.to_dict(orient="records")
        ]
    if isinstance(value, dict):
        return {
            key: _records(item)
            for key, item in value.items()
            if key not in {"row_nnz", "column_nnz"}
        }
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    return _json_value(value)


def run(root: Path, outdir: Path | None = None) -> Path:
    """Generate FoM tables and figures from outputs found below ``root``."""

    collector = FoMCollector(root)
    roots = _search_roots(root)
    if not roots:
        collector.warn("No data/results or characterization output directory was found.")
    config = _load_yaml_config(root, collector)
    geom = _geometry_from_config(config)
    target = outdir or root / "data" / "results" / "figures_of_merit"
    figures = target / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    geometry = _add_geometric_fom(roots, geom, collector, figures)
    robustness, tolerance_contributions = _add_robustness_fom(roots, collector, figures)
    report_path, phase_c_tables = _phase_c_report(roots, collector)
    sensitivity = _add_sensitivity_fom(report_path, phase_c_tables, collector, figures)
    forward = _add_forward_fom(roots, report_path, phase_c_tables, collector, figures)
    reconstruction = _add_reconstruction_fom(roots, collector, figures)
    collector.warn(
        "No dedicated FUENTES/source-model output is currently exposed; source-model "
        "comparisons were not aggregated."
    )
    _write_outputs(
        target,
        collector,
        {
            "geometric": geometry,
            "robustness": robustness,
            "tolerance_contributions": tolerance_contributions,
            "sensitivity": sensitivity,
            "forward_operator": forward,
            "reconstruction": reconstruction,
        },
    )
    LOGGER.info("Wrote figures of merit to %s", target)
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=_repo_root())
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
    run(args.root.resolve(), args.outdir.resolve() if args.outdir else None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
