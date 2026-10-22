from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import sparse

from characterization.figures_of_merit.run_fom import (
    compute_sparsity_stats,
    robustness_index,
    run,
    summarize_fractional_uncertainty,
)


def test_missing_files_write_csv_json_and_markdown(tmp_path: Path) -> None:
    outdir = run(tmp_path)
    csv_path = outdir / "fom_summary.csv"
    json_path = outdir / "fom_summary.json"
    report_path = outdir / "fom_report.md"

    assert csv_path.exists()
    assert json_path.exists()
    assert report_path.exists()
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["warnings"]
    report = report_path.read_text(encoding="utf-8")
    assert "## [SIMULATION]" in report
    assert "## [ILLUSTRATIVE]" in report
    assert "## [FUTURE WORK]" in report
    assert "shared-grid posterior covariance" in report


def test_compute_sparsity_stats() -> None:
    matrix = sparse.csr_matrix(
        np.asarray(
            [
                [1.0, 0.0, 2.0, 0.0],
                [0.0, 0.0, 3.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ]
        )
    )
    stats = compute_sparsity_stats(matrix)

    assert stats["shape"] == [3, 4]
    assert stats["nnz"] == 3
    assert stats["sparsity"] == 0.75
    assert stats["voxel_coverage_fraction"] == 0.5
    assert stats["track_coverage_fraction"] == 2.0 / 3.0


def test_fractional_uncertainty_percentiles() -> None:
    stats = summarize_fractional_uncertainty(
        np.asarray([10.0, 10.0, 10.0, 10.0]),
        np.asarray([1.0, 2.0, 3.0, 4.0]),
    )

    assert stats["median"] == 0.25
    assert np.isclose(stats["p90"], 0.37)
    assert np.isclose(stats["p95"], 0.385)
    assert stats["max"] == 0.4


def test_robustness_index_ignores_zeros_and_nan() -> None:
    value = robustness_index(
        np.asarray([10.0, 0.0, np.nan, 4.0, 3.0]),
        np.asarray([1.0, 1.0, 1.0, 2.0, 0.0]),
    )

    assert value == 6.0
