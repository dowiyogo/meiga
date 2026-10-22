from __future__ import annotations

import numpy as np

from characterization.geometry import DetectorGeometry
from characterization.solid_angle import omega_analytical
from characterization.tolerances.analytical import dOmega_dD
from characterization.tolerances.error_budget import (
    dominant_term,
    propagate_to_acceptance,
)
from characterization.tolerances.finite_diff import fd_dOmega


def test_domega_dD_agrees_with_central_finite_difference() -> None:
    geom = DetectorGeometry()
    pixel = geom.n_channels * (geom.n_channels // 2) + geom.n_channels // 2
    analytical = dOmega_dD(pixel, pixel, geom)
    fd = fd_dOmega(geom, "D", 1e-3, pixel_pairs=[(pixel, pixel)])[pixel, pixel]
    assert abs(analytical - fd) / abs(fd) < 1e-5


def test_uniform_translation_of_both_modules_is_invariant() -> None:
    geom = DetectorGeometry(n_channels=2, D=20.0)

    class ShiftBothGeometry(DetectorGeometry):
        def rectangle_for_pixel(
            self, pixel_id: int
        ) -> tuple[float, float, float, float]:
            rect = super().rectangle_for_pixel(pixel_id)
            return rect[0] + 7.0, rect[1] + 7.0, rect[2] - 3.0, rect[3] - 3.0

    shifted = ShiftBothGeometry(
        D=geom.D,
        bar_length=geom.bar_length,
        bar_width=geom.bar_width,
        n_channels=geom.n_channels,
        pixel_pitch=geom.pixel_pitch,
        sub_plane_offset=geom.sub_plane_offset,
    )
    assert abs(omega_analytical(0, 3, geom) - omega_analytical(0, 3, shifted)) < 1e-12


def test_positive_dD_decreases_center_omega() -> None:
    geom = DetectorGeometry()
    pixel = geom.n_channels * (geom.n_channels // 2) + geom.n_channels // 2
    assert dOmega_dD(pixel, pixel, geom) < 0.0


def test_error_budget_single_variable_reduction() -> None:
    geom = DetectorGeometry(n_channels=1, D=30.0)
    sigma = 0.5
    variance = propagate_to_acceptance(geom, {"D": sigma})
    expected = dOmega_dD(0, 0, geom) ** 2 * sigma * sigma
    assert abs(variance[0, 0] - expected) < 1e-18


def test_dominant_term_toy_manual(monkeypatch) -> None:
    geom = DetectorGeometry(n_channels=2, D=10.0)

    def fake_jacobian_matrix(_geom, fields, **_kwargs):
        del _geom
        return {
            "D": np.asarray(
                [
                    [4.0, 1.0, 1.0, 1.0],
                    [1.0, 4.0, 1.0, 1.0],
                    [1.0, 1.0, 4.0, 1.0],
                    [1.0, 1.0, 1.0, 4.0],
                ]
            ),
            "dx": np.asarray(
                [
                    [1.0, 5.0, 5.0, 5.0],
                    [5.0, 1.0, 5.0, 5.0],
                    [5.0, 5.0, 1.0, 5.0],
                    [5.0, 5.0, 5.0, 1.0],
                ]
            ),
        }

    monkeypatch.setattr(
        "characterization.tolerances.error_budget.jacobian_matrix", fake_jacobian_matrix
    )
    result = dominant_term(geom, {"D": 1.0, "dx": 1.0}, fields=("D", "dx"))
    assert result[0, 0] == "D"
    assert result[0, 1] == "dx"
