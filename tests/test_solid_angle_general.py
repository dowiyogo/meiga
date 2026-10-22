from __future__ import annotations

import numpy as np
import pytest

from characterization.geometry import DetectorGeometry
from characterization.solid_angle import omega_analytical
from characterization.solid_angle_general import Rectangle3D, omega_general
from characterization.tolerances.perturbation import apply_perturbation
from characterization.tolerances.tilt_validation import build_rectangles


def _central_pair(geom: DetectorGeometry) -> tuple[int, int]:
    pixel = geom.n_channels * (geom.n_channels // 2) + geom.n_channels // 2
    return pixel, pixel


def test_parallel_limit_matches_parallel_kernel() -> None:
    geom = DetectorGeometry()
    top, bot = _central_pair(geom)
    rect_top, rect_bot = build_rectangles(
        apply_perturbation(geom, {}),
        top % geom.n_channels,
        top // geom.n_channels,
        bot % geom.n_channels,
        bot // geom.n_channels,
    )
    omega_gl = omega_general(rect_top, rect_bot, method="gauss_legendre", n_quad=12)
    omega_parallel = omega_analytical(top, bot, geom)
    assert abs(omega_gl - omega_parallel) / omega_parallel < 1e-8


def test_gauss_legendre_convergence() -> None:
    geom = apply_perturbation(DetectorGeometry(), {"tilt_x": 0.012, "tilt_y": -0.007})
    top, bot = _central_pair(geom)
    rect_top, rect_bot = build_rectangles(
        geom,
        top % geom.n_channels,
        top // geom.n_channels,
        bot % geom.n_channels,
        bot // geom.n_channels,
    )
    omega_8 = omega_general(rect_top, rect_bot, method="gauss_legendre", n_quad=8)
    omega_16 = omega_general(rect_top, rect_bot, method="gauss_legendre", n_quad=16)
    assert abs(omega_16 - omega_8) < 1e-7


def test_monte_carlo_cross_check() -> None:
    geom = apply_perturbation(DetectorGeometry(), {"tilt_x": 0.01})
    top, bot = _central_pair(geom)
    rect_top, rect_bot = build_rectangles(
        geom,
        top % geom.n_channels,
        top // geom.n_channels,
        bot % geom.n_channels,
        bot // geom.n_channels,
    )
    omega_gl = omega_general(rect_top, rect_bot, method="gauss_legendre", n_quad=12)
    omega_mc = omega_general(
        rect_top, rect_bot, method="monte_carlo", n_samples=2_000_000, seed=2025
    )
    assert abs(omega_mc - omega_gl) / omega_gl < 0.003


def test_coplanar_rectangles_raise() -> None:
    rect_top = Rectangle3D([0, 0, 0], [0, 0, 1], [1, 0, 0], [0, 1, 0], 1, 1)
    rect_bot = Rectangle3D([10, 0, 0], [0, 0, -1], [1, 0, 0], [0, 1, 0], 1, 1)
    with pytest.raises(ValueError, match="Coplanar"):
        omega_general(rect_top, rect_bot)


def test_translation_invariance() -> None:
    geom = apply_perturbation(DetectorGeometry(), {"tilt_x": 0.006})
    top, bot = _central_pair(geom)
    rect_top, rect_bot = build_rectangles(
        geom,
        top % geom.n_channels,
        top // geom.n_channels,
        bot % geom.n_channels,
        bot // geom.n_channels,
    )
    shift = np.asarray([7.0, -3.0, 11.0])
    shifted_top = Rectangle3D(
        rect_top.center + shift,
        rect_top.normal,
        rect_top.u_axis,
        rect_top.v_axis,
        rect_top.half_extent_u,
        rect_top.half_extent_v,
    )
    shifted_bot = Rectangle3D(
        rect_bot.center + shift,
        rect_bot.normal,
        rect_bot.u_axis,
        rect_bot.v_axis,
        rect_bot.half_extent_u,
        rect_bot.half_extent_v,
    )
    original = omega_general(rect_top, rect_bot, method="gauss_legendre", n_quad=12)
    translated = omega_general(
        shifted_top, shifted_bot, method="gauss_legendre", n_quad=12
    )
    assert abs(original - translated) / original < 1e-14
