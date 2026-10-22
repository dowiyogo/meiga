from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from characterization.acceptance import (
    cos2_bin_prob,
    geometric_efficiency,
    pixel_acceptance,
)
from characterization.geometry import DetectorGeometry
from characterization.io import load_events
from characterization.solid_angle import (
    omega_analytical,
    omega_mc,
    synthetic_pair_events,
)


def test_solid_angle_coaxial_unit_squares() -> None:
    geom = DetectorGeometry(
        D=1.0,
        bar_length=1.0,
        bar_width=1.0,
        bar_thickness=1.0,
        n_channels=1,
        pixel_pitch=1.0,
        sub_plane_offset=0.0,
        module_size=(1.0, 1.0),
    )
    assert abs(omega_analytical(0, 0, geom) - 0.6277684243304004) < 1e-4


def test_omega_mc_converges_to_analytical() -> None:
    geom = DetectorGeometry(
        D=1.0,
        bar_length=1.0,
        bar_width=1.0,
        bar_thickness=1.0,
        n_channels=1,
        pixel_pitch=1.0,
        sub_plane_offset=0.0,
        module_size=(1.0, 1.0),
    )
    events = synthetic_pair_events(0, 0, geom, n_rays=250_000, seed=7)
    estimated = omega_mc(events, geom)[0, 0]
    assert (
        abs(estimated - omega_analytical(0, 0, geom)) / omega_analytical(0, 0, geom)
        < 0.03
    )


def test_load_synthetic_minimal_json(tmp_path: Path) -> None:
    payload = {
        "Geometry": {
            "D": 98.4,
            "bar_length": 75.0,
            "bar_width": 3.0,
            "bar_thickness": 1.0,
            "n_channels": 16,
            "pixel_pitch": 2.0,
            "sub_plane_offset": 2.0,
            "module_size": [75.0, 75.0],
        },
        "Output": {
            "Event_0": {
                "InputFlux": {"Momentum": [0.0, 0.0, -1.0]},
                "HodoscopeNA64Mini": {
                    "four_plane_coincidence": True,
                    "upper_pixel_id": 108,
                    "lower_pixel_id": 101,
                    "active_hits": [
                        {"module": "UPPER", "edep_MeV": 0.2, "first_time_ns": 1.0},
                        {"module": "LOWER", "edep_MeV": 0.3, "first_time_ns": 2.0},
                    ],
                },
            },
            "Event_1": {"HodoscopeNA64Mini": {"four_plane_coincidence": False}},
        },
    }
    path = tmp_path / "output.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    events = load_events(path)
    assert len(events) == 1
    assert list(events.columns) == [
        "event_id",
        "pixel_top_x",
        "pixel_top_y",
        "pixel_bot_x",
        "pixel_bot_y",
        "theta_rec",
        "phi_rec",
        "energy_top",
        "energy_bot",
        "t_top",
        "t_bot",
    ]
    assert pd.api.types.is_integer_dtype(events["event_id"])
    assert pd.api.types.is_float_dtype(events["theta_rec"])
    assert events.loc[0, "energy_top"] == 0.2
    assert events.loc[0, "energy_bot"] == 0.3


def test_pixel_acceptance_integrates_to_geometric_efficiency() -> None:
    geom = DetectorGeometry(n_channels=2, D=10.0, bar_width=1.0, pixel_pitch=1.0)
    events = pd.DataFrame(
        {
            "event_id": [0, 1, 2],
            "pixel_top_x": [0, 0, 1],
            "pixel_top_y": [0, 0, 1],
            "pixel_bot_x": [0, 1, 1],
            "pixel_bot_y": [0, 0, 1],
            "theta_rec": [0.1, 0.2, 0.6],
            "phi_rec": [0.0, 1.0, -1.0],
            "energy_top": [1.0, 1.0, 1.0],
            "energy_bot": [1.0, 1.0, 1.0],
            "t_top": [0.0, 0.0, 0.0],
            "t_bot": [1.0, 1.0, 1.0],
        }
    )
    events.attrs["n_incident_top"] = 1000
    theta_bins = np.array([0.0, 0.3, math.pi / 2.0])
    phi_bins = np.array([-math.pi, 0.0, math.pi])
    acceptance = pixel_acceptance(events, geom, theta_bins, phi_bins)
    bin_prob = cos2_bin_prob(theta_bins, phi_bins)
    integrated_acceptance = float(np.sum(acceptance * bin_prob[None, None, :, :]))
    epsilon = geometric_efficiency(events, geom, theta_bins)
    theta_prob = epsilon["n_incident_top"].to_numpy() / 1000.0
    integrated_epsilon = float(np.sum(epsilon["epsilon"].to_numpy() * theta_prob))
    assert abs(integrated_acceptance - integrated_epsilon) < 1e-12
