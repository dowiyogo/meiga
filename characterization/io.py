"""Input/output utilities for Meiga hodoscope characterization.

Minimal reproducible example:

```python
from characterization.io import load_events, load_geometry

geom = load_geometry("src/Applications/G4HodoscopeSimulator/G4HodoscopeSimulator.json")
events = load_events("src/Applications/G4HodoscopeSimulator/G4HodoscopeSimulator.json")
```

Passing a simulator config resolves ``Output.OutputFile`` relative to the
config file location. Passing an output JSON reads that file directly.
"""

from __future__ import annotations

import gzip
import json
import logging
import math
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from characterization.geometry import DetectorGeometry

LOGGER = logging.getLogger(__name__)


EVENT_COLUMNS = [
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


def _read_json(path: Path) -> dict[str, Any]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return json.load(handle)
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _unit_scale(unit: str | None) -> float:
    unit = (unit or "mm").strip().lower()
    scales = {"mm": 1.0, "millimeter": 1.0, "cm": 10.0, "m": 1000.0, "um": 1e-3}
    return scales.get(unit, 1.0)


def _float_child(node: ET.Element, name: str, default: float) -> float:
    child = node.find(name)
    if child is None or child.text is None:
        return default
    return float(child.text) * _unit_scale(child.attrib.get("unit"))


def _geometry_from_xml(path: Path, d_override: float | None = None) -> DetectorGeometry:
    root = ET.parse(path).getroot()
    detector = root.find(".//detector[@type='eHodoscope']")
    if detector is None:
        LOGGER.warning("No eHodoscope detector found in %s; using defaults", path)
        return DetectorGeometry(D=d_override or DetectorGeometry.D)

    n_channels = int(_float_child(detector, "numberOfBars", 16.0))
    bar_length = _float_child(detector, "barLength", 75.0)
    bar_width = _float_child(detector, "barWidth", 3.0)
    bar_thickness = _float_child(detector, "barThickness", 1.0)
    pixel_pitch = 2.0
    sub_plane_offset = 2.0
    module_size = (bar_length, bar_length)
    return DetectorGeometry(
        D=d_override or 98.4,
        bar_length=bar_length,
        bar_width=bar_width,
        bar_thickness=bar_thickness,
        n_channels=n_channels,
        pixel_pitch=pixel_pitch,
        sub_plane_offset=sub_plane_offset,
        module_size=module_size,
    )


def resolve_output_path(path: str | Path) -> Path:
    """Resolve a simulator output path from a config or output file.

    Args:
        path: Path to a Meiga config JSON, output JSON, or ``.json.gz`` file.

    Returns:
        Absolute path to the event output JSON. Units are not involved.

    Edge cases:
        If the file has no ``Output.OutputFile`` entry, the input path is
        treated as the output file itself.
    """

    candidate = Path(path).expanduser().resolve()
    if candidate.suffix not in {".json", ".gz"}:
        return candidate
    data = _read_json(candidate)
    output = data.get("Output", {})
    if isinstance(output, dict) and "OutputFile" in output:
        output_path = Path(str(output["OutputFile"]))
        if not output_path.is_absolute():
            output_path = candidate.parent / output_path
        return output_path.resolve()
    return candidate


def _resolve_detector_list(path: Path, data: dict[str, Any]) -> Path | None:
    detector_list = data.get("DetectorList")
    if isinstance(detector_list, str):
        xml_path = Path(detector_list)
        if not xml_path.is_absolute():
            xml_path = path.parent / xml_path
        if xml_path.exists():
            return xml_path.resolve()
    sibling = path.parent / "DetectorList.xml"
    if sibling.exists():
        return sibling.resolve()
    source_default = (
        Path.cwd()
        / "src"
        / "Applications"
        / "G4HodoscopeSimulator"
        / "DetectorList.xml"
    )
    if source_default.exists():
        return source_default.resolve()
    return None


def load_geometry(path: str | Path) -> DetectorGeometry:
    """Load detector geometry with all lengths in mm.

    Args:
        path: Simulator config JSON, output JSON/JSON.GZ, or DetectorList XML.

    Returns:
        ``DetectorGeometry`` with ``D``, bar dimensions, channel count, pitch,
        sub-plane offset, and module size in mm.

    Edge cases:
        Current Meiga event JSON has no geometry header. In that case this
        loader resolves a neighboring or configured ``DetectorList.xml`` and
        uses the NA64-mini nominal module separation ``D = 98.4 mm`` unless a
        synthetic ``Geometry``/``DetectorGeometry`` header provides ``D``.
    """

    source = Path(path).expanduser().resolve()
    if source.suffix.lower() == ".xml":
        return _geometry_from_xml(source)

    data = _read_json(source)
    header = data.get("Geometry") or data.get("DetectorGeometry") or {}
    d_override = header.get("D") if isinstance(header, dict) else None
    if d_override is not None:
        d_override = float(d_override)

    if isinstance(header, dict) and {"bar_length", "bar_width"}.issubset(header):
        return DetectorGeometry(
            D=float(header.get("D", 98.4)),
            bar_length=float(header.get("bar_length", 75.0)),
            bar_width=float(header.get("bar_width", 3.0)),
            bar_thickness=float(header.get("bar_thickness", 1.0)),
            n_channels=int(header.get("n_channels", 16)),
            pixel_pitch=float(header.get("pixel_pitch", 2.0)),
            sub_plane_offset=float(header.get("sub_plane_offset", 2.0)),
            module_size=tuple(header.get("module_size", (75.0, 75.0))),
        )

    xml_path = _resolve_detector_list(source, data)
    if xml_path is not None:
        return _geometry_from_xml(xml_path, d_override=d_override)

    LOGGER.warning("Could not resolve geometry metadata from %s; using defaults", path)
    return DetectorGeometry(D=d_override or 98.4)


def _event_id_from_name(name: str) -> int:
    try:
        return int(name.rsplit("_", 1)[-1])
    except ValueError:
        return -1


def _module_energy_time(hodo: dict[str, Any], module: str) -> tuple[float, float]:
    hits = hodo.get("active_hits", [])
    energies: list[float] = []
    times: list[float] = []
    for hit in hits:
        if hit.get("module") != module:
            continue
        energy = hit.get("edep_MeV")
        if energy is not None:
            energies.append(float(energy))
        time = hit.get("first_time_ns")
        if time is not None:
            times.append(float(time))
    return float(np.sum(energies)), float(np.min(times)) if times else math.nan


def _truth_angles(event: dict[str, Any]) -> tuple[float, float] | None:
    flux = event.get("InputFlux", {})
    momentum = flux.get("Momentum")
    if not isinstance(momentum, list) or len(momentum) != 3:
        return None
    px, py, pz = (float(momentum[0]), float(momentum[1]), float(momentum[2]))
    norm = math.sqrt(px * px + py * py + pz * pz)
    if norm == 0:
        return None
    theta = math.acos(max(-1.0, min(1.0, -pz / norm)))
    phi = math.atan2(py, px)
    return theta, phi


def load_events(path: str | Path) -> pd.DataFrame:
    """Load one row per valid four-plane hodoscope coincidence.

    Args:
        path: Simulator config JSON or event output JSON/JSON.GZ. A config path
            resolves ``Output.OutputFile`` relative to that config directory.

    Returns:
        ``pandas.DataFrame`` with columns ``event_id``, top/bottom pixel X-Y
        indices, reconstructed ``theta_rec`` and ``phi_rec`` in rad,
        ``energy_top``/``energy_bot`` in MeV, and ``t_top``/``t_bot`` in ns.

    Edge cases:
        Events missing ``HodoscopeNA64Mini`` or without a four-plane match are
        skipped. The returned frame carries ``attrs`` with ``n_total_events``,
        ``n_4plane``, ``n_incident_top`` and truth arrays when available.
    """

    input_path = Path(path).expanduser().resolve()
    output_path = resolve_output_path(input_path)
    data = _read_json(output_path)
    geom = load_geometry(input_path if input_path != output_path else output_path)
    centers = geom.pixel_centers()
    events = data.get("Output", {})
    rows: list[dict[str, float | int]] = []
    truth_theta: list[float] = []
    truth_phi: list[float] = []

    for event_name, event in sorted(
        events.items(), key=lambda item: _event_id_from_name(item[0])
    ):
        hodo = event.get("HodoscopeNA64Mini", {})
        if not hodo or not hodo.get("four_plane_coincidence"):
            continue
        top_pixel = hodo.get("upper_pixel_id")
        bot_pixel = hodo.get("lower_pixel_id")
        if top_pixel is None or bot_pixel is None:
            continue
        top_pixel = int(top_pixel)
        bot_pixel = int(bot_pixel)
        if not (0 <= top_pixel < geom.n_pixels and 0 <= bot_pixel < geom.n_pixels):
            continue

        top_x, top_y = centers[top_pixel]
        bot_x, bot_y = centers[bot_pixel]
        dx = bot_x - top_x
        dy = bot_y - top_y
        theta_rec = math.atan2(math.hypot(dx, dy), geom.D)
        phi_rec = math.atan2(dy, dx)
        energy_bot, t_bot = _module_energy_time(hodo, "LOWER")
        energy_top, t_top = _module_energy_time(hodo, "UPPER")
        rows.append(
            {
                "event_id": _event_id_from_name(event_name),
                "pixel_top_x": top_pixel % geom.n_channels,
                "pixel_top_y": top_pixel // geom.n_channels,
                "pixel_bot_x": bot_pixel % geom.n_channels,
                "pixel_bot_y": bot_pixel // geom.n_channels,
                "theta_rec": theta_rec,
                "phi_rec": phi_rec,
                "energy_top": energy_top,
                "energy_bot": energy_bot,
                "t_top": t_top,
                "t_bot": t_bot,
            }
        )
        truth = _truth_angles(event)
        if truth is not None:
            truth_theta.append(truth[0])
            truth_phi.append(truth[1])

    frame = pd.DataFrame(rows, columns=EVENT_COLUMNS)
    int_cols = ["event_id", "pixel_top_x", "pixel_top_y", "pixel_bot_x", "pixel_bot_y"]
    for col in int_cols:
        frame[col] = frame[col].astype("int64")
    for col in [c for c in EVENT_COLUMNS if c not in int_cols]:
        frame[col] = frame[col].astype("float64")

    frame.attrs["geometry"] = geom
    frame.attrs["n_total_events"] = len(events)
    frame.attrs["n_4plane"] = len(frame)
    frame.attrs["n_incident_top"] = len(events)
    if len(truth_theta) == len(frame):
        frame.attrs["theta_true"] = np.asarray(truth_theta, dtype=float)
        frame.attrs["phi_true"] = np.asarray(truth_phi, dtype=float)
    return frame
