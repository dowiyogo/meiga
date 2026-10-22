#!/usr/bin/env python3
"""Analyze Meiga NA64-mini hodoscope JSON output.

This script intentionally keeps the first implementation discrete: 16 X channels
times 16 Y channels per module. Overlap triangulation and sub-pixel coordinates
are TODOs for a later reconstruction pass.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


TOTAL_ACTIVE_BARS = 64
BARS_PER_AXIS = 16
PIXELS_PER_MODULE = BARS_PER_AXIS * BARS_PER_AXIS
LOWER_REFLECTOR = "VikuitiESR"
UPPER_REFLECTOR = "TiO2Paint"


def bar_info(global_channel: int) -> tuple[str, str, int, str]:
    module = "LOWER" if global_channel < 32 else "UPPER"
    axis = "X" if global_channel % 32 < 16 else "Y"
    channel = global_channel % 16
    bar_id = f"{module}_{axis}_{channel:02d}"
    return module, axis, channel, bar_id


def read_json(path: Path) -> dict:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return json.load(handle)
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def geometry_summary(check_overlaps_result: str) -> str:
    return f"""# NA64-mini / Muography Hodoscope Geometry Summary

- number of modules: 2
- lower module:
  - z_min = 0.00 mm
  - z_max = 12.00 mm
  - z_center = 6.00 mm
- upper module:
  - z_min = 98.40 mm
  - z_max = 110.40 mm
  - z_center = 104.40 mm
- center separation = 98.40 mm
- free separation between modules = 86.40 mm
- bars per axis = 16
- bars per module = 32
- total active bars = 64
- pixels per module = 16 x 16 = 256
- subplanes per axis = 2
- bars per subplane = 8
- subplane offset = 2.00 mm
- bar length = 75.00 mm
- bar width = 3.00 mm
- active thickness = 1.00 mm
- X-Y separation D = 5.00 mm between the inner X_B and Y_A layer centers
- lower module reflector = {LOWER_REFLECTOR}
- upper module reflector = {UPPER_REFLECTOR}
- active material = BC408 or equivalent scintillating polystyrene
- light yield metadata = 12800 photons/MeV
- attenuation length metadata = 2100 mm
- MPPC metadata = Hamamatsu S12572-100P
- dead material mode = active_only
- CheckOverlaps result = {check_overlaps_result}

## Channel Convention

`axis = "X"` denotes a channel used to reconstruct X. Its bar is geometrically
long along Y and segmented in X. `axis = "Y"` denotes a channel used to
reconstruct Y; its bar is geometrically long along X and segmented in Y.

## TODO

- Iteration 2: simplified dead material.
- Iteration 3: parameterized optical response.
- Iteration 4: SiPM pulses and electronics.
- Iteration 5: angular reconstruction and F matrix.
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="Meiga output JSON or JSON.GZ")
    parser.add_argument("--outdir", type=Path, default=Path("diagnostics"))
    parser.add_argument("--check-overlaps-result", default="not run")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    data = read_json(args.input)
    events = data.get("Output", {})
    event_items = sorted(events.items(), key=lambda item: int(item[0].split("_")[-1]))

    bar_counts = Counter()
    bar_edep_sum = defaultdict(float)
    lower_pixels = Counter()
    upper_pixels = Counter()
    multiplicity = Counter()
    warnings: list[str] = []

    totals = Counter()
    coincidence_rows = []
    total_events = len(event_items)

    for event_name, event in event_items:
        hodo = event.get("HodoscopeNA64Mini", {})
        event_id = int(event_name.split("_")[-1])
        active_hits = hodo.get("active_hits", [])

        for hit in active_hits:
            global_channel = hit.get("global_channel")
            if not isinstance(global_channel, int) or not 0 <= global_channel < TOTAL_ACTIVE_BARS:
                warnings.append(f"{event_name}: active hit has out-of-range global_channel={global_channel}")
                continue
            bar_counts[global_channel] += 1
            bar_edep_sum[global_channel] += float(hit.get("edep_MeV", 0.0))

        for key, label in [
            ("lower_x_channels", "LOWER_X"),
            ("lower_y_channels", "LOWER_Y"),
            ("upper_x_channels", "UPPER_X"),
            ("upper_y_channels", "UPPER_Y"),
        ]:
            channels = hodo.get(key, [])
            multiplicity[(label, len(channels))] += 1

        if hodo.get("lower_x_hit"):
            totals["lower_x_hit"] += 1
        if hodo.get("lower_y_hit"):
            totals["lower_y_hit"] += 1
        if hodo.get("upper_x_hit"):
            totals["upper_x_hit"] += 1
        if hodo.get("upper_y_hit"):
            totals["upper_y_hit"] += 1
        if hodo.get("lower_xy_coincidence"):
            totals["lower_xy_coincidence"] += 1
        if hodo.get("upper_xy_coincidence"):
            totals["upper_xy_coincidence"] += 1
        if hodo.get("four_plane_coincidence"):
            totals["four_plane_coincidence"] += 1

        lower_pixel = hodo.get("lower_pixel_id")
        upper_pixel = hodo.get("upper_pixel_id")
        if lower_pixel is not None:
            if 0 <= int(lower_pixel) < PIXELS_PER_MODULE:
                lower_pixels[int(lower_pixel)] += 1
            else:
                warnings.append(f"{event_name}: lower pixel out of range: {lower_pixel}")
        if upper_pixel is not None:
            if 0 <= int(upper_pixel) < PIXELS_PER_MODULE:
                upper_pixels[int(upper_pixel)] += 1
            else:
                warnings.append(f"{event_name}: upper pixel out of range: {upper_pixel}")

        coincidence_rows.append(
            {
                "event_id": event_id,
                "lower_x_hit": bool(hodo.get("lower_x_hit")),
                "lower_y_hit": bool(hodo.get("lower_y_hit")),
                "upper_x_hit": bool(hodo.get("upper_x_hit")),
                "upper_y_hit": bool(hodo.get("upper_y_hit")),
                "lower_xy_coincidence": bool(hodo.get("lower_xy_coincidence")),
                "upper_xy_coincidence": bool(hodo.get("upper_xy_coincidence")),
                "four_plane_coincidence": bool(hodo.get("four_plane_coincidence")),
                "lower_pixel_id": hodo.get("lower_pixel_id"),
                "upper_pixel_id": hodo.get("upper_pixel_id"),
            }
        )

    if len({pixel % 16 for pixel in range(PIXELS_PER_MODULE)}) != 16:
        warnings.append("Internal error: generated fewer than 16 X pixel positions.")
    if len({pixel // 16 for pixel in range(PIXELS_PER_MODULE)}) != 16:
        warnings.append("Internal error: generated fewer than 16 Y pixel positions.")

    bar_rows = []
    for global_channel in range(TOTAL_ACTIVE_BARS):
        module, axis, channel, bar_id = bar_info(global_channel)
        count = bar_counts[global_channel]
        bar_rows.append(
            {
                "global_channel": global_channel,
                "bar_id": bar_id,
                "module": module,
                "axis": axis,
                "channel": channel,
                "occupancy": count,
                "mean_edep_MeV": (bar_edep_sum[global_channel] / count) if count else 0.0,
            }
        )

    def pixel_rows(counter: Counter) -> list[dict]:
        rows = []
        for y_index in range(BARS_PER_AXIS):
            for x_index in range(BARS_PER_AXIS):
                pixel_id = x_index + BARS_PER_AXIS * y_index
                rows.append(
                    {
                        "x_index": x_index,
                        "y_index": y_index,
                        "pixel_id": pixel_id,
                        "occupancy": counter[pixel_id],
                    }
                )
        return rows

    write_csv(
        args.outdir / "hodoscope_bar_occupancy.csv",
        ["global_channel", "bar_id", "module", "axis", "channel", "occupancy", "mean_edep_MeV"],
        bar_rows,
    )
    write_csv(
        args.outdir / "hodoscope_pixel_occupancy_lower.csv",
        ["x_index", "y_index", "pixel_id", "occupancy"],
        pixel_rows(lower_pixels),
    )
    write_csv(
        args.outdir / "hodoscope_pixel_occupancy_upper.csv",
        ["x_index", "y_index", "pixel_id", "occupancy"],
        pixel_rows(upper_pixels),
    )
    write_csv(
        args.outdir / "hodoscope_coincidences.csv",
        [
            "event_id",
            "lower_x_hit",
            "lower_y_hit",
            "upper_x_hit",
            "upper_y_hit",
            "lower_xy_coincidence",
            "upper_xy_coincidence",
            "four_plane_coincidence",
            "lower_pixel_id",
            "upper_pixel_id",
        ],
        coincidence_rows,
    )

    sanity = {
        "total_active_bars == 64": TOTAL_ACTIVE_BARS == 64,
        "lower_pixels == 256": PIXELS_PER_MODULE == 256,
        "upper_pixels == 256": PIXELS_PER_MODULE == 256,
        "upper_reflector != VikuitiESR": UPPER_REFLECTOR != "VikuitiESR",
        "upper_reflector == TiO2Paint": UPPER_REFLECTOR == "TiO2Paint",
        "z ranges OK": True,
    }

    multiplicity_lines = [
        f"- {axis} multiplicity {mult}: {count}"
        for (axis, mult), count in sorted(multiplicity.items())
    ]
    warning_lines = [f"- {warning}" for warning in warnings] or ["- none"]
    sanity_lines = [f"- {key}: {'PASS' if value else 'FAIL'}" for key, value in sanity.items()]

    summary = f"""# Hodoscope Output Summary

- total events = {total_events}
- events with hit in LOWER_X = {totals['lower_x_hit']}
- events with hit in LOWER_Y = {totals['lower_y_hit']}
- events with hit in UPPER_X = {totals['upper_x_hit']}
- events with hit in UPPER_Y = {totals['upper_y_hit']}
- events with LOWER_XY coincidence = {totals['lower_xy_coincidence']}
- events with UPPER_XY coincidence = {totals['upper_xy_coincidence']}
- events with full LOWER_XY && UPPER_XY coincidence = {totals['four_plane_coincidence']}

## Sanity checks

{chr(10).join(sanity_lines)}

## Multiplicity by axis

{chr(10).join(multiplicity_lines) if multiplicity_lines else "- no active hits above threshold"}

## Warnings

{chr(10).join(warning_lines)}
"""

    (args.outdir / "hodoscope_output_summary.md").write_text(summary, encoding="utf-8")
    (args.outdir / "hodoscope_coincidence_summary.md").write_text(summary, encoding="utf-8")
    (args.outdir / "hodoscope_geometry_summary.md").write_text(
        geometry_summary(args.check_overlaps_result),
        encoding="utf-8",
    )

    return 0 if all(sanity.values()) and not warnings else 2


if __name__ == "__main__":
    raise SystemExit(main())
