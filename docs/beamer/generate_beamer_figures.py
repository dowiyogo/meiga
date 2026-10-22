#!/usr/bin/env python3
"""Generate Beamer figures for the NA64-mini hodoscope results."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


FIGURE_NAMES = [
    "hodoscope_z_layout.pdf",
    "module_xy_structure.pdf",
    "pixel_matrix_16x16.pdf",
    "coincidence_counts.pdf",
    "multiplicity_by_axis.pdf",
    "bar_occupancy.pdf",
    "pixel_occupancy_lower.pdf",
    "pixel_occupancy_upper.pdf",
]


def find_repo_root() -> Path:
    starts = [Path.cwd(), Path(__file__).resolve().parent]
    for start in starts:
        for candidate in [start, *start.parents]:
            if (candidate / "diagnostics").is_dir():
                return candidate
    return Path.cwd()


ROOT = find_repo_root()
OUTDIR = ROOT / "docs" / "beamer" / "figures"
DIAG = ROOT / "diagnostics"
LOG_LINES: list[str] = []


def log(message: str) -> None:
    LOG_LINES.append(message)
    print(message)


def save_log() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "generation_log.txt").write_text("\n".join(LOG_LINES) + "\n", encoding="utf-8")


def placeholder(name: str, message: str) -> None:
    log(f"[placeholder] {name}: {message}")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.axis("off")
    ax.text(
        0.5,
        0.58,
        "missing input file" if "missing" in message.lower() else "figure unavailable",
        ha="center",
        va="center",
        fontsize=18,
        weight="bold",
    )
    ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=10, wrap=True)
    fig.tight_layout()
    fig.savefig(OUTDIR / name)
    plt.close(fig)


def read_text(path: Path) -> str | None:
    if not path.exists():
        log(f"[missing] {path}")
        return None
    return path.read_text(encoding="utf-8")


def parse_summary() -> tuple[dict[str, int], dict[str, dict[int, int]]]:
    text = read_text(DIAG / "hodoscope_output_summary.md")
    if text is None:
        return {}, {}

    patterns = {
        "total events": r"total events = (\d+)",
        "LOWER_X": r"events with hit in LOWER_X = (\d+)",
        "LOWER_Y": r"events with hit in LOWER_Y = (\d+)",
        "UPPER_X": r"events with hit in UPPER_X = (\d+)",
        "UPPER_Y": r"events with hit in UPPER_Y = (\d+)",
        "LOWER_XY": r"events with LOWER_XY coincidence = (\d+)",
        "UPPER_XY": r"events with UPPER_XY coincidence = (\d+)",
        "FULL": r"events with full LOWER_XY && UPPER_XY coincidence = (\d+)",
    }
    counts = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            counts[key] = int(match.group(1))
        else:
            log(f"[warning] Could not parse summary count for {key}")

    multiplicity: dict[str, dict[int, int]] = {}
    for axis, mult, count in re.findall(r"- (LOWER_X|LOWER_Y|UPPER_X|UPPER_Y) multiplicity (\d+): (\d+)", text):
        multiplicity.setdefault(axis, {})[int(mult)] = int(count)

    return counts, multiplicity


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        log(f"[missing] {path}")
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if rows:
        log(f"[csv] {path.name}: columns={list(rows[0].keys())}")
    return rows


def infer_column(rows: list[dict[str, str]], candidates: list[str], role: str) -> str | None:
    if not rows:
        return None
    columns = list(rows[0].keys())
    for candidate in candidates:
        if candidate in columns:
            return candidate
    lower_map = {column.lower(): column for column in columns}
    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    log(f"[warning] Could not infer {role}. Available columns: {columns}")
    return None


def figure_z_layout() -> None:
    fig, ax = plt.subplots(figsize=(8.8, 3.2))
    ax.set_xlim(-8, 118)
    ax.set_ylim(-1.5, 1.5)
    ax.axhline(0, color="0.25", linewidth=1.2)
    modules = [("LOWER", 0.0, 12.0, "#3b6ea8"), ("UPPER", 98.4, 110.4, "#c45f3c")]
    for label, z0, z1, color in modules:
        ax.add_patch(Rectangle((z0, -0.42), z1 - z0, 0.84, facecolor=color, edgecolor="black", alpha=0.85))
        ax.text((z0 + z1) / 2, 0.65, f"{label}\n{z0:.1f}-{z1:.1f} mm", ha="center", va="bottom", fontsize=10)
        ax.plot([(z0 + z1) / 2, (z0 + z1) / 2], [-0.62, 0.62], color="black", linewidth=0.9, linestyle="--")
    ax.annotate("", xy=(12.0, -0.9), xytext=(98.4, -0.9), arrowprops=dict(arrowstyle="<->", lw=1.4))
    ax.text(55.2, -1.08, "separacion libre = 86.40 mm", ha="center", va="top", fontsize=10)
    ax.annotate("", xy=(6.0, 1.15), xytext=(104.4, 1.15), arrowprops=dict(arrowstyle="<->", lw=1.4))
    ax.text(55.2, 1.26, "separacion entre centros = 98.40 mm", ha="center", va="bottom", fontsize=10)
    ax.set_xlabel("z [mm]")
    ax.set_yticks([])
    ax.set_title("Disposicion global en Z")
    ax.spines[["left", "right", "top"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "hodoscope_z_layout.pdf")
    plt.close(fig)


def figure_module_structure() -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.set_xlim(-20, 20)
    ax.set_ylim(-21, 21)
    ax.set_aspect("equal")
    base = [i * 4.0 - 15.0 for i in range(8)]
    colors = {"X_A": "#2f6f9f", "X_B": "#75aadb", "Y_A": "#b7523b", "Y_B": "#e08a6a"}
    for i, x in enumerate(base):
        ax.add_patch(Rectangle((x - 1.5, -18), 3.0, 15.0, facecolor=colors["X_A"], edgecolor="black", alpha=0.65))
        ax.add_patch(Rectangle((x + 2.0 - 1.5, -18), 3.0, 15.0, facecolor=colors["X_B"], edgecolor="black", alpha=0.45))
        ax.add_patch(Rectangle((-18, x - 1.5), 15.0, 3.0, facecolor=colors["Y_A"], edgecolor="black", alpha=0.65))
        ax.add_patch(Rectangle((-18, x + 2.0 - 1.5), 15.0, 3.0, facecolor=colors["Y_B"], edgecolor="black", alpha=0.45))
    ax.annotate("offset 2 mm", xy=(base[4], -19.5), xytext=(base[4] + 2.0, -19.5),
                arrowprops=dict(arrowstyle="<->"), ha="center", va="top", fontsize=9)
    ax.text(9, -14, "canales X\nreconstruyen X", ha="center", va="center", fontsize=10, color="#173a55")
    ax.text(-13, 11, "canales Y\nreconstruyen Y", ha="center", va="center", fontsize=10, color="#6b271d")
    ax.text(0, 19, "Modulo X-Y: 2 ejes, 2 subplanos/eje, 8 barras/subplano", ha="center", fontsize=11)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.grid(True, color="0.9")
    fig.tight_layout()
    fig.savefig(OUTDIR / "module_xy_structure.pdf")
    plt.close(fig)


def figure_pixel_matrix() -> None:
    fig, ax = plt.subplots(figsize=(6.2, 5.7))
    for x in range(17):
        ax.plot([x, x], [0, 16], color="0.75", linewidth=0.5)
    for y in range(17):
        ax.plot([0, 16], [y, y], color="0.75", linewidth=0.5)
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 16)
    ax.set_aspect("equal")
    ax.set_xticks(range(16))
    ax.set_yticks(range(16))
    ax.set_xticklabels(range(16), fontsize=7)
    ax.set_yticklabels(range(16), fontsize=7)
    ax.set_xlabel("x_index")
    ax.set_ylabel("y_index")
    ax.set_title("Matriz discreta 16 x 16 = 256 pixeles")
    ax.text(8, -2.2, r"$pixel\_id = x + 16\,y$", ha="center", va="top", fontsize=12)
    fig.tight_layout()
    fig.savefig(OUTDIR / "pixel_matrix_16x16.pdf")
    plt.close(fig)


def figure_coincidences(counts: dict[str, int]) -> None:
    labels = ["LOWER_X", "LOWER_Y", "UPPER_X", "UPPER_Y", "LOWER_XY", "UPPER_XY", "FULL"]
    if not all(label in counts for label in labels):
        placeholder("coincidence_counts.pdf", "missing coincidence counts in diagnostics/hodoscope_output_summary.md")
        return
    values = [counts[label] for label in labels]
    fig, ax = plt.subplots(figsize=(8.4, 4.2))
    bars = ax.bar(labels, values, color=["#406c9c", "#406c9c", "#bf5f45", "#bf5f45", "#5c8d5a", "#b37b4c", "#303030"])
    ax.set_ylabel("eventos")
    ax.set_title(f"Conteos de hits y coincidencias ({counts.get('total events', 0)} eventos)")
    ax.tick_params(axis="x", labelrotation=25)
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_ylim(0, max(values) * 1.18)
    ax.grid(axis="y", color="0.9")
    fig.tight_layout()
    fig.savefig(OUTDIR / "coincidence_counts.pdf")
    plt.close(fig)


def figure_multiplicity(multiplicity: dict[str, dict[int, int]]) -> None:
    axes = ["LOWER_X", "LOWER_Y", "UPPER_X", "UPPER_Y"]
    if not all(axis in multiplicity for axis in axes):
        placeholder("multiplicity_by_axis.pdf", "missing multiplicity data in diagnostics/hodoscope_output_summary.md")
        return
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    x = list(range(len(axes)))
    width = 0.24
    colors = ["#bdbdbd", "#5b8cc0", "#c65f46"]
    for idx, mult in enumerate([0, 1, 2]):
        vals = [multiplicity.get(axis, {}).get(mult, 0) for axis in axes]
        bars = ax.bar([pos + (idx - 1) * width for pos in x], vals, width=width, label=f"mult. {mult}", color=colors[idx])
        ax.bar_label(bars, padding=2, fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(axes)
    ax.set_ylabel("eventos")
    ax.set_title("Multiplicidad por eje")
    ax.legend(frameon=False, ncol=3)
    ax.grid(axis="y", color="0.9")
    fig.tight_layout()
    fig.savefig(OUTDIR / "multiplicity_by_axis.pdf")
    plt.close(fig)


def figure_bar_occupancy() -> None:
    rows = read_csv_rows(DIAG / "hodoscope_bar_occupancy.csv")
    if not rows:
        placeholder("bar_occupancy.pdf", "missing diagnostics/hodoscope_bar_occupancy.csv")
        return
    ch_col = infer_column(rows, ["global_channel", "global_channel_id", "channel_id"], "global channel")
    occ_col = infer_column(rows, ["occupancy", "counts", "count"], "occupancy")
    if ch_col is None or occ_col is None:
        placeholder("bar_occupancy.pdf", "could not infer channel/occupancy columns")
        return
    pairs = []
    for row in rows:
        try:
            pairs.append((int(float(row[ch_col])), float(row[occ_col])))
        except (TypeError, ValueError):
            continue
    pairs.sort()
    if not pairs:
        placeholder("bar_occupancy.pdf", "bar occupancy CSV has no parseable rows")
        return
    channels = [p[0] for p in pairs]
    values = [p[1] for p in pairs]
    colors = []
    for channel in channels:
        if channel < 16:
            colors.append("#406c9c")
        elif channel < 32:
            colors.append("#75aadb")
        elif channel < 48:
            colors.append("#bf5f45")
        else:
            colors.append("#e08a6a")
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    ax.bar(channels, values, color=colors, width=0.85)
    for boundary in [15.5, 31.5, 47.5]:
        ax.axvline(boundary, color="black", linewidth=0.8, linestyle="--")
    for xmid, label in [(7.5, "LOWER_X"), (23.5, "LOWER_Y"), (39.5, "UPPER_X"), (55.5, "UPPER_Y")]:
        ax.text(xmid, max(values) * 1.04, label, ha="center", va="bottom", fontsize=9)
    ax.set_xlabel("canal global")
    ax.set_ylabel("ocupancia")
    ax.set_title("Ocupancia por barra activa")
    ax.set_xlim(min(channels) - 1, max(channels) + 1)
    ax.set_ylim(0, max(values) * 1.18)
    ax.grid(axis="y", color="0.9")
    fig.tight_layout()
    fig.savefig(OUTDIR / "bar_occupancy.pdf")
    plt.close(fig)


def figure_pixel_heatmap(csv_name: str, pdf_name: str, png_name: str, title: str) -> None:
    rows = read_csv_rows(DIAG / csv_name)
    if not rows:
        placeholder(pdf_name, f"missing diagnostics/{csv_name}")
        return
    x_col = infer_column(rows, ["x_index", "x", "ix"], "x index")
    y_col = infer_column(rows, ["y_index", "y", "iy"], "y index")
    occ_col = infer_column(rows, ["occupancy", "counts", "count"], "occupancy")
    if x_col is None or y_col is None or occ_col is None:
        placeholder(pdf_name, "could not infer pixel heatmap columns")
        return
    matrix = [[0.0 for _ in range(16)] for _ in range(16)]
    for row in rows:
        try:
            x = int(float(row[x_col]))
            y = int(float(row[y_col]))
            val = float(row[occ_col])
        except (TypeError, ValueError):
            continue
        if 0 <= x < 16 and 0 <= y < 16:
            matrix[y][x] = val
    fig, ax = plt.subplots(figsize=(5.6, 4.9))
    im = ax.imshow(matrix, origin="lower", cmap="viridis", interpolation="nearest")
    ax.set_xticks(range(16))
    ax.set_yticks(range(16))
    ax.set_xticklabels(range(16), fontsize=7)
    ax.set_yticklabels(range(16), fontsize=7)
    ax.set_xlabel("x_index")
    ax.set_ylabel("y_index")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("ocupancia")
    fig.tight_layout()
    fig.savefig(OUTDIR / pdf_name)
    fig.savefig(OUTDIR / png_name, dpi=180)
    plt.close(fig)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    log(f"repo_root={ROOT}")
    log(f"diagnostics_dir={DIAG}")
    counts, multiplicity = parse_summary()

    try:
        figure_z_layout()
    except Exception as exc:  # noqa: BLE001
        placeholder("hodoscope_z_layout.pdf", str(exc))
    try:
        figure_module_structure()
    except Exception as exc:  # noqa: BLE001
        placeholder("module_xy_structure.pdf", str(exc))
    try:
        figure_pixel_matrix()
    except Exception as exc:  # noqa: BLE001
        placeholder("pixel_matrix_16x16.pdf", str(exc))
    try:
        figure_coincidences(counts)
    except Exception as exc:  # noqa: BLE001
        placeholder("coincidence_counts.pdf", str(exc))
    try:
        figure_multiplicity(multiplicity)
    except Exception as exc:  # noqa: BLE001
        placeholder("multiplicity_by_axis.pdf", str(exc))
    try:
        figure_bar_occupancy()
    except Exception as exc:  # noqa: BLE001
        placeholder("bar_occupancy.pdf", str(exc))
    try:
        figure_pixel_heatmap(
            "hodoscope_pixel_occupancy_lower.csv",
            "pixel_occupancy_lower.pdf",
            "pixel_occupancy_lower.png",
            "Ocupancia de pixeles LOWER",
        )
    except Exception as exc:  # noqa: BLE001
        placeholder("pixel_occupancy_lower.pdf", str(exc))
    try:
        figure_pixel_heatmap(
            "hodoscope_pixel_occupancy_upper.csv",
            "pixel_occupancy_upper.pdf",
            "pixel_occupancy_upper.png",
            "Ocupancia de pixeles UPPER",
        )
    except Exception as exc:  # noqa: BLE001
        placeholder("pixel_occupancy_upper.pdf", str(exc))

    for name in FIGURE_NAMES:
        if not (OUTDIR / name).exists():
            placeholder(name, "figure was not produced")

    save_log()
    return 0


if __name__ == "__main__":
    sys.exit(main())
