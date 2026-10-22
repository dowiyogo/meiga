#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path


def build_viewpoint(name: str, position: str, orientation: str, field_of_view: str) -> str:
    return (
        f'DEF {name} Viewpoint {{\n'
        f'  description "{name}"\n'
        f'  position {position}\n'
        f'  orientation {orientation}\n'
        f'  fieldOfView {field_of_view}\n'
        f'}}\n'
    )


def insert_viewpoint(text: str, viewpoint_block: str) -> str:
    header = "#VRML V2.0 utf8"
    if not text.startswith(header):
        raise ValueError("input is not a VRML 2.0 file")

    lines = text.splitlines(keepends=True)
    insert_at = 1
    while insert_at < len(lines) and lines[insert_at].lstrip().startswith("#"):
        insert_at += 1

    while insert_at < len(lines) and lines[insert_at].strip() == "":
        insert_at += 1

    return "".join(lines[:insert_at]) + "\n" + viewpoint_block + "\n" + "".join(lines[insert_at:])


def main() -> int:
    parser = argparse.ArgumentParser(description="Insert a Viewpoint into a WRL file.")
    parser.add_argument("--input", required=True, help="Input WRL file")
    parser.add_argument("--output", required=True, help="Output WRL file")
    parser.add_argument("--name", required=True, help="Viewpoint DEF name")
    parser.add_argument("--position", required=True, help="Viewpoint position, e.g. '0 0 220'")
    parser.add_argument("--orientation", required=True, help="Viewpoint orientation, e.g. '0 0 1 0'")
    parser.add_argument("--field-of-view", required=True, help="Viewpoint fieldOfView value")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    input_text = input_path.read_text(encoding="utf-8")
    viewpoint_block = build_viewpoint(args.name, args.position, args.orientation, args.field_of_view)
    output_text = insert_viewpoint(input_text, viewpoint_block)
    output_path.write_text(output_text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())