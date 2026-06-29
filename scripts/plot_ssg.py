#!/usr/bin/env python3
# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Compute line-wise distance gain and net gain factors from the VTK field ``gain``.

Example usages:
    python3 scripts/plot_ssg.py \
        --input-dir example/python_example \
        --x 0 --y 0

    python3 scripts/plot_ssg.py \
        --input-dir example/python_example \
        --direction z --x 0 --y 0 \
        --back-reflection --reflectivity 1.0 \
        --output-prefix scripts/gain_z_origin
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


AXES = ("x", "y", "z")
AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
TIMESTEP_RE = re.compile(r"laserPumpCladding_(\d+)\.vtk$")


@dataclass(frozen=True)
class LineSelection:
    varying_axis: str
    fixed_axes: tuple[str, str]
    requested_coords: tuple[float, float]
    actual_coords: tuple[float, float]
    coordinate_groups: tuple[tuple[float, tuple[int, ...]], ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read laserPumpCladding_*.vtk files, integrate the local gain field "
            "'gain' along one coordinate axis, and plot net gain versus timestep."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing laserPumpCladding_*.vtk files.",
    )
    parser.add_argument(
        "--direction",
        choices=AXES,
        help=(
            "Integration axis. If omitted, the script infers it from the "
            "missing coordinate when exactly two of --x/--y/--z are provided."
        ),
    )
    parser.add_argument("--x", type=float, help="Requested x coordinate.")
    parser.add_argument("--y", type=float, help="Requested y coordinate.")
    parser.add_argument("--z", type=float, help="Requested z coordinate.")
    parser.add_argument(
        "--field",
        default="gain",
        help="POINT_DATA scalar field to integrate (default: gain).",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        help="Output prefix for the generated CSV and PNG files.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively after saving it.",
    )
    parser.add_argument(
        "--back-reflection",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Apply the legacy round-trip rule used by gain.py: "
            "net_gain = reflectivity * single_pass_gain^2. Enabled by default."
        ),
    )
    parser.add_argument(
        "--reflectivity",
        type=float,
        default=1.0,
        help="Surface reflectivity R used when --back-reflection is enabled.",
    )
    return parser.parse_args()


def resolve_direction_and_coords(args: argparse.Namespace) -> tuple[str, tuple[str, str], tuple[float, float]]:
    supplied = {axis for axis in AXES if getattr(args, axis) is not None}

    if args.direction is None:
        if len(supplied) != 2:
            raise SystemExit(
                "Provide exactly two coordinates among --x/--y/--z when --direction is omitted."
            )
        varying_axis = next(axis for axis in AXES if axis not in supplied)
    else:
        varying_axis = args.direction

    fixed_axes = tuple(axis for axis in AXES if axis != varying_axis)
    missing_fixed_axes = [axis for axis in fixed_axes if getattr(args, axis) is None]
    if missing_fixed_axes:
        missing = ", ".join(f"--{axis}" for axis in missing_fixed_axes)
        raise SystemExit(
            f"Integration along {varying_axis} requires fixed coordinates {missing}."
        )

    requested_coords = tuple(float(getattr(args, axis)) for axis in fixed_axes)
    return varying_axis, fixed_axes, requested_coords


def find_vtk_files(input_dir: Path) -> list[tuple[int, Path]]:
    matches: list[tuple[int, Path]] = []
    for path in sorted(input_dir.glob("laserPumpCladding_*.vtk")):
        match = TIMESTEP_RE.fullmatch(path.name)
        if match:
            matches.append((int(match.group(1)), path))

    if not matches:
        raise SystemExit(f"No laserPumpCladding_*.vtk files found in {input_dir}.")

    return sorted(matches, key=lambda item: item[0])


def parse_ascii_vtk_scalar(
    path: Path,
    scalar_name: str,
    *,
    include_points: bool,
) -> tuple[list[tuple[float, float, float]] | None, list[float], int]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 4 or lines[2].strip().upper() != "ASCII":
        raise ValueError(f"{path} is not an ASCII VTK file.")

    points: list[tuple[float, float, float]] | None = None
    point_count: int | None = None

    if include_points:
        for index, line in enumerate(lines):
            if line.startswith("POINTS "):
                parts = line.split()
                point_count = int(parts[1])
                flat_points = collect_float_values(lines, index + 1, point_count * 3)
                points = [
                    (
                        flat_points[offset],
                        flat_points[offset + 1],
                        flat_points[offset + 2],
                    )
                    for offset in range(0, len(flat_points), 3)
                ]
                break
        if points is None or point_count is None:
            raise ValueError(f"{path} does not contain a POINTS section.")

    point_data_count: int | None = None
    scalar_values: list[float] | None = None
    scalar_header = f"SCALARS {scalar_name} "

    for index, line in enumerate(lines):
        if line.startswith("POINT_DATA "):
            point_data_count = int(line.split()[1])
        if line.startswith(scalar_header):
            if point_data_count is None:
                raise ValueError(f"{path} declares scalar '{scalar_name}' before POINT_DATA.")
            scalar_values = collect_float_values(lines, index + 2, point_data_count)
            break

    if point_data_count is None:
        raise ValueError(f"{path} does not contain POINT_DATA.")
    if scalar_values is None:
        raise ValueError(f"{path} does not contain scalar '{scalar_name}'.")
    if include_points and point_count != point_data_count:
        raise ValueError(
            f"{path} POINTS count ({point_count}) does not match POINT_DATA count ({point_data_count})."
        )

    return points, scalar_values, point_data_count


def collect_float_values(lines: list[str], start_index: int, expected_count: int) -> list[float]:
    values: list[float] = []
    line_index = start_index
    while len(values) < expected_count and line_index < len(lines):
        stripped = lines[line_index].strip()
        if stripped:
            values.extend(float(token) for token in stripped.split())
        line_index += 1

    if len(values) < expected_count:
        raise ValueError(
            f"Expected {expected_count} float values but found only {len(values)}."
        )

    return values[:expected_count]


def select_line(
    points: Iterable[tuple[float, float, float]],
    varying_axis: str,
    fixed_axes: tuple[str, str],
    requested_coords: tuple[float, float],
) -> LineSelection:
    grouped_indices: dict[tuple[float, float], list[tuple[float, int]]] = defaultdict(list)
    varying_index = AXIS_INDEX[varying_axis]
    fixed_indices = tuple(AXIS_INDEX[axis] for axis in fixed_axes)

    for point_index, point in enumerate(points):
        fixed_key = (point[fixed_indices[0]], point[fixed_indices[1]])
        grouped_indices[fixed_key].append((point[varying_index], point_index))

    if not grouped_indices:
        raise ValueError("No candidate integration lines were found in the point set.")

    actual_coords = min(
        grouped_indices,
        key=lambda key: (key[0] - requested_coords[0]) ** 2 + (key[1] - requested_coords[1]) ** 2,
    )

    coordinate_buckets: dict[float, list[int]] = defaultdict(list)
    for coordinate, point_index in grouped_indices[actual_coords]:
        coordinate_buckets[coordinate].append(point_index)

    sorted_groups = tuple(
        (coordinate, tuple(sorted(indices)))
        for coordinate, indices in sorted(coordinate_buckets.items(), key=lambda item: item[0])
    )
    if len(sorted_groups) < 2:
        raise ValueError(
            f"The selected {varying_axis}-line at {fixed_axes[0]}={actual_coords[0]}, "
            f"{fixed_axes[1]}={actual_coords[1]} has fewer than two sample points."
        )

    return LineSelection(
        varying_axis=varying_axis,
        fixed_axes=fixed_axes,
        requested_coords=requested_coords,
        actual_coords=actual_coords,
        coordinate_groups=sorted_groups,
    )


def line_point_gains(selection: LineSelection, gains: list[float], np_module):
    axis_coordinates = np_module.asarray([group[0] for group in selection.coordinate_groups], dtype=float)
    point_gains = np_module.asarray(
        [
            sum(gains[index] for index in group[1]) / len(group[1])
            for group in selection.coordinate_groups
        ],
        dtype=float,
    )
    return axis_coordinates, point_gains


def compute_distance_gain(selection: LineSelection, gains: list[float], np_module) -> float:
    axis_coordinates, point_gains = line_point_gains(selection, gains, np_module)
    line_span = float(axis_coordinates[-1] - axis_coordinates[0])
    if line_span <= 0.0:
        raise ValueError(
            f"The selected {selection.varying_axis}-line must span a positive distance, got {line_span}."
        )
    return float(np_module.trapezoid(point_gains, axis_coordinates))


def compute_single_pass_gain_factor(distance_gain: float, np_module) -> float:
    return float(np_module.exp(distance_gain))


def compute_net_gain_factor(single_pass_gain: float, *, back_reflection: bool, reflectivity: float) -> float:
    if back_reflection:
        return float(reflectivity) * single_pass_gain * single_pass_gain
    return single_pass_gain


def build_output_prefix(args: argparse.Namespace, varying_axis: str, fixed_axes: tuple[str, str]) -> Path:
    if args.output_prefix is not None:
        return args.output_prefix.expanduser().resolve()

    base_dir = Path(__file__).resolve().parent
    labels = "_".join(
        f"{axis}_{format_coordinate_for_filename(float(getattr(args, axis)))}"
        for axis in fixed_axes
    )
    return base_dir / f"net_gain_factor_{varying_axis}_{labels}"


def format_coordinate_for_filename(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def ensure_parent_directory(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    varying_axis, fixed_axes, requested_coords = resolve_direction_and_coords(args)
    vtk_files = find_vtk_files(input_dir)

    try:
        import numpy as np
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "NumPy is required to run this script. Install project dependencies first."
        ) from exc

    if not hasattr(np, "trapezoid"):
        raise SystemExit("This script requires a NumPy version that provides np.trapezoid.")

    try:
        import matplotlib
        if not args.show and not os.environ.get("DISPLAY"):
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "matplotlib is required to generate the gain-versus-timestep plot."
        ) from exc

    reference_points, reference_gains, point_count = parse_ascii_vtk_scalar(
        vtk_files[0][1],
        args.field,
        include_points=True,
    )
    assert reference_points is not None

    selection = select_line(reference_points, varying_axis, fixed_axes, requested_coords)
    results: list[tuple[int, float, float, float]] = []
    first_distance_gain = compute_distance_gain(selection, reference_gains, np)
    first_single_pass = compute_single_pass_gain_factor(first_distance_gain, np)
    results.append(
        (
            vtk_files[0][0],
            first_distance_gain,
            first_single_pass,
            compute_net_gain_factor(
                first_single_pass,
                back_reflection=args.back_reflection,
                reflectivity=args.reflectivity,
            ),
        )
    )

    for timestep, path in vtk_files[1:]:
        _, gains, current_point_count = parse_ascii_vtk_scalar(path, args.field, include_points=False)
        if current_point_count != point_count:
            raise ValueError(
                f"{path} has {current_point_count} point values, expected {point_count}."
            )
        distance_gain = compute_distance_gain(selection, gains, np)
        single_pass_gain = compute_single_pass_gain_factor(distance_gain, np)
        results.append(
            (
                timestep,
                distance_gain,
                single_pass_gain,
                compute_net_gain_factor(
                    single_pass_gain,
                    back_reflection=args.back_reflection,
                    reflectivity=args.reflectivity,
                ),
            )
        )

    output_prefix = build_output_prefix(args, varying_axis, fixed_axes)
    csv_path = output_prefix.with_suffix(".csv")
    png_path = output_prefix.with_suffix(".png")
    ensure_parent_directory(csv_path)
    ensure_parent_directory(png_path)

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["timestep", "distance_gain", "single_pass_gain", "net_gain_factor"])
        writer.writerows(results)

    timesteps = np.asarray([row[0] for row in results], dtype=float)
    net_gain_factor = np.asarray([row[3] for row in results], dtype=float)

    plt.figure(figsize=(9, 5))
    plt.plot(timesteps, net_gain_factor, marker="o", markersize=3, linewidth=1.5)
    plt.xlabel("Timestep")
    plt.ylabel("Net gain factor")
    plt.title(
        f"Net gain factor along {varying_axis} at "
        f"{fixed_axes[0]}={selection.actual_coords[0]:.6g}, "
        f"{fixed_axes[1]}={selection.actual_coords[1]:.6g}"
    )
    plt.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig(png_path, dpi=200)
    if args.show:
        plt.show()
    plt.close()

    requested_text = ", ".join(
        f"{axis}={value:.6g}" for axis, value in zip(fixed_axes, selection.requested_coords)
    )
    actual_text = ", ".join(
        f"{axis}={value:.6g}" for axis, value in zip(fixed_axes, selection.actual_coords)
    )
    mode_text = (
        f"round-trip reflection enabled (R={args.reflectivity:.6g})"
        if args.back_reflection
        else "single-pass only"
    )
    print(f"Processed {len(results)} files from {input_dir}")
    print(f"Integrated local gain along {varying_axis} using requested coordinates: {requested_text}")
    print(f"Nearest available line used: {actual_text}")
    print(f"Samples along {varying_axis}: {len(selection.coordinate_groups)}")
    print(f"Mode: {mode_text}")
    print(f"Saved CSV: {csv_path}")
    print(f"Saved plot: {png_path}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
