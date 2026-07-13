#!/usr/bin/env python3
# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Convert ASCII VTK topology between legacy wedge and Tet4 cells."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pyInclude.geometry.vtk import _dtypeName, _parseVtk  # noqa: E402


VTK_TETRA = np.uint32(10)
VTK_WEDGE = np.uint32(13)


def _tet_volume(points):
    a, b, c, d = np.asarray(points, dtype=np.float64)
    return abs(float(np.dot(b - a, np.cross(c - a, d - a)))) / 6.0


def _tetrahedra_from_wedge(cell):
    # Use a global base-vertex order so neighboring prisms choose the same
    # diagonal on every shared quadrilateral face.
    vertices = np.asarray(cell, dtype=np.uint32)
    order = np.argsort(vertices[:3])
    a, b, c = [int(vertex) for vertex in vertices[:3][order]]
    d, e, f = [int(vertex) for vertex in vertices[3:][order]]
    return np.asarray(
        [
            [a, b, c, d],
            [b, e, c, d],
            [c, e, f, d],
        ],
        dtype=np.uint32,
    )


def _unique_in_order(values):
    seen = set()
    ordered = []
    for value in values:
        value = int(value)
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def _same_xy(a, b, *, atol=1e-12):
    return np.allclose(a[:2], b[:2], atol=atol, rtol=0.0)


def _wedge_from_tetrahedra(points, tet_group):
    flat = tet_group.reshape(-1)
    vertices = _unique_in_order(flat)
    if len(vertices) != 6:
        raise ValueError("three Tet4 cells must reconstruct exactly six wedge vertices")

    z_values = np.asarray([points[index, 2] for index in vertices], dtype=np.float64)
    z_min = float(z_values.min())
    z_max = float(z_values.max())
    lower = [index for index in vertices if np.isclose(points[index, 2], z_min)]
    upper_candidates = [index for index in vertices if np.isclose(points[index, 2], z_max)]
    if len(lower) != 3 or len(upper_candidates) != 3:
        raise ValueError("Tet4 group does not contain three lower and three upper wedge vertices")

    upper = []
    for lower_index in lower:
        matches = [index for index in upper_candidates if _same_xy(points[lower_index], points[index])]
        if len(matches) != 1:
            raise ValueError("could not match lower wedge vertex to an upper vertex with the same x/y")
        upper.append(matches[0])
    return np.asarray([*lower, *upper], dtype=np.uint32)


def _infer_structured_fields(points, fields):
    result = dict(fields)
    if {
        "structuredNumberOfPoints",
        "structuredNumberOfLevels",
        "structuredThickness",
    }.issubset(result):
        return result

    z_values = np.unique(np.asarray(points[:, 2], dtype=np.float64))
    z_values.sort()
    if z_values.size < 2 or points.shape[0] % z_values.size != 0:
        return result
    points_per_level = points.shape[0] // z_values.size
    for level, z_value in enumerate(z_values):
        start = level * points_per_level
        stop = start + points_per_level
        if not np.allclose(points[start:stop, 2], z_value, atol=1e-12, rtol=0.0):
            return result
    result.setdefault("structuredNumberOfPoints", np.asarray([points_per_level], dtype=np.uint32))
    result.setdefault("structuredNumberOfLevels", np.asarray([z_values.size], dtype=np.uint32))
    result.setdefault("structuredThickness", np.asarray([z_values[1] - z_values[0]], dtype=np.float64))
    return result


def _repeat_cell_array(values, cell_count):
    arr = np.asarray(values)
    if arr.ndim == 0:
        return arr
    if arr.shape[0] != cell_count:
        return arr
    return np.repeat(arr, 3, axis=0)


def _project_cell_array(values, tet_volumes):
    arr = np.asarray(values)
    if arr.ndim == 0:
        return arr
    if arr.shape[0] != tet_volumes.size:
        return arr
    grouped = arr.reshape((tet_volumes.size // 3, 3, *arr.shape[1:]))
    constant = np.all(grouped == grouped[:, :1, ...], axis=1)
    weights = tet_volumes.reshape((-1, 3))
    weighted = np.sum(grouped * weights.reshape((weights.shape[0], 3, *([1] * (arr.ndim - 1)))), axis=1) / weights.sum(axis=1).reshape((-1, *([1] * (arr.ndim - 1))))
    if np.all(constant):
        return grouped[:, 0, ...].astype(arr.dtype, copy=False)
    if np.issubdtype(arr.dtype, np.integer):
        return weighted
    weighted[constant] = grouped[:, 0, ...][constant]
    return weighted


def _write_ascii_vtk(path, *, title, points, cells, cell_types, fields=None, point_data=None, cell_data=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or {}
    point_data = point_data or {}
    cell_data = cell_data or {}

    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 2.0\n")
        handle.write(f"{title}\n")
        handle.write("ASCII\n")
        handle.write("DATASET UNSTRUCTURED_GRID\n")
        handle.write(f"POINTS {points.shape[0]} double\n")
        for x, y, z in np.asarray(points, dtype=np.float64):
            handle.write(f"{x:.17g} {y:.17g} {z:.17g}\n")

        width = np.asarray(cells).shape[1]
        handle.write(f"CELLS {len(cells)} {len(cells) * (width + 1)}\n")
        for cell in cells:
            handle.write(str(width) + " " + " ".join(str(int(vertex)) for vertex in cell) + "\n")

        handle.write(f"CELL_TYPES {len(cell_types)}\n")
        for cell_type in cell_types:
            handle.write(f"{int(cell_type)}\n")

        if fields:
            handle.write(f"FIELD HASEonGPU {len(fields)}\n")
            for name, values in fields.items():
                arr = np.asarray(values)
                if arr.ndim == 0:
                    arr = arr.reshape(1)
                components = 1 if arr.ndim == 1 else int(np.prod(arr.shape[1:]))
                tuples = arr.size // components
                handle.write(f"{name} {components} {tuples} {_dtypeName(arr)}\n")
                handle.write(" ".join(str(value) for value in arr.reshape(-1).tolist()) + "\n")

        for label, data in (("POINT_DATA", point_data), ("CELL_DATA", cell_data)):
            if not data:
                continue
            count = points.shape[0] if label == "POINT_DATA" else len(cells)
            handle.write(f"{label} {count}\n")
            for name, values in data.items():
                arr = np.asarray(values)
                if arr.shape[0] != count:
                    continue
                components = 1 if arr.ndim == 1 else int(np.prod(arr.shape[1:]))
                handle.write(f"SCALARS {name} {_dtypeName(arr)} {components}\n")
                handle.write("LOOKUP_TABLE default\n")
                handle.write("\n".join(str(value) for value in arr.reshape(-1).tolist()) + "\n")
    return path


def wedgeToTet4(inputPath, outputPath):
    points, cells, cell_types, point_data, cell_data, fields = _parseVtk(inputPath)
    cell_types = np.asarray(cell_types, dtype=np.uint32)
    if not np.all(cell_types == VTK_WEDGE):
        raise ValueError("wedge-to-tet4 expects only VTK wedge cells (cell type 13)")
    wedge_cells = np.asarray(cells, dtype=np.uint32)
    if wedge_cells.ndim != 2 or wedge_cells.shape[1] != 6:
        raise ValueError("wedge-to-tet4 expects six-node wedge connectivity")

    tet_cells = np.vstack([_tetrahedra_from_wedge(cell) for cell in wedge_cells])
    converted_fields = {
        name: _repeat_cell_array(values, wedge_cells.shape[0])
        for name, values in _infer_structured_fields(points, fields).items()
    }
    converted_cell_data = {
        name: _repeat_cell_array(values, wedge_cells.shape[0])
        for name, values in cell_data.items()
    }
    return _write_ascii_vtk(
        outputPath,
        title="HASEonGPU wedge converted to Tet4",
        points=np.asarray(points, dtype=np.float64),
        cells=tet_cells,
        cell_types=np.full(tet_cells.shape[0], VTK_TETRA, dtype=np.uint32),
        fields=converted_fields,
        point_data=point_data,
        cell_data=converted_cell_data,
    )


def tet4ToWedge(inputPath, outputPath):
    points, cells, cell_types, point_data, cell_data, fields = _parseVtk(inputPath)
    cell_types = np.asarray(cell_types, dtype=np.uint32)
    if not np.all(cell_types == VTK_TETRA):
        raise ValueError("tet4-to-wedge expects only VTK Tet4 cells (cell type 10)")
    tet_cells = np.asarray(cells, dtype=np.uint32)
    if tet_cells.ndim != 2 or tet_cells.shape[1] != 4:
        raise ValueError("tet4-to-wedge expects four-node Tet4 connectivity")
    if tet_cells.shape[0] % 3 != 0:
        raise ValueError("tet4-to-wedge expects three Tet4 cells per legacy wedge")

    wedge_cells = np.asarray(
        [_wedge_from_tetrahedra(points, group) for group in tet_cells.reshape((-1, 3, 4))],
        dtype=np.uint32,
    )
    tet_volumes = np.asarray(
        [_tet_volume(points[cell]) for cell in tet_cells],
        dtype=np.float64,
    )
    converted_fields = {
        name: _project_cell_array(values, tet_volumes)
        for name, values in fields.items()
    }
    converted_cell_data = {
        name: _project_cell_array(values, tet_volumes)
        for name, values in cell_data.items()
    }
    return _write_ascii_vtk(
        outputPath,
        title="HASEonGPU Tet4 converted to wedge",
        points=np.asarray(points, dtype=np.float64),
        cells=wedge_cells,
        cell_types=np.full(wedge_cells.shape[0], VTK_WEDGE, dtype=np.uint32),
        fields=converted_fields,
        point_data=point_data,
        cell_data=converted_cell_data,
    )


def defaultOutputPath(inputPath, direction):
    path = Path(inputPath)
    suffix = "_tet4" if direction == "wedge-to-tet4" else "_wedge"
    return path.with_name(f"{path.stem}{suffix}{path.suffix or '.vtk'}")


def convertVtk(inputPath, outputPath=None, *, direction):
    outputPath = defaultOutputPath(inputPath, direction) if outputPath is None else Path(outputPath)
    if direction == "wedge-to-tet4":
        return wedgeToTet4(inputPath, outputPath)
    if direction == "tet4-to-wedge":
        return tet4ToWedge(inputPath, outputPath)
    raise ValueError(f"unknown conversion direction: {direction}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input ASCII VTK file.")
    parser.add_argument("output", type=Path, nargs="?", help="Output ASCII VTK file. Defaults to a suffixed filename.")
    parser.add_argument("--direction", required=True, choices=("wedge-to-tet4", "tet4-to-wedge"))
    args = parser.parse_args(argv)

    output = convertVtk(args.input, args.output, direction=args.direction)
    print(output)


if __name__ == "__main__":
    main()
