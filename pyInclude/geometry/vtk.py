# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Read and write legacy ASCII VTK wedge meshes for gain media."""

from __future__ import annotations

from pathlib import Path

import numpy as np


_VTK_TO_DTYPE = {
    "float": np.float64,
    "double": np.float64,
    "int": np.int32,
    "unsigned_int": np.uint32,
    "uint": np.uint32,
}


def _dtypeName(values):
    dtype = np.asarray(values).dtype
    if np.issubdtype(dtype, np.unsignedinteger):
        return "unsigned_int"
    if np.issubdtype(dtype, np.integer):
        return "int"
    return "double"


def _tokens(path):
    return Path(path).read_text(encoding="utf-8").split()


class _TokenReader:
    """Small token stream used by the legacy VTK parser."""

    def __init__(self, tokens):
        self.tokens = tokens
        self.index = 0

    def done(self):
        return self.index >= len(self.tokens)

    def next(self):
        value = self.tokens[self.index]
        self.index += 1
        return value

    def take(self, count, dtype=float):
        values = self.tokens[self.index:self.index + count]
        self.index += count
        return np.asarray(values, dtype=dtype)


def _parseVtk(path):
    reader = _TokenReader(_tokens(path))
    points = None
    cells = []
    cellTypes = None
    pointData = {}
    cellData = {}
    fields = {}
    activeData = None
    activeCount = 0

    while not reader.done():
        token = reader.next()
        upper = token.upper()

        if upper == "POINTS":
            count = int(reader.next())
            reader.next()
            points = reader.take(count * 3, float).reshape(count, 3)
        elif upper == "CELLS":
            count = int(reader.next())
            reader.next()
            cells = []
            for _ in range(count):
                width = int(reader.next())
                cells.append(reader.take(width, int))
        elif upper == "CELL_TYPES":
            count = int(reader.next())
            cellTypes = reader.take(count, int)
        elif upper == "POINT_DATA":
            activeData = pointData
            activeCount = int(reader.next())
        elif upper == "CELL_DATA":
            activeData = cellData
            activeCount = int(reader.next())
        elif upper == "SCALARS":
            if activeData is None:
                raise ValueError("SCALARS section must follow POINT_DATA or CELL_DATA")
            name = reader.next()
            dtype = _VTK_TO_DTYPE.get(reader.next().lower(), np.float64)
            components = 1
            if not reader.done() and reader.tokens[reader.index].upper() != "LOOKUP_TABLE":
                components = int(reader.next())
            if not reader.done() and reader.tokens[reader.index].upper() == "LOOKUP_TABLE":
                reader.next()
                reader.next()
            values = reader.take(activeCount * components, dtype)
            activeData[name] = values.reshape(activeCount, components).squeeze()
        elif upper == "FIELD":
            reader.next()
            arrayCount = int(reader.next())
            for _ in range(arrayCount):
                name = reader.next()
                components = int(reader.next())
                tuples = int(reader.next())
                dtype = _VTK_TO_DTYPE.get(reader.next().lower(), np.float64)
                values = reader.take(components * tuples, dtype)
                fields[name] = values.reshape(tuples, components).squeeze()

    if points is None or not cells:
        raise ValueError(f"'{path}' does not contain POINTS/CELLS data")
    return points, cells, cellTypes, pointData, cellData, fields


def _topologyFromUnstructuredGrid(path):
    points3d, cells, cellTypes, pointData, cellData, fields = _parseVtk(path)
    if cellTypes is not None and not np.all(np.asarray(cellTypes) == 13):
        raise ValueError("only VTK wedge cells (cell type 13) are supported")

    zValues = np.unique(points3d[:, 2])
    zValues.sort()
    levels = int(zValues.size)
    if levels < 2:
        raise ValueError("VTK topology requires at least two z levels")
    thickness = float(zValues[1] - zValues[0])

    xy = points3d[:, :2]
    if points3d.shape[0] % levels != 0:
        raise ValueError("VTK point count must be divisible by the number of z levels")
    pointsPerLevel = points3d.shape[0] // levels
    points2d = xy[:pointsPerLevel]

    cellsPerLayer = len(cells) // (levels - 1)
    triangles = np.empty((cellsPerLayer, 3), dtype=np.uint32)
    for cellIndex, cell in enumerate(cells[:cellsPerLayer]):
        if len(cell) != 6:
            raise ValueError("only six-node wedge cells are supported")
        triangles[cellIndex, :] = np.asarray(cell[:3], dtype=np.uint32) % np.uint32(pointsPerLevel)

    physical = {}
    physical.update(pointData)
    physical.update(cellData)
    for name in (
        "claddingCellTypes",
        "refractiveIndices",
        "reflectivities",
        "nTot",
        "crystalTFluo",
        "claddingNumber",
        "claddingAbsorption",
    ):
        if name in fields:
            physical[name] = fields[name]

    return points2d, triangles, levels, thickness, physical


def topologyFromVtk(path, topologyCls):
    """Create a ``MeshTopology`` from legacy VTK wedge cells."""
    points, triangles, levels, thickness, _ = _topologyFromUnstructuredGrid(path)
    return topologyCls(
        points=points,
        trianglePointIndices=triangles,
        levels=levels,
        thickness=thickness,
        metadata={"source": str(path), "format": "vtk"},
    )


def gainMediumFromVtk(path, topologyCls, gainMediumCls, *, numberOfLevels=None, thickness=None):
    """Load topology plus beta/material arrays from a legacy VTK file."""
    topology = topologyFromVtk(path, topologyCls)
    if numberOfLevels is not None:
        topology.numberOfLevels(numberOfLevels)
    if thickness is not None:
        topology.withThickness(thickness)
    _, _, _, _, physical = _topologyFromUnstructuredGrid(path)
    medium = gainMediumCls(topology=topology)
    for name, value in physical.items():
        medium.set(name, value)
    return medium


def writeGainMediumVtk(path, gainMedium):
    """Write a ``GainMedium`` as legacy ASCII VTK wedge cells.

    ``betaCells`` is written as point data and ``betaVolume`` as cell data;
    scalar material properties are stored in a ``FIELD`` block.
    """
    topology = gainMedium.topology
    topology._require_levels()
    topology._require_thickness()

    points = np.asarray(topology.points, dtype=np.float64)
    triangles = np.asarray(topology.trianglePointIndices, dtype=np.uint32)
    levels = int(topology.levels)
    zValues = topology.levelCoordinates()
    pointCount = topology.numberOfPoints * levels
    cellCount = topology.numberOfTriangles * (levels - 1)

    fieldArrays = {
        "claddingCellTypes": gainMedium.get("claddingCellTypes").value,
        "refractiveIndices": gainMedium.get("refractiveIndices").value,
        "reflectivities": gainMedium.get("reflectivities").value,
        "nTot": np.asarray([gainMedium.get("nTot").value], dtype=np.float64),
        "crystalTFluo": np.asarray([gainMedium.get("crystalTFluo").value], dtype=np.float64),
        "claddingNumber": np.asarray([gainMedium.get("claddingNumber").value], dtype=np.uint32),
        "claddingAbsorption": np.asarray([gainMedium.get("claddingAbsorption").value], dtype=np.float64),
    }
    pointArrays = {"betaCells": gainMedium.get("betaCells").value}
    cellArrays = {"betaVolume": gainMedium.get("betaVolume").value}

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 2.0\n")
        handle.write("HASEonGPU gain medium input\n")
        handle.write("ASCII\n")
        handle.write("DATASET UNSTRUCTURED_GRID\n")
        handle.write(f"POINTS {pointCount} double\n")
        for z in zValues:
            for x, y in points:
                handle.write(f"{x:.17g} {y:.17g} {z:.17g}\n")

        handle.write(f"CELLS {cellCount} {cellCount * 7}\n")
        for level in range(levels - 1):
            lower = level * topology.numberOfPoints
            upper = (level + 1) * topology.numberOfPoints
            for tri in triangles:
                ids = [int(vertex) for vertex in tri]
                handle.write(
                    "6 "
                    f"{ids[0] + lower} {ids[1] + lower} {ids[2] + lower} "
                    f"{ids[0] + upper} {ids[1] + upper} {ids[2] + upper}\n"
                )

        handle.write(f"CELL_TYPES {cellCount}\n")
        for _ in range(cellCount):
            handle.write("13\n")

        handle.write(f"FIELD HASEonGPU {len(fieldArrays)}\n")
        for name, values in fieldArrays.items():
            arr = np.asarray(values).reshape(-1, order="F")
            handle.write(f"{name} 1 {arr.size} {_dtypeName(arr)}\n")
            handle.write(" ".join(str(value) for value in arr.tolist()) + "\n")

        for section, arrays, count in (
            ("POINT_DATA", pointArrays, pointCount),
            ("CELL_DATA", cellArrays, cellCount),
        ):
            handle.write(f"{section} {count}\n")
            for name, values in arrays.items():
                arr = np.asarray(values, dtype=np.float64).reshape(-1, order="F")
                if arr.size != count:
                    raise ValueError(f"{name} has {arr.size} values, expected {count}")
                handle.write(f"SCALARS {name} double 1\n")
                handle.write("LOOKUP_TABLE default\n")
                for value in arr:
                    handle.write(f"{float(value):.17g}\n")
    return path
