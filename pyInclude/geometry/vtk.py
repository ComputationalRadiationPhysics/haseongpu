# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Read and write ASCII VTK Tet4 meshes for gain media."""

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
    """Small token stream used by the ASCII VTK parser."""

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


def _tet4TopologyFromUnstructuredGrid(path):
    points, cells, cellTypes, pointData, cellData, fields = _parseVtk(path)
    if cellTypes is None:
        raise ValueError("VTK Tet4 topology requires CELL_TYPES")
    cellTypes = np.asarray(cellTypes, dtype=np.uint32)
    if not np.all(cellTypes == 10):
        raise ValueError("only VTK Tet4 cells (cell type 10) are supported")
    if any(len(cell) != 4 for cell in cells):
        raise ValueError("only four-node tetrahedral VTK cells are supported")

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
        "betaCells",
        "betaVolume",
        "cellDomains",
        "faceBoundaries",
        "structuredNumberOfPoints",
        "structuredNumberOfLevels",
        "structuredThickness",
    ):
        if name in fields:
            physical[name] = fields[name]

    return points, np.asarray(cells, dtype=np.uint32), cellTypes, physical


def _structured_metadata_from_vtk(points, physical):
    structured = {}
    if "structuredNumberOfPoints" in physical:
        structured["numberOfPoints"] = int(np.asarray(physical["structuredNumberOfPoints"]).reshape(-1)[0])
    if "structuredNumberOfLevels" in physical:
        structured["numberOfLevels"] = int(np.asarray(physical["structuredNumberOfLevels"]).reshape(-1)[0])
    if "structuredThickness" in physical:
        structured["thickness"] = float(np.asarray(physical["structuredThickness"]).reshape(-1)[0])

    if not structured:
        points_array = np.asarray(points, dtype=np.float64)
        z_values = np.unique(points_array[:, 2])
        z_values.sort()
        if z_values.size >= 2 and points_array.shape[0] % z_values.size == 0:
            per_level = points_array.shape[0] // z_values.size
            xy_reference = None
            layered = True
            for z in z_values:
                level_xy = points_array[np.isclose(points_array[:, 2], z, atol=1e-12, rtol=0.0), :2]
                if level_xy.shape[0] != per_level:
                    layered = False
                    break
                level_xy = level_xy[np.lexsort((level_xy[:, 1], level_xy[:, 0]))]
                if xy_reference is None:
                    xy_reference = level_xy
                elif not np.allclose(level_xy, xy_reference, atol=1e-12, rtol=0.0):
                    layered = False
                    break
            if layered:
                structured["numberOfPoints"] = int(per_level)
                structured["numberOfLevels"] = int(z_values.size)
                structured["thickness"] = float(z_values[1] - z_values[0])
    return {"structured": structured} if structured else {}


def volumeTopologyFromVtk(path, topologyCls):
    """Create a Tet4 ``VolumeTopology`` from an ASCII VTK unstructured grid."""
    points, cells, cellTypes, physical = _tet4TopologyFromUnstructuredGrid(path)
    domains = physical.get("cellDomains")
    if domains is not None:
        domains = np.asarray(domains, dtype=np.int32).reshape(cells.shape[0])
    boundaries = physical.get("faceBoundaries")
    if boundaries is not None:
        boundaries = np.asarray(boundaries, dtype=np.int32).reshape((cells.shape[0], 4))
    topology = topologyCls(
        points=points,
        cellPointIndices=cells,
        cellTypes=cellTypes,
        cellDomains=domains,
        faceBoundaries=boundaries,
        metadata={"source": str(path), "format": "vtk", **_structured_metadata_from_vtk(points, physical)},
    )
    beta_cells = physical.get("betaCells")
    if beta_cells is not None and np.asarray(beta_cells).size == points.shape[0]:
        topology.samplePoints = np.asarray(points, dtype=np.float64)
    return topology


def topologyFromVtk(path, topologyCls):
    """Create a Tet4 ``VolumeTopology`` from an ASCII VTK unstructured grid."""
    return volumeTopologyFromVtk(path, topologyCls)


def gainMediumFromVtk(path, topologyCls, gainMediumCls, *, numberOfLevels=None, thickness=None):
    """Load Tet4 topology plus material arrays from a VTK file."""
    from ..openpmd import backendFlat

    if numberOfLevels is not None or thickness is not None:
        raise ValueError("Tet4 VTK gain-medium import does not use numberOfLevels or thickness")
    points, cells, cellTypes, physical = _tet4TopologyFromUnstructuredGrid(path)
    topology = topologyCls(
        points=points,
        cellPointIndices=cells,
        cellTypes=cellTypes,
        cellDomains=(
            None
            if physical.get("cellDomains") is None
            else np.asarray(physical["cellDomains"], dtype=np.int32).reshape(cells.shape[0])
        ),
        faceBoundaries=(
            None
            if physical.get("faceBoundaries") is None
            else np.asarray(physical["faceBoundaries"], dtype=np.int32).reshape((cells.shape[0], 4))
        ),
        metadata={"source": str(path), "format": "vtk", **_structured_metadata_from_vtk(points, physical)},
    )
    beta_cells = physical.get("betaCells")
    if beta_cells is not None and np.asarray(beta_cells).size == points.shape[0]:
        topology.samplePoints = np.asarray(points, dtype=np.float64)
    medium = gainMediumCls(topology=topology)
    for name, value in physical.items():
        if name in {"cellDomains", "faceBoundaries", "structuredNumberOfPoints", "structuredNumberOfLevels", "structuredThickness"}:
            continue
        try:
            arr = np.asarray(value)
            if medium.get(name).expectedShape == ():
                medium.set(name, arr.reshape(-1)[0])
            else:
                medium.set(name, backendFlat(arr.reshape(-1)) if arr.ndim <= 1 else value)
        except KeyError:
            continue
    return medium


def writeGainMediumVtk(path, gainMedium):
    """Write a ``GainMedium`` as an ASCII VTK Tet4 unstructured grid."""
    topology = gainMedium.topology
    if not hasattr(topology, "cellPointIndices"):
        raise TypeError("writeGainMediumVtk requires a Tet4 VolumeTopology")

    points = np.asarray(topology.points, dtype=np.float64)
    cells = np.asarray(topology.cellPointIndices, dtype=np.uint32)
    cellTypes = np.asarray(topology.cellTypes, dtype=np.uint32)
    if not np.all(cellTypes == 10):
        raise ValueError("writeGainMediumVtk supports Tet4 cells only")

    structured = topology.metadata.get("structured", {}) if isinstance(topology.metadata, dict) else {}
    fieldArrays = {
        "cellDomains": topology.cellDomains,
        "faceBoundaries": topology.faceBoundaries,
        "structuredNumberOfPoints": np.asarray([structured.get("numberOfPoints", topology.numberOfSamplePoints)], dtype=np.uint32),
        "structuredNumberOfLevels": np.asarray([structured.get("numberOfLevels", 1)], dtype=np.uint32),
        "structuredThickness": np.asarray([structured.get("thickness", 0.0)], dtype=np.float64),
        "claddingCellTypes": gainMedium.get("claddingCellTypes").value,
        "refractiveIndices": gainMedium.get("refractiveIndices").value,
        "reflectivities": gainMedium.get("reflectivities").value,
        "nTot": np.asarray([gainMedium.get("nTot").value], dtype=np.float64),
        "crystalTFluo": np.asarray([gainMedium.get("crystalTFluo").value], dtype=np.float64),
        "claddingNumber": np.asarray([gainMedium.get("claddingNumber").value], dtype=np.uint32),
        "claddingAbsorption": np.asarray([gainMedium.get("claddingAbsorption").value], dtype=np.float64),
    }
    pointArrays = {}
    cellArrays = {"betaVolume": gainMedium.get("betaVolume").value}
    betaCells = gainMedium.get("betaCells").value
    if betaCells.size == topology.numberOfPoints:
        pointArrays["betaCells"] = betaCells
    elif betaCells.size == topology.numberOfCells:
        cellArrays["betaCells"] = betaCells

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 2.0\n")
        handle.write("HASEonGPU Tet4 gain medium input\n")
        handle.write("ASCII\n")
        handle.write("DATASET UNSTRUCTURED_GRID\n")
        handle.write(f"POINTS {points.shape[0]} double\n")
        for x, y, z in points:
            handle.write(f"{x:.17g} {y:.17g} {z:.17g}\n")

        handle.write(f"CELLS {cells.shape[0]} {cells.shape[0] * 5}\n")
        for cell in cells:
            handle.write("4 " + " ".join(str(int(vertex)) for vertex in cell) + "\n")

        handle.write(f"CELL_TYPES {cells.shape[0]}\n")
        for _ in range(cells.shape[0]):
            handle.write("10\n")

        handle.write(f"FIELD HASEonGPU {len(fieldArrays)}\n")
        for name, values in fieldArrays.items():
            arr = np.asarray(values).reshape(-1)
            handle.write(f"{name} 1 {arr.size} {_dtypeName(arr)}\n")
            handle.write(" ".join(str(value) for value in arr.tolist()) + "\n")

        if pointArrays:
            handle.write(f"POINT_DATA {points.shape[0]}\n")
            for name, values in pointArrays.items():
                arr = np.asarray(values).reshape(-1)
                if arr.size != points.shape[0]:
                    continue
                handle.write(f"SCALARS {name} {_dtypeName(arr)} 1\n")
                handle.write("LOOKUP_TABLE default\n")
                handle.write("\n".join(str(value) for value in arr.tolist()) + "\n")

        handle.write(f"CELL_DATA {cells.shape[0]}\n")
        for name, values in cellArrays.items():
            arr = np.asarray(values).reshape(-1)
            if arr.size != cells.shape[0]:
                continue
            handle.write(f"SCALARS {name} {_dtypeName(arr)} 1\n")
            handle.write("LOOKUP_TABLE default\n")
            handle.write("\n".join(str(value) for value in arr.tolist()) + "\n")
    return path
