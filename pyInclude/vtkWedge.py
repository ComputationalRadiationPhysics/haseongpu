# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from __future__ import annotations

from pathlib import Path

import numpy as np

from .geometry import GainMedium, MeshTopology


def _topology_from_geometry(geometry):
    if isinstance(geometry, GainMedium):
        return geometry.topology
    if isinstance(geometry, MeshTopology):
        return geometry
    raise TypeError("geometry must be a GainMedium or MeshTopology")


def _vtk_filename(file_name, state=None, field="phiAse"):
    text = str(file_name)
    if state is not None:
        text = text.format(
            step=getattr(state, "step", 0),
            time=getattr(state, "time", 0.0),
            field=field,
        )
    if not text.endswith(".vtk"):
        text += ".vtk"
    return Path(text)


def _as_named_array(data, field, scalar_name):
    if hasattr(data, field):
        values = getattr(data, field)
        name = scalar_name or field
    else:
        values = data
        name = scalar_name or "scalars"
    if values is None:
        raise ValueError(f"no data available for VTK field '{field}'")
    return np.asarray(values, dtype=np.float64), name


def _resolve_data_shape(values, topology):
    point_shape = (topology.numberOfPoints, int(topology.levels))
    cell_shape = (topology.numberOfTriangles, int(topology.levels) - 1)
    if values.shape == point_shape:
        return "POINT_DATA", point_shape
    if values.shape == cell_shape:
        return "CELL_DATA", cell_shape
    if values.ndim == 1 and values.size == int(np.prod(point_shape)):
        return "POINT_DATA", point_shape
    if values.ndim == 1 and values.size == int(np.prod(cell_shape)):
        return "CELL_DATA", cell_shape
    raise ValueError(
        "VTK wedge data must match point data shape "
        f"{point_shape} or cell data shape {cell_shape}, got {values.shape}"
    )


def _write_ascii_wedge(file_name, values, topology, scalar_name="scalars"):
    topology._require_levels()
    topology._require_thickness()

    values = np.asarray(values, dtype=np.float64)
    data_kind, data_shape = _resolve_data_shape(values, topology)
    values = values.reshape(data_shape, order="F")

    points = np.asarray(topology.points, dtype=np.float64)
    triangles = np.asarray(topology.trianglePointIndices, dtype=np.uint32)
    levels = int(topology.levels)
    z_values = topology.levelCoordinates()

    number_of_points = topology.numberOfPoints * levels
    number_of_cells = topology.numberOfTriangles * (levels - 1)

    path = Path(file_name)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# vtk DataFile Version 2.0\n")
        handle.write("HASEonGPU wedge output\n")
        handle.write("ASCII\n")
        handle.write("DATASET UNSTRUCTURED_GRID\n")
        handle.write(f"POINTS {number_of_points} float\n")
        for z in z_values:
            for x, y in points:
                handle.write(f"{x:.17g} {y:.17g} {z:.17g}\n")

        handle.write(f"CELLS {number_of_cells} {number_of_cells * 7}\n")
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

        handle.write(f"CELL_TYPES {number_of_cells}\n")
        for _ in range(number_of_cells):
            handle.write("13\n")

        handle.write(f"{data_kind} {values.size}\n")
        handle.write(f"SCALARS {scalar_name} float 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for value in values.reshape(-1, order="F"):
            handle.write(f"{float(value):.17g}\n")

    return path


def _legacy_vtk_wedge(file_name, values, points, triangle_point_indices, mesh_z, z_mesh, scalar_name):
    topology = MeshTopology(
        points=np.asarray(points, dtype=np.float64),
        trianglePointIndices=np.asarray(triangle_point_indices, dtype=np.uint32),
        levels=int(mesh_z),
        thickness=float(z_mesh),
    )
    return _write_ascii_wedge(file_name, values, topology, scalar_name=scalar_name)


def vtkWedge(file_name, data=None, geometry=None, field="phiAse", scalar_name=None, every=1):
    """
        Write HASEonGPU wedge data to a legacy ASCII VTK file.
    """

    if isinstance(data, (GainMedium, MeshTopology)) and geometry is None:
        topology = _topology_from_geometry(data)

        def write_state(state):
            if int(every) <= 0:
                raise ValueError("every must be a positive integer")
            if getattr(state, "step", 0) % int(every) != 0:
                return None
            values, name = _as_named_array(state, field, scalar_name)
            return _write_ascii_wedge(_vtk_filename(file_name, state, field), values, topology, name)

        return write_state

    if geometry is None:
        raise TypeError("new-interface vtkWedge requires a GainMedium or MeshTopology")

    topology = _topology_from_geometry(geometry)
    values, name = _as_named_array(data, field, scalar_name)
    return _write_ascii_wedge(_vtk_filename(file_name, data, field), values, topology, name)
