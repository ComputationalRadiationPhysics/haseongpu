# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Legacy ASCII VTK writer for wedge-prism HASEonGPU data."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from .geometry import GainMedium, MeshTopology


def _topology_from_geometry(geometry):
    """Extract ``MeshTopology`` from a topology or owning ``GainMedium``."""
    if isinstance(geometry, GainMedium):
        return geometry.topology
    if isinstance(geometry, MeshTopology):
        return geometry
    raise TypeError("geometry must be a GainMedium or MeshTopology")


def _field_label(field):
    """Create a stable filename label for one or more VTK fields."""
    if isinstance(field, Mapping):
        names = field.keys()
    elif isinstance(field, Sequence) and not isinstance(field, (str, bytes)):
        names = field
    else:
        names = (field,)
    return "_".join(str(name) for name in names)


def _vtk_filename(file_name, state=None, field="phiAse"):
    """Resolve ``{step}``, ``{time}``, and ``{field}`` filename placeholders."""
    text = str(file_name)
    if state is not None:
        text = text.format(
            step=getattr(state, "step", 0),
            time=getattr(state, "time", 0.0),
            field=_field_label(field),
        )
    if not text.endswith(".vtk"):
        text += ".vtk"
    return Path(text)


def _as_named_array(data, field, scalar_name):
    """Select a named scalar array from a state object or raw array."""
    if isinstance(data, Mapping) and field in data:
        values = data[field]
        name = scalar_name or field
    elif hasattr(data, field):
        values = getattr(data, field)
        name = scalar_name or field
    else:
        values = data
        name = scalar_name or "scalars"
    if values is None:
        raise ValueError(f"no data available for VTK field '{field}'")
    return name, np.asarray(values, dtype=np.float64)


def _scalar_name_for_field(scalar_name, field, index):
    if scalar_name is None:
        return None
    if isinstance(scalar_name, Mapping):
        return scalar_name.get(field)
    if isinstance(scalar_name, Sequence) and not isinstance(scalar_name, (str, bytes)):
        return scalar_name[index]
    return scalar_name


def _as_named_arrays(data, field, scalar_name=None, fields=None):
    """Resolve one or more named scalar arrays from object attributes or mappings."""
    if fields is not None:
        if not isinstance(fields, Mapping):
            raise TypeError("fields must be a mapping of scalar names to arrays or attribute names")
        arrays = []
        for name, source in fields.items():
            if isinstance(source, str) and hasattr(data, source):
                values = getattr(data, source)
            elif isinstance(data, Mapping) and isinstance(source, str) and source in data:
                values = data[source]
            else:
                values = source
            if values is None:
                raise ValueError(f"no data available for VTK field '{name}'")
            arrays.append((str(name), np.asarray(values, dtype=np.float64)))
        return arrays

    if isinstance(field, Mapping):
        arrays = []
        for index, (name, source) in enumerate(field.items()):
            alias = _scalar_name_for_field(scalar_name, name, index) or name
            if isinstance(source, str):
                _, values = _as_named_array(data, source, alias)
            else:
                values = np.asarray(source, dtype=np.float64)
            arrays.append((str(alias), values))
        return arrays

    if isinstance(field, Sequence) and not isinstance(field, (str, bytes)):
        return [
            _as_named_array(data, item, _scalar_name_for_field(scalar_name, item, index))
            for index, item in enumerate(field)
        ]

    return [_as_named_array(data, field, scalar_name)]


def _resolve_data_shape(values, topology):
    """Classify arrays as point-level or prism-level VTK scalar data."""
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


def _write_ascii_wedge(file_name, arrays, topology, scalar_name="scalars"):
    """Write scalar data on the topology's extruded wedge-prism mesh."""
    topology._require_levels()
    topology._require_thickness()

    if isinstance(arrays, tuple) and len(arrays) == 2 and isinstance(arrays[0], str):
        named_arrays = [arrays]
    elif isinstance(arrays, Sequence) and arrays and isinstance(arrays[0], tuple):
        named_arrays = list(arrays)
    else:
        named_arrays = [(scalar_name, arrays)]

    data_arrays = {"POINT_DATA": [], "CELL_DATA": []}
    data_counts = {}
    for name, values in named_arrays:
        values = np.asarray(values, dtype=np.float64)
        data_kind, data_shape = _resolve_data_shape(values, topology)
        values = values.reshape(data_shape, order="F")
        data_arrays[data_kind].append((str(name), values))
        data_counts[data_kind] = values.size

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

        for data_kind in ("POINT_DATA", "CELL_DATA"):
            arrays_for_kind = data_arrays[data_kind]
            if not arrays_for_kind:
                continue
            handle.write(f"{data_kind} {data_counts[data_kind]}\n")
            for name, values in arrays_for_kind:
                handle.write(f"SCALARS {name} float 1\n")
                handle.write("LOOKUP_TABLE default\n")
                for value in values.reshape(-1, order="F"):
                    handle.write(f"{float(value):.17g}\n")

    return path


def _legacy_vtk_wedge(file_name, values, points, triangle_point_indices, mesh_z, z_mesh, scalar_name):
    """Compatibility wrapper for the historical points/triangles signature."""
    topology = MeshTopology(
        points=np.asarray(points, dtype=np.float64),
        trianglePointIndices=np.asarray(triangle_point_indices, dtype=np.uint32),
        levels=int(mesh_z),
        thickness=float(z_mesh),
    )
    return _write_ascii_wedge(file_name, values, topology, scalar_name=scalar_name)


def vtkWedge(file_name, data=None, geometry=None, field="phiAse", scalar_name=None, every=1, fields=None):
    """Write scalar data on a wedge mesh to a legacy ASCII VTK file.

    The preferred simulation callback form is ``vtkWedge(file_name, state,
    fields=...)`` inside a user ``onStep`` function. ``state`` must be a
    ``TimeStepState`` produced by ``Simulation``; it carries both the dynamic
    arrays and the static topology needed to write VTK geometry. ``file_name``
    may contain ``{step}``, ``{time}``, and ``{field}`` placeholders.

    ``field`` selects attributes or mapping keys from ``data``. A string writes
    one scalar array, a sequence writes several arrays with their original
    names, and a mapping writes aliased names, for example
    ``field={"phi": "phiAse"}``. ``fields`` is for explicit arrays and maps
    VTK scalar names to arrays, for example
    ``fields={"phiASE": state.phiAse, "cladAbs": state.phiAse * 5.5}``.
    ``fields`` takes precedence over ``field``.

    Point-shaped arrays ``(numberOfPoints, levels)`` are written as
    ``POINT_DATA``. Prism-shaped arrays ``(numberOfTriangles, levels - 1)`` are
    written as ``CELL_DATA``. Mixed point and cell arrays can be written to the
    same file. For standalone array exports outside a ``TimeStepState``, pass
    ``geometry`` as a ``GainMedium`` or ``MeshTopology``.

    The legacy callback-factory form ``simulation.onStep(vtkWedge(path,
    geometry, ...))`` is still accepted. In that form, ``every`` controls the
    output period and must be a positive integer.
    """

    if isinstance(data, (GainMedium, MeshTopology)) and geometry is None:
        topology = _topology_from_geometry(data)
        def write_state(state):
            if int(every) <= 0:
                raise ValueError("every must be a positive integer")
            if getattr(state, "step", 0) % int(every) != 0:
                return None
            arrays = _as_named_arrays(state, field, scalar_name, fields=fields)
            return _write_ascii_wedge(_vtk_filename(file_name, state, field), arrays, topology)

        return write_state

    if geometry is None:
        topology = getattr(data, "topology", None)
        if topology is None:
            raise TypeError("vtkWedge requires a TimeStepState with topology or an explicit GainMedium/MeshTopology")
    else:
        topology = _topology_from_geometry(geometry)

    arrays = _as_named_arrays(data, field, scalar_name, fields=fields)
    return _write_ascii_wedge(_vtk_filename(file_name, data, field), arrays, topology)
