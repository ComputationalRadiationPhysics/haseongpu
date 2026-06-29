# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Geometry, topology, and gain-medium state containers for Python users."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
import warnings


import numpy as np
from ..openpmd import BackendFlatArray, FieldSpec, GroupFieldSpec, PrimitiveSchema, PrimitiveSchemaDefinition, schemaFields
from .msh import Gmsh
from .stl import from_stl
from .vtk import gainMediumFromVtk, topologyFromVtk, writeGainMediumVtk
try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        def decorator(fn):
            return fn
        return decorator


@dataclass(frozen=True)
class PhysicalPropertySpec:
    """Metadata for a physical array or scalar stored on a ``GainMedium``.

    ``shape`` is evaluated against the current ``MeshTopology`` so properties
    can describe point, prism, triangle, surface, or scalar data.
    """

    name: str
    dtype: object
    shape: callable
    description: str
    default: callable | None = None

    def expectedShape(self, topology):
        return self.shape(topology)

    def defaultValue(self, topology):
        if self.default is None:
            return None
        return self.default(topology)



_FIELD_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_ENTITY_ALIASES = {
    "triangle": "cell",
    "triangles": "cell",
    "prism": "cell_layer",
    "prisms": "cell_layer",
    "point_level": "point_level",
    "cell_layer": "cell_layer",
}
_ALLOWED_ENTITY_AXES = {
    "point",
    "cell",
    "layer",
    "level",
    "interface",
    "wavelength",
    "local_vertex",
    "local_side",
    "domain",
}
_CUSTOM_FIELD_RESERVED_NAMES = {
    "coordinates",
    "pointIndices",
    "betaCells",
    "dntdAse",
    "betaVolume",
    "claddingCellTypes",
    "reflectivities",
}


def _camel_to_snake(name):
    out = []
    for index, char in enumerate(name):
        if char.isupper() and index > 0 and (not name[index - 1].isupper()):
            out.append("_")
        out.append(char.lower())
    return "".join(out)


def _normalize_custom_entity(entity):
    if isinstance(entity, str):
        entity = _ENTITY_ALIASES.get(entity, entity)
        axes = tuple(part for part in entity.split("_") if part)
    else:
        axes = tuple(str(part) for part in entity)
        if len(axes) == 1:
            alias = _ENTITY_ALIASES.get(axes[0])
            if alias is not None:
                axes = tuple(alias.split("_"))
    if not axes:
        raise ValueError("field entity must not be empty")
    unknown = [axis for axis in axes if axis not in _ALLOWED_ENTITY_AXES]
    if unknown:
        allowed = ", ".join(sorted(_ALLOWED_ENTITY_AXES | set(_ENTITY_ALIASES)))
        raise ValueError(f"unknown field entity axis '{unknown[0]}'. Allowed axes: {allowed}")
    if axes == ("cell", "layer") or axes == ("point", "level"):
        return axes
    if len(axes) == 1 and axes[0] in {"point", "cell", "interface", "wavelength", "domain"}:
        return axes
    if axes in {("cell", "local_vertex"), ("cell", "local_side")}:
        return axes
    raise ValueError(
        "unsupported field entity " + repr(axes) + "; supported entities include "
        "'point', 'cell'/'triangle', 'prism' or ('cell', 'layer'), "
        "('point', 'level'), 'interface', 'wavelength', and 'domain'"
    )


def _shape_for_custom_entity(topology, axes, values=None):
    if axes == ("point",):
        return (topology.numberOfPoints,)
    if axes == ("cell",):
        return (topology.numberOfTriangles,)
    if axes == ("cell", "layer"):
        topology._require_levels()
        return (topology.numberOfTriangles, topology.levels - 1)
    if axes == ("point", "level"):
        topology._require_levels()
        return (topology.numberOfPoints, topology.levels)
    if axes == ("cell", "local_vertex") or axes == ("cell", "local_side"):
        return (topology.numberOfTriangles, 3)
    if axes == ("interface",):
        return (4,)
    if axes == ("wavelength",):
        if values is None:
            raise ValueError("wavelength field values are required to determine shape")
        return (np.asarray(values).size,)
    if axes == ("domain",):
        number_of_domains = topology.metadata.get("numberOfDomains", topology.metadata.get("number_of_domains"))
        if number_of_domains is not None:
            return (int(number_of_domains),)
        if values is None:
            raise ValueError("domain field values are required when topology has no numberOfDomains metadata")
        arr = np.asarray(values)
        if arr.ndim != 1:
            raise ValueError(f"domain field values must be 1-D, got shape {arr.shape}")
        return (arr.size,)
    raise ValueError(f"unsupported field entity {axes}")


@dataclass(frozen=True)
class CustomField:
    spec: FieldSpec
    value: np.ndarray
    primitiveShape: tuple[int, ...]

    @property
    def name(self):
        return self.spec.name

    @property
    def entity(self):
        return self.spec.entity

    @property
    def axes(self):
        return self.spec.axes

    @property
    def dtype(self):
        return self.spec.dtypeObject

    @property
    def dynamic(self):
        return self.spec.dynamic

    @property
    def unit(self):
        return self.spec.unit

    def primitiveView(self):
        return self.value.reshape(self.primitiveShape, order="F")

    def meta(self):
        return {
            "name": self.name,
            "entity": self.entity,
            "axes": self.axes,
            "dtype": str(self.dtype),
            "unit": self.unit,
            "dynamic": self.dynamic,
            "backendRequired": self.spec.backendRequired,
            "expectedShape": self.primitiveShape,
            "isSet": True,
        }


@dataclass(frozen=True)
class OpenPmdScalarField:
    name: str
    values: object
    context: object | None = None
    prefix: str = "core_"
    spec: object | None = None


@dataclass(frozen=True)
class OpenPmdComponentField:
    name: str
    components: dict[str, object]
    axisLabels: list[str]
    context: object | None = None
    prefix: str = "core_"
    recordName: str | None = None


class PrimitiveElement:
    def __init__(self, view, index):
        object.__setattr__(self, "_view", view)
        object.__setattr__(self, "index", index if isinstance(index, tuple) else (index,))

    def _value(self, name):
        value = self._view[name]
        if hasattr(value, "shape") and value.shape[:len(self.index)] == self._view.shape[:len(self.index)]:
            result = value[self.index]
            if hasattr(result, "item") and getattr(result, "shape", ()) == ():
                return result.item()
            return result
        return value

    def _set_value(self, name, value):
        field = self._view[name]
        if hasattr(field, "shape") and field.shape[:len(self.index)] == self._view.shape[:len(self.index)]:
            try:
                field[self.index] = value
            except (TypeError, ValueError) as exc:
                expected = field[self.index].shape if hasattr(field[self.index], "shape") else ()
                warnings.warn(
                    f"could not assign value to field '{name}' on {self._view.name} element "
                    f"{self.index}; expected value compatible with shape {expected}: {exc}",
                    RuntimeWarning,
                    stacklevel=2,
                )
                raise
            return
        raise TypeError(f"field '{name}' is not indexed by {self._view.name} elements")

    def __getitem__(self, name):
        return self._value(name)

    def __setitem__(self, name, value):
        self._set_value(name, value)

    def __getattr__(self, name):
        if name in self._view.keys():
            return self._value(name)
        raise AttributeError(name)

    def __setattr__(self, name, value):
        if name.startswith("_") or name == "index":
            object.__setattr__(self, name, value)
            return
        if name in self._view.keys():
            self._set_value(name, value)
            return
        raise AttributeError(name)

    def getFields(self):
        return [PrimitiveElementField(self, name) for name in self._view.keys()]

    def _applyGroup(self, group):
        self._view._applyGroup(self.index, group)


_MISSING_FIELD_VALUE = object()


class PrimitiveElementField:
    def __init__(self, element, name):
        self._element = element
        self.name = name

    def value(self, newValue=_MISSING_FIELD_VALUE):
        if newValue is _MISSING_FIELD_VALUE:
            return self._element[self.name]
        self._element[self.name] = newValue
        return self

    def meta(self):
        return self._element._view.fieldMeta(self.name)

    def __repr__(self):
        meta = self.meta()
        return (
            f"PrimitiveElementField(name={self.name!r}, dtype={meta['dtype']!r}, "
            f"unit={meta.get('unit', '1')!r}, shape={meta['shape']!r})"
        )


class PrimitiveView:
    def __init__(self, name, shape, fields, metadata=None):
        self.name = name
        self.shape = tuple(shape)
        self._fields = dict(fields)
        self._metadata = dict(metadata or {})
        self._groupAssignments = {}

    def keys(self):
        return self._fields.keys()

    def items(self):
        return self._fields.items()

    def __iter__(self):
        for index in np.ndindex(self.shape):
            yield PrimitiveElement(self, index)

    def get(self, name):
        return self._fields[name]


    def _applyGroup(self, index, group):
        for spec, value in group.fieldItems():
            self._assignGroupField(index, group, spec, value)

    def _assignGroupField(self, index, group, spec: GroupFieldSpec, value):
        name = spec.target
        assignments = self._groupAssignments.setdefault(name, {})
        existing = assignments.get(tuple(index))
        if existing is not None and existing is not group._groupToken:
            raise ValueError(
                f"{self.name} element {tuple(index)} already belongs to a group assigning field '{name}'"
            )

        arr_value = np.asarray(value, dtype=spec.dtype)
        if name not in self._fields:
            field_shape = self.shape + tuple(arr_value.shape)
            self._fields[name] = np.zeros(field_shape, dtype=np.dtype(spec.dtype))
            self._metadata[name] = {
                "name": name,
                "entity": self.name,
                "axes": (self.name,),
                "dtype": str(np.dtype(spec.dtype)),
                "unit": spec.unit,
                "shape": field_shape,
                "isSet": True,
            }

        field = self._fields[name]
        expected = field[tuple(index)].shape if hasattr(field[tuple(index)], "shape") else ()
        if tuple(arr_value.shape) != tuple(expected):
            raise ValueError(
                f"group field '{name}' value has shape {arr_value.shape}, expected {expected} for {self.name} element"
            )
        field[tuple(index)] = arr_value
        assignments[tuple(index)] = group._groupToken

    def asDict(self):
        return dict(self._fields)

    def __getitem__(self, name):
        return self.get(name)

    def getFields(self):
        return [self.fieldMeta(name) for name in self._fields]

    def fieldMeta(self, name):
        if name not in self._fields:
            raise KeyError(name)
        if name in self._metadata:
            return dict(self._metadata[name])
        value = np.asarray(self._fields[name])
        return {
            "name": name,
            "entity": self.name,
            "axes": (),
            "dtype": str(value.dtype),
            "unit": "1",
            "shape": value.shape,
            "isSet": True,
        }


class GainMediumProperty:
    """Handle returned by ``GainMedium.get(...)``.

    It carries the property's physical description, dtype, expected shape, and
    a ``value`` setter that validates arrays before storing them on the medium.
    """

    def __init__(self, medium, spec):
        self._medium = medium
        self.spec = spec

    @property
    def name(self):
        return self.spec.name

    @property
    def description(self):
        return self.spec.description

    @property
    def dtype(self):
        return np.dtype(self.spec.dtype)

    @property
    def expectedShape(self):
        return self.spec.expectedShape(self._medium.topology)

    @property
    def value(self):
        return self._medium.physical.get(self.name, self.spec.defaultValue(self._medium.topology))

    @value.setter
    def value(self, values):
        self._medium.set(self.name, values)

    def meta(self):
        """Return serializable metadata for docs, inspection, or UI tooling."""
        return {
            "name": self.name,
            "description": self.description,
            "dtype": str(self.dtype),
            "expectedShape": self.expectedShape,
            "isSet": self.name in self._medium.physical,
        }


def _shape_prisms(topology):
    topology._require_levels()
    return (topology.numberOfTriangles, topology.levels - 1)


def _shape_cells(topology):
    topology._require_levels()
    return (topology.numberOfPoints, topology.levels)


def _shape_triangles(topology):
    return (topology.numberOfTriangles,)


def _shape_reflectivities(topology):
    return (topology.numberOfTriangles, 2)


def _shape_refractive(_):
    return (4,)


def _shape_scalar(_):
    return ()


PHYSICAL_PROPERTY_SPECS = {
    "betaVolume": PhysicalPropertySpec(
        name="betaVolume",
        dtype=np.float64,
        shape=_shape_prisms,
        description=(
            "Prism-centered excited-state fraction beta_j used by ASE. "
            "Matrix layout is (numberOfTriangles, numberOfLevels - 1)."
        ),
        default=lambda topology: np.zeros(topology.numberOfPrisms, dtype=np.float64),
    ),
    "betaCells": PhysicalPropertySpec(
        name="betaCells",
        dtype=np.float64,
        shape=_shape_cells,
        description=(
            "Point-centered excited-state fraction beta_i. Matrix layout is "
            "(numberOfPoints, numberOfLevels)."
        ),
        default=lambda topology: np.zeros(topology.numberOfPoints * topology.levels, dtype=np.float64),
    ),
    "dntdAse": PhysicalPropertySpec(
        name="dntdAse",
        dtype=np.float64,
        shape=_shape_cells,
        description="ASE contribution to d beta / dt at sample points.",
        default=lambda topology: np.zeros(topology.numberOfPoints * topology.levels, dtype=np.float64),
    ),
    "claddingCellTypes": PhysicalPropertySpec(
        name="claddingCellTypes",
        dtype=np.uint32,
        shape=_shape_triangles,
        description="Cladding type index for each triangle.",
        default=lambda topology: np.zeros(topology.numberOfTriangles, dtype=np.uint32),
    ),
    "refractiveIndices": PhysicalPropertySpec(
        name="refractiveIndices",
        dtype=np.float32,
        shape=_shape_refractive,
        description="[bottomInside, bottomOutside, topInside, topOutside].",
        default=lambda _: np.ones(4, dtype=np.float32),
    ),
    "reflectivities": PhysicalPropertySpec(
        name="reflectivities",
        dtype=np.float32,
        shape=_shape_reflectivities,
        description=(
            "Surface reflectivity per triangle. Matrix layout is "
            "(numberOfTriangles, 2): column 0 bottom, column 1 top."
        ),
        default=lambda topology: np.ones(2 * topology.numberOfTriangles, dtype=np.float32),
    ),
    "nTot": PhysicalPropertySpec(
        name="nTot",
        dtype=np.float64,
        shape=_shape_scalar,
        description=(
            "Total active-ion concentration of the gain medium in cm^-3."
        ),
        default=lambda _: 0.0,
    ),
    "crystalTFluo": PhysicalPropertySpec(
        name="crystalTFluo",
        dtype=np.float64,
        shape=_shape_scalar,
        description="Fluorescence lifetime tau of the gain medium in seconds.",
        default=lambda _: 0.0,
    ),
    "claddingNumber": PhysicalPropertySpec(
        name="claddingNumber",
        dtype=np.uint32,
        shape=_shape_scalar,
        description="Corresponding cladding cells are selected for cladding absorption handling.",
        default=lambda _: 1,
    ),
    "claddingAbsorption": PhysicalPropertySpec(
        name="claddingAbsorption",
        dtype=np.float64,
        shape=_shape_scalar,
        description="Absorption coefficient of the cladding.",
        default=lambda _: 0.0,
    ),
}

PROPERTY_ALIASES = {
    "refractiveIndicies": "refractiveIndices",
    "refractive_indices": "refractiveIndices",
    "cladding_cell_types": "claddingCellTypes",
    "cladding_absorption": "claddingAbsorption",
    "cladding_number": "claddingNumber",
    "crystal_t_fluo": "crystalTFluo",
    "beta_volume": "betaVolume",
    "beta_cells": "betaCells",
    "dndtAse": "dntdAse",
    "dndt_ASE": "dntdAse",
    "phi_ASE": "dntdAse",
}


def _as_points(points):
    arr = np.asarray(points, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"points must have shape (N, 2), got {arr.shape}")
    if arr.shape[0] < 3:
        raise ValueError("at least three points are required")
    _, unique_idx = np.unique(arr, axis=0, return_index=True)
    return arr[np.sort(unique_idx)]


def _cross2d(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _circumcircle_contains(points, tri, point):
    a, b, c = points[list(tri)]
    ax, ay = a - point
    bx, by = b - point
    cx, cy = c - point
    det = (
        (ax * ax + ay * ay) * (bx * cy - cx * by)
        - (bx * bx + by * by) * (ax * cy - cx * ay)
        + (cx * cx + cy * cy) * (ax * by - bx * ay)
    )
    orient = _cross2d(b - a, c - a)
    return det > 1e-12 if orient > 0 else det < -1e-12


def _delaunay(points):
    points = _as_points(points)
    min_xy = points.min(axis=0)
    max_xy = points.max(axis=0)
    delta = max(max_xy - min_xy)
    if delta <= 0.0:
        raise ValueError("points must span a two-dimensional area")

    center = (min_xy + max_xy) / 2.0
    super_points = np.array(
        [
            [center[0] - 20.0 * delta, center[1] - delta],
            [center[0], center[1] + 20.0 * delta],
            [center[0] + 20.0 * delta, center[1] - delta],
        ],
        dtype=np.float64,
    )
    all_points = np.vstack([points, super_points])
    super_ids = set(range(len(points), len(points) + 3))
    triangles = [tuple(super_ids)]

    for point_id in range(len(points)):
        point = all_points[point_id]
        bad = [tri for tri in triangles if _circumcircle_contains(all_points, tri, point)]
        boundary = {}
        for tri in bad:
            for edge in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                key = tuple(sorted(edge))
                boundary[key] = boundary.get(key, 0) + 1
        triangles = [tri for tri in triangles if tri not in bad]
        for edge, count in boundary.items():
            if count == 1:
                tri = (edge[0], edge[1], point_id)
                a, b, c = all_points[list(tri)]
                if _cross2d(b - a, c - a) < 0:
                    tri = (edge[1], edge[0], point_id)
                triangles.append(tri)

    result = [
        tri for tri in triangles
        if not any(vertex in super_ids for vertex in tri)
    ]
    if not result:
        raise ValueError("could not construct triangles; points may be collinear")
    return points, np.asarray(result, dtype=np.uint32)


def _grid_points_and_triangles(grid):
    ox, oy = grid.origin
    x_values = ox + Grid._axisValues(grid.xExtent, grid.tileSizeX)
    y_values = oy + Grid._axisValues(grid.yExtent, grid.tileSizeY)
    nx = x_values.size
    ny = y_values.size
    if nx < 2 or ny < 2:
        raise ValueError("grid must contain at least two points along x and y")

    points = np.empty((nx * ny, 2), dtype=np.float64)
    points[:, 0] = np.tile(x_values, ny)
    points[:, 1] = np.repeat(y_values, nx)

    cells_x = nx - 1
    cells_y = ny - 1
    cells = np.arange(cells_x * cells_y, dtype=np.uint32)
    i = cells % cells_x
    j = cells // cells_x
    p00 = j * nx + i
    p10 = p00 + 1
    p01 = p00 + nx
    p11 = p01 + 1

    first = np.empty((cells_y, cells_x, 3), dtype=np.uint32)
    first[:, :, 0] = p00.reshape(cells_y, cells_x)
    first[:, :, 1] = p10.reshape(cells_y, cells_x)
    first[:, :, 2] = p01.reshape(cells_y, cells_x)

    second = np.empty((cells_y, cells_x, 3), dtype=np.uint32)
    second[:, :, 0] = p01.reshape(cells_y, cells_x)
    second[:, :, 1] = p10.reshape(cells_y, cells_x)
    second[:, :, 2] = p11.reshape(cells_y, cells_x)

    triangles_by_row = np.empty((cells_y, 2 * cells_x, 3), dtype=np.uint32)
    triangles_by_row[:, :cells_x, :] = first
    triangles_by_row[:, cells_x:, :] = second
    triangles = triangles_by_row.reshape(-1, 3)
    return points, triangles


@njit(cache=True, fastmath=True)
def _triangle_geometry_kernel(points, triangles):
    number_of_triangles = triangles.shape[0]
    centers_x = np.empty(number_of_triangles, dtype=np.float64)
    centers_y = np.empty(number_of_triangles, dtype=np.float64)
    surfaces = np.empty(number_of_triangles, dtype=np.float32)
    normals_x = np.empty((number_of_triangles, 3), dtype=np.float64)
    normals_y = np.empty((number_of_triangles, 3), dtype=np.float64)
    normal_points = np.empty((number_of_triangles, 3), dtype=np.uint32)
    edge_keys = np.empty((number_of_triangles * 3, 2), dtype=np.uint32)
    edge_triangles = np.empty(number_of_triangles * 3, dtype=np.int32)
    edge_positions = np.empty(number_of_triangles * 3, dtype=np.int32)

    for tri_i in range(number_of_triangles):
        p0 = triangles[tri_i, 0]
        p1 = triangles[tri_i, 1]
        p2 = triangles[tri_i, 2]

        x0 = points[p0, 0]
        y0 = points[p0, 1]
        x1 = points[p1, 0]
        y1 = points[p1, 1]
        x2 = points[p2, 0]
        y2 = points[p2, 1]

        centers_x[tri_i] = (x0 + x1 + x2) / 3.0
        centers_y[tri_i] = (y0 + y1 + y2) / 3.0
        cross = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
        if cross < 0.0:
            cross = -cross
        surfaces[tri_i] = 0.5 * cross

        for edge_i in range(3):
            if edge_i == 0:
                a = p0
                b = p1
            elif edge_i == 1:
                a = p1
                b = p2
            else:
                a = p2
                b = p0

            edge_x = points[b, 0] - points[a, 0]
            edge_y = points[b, 1] - points[a, 1]
            inv_norm = 1.0 / np.sqrt(edge_x * edge_x + edge_y * edge_y)
            normals_x[tri_i, edge_i] = edge_y * inv_norm
            normals_y[tri_i, edge_i] = -edge_x * inv_norm
            normal_points[tri_i, edge_i] = a

            edge_index = tri_i * 3 + edge_i
            if a < b:
                edge_keys[edge_index, 0] = a
                edge_keys[edge_index, 1] = b
            else:
                edge_keys[edge_index, 0] = b
                edge_keys[edge_index, 1] = a
            edge_triangles[edge_index] = tri_i
            edge_positions[edge_index] = edge_i

    return centers_x, centers_y, surfaces, normals_x, normals_y, normal_points, edge_keys, edge_triangles, edge_positions


def _derived_topology(points, triangles):
    (
        centers_x,
        centers_y,
        surfaces,
        normals_x,
        normals_y,
        normal_points,
        edge_keys,
        edge_triangles,
        edge_positions,
    ) = _triangle_geometry_kernel(points, triangles)

    number_of_triangles = triangles.shape[0]
    neighbors = np.full((number_of_triangles, 3), -1, dtype=np.int32)
    forbidden = np.full((number_of_triangles, 3), -1, dtype=np.int32)
    order = np.lexsort((edge_keys[:, 1], edge_keys[:, 0]))
    sorted_keys = edge_keys[order]
    matches = np.flatnonzero(np.all(sorted_keys[1:] == sorted_keys[:-1], axis=1))

    if matches.size:
        left_order = order[matches]
        right_order = order[matches + 1]
        left_triangles = edge_triangles[left_order]
        left_positions = edge_positions[left_order]
        right_triangles = edge_triangles[right_order]
        right_positions = edge_positions[right_order]
        neighbors[left_triangles, left_positions] = right_triangles
        forbidden[left_triangles, left_positions] = right_positions
        neighbors[right_triangles, right_positions] = left_triangles
        forbidden[right_triangles, right_positions] = left_positions

    return {
        "triangleCenterX": centers_x,
        "triangleCenterY": centers_y,
        "triangleSurfaces": surfaces,
        "triangleNormalsX": normals_x.reshape(-1, order="F"),
        "triangleNormalsY": normals_y.reshape(-1, order="F"),
        "triangleNormalPoint": normal_points.reshape(-1, order="F"),
        "triangleNeighbors": neighbors.reshape(-1, order="F"),
        "forbiddenEdge": forbidden.reshape(-1, order="F"),
    }


def _flat(arr, width, dtype, name):
    arr = np.asarray(arr, dtype=dtype)
    if arr.ndim == 1:
        if width and arr.size % width != 0:
            raise ValueError(f"{name} length must be divisible by {width}")
        return arr
    if arr.ndim == 2 and (width is None or arr.shape[1] == width):
        return arr.reshape(-1, order="F")
    raise ValueError(f"{name} must be flat or have width {width}")


class Grid:
    """Regular rectangular helper that generates a 2D triangular base mesh.

    ``xExtent`` and ``yExtent`` describe the transverse footprint. ``zExtent``
    and ``tileSizeZ`` define the sampled crystal levels used for beta arrays.
    """

    def __init__(
        self,
        xExtent=None,
        yExtent=None,
        zExtent=None,
        tileSizeX=1.0,
        tileSizeY=None,
        tileSizeZ=None,
        origin=(0.0, 0.0),
        **aliases,
    ):
        """Create a grid from extents, tile sizes, and optional ``(x, y)`` origin."""
        if xExtent is None and "x" in aliases:
            xExtent = aliases.pop("x")
        if yExtent is None and "y" in aliases:
            yExtent = aliases.pop("y")
        if zExtent is None and "z" in aliases:
            zExtent = aliases.pop("z")
        if "tilesSizeX" in aliases:
            tileSizeX = aliases.pop("tilesSizeX")
        if "tilesSizeY" in aliases:
            tileSizeY = aliases.pop("tilesSizeY")
        if "tilesSizeZ" in aliases:
            tileSizeZ = aliases.pop("tilesSizeZ")
        if aliases:
            unknown = ", ".join(aliases)
            raise TypeError(f"unknown Grid arguments: {unknown}")
        if xExtent is None or yExtent is None or zExtent is None:
            raise TypeError("Grid requires xExtent, yExtent and zExtent")

        self.xExtent = float(xExtent)
        self.yExtent = float(yExtent)
        self.zExtent = float(zExtent)
        self.tileSizeX = float(tileSizeX)
        self.tileSizeY = float(self.tileSizeX if tileSizeY is None else tileSizeY)
        self.tileSizeZ = float(self.tileSizeX if tileSizeZ is None else tileSizeZ)
        self.origin = origin
        if self.xExtent <= 0.0 or self.yExtent <= 0.0 or self.zExtent <= 0.0:
            raise ValueError("grid dimensions x, y and z must be positive")
        if self.tileSizeX <= 0.0 or self.tileSizeY <= 0.0 or self.tileSizeZ <= 0.0:
            raise ValueError("tile sizes must be positive")

    @property
    def numberOfLevels(self):
        """Number of z sample planes generated from ``zExtent`` and ``tileSizeZ``."""
        return len(self._axisValues(self.zExtent, self.tileSizeZ))

    @property
    def thickness(self):
        """Distance between adjacent z levels."""
        return self.tileSizeZ

    def constructPoints(self):
        """Return transverse topology points with shape ``(numberOfPoints, 2)``."""
        ox, oy = self.origin
        x_values = ox + self._axisValues(self.xExtent, self.tileSizeX)
        y_values = oy + self._axisValues(self.yExtent, self.tileSizeY)

        points = np.empty((x_values.size * y_values.size, 2), dtype=np.float64)
        points[:, 0] = np.tile(x_values, y_values.size)
        points[:, 1] = np.repeat(y_values, x_values.size)
        return points

    @staticmethod
    def _axisValues(extent, tileSize):
        # Add half a step so exactly divisible extents include the endpoint.
        values = np.arange(0.0, extent + tileSize * 0.5, tileSize)
        return values[values <= extent + 1e-12]


@dataclass
class MeshTopology:
    """Triangular base mesh extruded through z levels into wedge prisms.

    ``points`` are transverse ``(x, y)`` coordinates. ``trianglePointIndices``
    connect those points into the 2D base mesh. ``levels`` and ``thickness``
    describe the z sampling used by ``betaCells`` and ``betaVolume``.
    """

    points: np.ndarray
    """Transverse coordinates with shape ``(numberOfPoints, 2)``."""
    trianglePointIndices: np.ndarray
    """Triangle vertex indices with shape ``(numberOfTriangles, 3)``."""
    levels: int | None = None
    """Number of z sample planes; at least two when simulation data is used."""
    thickness: float | None = None
    """Distance between adjacent z sample planes."""
    metadata: dict = field(default_factory=dict)
    """Importer or construction metadata kept for tooling and round-trips."""

    @classmethod
    def fromPoints(cls, points, numberOfLevels=None):
        """Delaunay-triangulate transverse points into a ``MeshTopology``."""
        parsed_points, triangles = _delaunay(points)
        return cls(parsed_points, triangles, levels=numberOfLevels)

    @classmethod
    def fromGrid(cls, grid):
        """Construct topology from a regular ``Grid`` helper."""
        if not isinstance(grid, Grid):
            raise TypeError("fromGrid expects a Grid instance")
        points, triangles = _grid_points_and_triangles(grid)
        topology = cls(points, triangles, levels=grid.numberOfLevels)
        topology.thickness = grid.thickness
        topology.metadata.update(
            {
                "source": "grid",
                "grid": {
                    "x": grid.xExtent,
                    "y": grid.yExtent,
                    "z": grid.zExtent,
                    "tileSizeX": grid.tileSizeX,
                    "tileSizeY": grid.tileSizeY,
                    "tileSizeZ": grid.tileSizeZ,
                },
            }
        )
        return topology

    @classmethod
    def fromGmsh(cls, gmsh, *, numberOfLevels=None, thickness=None):
        """Import a planar gmsh triangle mesh and attach z sampling."""
        if isinstance(gmsh, (str, Path)):
            gmsh = Gmsh.fromFile(gmsh)
        if not isinstance(gmsh, Gmsh):
            raise TypeError("fromGmsh expects a Gmsh instance or a gmsh .msh filename")
        return gmsh.topology(numberOfLevels=numberOfLevels, thickness=thickness)

    @classmethod
    def fromVtk(cls, filename, *, numberOfLevels=None, thickness=None):
        """Load topology from a legacy VTK wedge file."""
        topology = topologyFromVtk(filename, cls)
        if numberOfLevels is not None:
            topology.numberOfLevels(numberOfLevels)
        if thickness is not None:
            topology.withThickness(thickness)
        return topology

    @classmethod
    def fromFile(cls, filename, format=None, numberOfLevels=None, thickness=None):
        """Load supported mesh formats: legacy VTK, planar STL, or gmsh."""
        path = Path(filename)
        mesh_format = (format or path.suffix.lstrip(".")).lower()
        if mesh_format in {"vtk", "legacy-vtk"}:
            return cls.fromVtk(path, numberOfLevels=numberOfLevels, thickness=thickness)
        if mesh_format in {"stl", "ascii-stl", "binary-stl", "dea/stl", "dae/stl"}:
            points, triangles = from_stl(path)
            return cls(
                points,
                triangles,
                levels=numberOfLevels,
                metadata={"source": str(path), "format": mesh_format},
            )
        if mesh_format in {"msh", "gmsh"}:
            gmsh = Gmsh.fromFile(path)
            return cls.fromGmsh(gmsh, numberOfLevels=numberOfLevels, thickness=thickness)
        raise NotImplementedError(
            f"mesh format '{mesh_format}' is not supported yet; supported formats: vtk, stl, gmsh"
        )

    def numberOfLevels(self, levels):
        """Set the number of z sample planes and return ``self``."""
        levels = int(levels)
        if levels < 2:
            raise ValueError("numberOfLevels must be at least 2")
        self.levels = levels
        return self

    def withThickness(self, thickness):
        """Set the spacing between adjacent z levels and return ``self``."""
        thickness = float(thickness)
        if thickness <= 0.0:
            raise ValueError("thickness must be positive")
        self.thickness = thickness
        return self

    def levelCoordinates(self):
        """Return z coordinates for every level as ``0, thickness, ...``."""
        self._require_levels()
        self._require_thickness()
        return np.arange(int(self.levels), dtype=np.float64) * float(self.thickness)

    def pointIndexAt(self, x, y, *, tol=1e-12):
        """Return the topology point index for a transverse coordinate."""
        query = np.asarray([x, y], dtype=np.float64)
        matches = np.flatnonzero(np.all(np.isclose(self.points, query, atol=tol, rtol=0.0), axis=1))
        if matches.size:
            return int(matches[0])
        distances = np.linalg.norm(self.points - query, axis=1)
        return int(np.argmin(distances))

    def levelIndexAt(self, z, *, tol=1e-12):
        """Return the z-level index for a physical z coordinate."""
        levels = self.levelCoordinates()
        matches = np.flatnonzero(np.isclose(levels, float(z), atol=tol, rtol=0.0))
        if matches.size == 0:
            raise ValueError(f"no topology level found at z={z}")
        return int(matches[0])

    def betaCellIndexAt(self, x, y, z, *, tol=1e-12, flat=False):
        """Return the ``betaCells`` index for a physical ``(x, y, z)`` point."""
        point_index = self.pointIndexAt(x, y, tol=tol)
        level_index = self.levelIndexAt(z, tol=tol)
        if flat:
            return point_index + level_index * self.numberOfPoints
        return point_index, level_index

    @property
    def numberOfTriangles(self):
        """Number of base triangles, and therefore prisms per z interval."""
        return int(self.trianglePointIndices.shape[0])

    @property
    def numberOfPoints(self):
        """Number of transverse sample points per z level."""
        return int(self.points.shape[0])

    @property
    def numberOfPrisms(self):
        """Total number of wedge prisms in the extruded mesh."""
        self._require_levels()
        return self.numberOfTriangles * (self.levels - 1)

    def getPoints(self):
        return PrimitiveView(
            "point",
            (self.numberOfPoints,),
            {"coordinates": np.asarray(self.points, dtype=np.float64)},
        )

    def getTriangles(self):
        return PrimitiveView(
            "triangle",
            (self.numberOfTriangles,),
            {"pointIndices": np.asarray(self.trianglePointIndices, dtype=np.uint32)},
        )

    def openPmdAttributes(self, context=None):
        self._require_levels()
        context = self._fieldContext() if context is None else context
        return {
            "numberOfPoints": context.numberOfPoints,
            "numberOfTriangles": context.numberOfTriangles,
            "numberOfLevels": context.numberOfLevels,
            "thickness": self.thickness,
        }

    def openPmdFields(self, context=None):
        self._require_levels()
        context = self._fieldContext() if context is None else context
        derived = self._topology()
        points = _flat(self.points, 2, np.float64, "points").reshape((2, context.numberOfPoints))
        yield OpenPmdComponentField(
            "points",
            components={"x": points[0], "y": points[1]},
            axisLabels=["point"],
            context=context,
        )
        yield OpenPmdScalarField("connectivity", _flat(self.trianglePointIndices, 3, np.uint32, "trianglePointIndices"), context)
        yield OpenPmdScalarField("neighbors", derived["triangleNeighbors"], context)
        yield OpenPmdScalarField("forbiddenEdges", derived["forbiddenEdge"], context)
        yield OpenPmdScalarField("normalPoints", derived["triangleNormalPoint"], context)
        yield OpenPmdComponentField(
            "cellCenterX",
            components={"x": derived["triangleCenterX"], "y": derived["triangleCenterY"]},
            axisLabels=["cell"],
            context=context,
            recordName="cell_center",
        )
        yield OpenPmdScalarField("cellNormalX", derived["triangleNormalsX"], context)
        yield OpenPmdScalarField("cellNormalY", derived["triangleNormalsY"], context)
        yield OpenPmdScalarField("surface", derived["triangleSurfaces"], context)

    def _fieldContext(self):
        self._require_levels()
        return type(
            "MeshTopologyOpenPmdContext",
            (),
            {
                "numberOfPoints": self.numberOfPoints,
                "numberOfTriangles": self.numberOfTriangles,
                "numberOfLevels": int(self.levels),
            },
        )()

    def asGainMedium(self):
        """Wrap this topology in an empty ``GainMedium``."""
        return GainMedium(topology=self)

    def _require_levels(self):
        if self.levels is None:
            raise ValueError("numberOfLevels is required")
        if int(self.levels) < 2:
            raise ValueError("numberOfLevels must be at least 2")

    def _require_thickness(self):
        if self.thickness is None:
            raise ValueError("thickness is required")
        if float(self.thickness) <= 0.0:
            raise ValueError("thickness must be positive")

    def _topology(self):
        return _derived_topology(
            np.asarray(self.points, dtype=np.float64),
            np.asarray(self.trianglePointIndices, dtype=np.uint32),
        )


@dataclass
class GainMedium:
    """Mesh plus material/state properties required by ASE and pump models.

    Use ``get(name)`` to inspect a property handle, including expected array
    shape, then set ``prop.value`` or call ``set(name, value)``. Arrays may be
    supplied in matrix shape and are stored internally in Fortran order.
    """

    topology: MeshTopology
    """Geometry and z-level sampling for this medium."""
    physical: dict = field(default_factory=dict)
    """Canonical physical arrays and scalars stored by property name."""
    customFields: dict = field(default_factory=dict)

    def __post_init__(self):
        gmsh = self.topology.metadata.get("gmsh")
        if isinstance(gmsh, Gmsh) and "claddingCellTypes" not in self.physical:
            cladding = gmsh.claddingCellTypes(self.topology)
            if np.any(cladding):
                self.physical["claddingCellTypes"] = cladding

    @classmethod
    def fromVtk(cls, filename, *, numberOfLevels=None, thickness=None):
        """Load topology and physical fields from a legacy VTK wedge file."""
        return gainMediumFromVtk(
            filename,
            MeshTopology,
            cls,
            numberOfLevels=numberOfLevels,
            thickness=thickness,
        )

    @classmethod
    def fromFile(cls, filename, format=None, *, numberOfLevels=None, thickness=None):
        """Load a gain medium; currently supports legacy VTK files."""
        path = Path(filename)
        mesh_format = (format or path.suffix.lstrip(".")).lower()
        if mesh_format in {"vtk", "legacy-vtk"}:
            return cls.fromVtk(path, numberOfLevels=numberOfLevels, thickness=thickness)
        raise NotImplementedError(
            f"gain medium format '{mesh_format}' is not supported yet; supported formats: vtk"
        )

    def withPhysicalProperties(self, **properties):
        """Set several physical properties and return ``self`` for chaining."""
        for name, value in properties.items():
            self.set(name, value)
        return self

    def toVtk(self, filename):
        """Write this gain medium to a legacy ASCII VTK wedge file."""
        return writeGainMediumVtk(filename, self)

    def defineField(
        self,
        name,
        *,
        entity,
        values,
        dtype=None,
        unit="1",
        unitSI=1.0,
        unitDimension=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        dynamic=False,
        backendRequired=False,
    ):
        if not isinstance(name, str) or not _FIELD_NAME_RE.match(name):
            raise ValueError(
                "field name must start with a letter and contain only letters, digits, and underscores"
            )
        canonical = PROPERTY_ALIASES.get(name, name)
        if canonical in PHYSICAL_PROPERTY_SPECS or name in _CUSTOM_FIELD_RESERVED_NAMES:
            raise ValueError(f"field name '{name}' is reserved by the built-in gain-medium field set")
        if not isinstance(unit, str) or not unit.strip():
            raise ValueError("field unit must be a non-empty string")
        try:
            unit_dimension = tuple(float(value) for value in unitDimension)
        except TypeError as exc:
            raise ValueError("field unitDimension must contain seven numeric SI exponents") from exc
        if len(unit_dimension) != 7:
            raise ValueError("field unitDimension must contain seven numeric SI exponents")

        axes = _normalize_custom_entity(entity)
        primitive_shape = _shape_for_custom_entity(self.topology, axes, values)
        raw_value = values.values if isinstance(values, BackendFlatArray) else values
        inferred_dtype = np.asarray(raw_value).dtype if dtype is None else dtype
        dtype_object = np.dtype(inferred_dtype)
        if dtype_object.kind in {"O", "S", "U", "V"}:
            raise ValueError(f"field dtype {dtype_object} is not supported for openPMD transport")

        spec = FieldSpec(
            name=name,
            recordName="custom_" + _camel_to_snake(name),
            axes=axes,
            dtype=dtype_object,
            shape=lambda _context, shape=primitive_shape: shape,
            unit=unit.strip(),
            unitSI=float(unitSI),
            unitDimension=unit_dimension,
            dynamic=bool(dynamic),
            backendRequired=bool(backendRequired),
            userDefined=True,
        )
        return self._storeCustomField(spec, values)

    def withPrimitiveSchema(self, schema, **values):
        for spec in self._schemaFieldSpecs(schema):
            if spec.name in schemaFields and not spec.userDefined:
                continue
            field_spec = self._asUserFieldSpec(spec)
            value = values.get(spec.name)
            if value is None:
                value = np.zeros(field_spec.expectedShape(self._fieldContext()), dtype=field_spec.dtypeObject)
            self._storeCustomField(field_spec, value)
        return self

    def _schemaFieldSpecs(self, schema):
        if isinstance(schema, PrimitiveSchema):
            return schema.fieldSpecs()
        if isinstance(schema, type) and issubclass(schema, PrimitiveSchemaDefinition):
            return schema.fieldSpecs()
        if hasattr(schema, "fieldSpecs"):
            return schema.fieldSpecs()
        raise TypeError("schema must be a Python PrimitiveSchema or PrimitiveSchemaDefinition subclass")

    def _asUserFieldSpec(self, spec):
        primitive_shape = spec.expectedShape(self._fieldContext())
        return FieldSpec(
            name=spec.name,
            recordName=spec.recordName,
            axes=spec.axes,
            dtype=spec.dtypeObject,
            shape=lambda _context, shape=primitive_shape: shape,
            unit=spec.unit,
            unitSI=spec.unitSI,
            unitDimension=spec.unitDimension,
            dynamic=spec.dynamic,
            backendRequired=spec.backendRequired,
            userDefined=True,
            schemaRole=spec.schemaRole,
        )

    def _fieldContext(self):
        self.topology._require_levels()
        return type(
            "GainMediumFieldContext",
            (),
            {
                "numberOfPoints": self.topology.numberOfPoints,
                "numberOfTriangles": self.topology.numberOfTriangles,
                "numberOfLevels": int(self.topology.levels),
            },
        )()

    def _storeCustomField(self, spec, values):
        primitive_shape = spec.expectedShape(self._fieldContext())
        raw_value = values.values if isinstance(values, BackendFlatArray) else values
        arr = np.asarray(raw_value, dtype=spec.dtypeObject)
        if spec.dtypeObject.kind in {"O", "S", "U", "V"}:
            raise ValueError(f"field dtype {spec.dtypeObject} is not supported for openPMD transport")

        expected_size = int(np.prod(primitive_shape, dtype=np.int64))
        if arr.size != expected_size:
            raise ValueError(
                f"{spec.name} expects {expected_size} values for entity {spec.entity} "
                f"with primitive shape {primitive_shape}, got shape {arr.shape}"
            )

        if isinstance(values, BackendFlatArray):
            if arr.ndim != 1:
                raise ValueError(
                    f"{spec.name} marked backend-flat must be a 1-D array with {expected_size} values, got {arr.shape}"
                )
            flat = arr
        elif arr.shape == primitive_shape:
            flat = arr.reshape(-1, order="F")
        elif arr.ndim == 1:
            raise ValueError(
                f"{spec.name} got ambiguous flat array with {arr.size} values for entity {spec.entity}; "
                "pass backendFlat(values) to declare canonical backend-flat order, "
                f"or pass primitive shape {primitive_shape}"
            )
        else:
            raise ValueError(f"{spec.name} expects primitive shape {primitive_shape}, got shape {arr.shape}")

        self.customFields[spec.name] = CustomField(spec=spec, value=flat, primitiveShape=primitive_shape)
        return self

    def getField(self, name):
        if name in self.customFields:
            return self.customFields[name]
        return self.get(name)

    def listFields(self):
        return [field.meta() for field in self.customFields.values()]

    def _customFieldViews(self, *entities):
        wanted = {tuple(entity) for entity in entities}
        return {
            name: field.primitiveView()
            for name, field in self.customFields.items()
            if tuple(field.axes) in wanted
        }

    def emptyBetaCells(self, fill=0.0):
        """Create a correctly shaped point-level beta array filled with ``fill``."""
        return np.full(self.get("betaCells").expectedShape, fill, dtype=np.float64)

    def betaCellIndexAt(self, x, y, z, *, tol=1e-12, flat=False):
        """Forward coordinate lookup for a ``betaCells`` entry."""
        return self.topology.betaCellIndexAt(x, y, z, tol=tol, flat=flat)

    def _fieldView(self, name):
        prop = self.get(name)
        if prop.name not in self.physical:
            value = prop.spec.defaultValue(self.topology)
            if prop.expectedShape == ():
                self.physical[prop.name] = prop.dtype.type(value).item()
            else:
                self.physical[prop.name] = np.asarray(value, dtype=prop.dtype).reshape(-1, order="F")
        return np.asarray(self.physical[prop.name], dtype=prop.dtype).reshape(prop.expectedShape, order="F")

    def _propertyFieldMeta(self, name, *, entity, axes=None):
        prop = self.get(name)
        return {
            "name": prop.name,
            "entity": entity,
            "axes": tuple(axes or ()),
            "dtype": str(prop.dtype),
            "unit": "1",
            "shape": prop.expectedShape,
            "isSet": prop.name in self.physical,
        }

    def _customFieldMetadata(self, *entities):
        wanted = {tuple(entity) for entity in entities}
        return {
            name: field.meta()
            for name, field in self.customFields.items()
            if tuple(field.axes) in wanted
        }

    def getPoints(self):
        position = np.asarray(self.topology.points, dtype=np.float64)
        fields = {
            "position": position,
            "coordinates": position,
            "betaCells": self._fieldView("betaCells"),
            "dntdAse": self._fieldView("dntdAse"),
        }
        custom = self._customFieldViews(("point",), ("point", "level"))
        fields.update(custom)
        position_meta = {
            "name": "position",
            "entity": "point_coordinate",
            "axes": ("point", "coordinate"),
            "dtype": str(np.dtype(np.float64)),
            "unit": "m",
            "shape": position.shape,
            "isSet": True,
        }
        metadata = {
            "position": position_meta,
            "coordinates": {**position_meta, "name": "coordinates"},
            "betaCells": self._propertyFieldMeta("betaCells", entity="point_level", axes=("point", "level")),
            "dntdAse": self._propertyFieldMeta("dntdAse", entity="point_level", axes=("point", "level")),
        }
        metadata.update(self._customFieldMetadata(("point",), ("point", "level")))
        return PrimitiveView("point", (self.topology.numberOfPoints,), fields, metadata)

    def getTriangles(self):
        derived = self.topology._topology()
        center = np.column_stack((derived["triangleCenterX"], derived["triangleCenterY"])).astype(np.float64, copy=False)
        normal = np.stack(
            (
                np.asarray(derived["triangleNormalsX"], dtype=np.float64).reshape((self.topology.numberOfTriangles, 3), order="F"),
                np.asarray(derived["triangleNormalsY"], dtype=np.float64).reshape((self.topology.numberOfTriangles, 3), order="F"),
            ),
            axis=-1,
        )
        connectivity = np.asarray(self.topology.trianglePointIndices, dtype=np.uint32)
        fields = {
            "connectivity": connectivity,
            "pointIndices": connectivity,
            "neighbors": np.asarray(derived["triangleNeighbors"], dtype=np.int32).reshape((self.topology.numberOfTriangles, 3), order="F"),
            "forbiddenEdges": np.asarray(derived["forbiddenEdge"], dtype=np.int32).reshape((self.topology.numberOfTriangles, 3), order="F"),
            "normalPoints": np.asarray(derived["triangleNormalPoint"], dtype=np.uint32).reshape((self.topology.numberOfTriangles, 3), order="F"),
            "center": center,
            "normal": normal,
            "surface": np.asarray(derived["triangleSurfaces"], dtype=np.float32),
            "claddingGroup": self._fieldView("claddingCellTypes"),
            "claddingCellTypes": self._fieldView("claddingCellTypes"),
            "reflectivities": self._fieldView("reflectivities"),
        }
        fields.update(self._customFieldViews(("cell",), ("cell", "local_vertex"), ("cell", "local_side")))
        connectivity_meta = {
            "name": "connectivity",
            "entity": "cell_local_vertex",
            "axes": ("cell", "local_vertex"),
            "dtype": str(np.dtype(np.uint32)),
            "unit": "1",
            "shape": connectivity.shape,
            "isSet": True,
        }
        metadata = {
            "connectivity": connectivity_meta,
            "pointIndices": {**connectivity_meta, "name": "pointIndices"},
            "neighbors": {"name": "neighbors", "entity": "cell_local_side", "axes": ("cell", "local_side"), "dtype": str(np.dtype(np.int32)), "unit": "1", "shape": fields["neighbors"].shape, "isSet": True},
            "forbiddenEdges": {"name": "forbiddenEdges", "entity": "cell_local_side", "axes": ("cell", "local_side"), "dtype": str(np.dtype(np.int32)), "unit": "1", "shape": fields["forbiddenEdges"].shape, "isSet": True},
            "normalPoints": {"name": "normalPoints", "entity": "cell_local_side", "axes": ("cell", "local_side"), "dtype": str(np.dtype(np.uint32)), "unit": "1", "shape": fields["normalPoints"].shape, "isSet": True},
            "center": {"name": "center", "entity": "cell_coordinate", "axes": ("cell", "coordinate"), "dtype": str(np.dtype(np.float64)), "unit": "m", "shape": center.shape, "isSet": True},
            "normal": {"name": "normal", "entity": "cell_local_side_coordinate", "axes": ("cell", "local_side", "coordinate"), "dtype": str(np.dtype(np.float64)), "unit": "1", "shape": normal.shape, "isSet": True},
            "surface": {"name": "surface", "entity": "cell", "axes": ("cell",), "dtype": str(np.dtype(np.float32)), "unit": "m^2", "shape": fields["surface"].shape, "isSet": True},
            "claddingGroup": {**self._propertyFieldMeta("claddingCellTypes", entity="cell", axes=("cell",)), "name": "claddingGroup"},
            "claddingCellTypes": self._propertyFieldMeta("claddingCellTypes", entity="cell", axes=("cell",)),
            "reflectivities": self._propertyFieldMeta("reflectivities", entity="cell_interface", axes=("cell", "interface")),
        }
        metadata.update(self._customFieldMetadata(("cell",), ("cell", "local_vertex"), ("cell", "local_side")))
        return PrimitiveView("triangle", (self.topology.numberOfTriangles,), fields, metadata)

    def getPrisms(self):
        fields = {"betaVolume": self._fieldView("betaVolume")}
        fields.update(self._customFieldViews(("cell", "layer")))
        metadata = {
            "betaVolume": self._propertyFieldMeta("betaVolume", entity="cell_layer", axes=("cell", "layer")),
        }
        metadata.update(self._customFieldMetadata(("cell", "layer")))
        return PrimitiveView("prism", (self.topology.numberOfTriangles, self.topology.levels - 1), fields, metadata)

    def getDomains(self):
        fields = self._customFieldViews(("domain",))
        size = 0
        if fields:
            size = next(iter(fields.values())).shape[0]
        return PrimitiveView("domain", (size,), fields, self._customFieldMetadata(("domain",)))

    def openPmdAttributes(self, context=None):
        return {
            "nTot": self.get("nTot").value,
            "crystalTFluo": self.get("crystalTFluo").value,
            "claddingNumber": self.get("claddingNumber").value,
            "claddingAbsorption": self.get("claddingAbsorption").value,
        }

    def openPmdFields(self, context=None):
        context = self._fieldContext() if context is None else context
        topology = self.topology
        yield OpenPmdScalarField("betaVolume", _flat(self.get("betaVolume").value, topology.levels - 1, np.float64, "betaVolume"), context)
        yield OpenPmdScalarField("pointBeta", _flat(self.get("betaCells").value, topology.levels, np.float64, "betaCells"), context)
        yield OpenPmdScalarField("claddingCellType", _flat(self.get("claddingCellTypes").value, None, np.uint32, "claddingCellTypes"), context)
        yield OpenPmdScalarField("refractiveIndex", _flat(self.get("refractiveIndices").value, None, np.float32, "refractiveIndices"), context)
        yield OpenPmdScalarField("reflectivity", _flat(self.get("reflectivities").value, None, np.float32, "reflectivities"), context)
        for field in self.customFields.values():
            yield OpenPmdScalarField(field.spec.name, field.value, context, prefix="", spec=field.spec)

    def listProperties(self):
        """Return metadata for all known gain-medium properties."""
        return [self.get(name).meta() for name in PHYSICAL_PROPERTY_SPECS]

    def get(self, name):
        """Return a ``GainMediumProperty`` handle by canonical name or alias."""
        canonical = PROPERTY_ALIASES.get(name, name)
        if canonical not in PHYSICAL_PROPERTY_SPECS:
            known = ", ".join(PHYSICAL_PROPERTY_SPECS)
            raise KeyError(f"unknown gain medium property '{name}'. Known properties: {known}")
        return GainMediumProperty(self, PHYSICAL_PROPERTY_SPECS[canonical])

    def set(self, name, value):
        """Validate and store one physical property."""
        prop = self.get(name)
        if prop.expectedShape == ():
            self.physical[prop.name] = prop.dtype.type(value).item()
            return self

        is_backend_flat = isinstance(value, BackendFlatArray)
        raw_value = value.values if is_backend_flat else value
        arr = np.asarray(raw_value, dtype=prop.dtype)
        expected = prop.expectedShape
        expected_size = int(np.prod(expected))
        if arr.size != expected_size:
            raise ValueError(
                f"{prop.name} expects {expected_size} values for shape {expected}, got shape {arr.shape}"
            )

        if is_backend_flat:
            if arr.ndim != 1:
                raise ValueError(
                    f"{prop.name} marked backend-flat must be a 1-D array with {expected_size} values, got {arr.shape}"
                )
            self.physical[prop.name] = arr
            return self

        if arr.shape == expected:
            self.physical[prop.name] = arr.reshape(-1, order="F")
            return self

        if arr.ndim == 1:
            raise ValueError(
                f"{prop.name} got ambiguous flat array with {arr.size} values; pass backendFlat(values) "
                f"to declare canonical backend-flat order, or pass primitive shape {expected}"
            )

        raise ValueError(f"{prop.name} expects primitive shape {expected}, got {arr.shape}")

    @property
    def dntdAse(self):
        """ASE contribution to ``d beta / dt`` at the sample points."""
        return self.get("dntdAse").value

    @property
    def numberOfTriangles(self):
        return self.topology.numberOfTriangles

    @property
    def numberOfPoints(self):
        return self.topology.numberOfPoints

    @property
    def numberOfPrisms(self):
        return self.topology.numberOfPrisms

    @property
    def numberOfLevels(self):
        """Number of z sample planes in the attached topology."""
        self.topology._require_levels()
        return int(self.topology.levels)


GainMediumGeometry = MeshTopology
