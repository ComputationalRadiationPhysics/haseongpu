# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_GMSH_TRIANGLE = 2


@dataclass(frozen=True)
class GmshElement:
    element_id: int
    element_type: int
    node_ids: tuple[int, ...]
    physical_tag: int | None = None


@dataclass
class Gmsh:
    """gmsh-backed mesh data accepted by MeshTopology."""

    nodes: dict[int, tuple[float, float, float]]
    elements: list[GmshElement]
    physical_names: dict[int, dict[int, str]] = field(default_factory=dict)
    source: str | None = None

    @classmethod
    def fromFile(cls, filename):
        return _read_gmsh(filename)

    def topology(self, *, numberOfLevels=None, thickness=None):
        return _from_gmsh(self, numberOfLevels=numberOfLevels, thickness=thickness)

    def claddingCellTypes(self, topology):
        if topology.metadata.get("gmsh") is not self:
            raise ValueError("claddingCellTypes requires the MeshTopology created from this Gmsh object")
        return _gmsh_cladding_cell_types(self, topology)


def _gmsh_module():
    try:
        import gmsh as gmsh_api
    except ImportError as exc:
        raise ImportError("gmsh mesh loading requires the 'gmsh' Python package") from exc
    return gmsh_api


def _read_gmsh(path):
    gmsh_api = _gmsh_module()
    path = Path(path)
    owned_session = not gmsh_api.isInitialized()
    if owned_session:
        gmsh_api.initialize()
    try:
        gmsh_api.option.setNumber("General.Terminal", 0)
        gmsh_api.clear()
        gmsh_api.open(str(path))
        nodes = _read_nodes(gmsh_api)
        physical_names = _read_physical_names(gmsh_api)
        elements = _read_elements(gmsh_api)
    finally:
        if owned_session:
            gmsh_api.finalize()
    return Gmsh(nodes=nodes, elements=elements, physical_names=physical_names, source=str(path))


def _read_nodes(gmsh_api):
    node_tags, coords, _ = gmsh_api.model.mesh.getNodes()
    coords = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
    return {int(tag): tuple(coords[i]) for i, tag in enumerate(node_tags)}


def _read_physical_names(gmsh_api):
    names = {}
    for dim, tag in gmsh_api.model.getPhysicalGroups():
        name = gmsh_api.model.getPhysicalName(dim, tag)
        names.setdefault(int(dim), {})[int(tag)] = name
    return names


def _read_elements(gmsh_api):
    physical_by_element = {}
    for _, entity_tag in gmsh_api.model.getEntities(2):
        physical = gmsh_api.model.getPhysicalGroupsForEntity(2, entity_tag)
        physical_tag = int(physical[0]) if len(physical) else None
        _, tags, _ = gmsh_api.model.mesh.getElements(2, entity_tag)
        for type_tags in tags:
            for element_id in type_tags:
                physical_by_element[int(element_id)] = physical_tag

    elements = []
    types, tags, node_tags = gmsh_api.model.mesh.getElements()
    for element_type, element_ids, flattened_nodes in zip(types, tags, node_tags):
        node_count = len(flattened_nodes) // len(element_ids)
        grouped_nodes = np.asarray(flattened_nodes, dtype=np.uint64).reshape(-1, node_count)
        for element_id, nodes in zip(element_ids, grouped_nodes):
            element_id = int(element_id)
            elements.append(
                GmshElement(
                    element_id=element_id,
                    element_type=int(element_type),
                    node_ids=tuple(int(node_id) for node_id in nodes),
                    physical_tag=physical_by_element.get(element_id),
                )
            )
    return elements


def _from_gmsh(gmsh, *, numberOfLevels=None, thickness=None):
    triangles = [e for e in gmsh.elements if e.element_type == _GMSH_TRIANGLE]
    if not triangles:
        raise ValueError("gmsh import supports two-dimensional triangle meshes only")
    if numberOfLevels is None or thickness is None:
        raise ValueError("2D gmsh triangle meshes require numberOfLevels and thickness")
    levels = int(numberOfLevels)
    dz = float(thickness)
    if levels < 2:
        raise ValueError("numberOfLevels must be at least 2")
    if dz <= 0.0:
        raise ValueError("thickness must be positive")
    points, triangle_indices = _gmsh_triangles_to_2d(gmsh, triangles)
    return _mesh_topology()(
        points,
        triangle_indices,
        levels=levels,
        thickness=dz,
        metadata={"source": gmsh.source, "format": "gmsh", "gmsh": gmsh, "dimension": 2},
    )


def _mesh_topology():
    from .core import MeshTopology

    return MeshTopology


def _gmsh_triangles_to_2d(gmsh, triangles):
    z_values = [gmsh.nodes[node_id][2] for tri in triangles for node_id in tri.node_ids[:3]]
    if not np.allclose(z_values, z_values[0]):
        raise ValueError("2D gmsh triangle meshes must be planar")
    return _base_triangles(gmsh, triangles)


def _base_triangles(gmsh, elements):
    index = {}
    points = []
    triangles = []
    for element in elements:
        ids = []
        for node_id in element.node_ids[:3]:
            x, y, _ = gmsh.nodes[node_id]
            key = (float(x), float(y))
            if key not in index:
                index[key] = len(points)
                points.append(key)
            ids.append(index[key])
        triangles.append(ids)
    return np.asarray(points, dtype=np.float64), np.asarray(triangles, dtype=np.uint32)


def _gmsh_cladding_cell_types(gmsh, topology):
    cladding_tags = {
        tag for tag, name in gmsh.physical_names.get(2, {}).items()
        if "cladding" in name.lower()
    }
    values = np.zeros(topology.numberOfTriangles, dtype=np.uint32)
    triangles = [e for e in gmsh.elements if e.element_type == _GMSH_TRIANGLE]
    if len(triangles) != topology.numberOfTriangles:
        raise ValueError("gmsh triangle count must match topology.numberOfTriangles to map cladding")
    for triangle_index, triangle in enumerate(triangles):
        if triangle.physical_tag in cladding_tags:
            values[triangle_index] = np.uint32(triangle.physical_tag)
    return values
