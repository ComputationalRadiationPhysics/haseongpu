# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Explicit 3D volume-cell topology helpers.

The legacy frontend represents wedge prisms implicitly as ``triangle × z-layer``.
This module provides the first explicit cell/face topology used by the migration
away from layout-derived neighbors.  The first supported runtime cell is a
linear six-node prism (VTK_WEDGE / gmsh prism type 6).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


VTK_WEDGE = np.uint32(13)
GMSH_PRISM6 = 6
GMSH_TET4 = 4
GMSH_TET10 = 11
GMSH_TRI3 = 2
GMSH_QUAD4 = 3

# Local faces for a six-node prism.  The first two are triangular bases; the
# remaining faces are quads.  A value of -1 pads triangular faces to width four.
PRISM6_FACE_VERTICES = np.asarray(
    [
        [0, 2, 1, -1],
        [3, 4, 5, -1],
        [0, 1, 4, 3],
        [1, 2, 5, 4],
        [2, 0, 3, 5],
    ],
    dtype=np.int8,
)

BOUND_STOP = np.int32(-1)
BOUND_INTERNAL = np.int32(0)


def _as_points3(points):
    arr = np.asarray(points, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"volume topology points must have shape (N, 3), got {arr.shape}")
    return np.ascontiguousarray(arr)


def _as_prisms(prism_point_indices):
    arr = np.asarray(prism_point_indices, dtype=np.uint32)
    if arr.ndim != 2 or arr.shape[1] != 6:
        raise ValueError(f"Prism6 connectivity must have shape (N, 6), got {arr.shape}")
    return np.ascontiguousarray(arr)


def _face_nodes(cell_nodes, local_face):
    local = PRISM6_FACE_VERTICES[int(local_face)]
    return tuple(int(cell_nodes[int(i)]) for i in local if i >= 0)


def _face_key(nodes):
    return tuple(sorted(int(node) for node in nodes))


def _triangle_area_normal(a, b, c):
    raw = np.cross(b - a, c - a)
    norm = float(np.linalg.norm(raw))
    if norm == 0.0:
        return 0.0, np.zeros(3, dtype=np.float64)
    return 0.5 * norm, raw / norm


def _face_geometry(points, nodes, cell_center):
    coords = points[np.asarray(nodes, dtype=np.uint32)]
    center = np.mean(coords, axis=0)
    if len(nodes) == 3:
        area, normal = _triangle_area_normal(coords[0], coords[1], coords[2])
    elif len(nodes) == 4:
        area_a, normal_a = _triangle_area_normal(coords[0], coords[1], coords[2])
        area_b, normal_b = _triangle_area_normal(coords[0], coords[2], coords[3])
        area = area_a + area_b
        raw = normal_a * area_a + normal_b * area_b
        norm = float(np.linalg.norm(raw))
        normal = raw / norm if norm else np.zeros(3, dtype=np.float64)
    else:
        raise ValueError("Prism6 faces must have three or four nodes")
    if np.dot(normal, cell_center - center) > 0.0:
        normal = -normal
    return center, np.float32(area), normal.astype(np.float32)


def _tet_volume(a, b, c, d):
    return abs(float(np.dot(np.cross(b - a, c - a), d - a))) / 6.0


def _prism_volume(points, nodes):
    p = points[np.asarray(nodes, dtype=np.uint32)]
    return np.float32(
        _tet_volume(p[0], p[1], p[2], p[3])
        + _tet_volume(p[1], p[2], p[4], p[3])
        + _tet_volume(p[2], p[4], p[5], p[3])
    )


@dataclass
class VolumeTopology:
    """Explicit 3D volume-cell topology.

    The initial implementation supports linear Prism6/WEDGE cells only.  The
    arrays are deliberately explicit so neighbor traversal does not depend on a
    ``triangle × level`` layout.
    """

    points: np.ndarray
    cellPointIndices: np.ndarray
    cellTypes: np.ndarray | None = None
    cellDomains: np.ndarray | None = None
    faceBoundaries: np.ndarray | None = None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.points = _as_points3(self.points)
        self.cellPointIndices = _as_prisms(self.cellPointIndices)
        n_cells = self.cellPointIndices.shape[0]
        if self.cellTypes is None:
            self.cellTypes = np.full(n_cells, VTK_WEDGE, dtype=np.uint32)
        else:
            self.cellTypes = np.asarray(self.cellTypes, dtype=np.uint32)
        if self.cellTypes.shape != (n_cells,):
            raise ValueError("cellTypes must have shape (numberOfCells,)")
        if np.any(self.cellTypes != VTK_WEDGE):
            raise NotImplementedError("VolumeTopology currently supports only Prism6 / VTK_WEDGE cells")
        if self.cellDomains is None:
            self.cellDomains = np.ones(n_cells, dtype=np.int32)
        else:
            self.cellDomains = np.asarray(self.cellDomains, dtype=np.int32)
        if self.cellDomains.shape != (n_cells,):
            raise ValueError("cellDomains must have shape (numberOfCells,)")
        derived = _derive_prism6_topology(self.points, self.cellPointIndices, self.faceBoundaries)
        self.facePointIndices = derived["facePointIndices"]
        self.neighborCells = derived["neighborCells"]
        self.neighborLocalFaces = derived["neighborLocalFaces"]
        self.faceBoundaries = derived["faceBoundaries"]
        self.faceCenters = derived["faceCenters"]
        self.faceNormals = derived["faceNormals"]
        self.faceAreas = derived["faceAreas"]
        self.cellCenters = derived["cellCenters"]
        self.cellVolumes = derived["cellVolumes"]
        self.samplePoints = np.asarray(self.cellCenters, dtype=np.float64)

    @classmethod
    def fromPrisms(
        cls,
        points,
        prismPointIndices,
        *,
        cellDomains=None,
        faceBoundaries=None,
        metadata=None,
    ):
        return cls(points, prismPointIndices, cellDomains=cellDomains, faceBoundaries=faceBoundaries, metadata=metadata or {})

    @classmethod
    def fromGmsh(cls, gmsh, *, boundaryDefault=BOUND_STOP):
        from .msh import Gmsh

        if isinstance(gmsh, (str, Path)):
            gmsh = Gmsh.fromFile(gmsh)
        if not isinstance(gmsh, Gmsh):
            raise TypeError("fromGmsh expects a Gmsh instance or a gmsh .msh filename")
        return _from_gmsh_volume(gmsh, boundaryDefault=boundaryDefault)

    @classmethod
    def fromFile(cls, filename, format=None, **kwargs):
        path = Path(filename)
        mesh_format = (format or path.suffix.lstrip(".")).lower()
        if mesh_format in {"msh", "gmsh"}:
            return cls.fromGmsh(path, **kwargs)
        raise NotImplementedError(f"volume mesh format '{mesh_format}' is not supported yet; supported formats: gmsh")

    @classmethod
    def fromExtrudedTopology(cls, topology):
        """Convert a legacy triangle×level topology into explicit Prism6 cells."""
        topology._require_levels()
        points2 = np.asarray(topology.points, dtype=np.float64)
        z_values = topology.levelCoordinates()
        points3 = np.empty((topology.numberOfPoints * int(topology.levels), 3), dtype=np.float64)
        for level, z in enumerate(z_values):
            start = level * topology.numberOfPoints
            stop = start + topology.numberOfPoints
            points3[start:stop, :2] = points2
            points3[start:stop, 2] = z
        rows = []
        for level in range(int(topology.levels) - 1):
            lower = level * topology.numberOfPoints
            upper = (level + 1) * topology.numberOfPoints
            for tri in np.asarray(topology.trianglePointIndices, dtype=np.uint32):
                ids = [int(v) for v in tri]
                rows.append([ids[0] + lower, ids[1] + lower, ids[2] + lower, ids[0] + upper, ids[1] + upper, ids[2] + upper])
        return cls(points3, np.asarray(rows, dtype=np.uint32), metadata={"source": "legacy-extruded", "legacyTopology": topology})

    @property
    def numberOfPoints(self):
        return int(self.points.shape[0])

    @property
    def numberOfCells(self):
        return int(self.cellPointIndices.shape[0])

    @property
    def numberOfFacesPerCell(self):
        return int(PRISM6_FACE_VERTICES.shape[0])

    @property
    def numberOfSamplePoints(self):
        return int(self.samplePoints.shape[0])

    @property
    def numberOfPrisms(self):
        return self.numberOfCells

    def openPmdAttributes(self, context=None):
        context = context or self
        return {
            "numberOfPoints": int(context.numberOfPoints),
            "numberOfTriangles": int(getattr(context, "numberOfCells", self.numberOfCells)),
            "numberOfLevels": 1,
            "thickness": 0.0,
        }

    def cellsConnectivityFlat(self):
        return np.asarray(self.cellPointIndices, dtype=np.uint32).reshape(-1)

    def cellsOffsets(self):
        return np.arange(self.numberOfCells + 1, dtype=np.uint32) * np.uint32(6)

    def faceConnectivityFlat(self):
        return np.asarray(self.facePointIndices, dtype=np.int32).reshape(-1)


def _derive_prism6_topology(points, cells, face_boundaries):
    n_cells = cells.shape[0]
    n_faces = PRISM6_FACE_VERTICES.shape[0]
    face_nodes = np.full((n_cells, n_faces, 4), -1, dtype=np.int32)
    neighbor_cells = np.full((n_cells, n_faces), -1, dtype=np.int32)
    neighbor_faces = np.full((n_cells, n_faces), -1, dtype=np.int32)
    boundaries = np.zeros((n_cells, n_faces), dtype=np.int32) if face_boundaries is None else np.asarray(face_boundaries, dtype=np.int32).copy()
    if boundaries.shape != (n_cells, n_faces):
        raise ValueError(f"faceBoundaries must have shape {(n_cells, n_faces)}, got {boundaries.shape}")

    owners = {}
    for cell_i in range(n_cells):
        for face_i in range(n_faces):
            nodes = _face_nodes(cells[cell_i], face_i)
            face_nodes[cell_i, face_i, : len(nodes)] = nodes
            key = _face_key(nodes)
            owner = owners.get(key)
            if owner is None:
                owners[key] = (cell_i, face_i)
            else:
                other_cell, other_face = owner
                neighbor_cells[cell_i, face_i] = other_cell
                neighbor_faces[cell_i, face_i] = other_face
                neighbor_cells[other_cell, other_face] = cell_i
                neighbor_faces[other_cell, other_face] = face_i
                if boundaries[cell_i, face_i] == BOUND_INTERNAL and boundaries[other_cell, other_face] != BOUND_INTERNAL:
                    boundaries[cell_i, face_i] = boundaries[other_cell, other_face]
                elif boundaries[other_cell, other_face] == BOUND_INTERNAL and boundaries[cell_i, face_i] != BOUND_INTERNAL:
                    boundaries[other_cell, other_face] = boundaries[cell_i, face_i]

    boundaries[neighbor_cells < 0] = np.where(boundaries[neighbor_cells < 0] == BOUND_INTERNAL, BOUND_STOP, boundaries[neighbor_cells < 0])

    cell_centers = np.empty((n_cells, 3), dtype=np.float64)
    cell_volumes = np.empty(n_cells, dtype=np.float32)
    face_centers = np.empty((n_cells, n_faces, 3), dtype=np.float32)
    face_normals = np.empty((n_cells, n_faces, 3), dtype=np.float32)
    face_areas = np.empty((n_cells, n_faces), dtype=np.float32)

    for cell_i in range(n_cells):
        nodes = cells[cell_i]
        cell_centers[cell_i] = np.mean(points[nodes], axis=0)
        cell_volumes[cell_i] = _prism_volume(points, nodes)
        for face_i in range(n_faces):
            valid = [int(n) for n in face_nodes[cell_i, face_i] if n >= 0]
            center, area, normal = _face_geometry(points, valid, cell_centers[cell_i])
            face_centers[cell_i, face_i] = center.astype(np.float32)
            face_areas[cell_i, face_i] = area
            face_normals[cell_i, face_i] = normal

    return {
        "facePointIndices": face_nodes,
        "neighborCells": neighbor_cells,
        "neighborLocalFaces": neighbor_faces,
        "faceBoundaries": boundaries,
        "faceCenters": face_centers,
        "faceNormals": face_normals,
        "faceAreas": face_areas,
        "cellCenters": cell_centers,
        "cellVolumes": cell_volumes,
    }


def _from_gmsh_volume(gmsh, *, boundaryDefault):
    prism_elements = [element for element in gmsh.elements if element.element_type == GMSH_PRISM6]
    unsupported = sorted({element.element_type for element in gmsh.elements if element.element_type in {GMSH_TET4, GMSH_TET10}})
    if unsupported:
        raise NotImplementedError("gmsh volume import M1 supports Prism6 only; TET4/TET10 are reserved for later milestones")
    if not prism_elements:
        raise ValueError("gmsh volume import requires at least one Prism6 element")

    used_node_ids = sorted({node_id for element in gmsh.elements for node_id in element.node_ids})
    node_to_index = {node_id: index for index, node_id in enumerate(used_node_ids)}
    points = np.asarray([gmsh.nodes[node_id] for node_id in used_node_ids], dtype=np.float64)

    cells = np.asarray([[node_to_index[node_id] for node_id in element.node_ids[:6]] for element in prism_elements], dtype=np.uint32)
    domains = np.asarray([element.physical_tag or 1 for element in prism_elements], dtype=np.int32)

    surface_boundaries = {}
    for element in gmsh.elements:
        if element.element_type not in {GMSH_TRI3, GMSH_QUAD4}:
            continue
        nodes = tuple(node_to_index[node_id] for node_id in element.node_ids[: (3 if element.element_type == GMSH_TRI3 else 4)] if node_id in node_to_index)
        if len(nodes) not in {3, 4}:
            continue
        surface_boundaries[_face_key(nodes)] = np.int32(element.physical_tag if element.physical_tag is not None else boundaryDefault)

    face_boundaries = np.zeros((cells.shape[0], PRISM6_FACE_VERTICES.shape[0]), dtype=np.int32)
    for cell_i, cell in enumerate(cells):
        for face_i in range(PRISM6_FACE_VERTICES.shape[0]):
            key = _face_key(_face_nodes(cell, face_i))
            if key in surface_boundaries:
                face_boundaries[cell_i, face_i] = surface_boundaries[key]

    return VolumeTopology(
        points,
        cells,
        cellDomains=domains,
        faceBoundaries=face_boundaries,
        metadata={"source": gmsh.source, "format": "gmsh", "dimension": 3, "gmsh": gmsh},
    )
