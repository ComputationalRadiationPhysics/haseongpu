# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Explicit Tet4 volume-cell topology helpers."""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from pathlib import Path
import math
import warnings

import numpy as np

try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        def decorator(fn):
            return fn
        return decorator

from .domains import DomainMap


VTK_TETRA = np.uint32(10)
GMSH_TET4 = 4
GMSH_TRI3 = 2
GMSH_TET10 = 11
GMSH_UNSUPPORTED_VOLUME_TYPES = frozenset({6, GMSH_TET10})

# Local faces for a four-node tetrahedron.  The orientation is chosen to be
# outward for a positively-oriented tetrahedron, but downstream code only
# requires consistent local lookup tables and vertex membership.
TET4_FACE_VERTICES = np.asarray(
    [
        [0, 2, 1],
        [0, 1, 3],
        [1, 2, 3],
        [2, 0, 3],
    ],
    dtype=np.int8,
)

BOUND_STOP = np.int32(-1)
BOUND_INTERNAL = np.int32(0)


def _asPoints3(points):
    arr = np.asarray(points, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"volume topology points must have shape (N, 3), got {arr.shape}")
    return np.ascontiguousarray(arr)


def _asTetrahedra(cellPointIndices):
    arr = np.asarray(cellPointIndices, dtype=np.uint32)
    if arr.ndim != 2 or arr.shape[1] != 4:
        raise ValueError(f"Tet4 connectivity must have shape (N, 4), got {arr.shape}")
    return np.ascontiguousarray(arr)


def _faceNodes(cellNodes, localFace):
    local = TET4_FACE_VERTICES[int(localFace)]
    return tuple(int(cellNodes[int(i)]) for i in local)


def _faceKey(nodes):
    return tuple(sorted(int(node) for node in nodes))


_TRIANGLE_KEY_DTYPE = np.dtype(
    [("node0", np.uint32), ("node1", np.uint32), ("node2", np.uint32)]
)


def _triangleKeys(triangles):
    rows = np.ascontiguousarray(
        np.sort(np.asarray(triangles, dtype=np.uint32), axis=1)
    )
    return rows.view(_TRIANGLE_KEY_DTYPE).reshape(-1)


def _surfaceBoundariesFromTriangles(cells, surfaceNodes, surfaceTags):
    """Map tagged Gmsh triangles to every matching tetra-local face."""
    flatBoundaries = np.zeros(cells.shape[0] * TET4_FACE_VERTICES.shape[0], dtype=np.int32)
    if len(surfaceNodes) == 0:
        return flatBoundaries.reshape((-1, TET4_FACE_VERTICES.shape[0]))

    surfaceKeys = _triangleKeys(surfaceNodes)
    surfaceTags = np.asarray(surfaceTags, dtype=np.int32)
    # Match the former dictionary behavior: the last duplicate TRI3 tag wins.
    _, reverseFirst = np.unique(surfaceKeys[::-1], return_index=True)
    keep = surfaceKeys.size - 1 - reverseFirst
    surfaceKeys = surfaceKeys[keep]
    surfaceTags = surfaceTags[keep]

    faceNodes = cells[:, TET4_FACE_VERTICES].reshape(-1, 3)
    faceKeys = _triangleKeys(faceNodes)
    faceOrder = np.argsort(faceKeys, kind="stable")
    sortedFaceKeys = faceKeys[faceOrder]
    lower = np.searchsorted(sortedFaceKeys, surfaceKeys, side="left")
    upper = np.searchsorted(sortedFaceKeys, surfaceKeys, side="right")
    counts = upper - lower
    totalMatches = int(np.sum(counts))
    if totalMatches == 0:
        return flatBoundaries.reshape((-1, TET4_FACE_VERTICES.shape[0]))

    repeatedLower = np.repeat(lower, counts)
    repeatedOffsets = np.repeat(np.cumsum(counts) - counts, counts)
    sortedPositions = repeatedLower + np.arange(totalMatches) - repeatedOffsets
    flatBoundaries[faceOrder[sortedPositions]] = np.repeat(surfaceTags, counts)
    return flatBoundaries.reshape((-1, TET4_FACE_VERTICES.shape[0]))


def _triangleAreaNormal(a, b, c):
    raw = np.cross(b - a, c - a)
    norm = float(np.linalg.norm(raw))
    if norm == 0.0:
        return 0.0, np.zeros(3, dtype=np.float64)
    return 0.5 * norm, raw / norm


def _faceGeometry(points, nodes, cellCenter):
    coords = points[np.asarray(nodes, dtype=np.uint32)]
    center = np.mean(coords, axis=0)
    area, normal = _triangleAreaNormal(coords[0], coords[1], coords[2])
    if np.dot(normal, cellCenter - center) > 0.0:
        normal = -normal
    return center, np.float32(area), normal.astype(np.float32)


def _tetVolume(a, b, c, d):
    return abs(float(np.dot(np.cross(b - a, c - a), d - a))) / 6.0


@njit(cache=True)
def _tet4GeometryKernel(points, cells, localFaces):
    """Build per-cell and per-face geometry without Python object overhead."""
    numberOfCells = cells.shape[0]
    numberOfFaces = localFaces.shape[0]
    faceNodes = np.empty((numberOfCells, numberOfFaces, 3), dtype=np.int32)
    cellCenters = np.empty((numberOfCells, 3), dtype=np.float64)
    cellVolumes = np.empty(numberOfCells, dtype=np.float32)
    faceCenters = np.empty((numberOfCells, numberOfFaces, 3), dtype=np.float32)
    faceNormals = np.empty((numberOfCells, numberOfFaces, 3), dtype=np.float32)
    faceAreas = np.empty((numberOfCells, numberOfFaces), dtype=np.float32)

    for cellIndex in range(numberOfCells):
        p0 = cells[cellIndex, 0]
        p1 = cells[cellIndex, 1]
        p2 = cells[cellIndex, 2]
        p3 = cells[cellIndex, 3]
        centerX = (points[p0, 0] + points[p1, 0] + points[p2, 0] + points[p3, 0]) * 0.25
        centerY = (points[p0, 1] + points[p1, 1] + points[p2, 1] + points[p3, 1]) * 0.25
        centerZ = (points[p0, 2] + points[p1, 2] + points[p2, 2] + points[p3, 2]) * 0.25
        cellCenters[cellIndex, 0] = centerX
        cellCenters[cellIndex, 1] = centerY
        cellCenters[cellIndex, 2] = centerZ

        abX = points[p1, 0] - points[p0, 0]
        abY = points[p1, 1] - points[p0, 1]
        abZ = points[p1, 2] - points[p0, 2]
        acX = points[p2, 0] - points[p0, 0]
        acY = points[p2, 1] - points[p0, 1]
        acZ = points[p2, 2] - points[p0, 2]
        adX = points[p3, 0] - points[p0, 0]
        adY = points[p3, 1] - points[p0, 1]
        adZ = points[p3, 2] - points[p0, 2]
        crossX = abY * acZ - abZ * acY
        crossY = abZ * acX - abX * acZ
        crossZ = abX * acY - abY * acX
        signedVolume = crossX * adX + crossY * adY + crossZ * adZ
        cellVolumes[cellIndex] = abs(signedVolume) / 6.0

        for faceIndex in range(numberOfFaces):
            node0 = cells[cellIndex, localFaces[faceIndex, 0]]
            node1 = cells[cellIndex, localFaces[faceIndex, 1]]
            node2 = cells[cellIndex, localFaces[faceIndex, 2]]
            faceNodes[cellIndex, faceIndex, 0] = node0
            faceNodes[cellIndex, faceIndex, 1] = node1
            faceNodes[cellIndex, faceIndex, 2] = node2

            faceCenterX = (points[node0, 0] + points[node1, 0] + points[node2, 0]) / 3.0
            faceCenterY = (points[node0, 1] + points[node1, 1] + points[node2, 1]) / 3.0
            faceCenterZ = (points[node0, 2] + points[node1, 2] + points[node2, 2]) / 3.0
            faceCenters[cellIndex, faceIndex, 0] = faceCenterX
            faceCenters[cellIndex, faceIndex, 1] = faceCenterY
            faceCenters[cellIndex, faceIndex, 2] = faceCenterZ

            edge0X = points[node1, 0] - points[node0, 0]
            edge0Y = points[node1, 1] - points[node0, 1]
            edge0Z = points[node1, 2] - points[node0, 2]
            edge1X = points[node2, 0] - points[node0, 0]
            edge1Y = points[node2, 1] - points[node0, 1]
            edge1Z = points[node2, 2] - points[node0, 2]
            normalX = edge0Y * edge1Z - edge0Z * edge1Y
            normalY = edge0Z * edge1X - edge0X * edge1Z
            normalZ = edge0X * edge1Y - edge0Y * edge1X
            normalLength = np.sqrt(normalX * normalX + normalY * normalY + normalZ * normalZ)
            faceAreas[cellIndex, faceIndex] = 0.5 * normalLength
            if normalLength > 0.0:
                normalX /= normalLength
                normalY /= normalLength
                normalZ /= normalLength
            else:
                normalX = 0.0
                normalY = 0.0
                normalZ = 0.0
            pointsInward = (
                normalX * (centerX - faceCenterX)
                + normalY * (centerY - faceCenterY)
                + normalZ * (centerZ - faceCenterZ)
            )
            if pointsInward > 0.0:
                normalX = -normalX
                normalY = -normalY
                normalZ = -normalZ
            faceNormals[cellIndex, faceIndex, 0] = normalX
            faceNormals[cellIndex, faceIndex, 1] = normalY
            faceNormals[cellIndex, faceIndex, 2] = normalZ

    return faceNodes, cellCenters, cellVolumes, faceCenters, faceNormals, faceAreas


@dataclass
class VolumeTopology:
    """Explicit Tet4 volume-cell topology."""

    points: np.ndarray
    cellPointIndices: np.ndarray
    cellTypes: np.ndarray | None = None
    cellDomains: np.ndarray | None = None
    faceBoundaries: np.ndarray | None = None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.points = _asPoints3(self.points)
        self.cellPointIndices = _asTetrahedra(self.cellPointIndices)
        numberOfCells = self.cellPointIndices.shape[0]
        if self.cellTypes is None:
            self.cellTypes = np.full(numberOfCells, VTK_TETRA, dtype=np.uint32)
        else:
            self.cellTypes = np.asarray(self.cellTypes, dtype=np.uint32)
        if self.cellTypes.shape != (numberOfCells,):
            raise ValueError("cellTypes must have shape (numberOfCells,)")
        if np.any(self.cellTypes != VTK_TETRA):
            raise NotImplementedError("VolumeTopology supports Tet4 / VTK_TETRA cells only")
        if self.cellDomains is None:
            self.cellDomains = np.ones(numberOfCells, dtype=np.int32)
        else:
            self.cellDomains = np.asarray(self.cellDomains, dtype=np.int32)
        if self.cellDomains.shape != (numberOfCells,):
            raise ValueError("cellDomains must have shape (numberOfCells,)")
        derived = _deriveTet4Topology(self.points, self.cellPointIndices, self.faceBoundaries)
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
    def fromTetrahedra(
        cls,
        points,
        cellPointIndices,
        *,
        cellDomains=None,
        faceBoundaries=None,
        metadata=None,
    ):
        return cls(
            points,
            cellPointIndices,
            cellDomains=cellDomains,
            faceBoundaries=faceBoundaries,
            metadata=metadata or {},
        )

    @classmethod
    def fromGmsh(cls, gmsh, *, boundaryDefault=BOUND_STOP):
        from .msh import Gmsh

        if isinstance(gmsh, (str, Path)):
            gmsh = Gmsh.fromFile(gmsh)
        if not isinstance(gmsh, Gmsh):
            raise TypeError("fromGmsh expects a Gmsh instance or a gmsh .msh filename")
        return _fromGmshVolume(gmsh, boundaryDefault=boundaryDefault)

    @classmethod
    def fromVtk(cls, filename):
        from .vtk import volumeTopologyFromVtk

        return volumeTopologyFromVtk(filename, cls)

    @classmethod
    def fromStl(cls, filename, *, boundaryDefault=BOUND_STOP, meshSize=None):
        return _fromStlVolume(filename, boundaryDefault=boundaryDefault, meshSize=meshSize)

    @classmethod
    def fromFile(cls, filename, format=None, **kwargs):
        path = Path(filename)
        meshFormat = (format or path.suffix.lstrip(".")).lower()
        if meshFormat in {"msh", "gmsh"}:
            return cls.fromGmsh(path, **kwargs)
        if meshFormat in {"vtk"}:
            return cls.fromVtk(path)
        if meshFormat in {"stl", "ascii-stl", "binary-stl", "dae/stl", "dea/stl"}:
            return cls.fromStl(path, **kwargs)
        raise NotImplementedError(f"volume mesh format '{meshFormat}' is not supported; supported formats: gmsh, vtk, stl")

    @property
    def numberOfPoints(self):
        return int(self.points.shape[0])

    @property
    def numberOfCells(self):
        return int(self.cellPointIndices.shape[0])

    @property
    def numberOfFacesPerCell(self):
        return int(TET4_FACE_VERTICES.shape[0])

    @property
    def numberOfSamplePoints(self):
        return int(self.samplePoints.shape[0])

    @property
    def structuredNumberOfPoints(self):
        structured = self.metadata.get("structured", {}) if isinstance(self.metadata, dict) else {}
        return int(structured.get("numberOfPoints", self.numberOfSamplePoints))

    @property
    def structuredNumberOfLevels(self):
        structured = self.metadata.get("structured", {}) if isinstance(self.metadata, dict) else {}
        return int(structured.get("numberOfLevels", 1))

    @property
    def structuredThickness(self):
        structured = self.metadata.get("structured", {}) if isinstance(self.metadata, dict) else {}
        return float(structured.get("thickness", 0.0))

    @property
    def numberOfPrisms(self):
        return self.numberOfCells

    @property
    def cellDomainNames(self):
        return dict(self.metadata.get("cellDomainNames", {})) if isinstance(self.metadata, dict) else {}

    @property
    def surfaceDomainNames(self):
        return dict(self.metadata.get("surfaceDomainNames", {})) if isinstance(self.metadata, dict) else {}

    def cellDomainMap(self):
        return DomainMap(self.cellDomainNames)

    def surfaceDomainMap(self):
        return DomainMap(self.surfaceDomainNames)

    def withDomains(self, *, cellDomains=None, surfaceDomains=None):
        topology = self
        if cellDomains is not None:
            topology = topology.withCellDomains(cellDomains)
        if surfaceDomains is not None:
            topology = topology.withSurfaceDomains(surfaceDomains)
        return topology

    def withCellDomains(self, assignments=None, **assignment):
        """Return a copy with selected cells assigned to positive domain ids."""
        cellDomains = np.asarray(self.cellDomains, dtype=np.int32).copy()
        names = self.cellDomainNames
        domainMap = DomainMap(names)
        for spec in _normalizeAssignments(assignments, assignment):
            mask = _selectCells(self, cellDomains, spec)
            target = _targetDomain(domainMap, spec)
            if not np.any(mask):
                raise ValueError(f"cell domain assignment selected no cells: {spec}")
            cellDomains[mask] = np.int32(target)
            _updateDomainNames(names, domainMap, spec, target)
        metadata = {**self.metadata, "cellDomainNames": names}
        return self._copyWith(cellDomains=cellDomains, metadata=metadata)

    def withSurfaceDomains(self, assignments=None, *, allowInternal=False, **assignment):
        """Return a copy with selected faces assigned to positive surface-domain ids."""
        faceBoundaries = np.asarray(self.faceBoundaries, dtype=np.int32).copy()
        names = self.surfaceDomainNames
        domainMap = DomainMap(names)
        for spec in _normalizeAssignments(assignments, assignment):
            mask = _selectFaces(self, faceBoundaries, spec)
            if not bool(spec.get("allowInternal", allowInternal)) and np.any(mask & (self.neighborCells >= 0)):
                raise ValueError("surface domain assignment selected internal faces; pass allowInternal=True to permit this")
            target = _targetDomain(domainMap, spec)
            if not np.any(mask):
                raise ValueError(f"surface domain assignment selected no faces: {spec}")
            faceBoundaries[mask] = np.int32(target)
            _updateDomainNames(names, domainMap, spec, target)
        metadata = {**self.metadata, "surfaceDomainNames": names}
        return self._copyWith(faceBoundaries=faceBoundaries, metadata=metadata)

    def _copyWith(self, *, cellDomains=None, faceBoundaries=None, metadata=None):
        # Domain reassignment does not alter connectivity or geometry.  Keep
        # the value semantics of the previous constructor-based copy without
        # rebuilding all Tet4 faces, neighbors, normals, and volumes.
        topology = copy.copy(self)
        topology.points = self.points.copy()
        topology.cellPointIndices = self.cellPointIndices.copy()
        topology.cellTypes = self.cellTypes.copy()
        topology.cellDomains = np.asarray(
            self.cellDomains if cellDomains is None else cellDomains,
            dtype=np.int32,
        ).copy()
        topology.faceBoundaries = np.asarray(
            self.faceBoundaries if faceBoundaries is None else faceBoundaries,
            dtype=np.int32,
        ).copy()
        topology.metadata = dict(self.metadata if metadata is None else metadata)
        for name in (
            "facePointIndices",
            "neighborCells",
            "neighborLocalFaces",
            "faceCenters",
            "faceNormals",
            "faceAreas",
            "cellCenters",
            "cellVolumes",
        ):
            setattr(topology, name, np.asarray(getattr(self, name)).copy())
        topology.samplePoints = np.asarray(self.samplePoints, dtype=np.float64).copy()
        return topology

    def openPmdAttributes(self, context=None):
        context = context or self
        return {
            "numberOfPoints": int(context.numberOfPoints),
            "numberOfTriangles": int(getattr(context, "numberOfCells", self.numberOfCells)),
            "numberOfLevels": int(getattr(context, "numberOfLevels", self.structuredNumberOfLevels)),
            "thickness": float(getattr(context, "thickness", self.structuredThickness)),
        }

    def cellsConnectivityFlat(self):
        return np.asarray(self.cellPointIndices, dtype=np.uint32).reshape(-1)

    def cellsOffsets(self):
        return np.arange(self.numberOfCells + 1, dtype=np.uint32) * np.uint32(4)

    def faceConnectivityFlat(self):
        return np.asarray(self.facePointIndices, dtype=np.int32).reshape(-1)


def _deriveTet4Topology(points, cells, faceBoundariesArg):
    numberOfCells = cells.shape[0]
    numberOfFaces = TET4_FACE_VERTICES.shape[0]
    faceNodes, cellCenters, cellVolumes, faceCenters, faceNormals, faceAreas = _tet4GeometryKernel(
        points,
        cells,
        TET4_FACE_VERTICES,
    )
    neighborCells = np.full((numberOfCells, numberOfFaces), -1, dtype=np.int32)
    neighborFaces = np.full((numberOfCells, numberOfFaces), -1, dtype=np.int32)
    boundaries = np.zeros((numberOfCells, numberOfFaces), dtype=np.int32) if faceBoundariesArg is None else np.asarray(faceBoundariesArg, dtype=np.int32).copy()
    if boundaries.shape != (numberOfCells, numberOfFaces):
        raise ValueError(f"faceBoundaries must have shape {(numberOfCells, numberOfFaces)}, got {boundaries.shape}")

    flatFaceKeys = np.sort(faceNodes.reshape(-1, 3), axis=1)
    order = np.lexsort((flatFaceKeys[:, 2], flatFaceKeys[:, 1], flatFaceKeys[:, 0]))
    sortedKeys = flatFaceKeys[order]
    matches = np.flatnonzero(np.all(sortedKeys[1:] == sortedKeys[:-1], axis=1))
    if matches.size:
        left = order[matches]
        right = order[matches + 1]
        flatNeighborCells = neighborCells.reshape(-1)
        flatNeighborFaces = neighborFaces.reshape(-1)
        flatNeighborCells[left] = right // numberOfFaces
        flatNeighborCells[right] = left // numberOfFaces
        flatNeighborFaces[left] = right % numberOfFaces
        flatNeighborFaces[right] = left % numberOfFaces

        flatBoundaries = boundaries.reshape(-1)
        leftBoundaries = flatBoundaries[left].copy()
        rightBoundaries = flatBoundaries[right].copy()
        takeRight = (leftBoundaries == BOUND_INTERNAL) & (rightBoundaries != BOUND_INTERNAL)
        takeLeft = (rightBoundaries == BOUND_INTERNAL) & (leftBoundaries != BOUND_INTERNAL)
        flatBoundaries[left[takeRight]] = rightBoundaries[takeRight]
        flatBoundaries[right[takeLeft]] = leftBoundaries[takeLeft]

    boundaryMask = neighborCells < 0
    boundaries[boundaryMask] = np.where(boundaries[boundaryMask] == BOUND_INTERNAL, BOUND_STOP, boundaries[boundaryMask])

    return {
        "facePointIndices": faceNodes,
        "neighborCells": neighborCells,
        "neighborLocalFaces": neighborFaces,
        "faceBoundaries": boundaries,
        "faceCenters": faceCenters,
        "faceNormals": faceNormals,
        "faceAreas": faceAreas,
        "cellCenters": cellCenters,
        "cellVolumes": cellVolumes,
    }


def _normalizeAssignments(assignments, assignment):
    if assignments is None:
        if not assignment:
            raise ValueError("at least one domain assignment is required")
        return [assignment]
    if assignment:
        raise ValueError("pass either assignments or keyword assignment arguments, not both")
    if isinstance(assignments, dict):
        return [assignments]
    return list(assignments)


def _targetDomain(domainMap, spec):
    if "domain" in spec:
        return domainMap.resolve(spec["domain"])
    if "gmshName" in spec:
        return domainMap.resolve(spec["gmshName"])
    if "gmshTag" in spec:
        return domainMap.resolve(spec["gmshTag"])
    raise ValueError("domain assignment requires domain, gmshName, or gmshTag")


def _updateDomainNames(names, domainMap, spec, target):
    name = spec.get("name")
    if name is None and isinstance(spec.get("domain"), str):
        name = spec["domain"]
    if name is None and "gmshName" in spec:
        name = spec["gmshName"]
    if name is None and "gmshTag" in spec:
        name = domainMap.names.get(int(spec["gmshTag"]))
    if name is not None:
        for domain, existingName in list(names.items()):
            if domain != int(target) and existingName == name:
                del names[domain]
        names[int(target)] = str(name)


def _selectCells(topology, cellDomains, spec):
    if "cellIndices" in spec:
        mask = np.zeros(topology.numberOfCells, dtype=bool)
        mask[np.asarray(spec["cellIndices"], dtype=np.int64)] = True
        return mask
    if spec.get("where") == "all":
        return np.ones(topology.numberOfCells, dtype=bool)
    if "gmshName" in spec:
        source = topology.cellDomainMap().resolve(spec["gmshName"])
        return cellDomains == source
    if "gmshTag" in spec:
        return cellDomains == int(spec["gmshTag"])
    raise ValueError("cell domain assignment requires where='all', cellIndices, gmshName, or gmshTag")


def _selectFaces(topology, faceBoundaries, spec):
    if "faceIndices" in spec:
        mask = np.zeros_like(faceBoundaries, dtype=bool)
        indices = np.asarray(spec["faceIndices"], dtype=np.int64)
        if indices.size == 0:
            return mask
        if indices.ndim != 2 or indices.shape[1] != 2:
            raise ValueError("faceIndices must have shape (N, 2)")
        mask[indices[:, 0], indices[:, 1]] = True
        return mask
    where = spec.get("where")
    if where in {"z_min", "z_max"}:
        z = np.asarray(topology.points, dtype=np.float64)[:, 2]
        target = np.min(z) if where == "z_min" else np.max(z)
        faceZ = z[np.asarray(topology.facePointIndices, dtype=np.uint32)]
        return (topology.neighborCells < 0) & np.all(np.isclose(faceZ, target), axis=2)
    if where == "all_exterior":
        return topology.neighborCells < 0
    if "gmshName" in spec:
        source = topology.surfaceDomainMap().resolve(spec["gmshName"])
        return faceBoundaries == source
    if "gmshTag" in spec:
        return faceBoundaries == int(spec["gmshTag"])
    raise ValueError("surface domain assignment requires where, faceIndices, gmshName, or gmshTag")


def _fromGmshVolume(gmsh, *, boundaryDefault):
    tetElements = [element for element in gmsh.elements if element.element_type == GMSH_TET4]
    unsupported = sorted({element.element_type for element in gmsh.elements if element.element_type in GMSH_UNSUPPORTED_VOLUME_TYPES})
    if unsupported:
        raise NotImplementedError("gmsh volume import supports Tet4 only")
    if not tetElements:
        raise ValueError("gmsh volume import requires at least one Tet4 element")

    surfaceElements = [element for element in gmsh.elements if element.element_type == GMSH_TRI3]
    usedNodeIds = sorted({nodeId for element in tetElements for nodeId in element.node_ids[:4]})
    nodeToIndex = {nodeId: index for index, nodeId in enumerate(usedNodeIds)}
    points = np.asarray([gmsh.nodes[nodeId] for nodeId in usedNodeIds], dtype=np.float64)

    cells = np.asarray([[nodeToIndex[nodeId] for nodeId in element.node_ids[:4]] for element in tetElements], dtype=np.uint32)
    domains = np.asarray([element.physical_tag or 1 for element in tetElements], dtype=np.int32)

    surfaceNodes = []
    surfaceTags = []
    for element in surfaceElements:
        if any(nodeId not in nodeToIndex for nodeId in element.node_ids[:3]):
            continue
        surfaceNodes.append([nodeToIndex[nodeId] for nodeId in element.node_ids[:3]])
        surfaceTags.append(
            element.physical_tag if element.physical_tag is not None else boundaryDefault
        )

    faceBoundariesArg = _surfaceBoundariesFromTriangles(cells, surfaceNodes, surfaceTags)

    return VolumeTopology(
        points,
        cells,
        cellDomains=domains,
        faceBoundaries=faceBoundariesArg,
        metadata={
            "source": gmsh.source,
            "format": "gmsh",
            "dimension": 3,
            "gmsh": gmsh,
            "cellDomainNames": {int(tag): str(name) for tag, name in gmsh.physical_names.get(3, {}).items()},
            "surfaceDomainNames": {int(tag): str(name) for tag, name in gmsh.physical_names.get(2, {}).items()},
        },
    )


def _fromStlVolume(filename, *, boundaryDefault, meshSize):
    try:
        import gmsh as gmshApi
    except (ImportError, OSError) as exc:
        raise ImportError("STL volume import requires a working 'gmsh' Python package for tetrahedralization") from exc

    warnings.warn(
        "STL volume import assumes a closed 3D surface suitable for Tet4 volume meshing; "
        "HASEonGPU does not run a full tetrahedral mesh validation pass.",
        RuntimeWarning,
        stacklevel=3,
    )

    from .msh import Gmsh, _read_elements as readElements, _read_nodes as readNodes, _read_physical_names as readPhysicalNames

    path = Path(filename)
    ownedSession = not gmshApi.isInitialized()
    if ownedSession:
        gmshApi.initialize()
    try:
        gmshApi.option.setNumber("General.Terminal", 0)
        gmshApi.clear()
        gmshApi.model.add("stlTet4Volume")
        gmshApi.merge(str(path))
        angle = 40.0 * math.pi / 180.0
        gmshApi.model.mesh.classifySurfaces(angle, True, True, math.pi)
        gmshApi.model.mesh.createGeometry()
        surfaces = [tag for dim, tag in gmshApi.model.getEntities(2)]
        if not surfaces:
            raise ValueError("STL tetrahedralization requires at least one closed surface")
        surfaceLoop = gmshApi.model.geo.addSurfaceLoop(surfaces)
        volume = gmshApi.model.geo.addVolume([surfaceLoop])
        gmshApi.model.geo.synchronize()
        gmshApi.model.addPhysicalGroup(3, [volume], 1)
        gmshApi.model.setPhysicalName(3, 1, "volume")
        gmshApi.model.addPhysicalGroup(2, surfaces, int(boundaryDefault) if int(boundaryDefault) > 0 else 1)
        if meshSize is not None:
            gmshApi.option.setNumber("Mesh.CharacteristicLengthMin", float(meshSize))
            gmshApi.option.setNumber("Mesh.CharacteristicLengthMax", float(meshSize))
        gmshApi.model.mesh.generate(3)
        gmsh = Gmsh(
            nodes=readNodes(gmshApi),
            elements=readElements(gmshApi),
            physical_names=readPhysicalNames(gmshApi),
            source=str(path),
        )
    finally:
        if ownedSession:
            gmshApi.finalize()
    return _fromGmshVolume(gmsh, boundaryDefault=boundaryDefault)
