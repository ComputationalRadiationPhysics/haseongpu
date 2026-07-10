# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Explicit Tet4 volume-cell topology helpers."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import math
import warnings

import numpy as np

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
        return type(self)(
            self.points.copy(),
            self.cellPointIndices.copy(),
            cellTypes=self.cellTypes.copy(),
            cellDomains=np.asarray(self.cellDomains if cellDomains is None else cellDomains, dtype=np.int32).copy(),
            faceBoundaries=np.asarray(self.faceBoundaries if faceBoundaries is None else faceBoundaries, dtype=np.int32).copy(),
            metadata=dict(self.metadata if metadata is None else metadata),
        )

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
    faceNodes = np.full((numberOfCells, numberOfFaces, 3), -1, dtype=np.int32)
    neighborCells = np.full((numberOfCells, numberOfFaces), -1, dtype=np.int32)
    neighborFaces = np.full((numberOfCells, numberOfFaces), -1, dtype=np.int32)
    boundaries = np.zeros((numberOfCells, numberOfFaces), dtype=np.int32) if faceBoundariesArg is None else np.asarray(faceBoundariesArg, dtype=np.int32).copy()
    if boundaries.shape != (numberOfCells, numberOfFaces):
        raise ValueError(f"faceBoundaries must have shape {(numberOfCells, numberOfFaces)}, got {boundaries.shape}")

    owners = {}
    for cellIndex in range(numberOfCells):
        for faceIndex in range(numberOfFaces):
            nodes = _faceNodes(cells[cellIndex], faceIndex)
            faceNodes[cellIndex, faceIndex] = nodes
            key = _faceKey(nodes)
            owner = owners.get(key)
            if owner is None:
                owners[key] = (cellIndex, faceIndex)
            else:
                otherCell, otherFace = owner
                neighborCells[cellIndex, faceIndex] = otherCell
                neighborFaces[cellIndex, faceIndex] = otherFace
                neighborCells[otherCell, otherFace] = cellIndex
                neighborFaces[otherCell, otherFace] = faceIndex
                if boundaries[cellIndex, faceIndex] == BOUND_INTERNAL and boundaries[otherCell, otherFace] != BOUND_INTERNAL:
                    boundaries[cellIndex, faceIndex] = boundaries[otherCell, otherFace]
                elif boundaries[otherCell, otherFace] == BOUND_INTERNAL and boundaries[cellIndex, faceIndex] != BOUND_INTERNAL:
                    boundaries[otherCell, otherFace] = boundaries[cellIndex, faceIndex]

    boundaryMask = neighborCells < 0
    boundaries[boundaryMask] = np.where(boundaries[boundaryMask] == BOUND_INTERNAL, BOUND_STOP, boundaries[boundaryMask])

    cellCenters = np.empty((numberOfCells, 3), dtype=np.float64)
    cellVolumes = np.empty(numberOfCells, dtype=np.float32)
    faceCenters = np.empty((numberOfCells, numberOfFaces, 3), dtype=np.float32)
    faceNormals = np.empty((numberOfCells, numberOfFaces, 3), dtype=np.float32)
    faceAreas = np.empty((numberOfCells, numberOfFaces), dtype=np.float32)

    for cellIndex in range(numberOfCells):
        nodes = cells[cellIndex]
        coordinates = points[nodes]
        cellCenters[cellIndex] = np.mean(coordinates, axis=0)
        cellVolumes[cellIndex] = np.float32(_tetVolume(coordinates[0], coordinates[1], coordinates[2], coordinates[3]))
        for faceIndex in range(numberOfFaces):
            center, area, normal = _faceGeometry(points, faceNodes[cellIndex, faceIndex], cellCenters[cellIndex])
            faceCenters[cellIndex, faceIndex] = center.astype(np.float32)
            faceAreas[cellIndex, faceIndex] = area
            faceNormals[cellIndex, faceIndex] = normal

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
        for cell, face in spec["faceIndices"]:
            mask[int(cell), int(face)] = True
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

    surfaceBoundariesArg = {}
    for element in surfaceElements:
        if any(nodeId not in nodeToIndex for nodeId in element.node_ids[:3]):
            continue
        nodes = tuple(nodeToIndex[nodeId] for nodeId in element.node_ids[:3])
        surfaceBoundariesArg[_faceKey(nodes)] = np.int32(element.physical_tag if element.physical_tag is not None else boundaryDefault)

    faceBoundariesArg = np.zeros((cells.shape[0], TET4_FACE_VERTICES.shape[0]), dtype=np.int32)
    for cellIndex, cell in enumerate(cells):
        for faceIndex in range(TET4_FACE_VERTICES.shape[0]):
            key = _faceKey(_faceNodes(cell, faceIndex))
            if key in surfaceBoundariesArg:
                faceBoundariesArg[cellIndex, faceIndex] = surfaceBoundariesArg[key]

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
