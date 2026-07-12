#!/usr/bin/env python3
# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Project converted Tet4 compareSerial fixtures back to legacy wedge cells."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import os
from pathlib import Path
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from HASEonGPU import GainMedium, MeshTopology, vtkWedge  # noqa: E402


DATASETS = ("cuboid", "cylindrical")
REFERENCE_PATH = REPO_ROOT / "tests" / "data" / "compareSerial" / "phiase_reference.npz"


@dataclass(frozen=True)
class WedgeProjection:
    topology: MeshTopology
    prismVolumes: np.ndarray
    originalBetaVolume: np.ndarray
    projectedBetaVolume: np.ndarray
    betaCells: np.ndarray


def loadCompareSerialReference(path: Path = REFERENCE_PATH):
    path = Path(path)
    if not path.is_file():
        return {}
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata"].item()))
        references = {}
        for name in metadata["datasets"]:
            references[name] = {
                "name": name,
                "phiASE": np.asarray(data[f"{name}_phiASE"], dtype=np.float64),
                "dndtAse": np.asarray(data[f"{name}_dndtAse"], dtype=np.float64),
                "metadata": metadata,
            }
        return references


def serialPhiAsePointFields(datasetReference, topology: MeshTopology):
    """Return full point fields for recorded compareSerial phiASE samples."""
    values = np.asarray(datasetReference.get("phiASE", []), dtype=np.float64).reshape(-1)
    point_count = int(topology.numberOfPoints) * int(topology.levels)
    if values.size != point_count:
        raise ValueError("compareSerial phiASE reference must cover every wedge point")
    return {"serialPhiASE": values}


def projectTet4CellValuesToLegacyPrisms(values, tetVolumes):
    """Volume-project values from three Tet4 cells back to one legacy prism."""
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    tetVolumes = np.asarray(tetVolumes, dtype=np.float64).reshape(-1)
    if values.shape != tetVolumes.shape:
        raise ValueError("values and tetVolumes must have the same flat shape")
    if values.size % 3 != 0:
        raise ValueError("converted legacy wedge fixtures must have three Tet4 cells per prism")
    prism_values = values.reshape((-1, 3))
    constant_prisms = np.logical_and(
        prism_values[:, 0] == prism_values[:, 1],
        prism_values[:, 0] == prism_values[:, 2],
    )
    prism_volumes = tetVolumes.reshape((-1, 3))
    projected = np.sum(prism_values * prism_volumes, axis=1) / np.sum(prism_volumes, axis=1)
    projected[constant_prisms] = prism_values[constant_prisms, 0]
    return projected


def tet4ResultToLegacyPointValues(medium: GainMedium, values):
    """Map current Tet4 PhiASE output to the legacy wedge point buffer order.

    Point-shaped results already use the legacy structured point order and are
    returned directly. Cell-shaped results are volume-averaged to incident
    vertices; vertices with no incident Tet4 cell are left as NaN because no
    cell-centered interpolation can define a pointwise value there.
    """
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    projection = projectionFromTet4Medium(medium)
    point_count = int(projection.topology.numberOfPoints) * int(projection.topology.levels)
    if values.size == point_count:
        return values
    if values.size != medium.topology.numberOfCells:
        raise ValueError(
            "Tet4 result must contain either one value per legacy point "
            f"({point_count}) or one value per Tet4 cell ({medium.topology.numberOfCells}), got {values.size}"
        )

    point_values = np.zeros(point_count, dtype=np.float64)
    point_weights = np.zeros(point_count, dtype=np.float64)
    cell_volumes = np.asarray(medium.topology.cellVolumes, dtype=np.float64).reshape(-1)
    cells = np.asarray(medium.topology.cellPointIndices, dtype=np.uint32)
    for cell_index, vertices in enumerate(cells):
        weight = cell_volumes[cell_index]
        for vertex in vertices:
            point_values[int(vertex)] += values[cell_index] * weight
            point_weights[int(vertex)] += weight
    result = np.full(point_count, np.nan, dtype=np.float64)
    valid = point_weights > 0.0
    result[valid] = point_values[valid] / point_weights[valid]
    return result


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


def _structured_metadata(topology):
    levels = int(getattr(topology, "structuredNumberOfLevels", 1))
    points_per_level = int(getattr(topology, "structuredNumberOfPoints", topology.numberOfSamplePoints))
    thickness = float(getattr(topology, "structuredThickness", 0.0))
    if levels < 2:
        raise ValueError("legacy wedge reconstruction requires at least two z levels")
    if topology.numberOfPoints != points_per_level * levels:
        raise ValueError("Tet4 points are not a structured legacy point stack")
    if topology.numberOfCells % (levels - 1) != 0:
        raise ValueError("Tet4 cell count is not divisible by legacy z intervals")
    return points_per_level, levels, thickness


def _triangle_area(points):
    a, b, c = np.asarray(points, dtype=np.float64)
    ab = b - a
    ac = c - a
    return 0.5 * abs(float(ab[0] * ac[1] - ab[1] * ac[0]))


def _extruded_prism_volume(points, vertex_indices):
    vertices = np.asarray(points, dtype=np.float64)[np.asarray(vertex_indices, dtype=np.uint32)]
    z_values = np.unique(vertices[:, 2])
    if z_values.size != 2:
        raise ValueError("legacy prism inverse must contain exactly two z levels")
    lower = vertices[np.isclose(vertices[:, 2], z_values[0])]
    upper = vertices[np.isclose(vertices[:, 2], z_values[1])]
    if lower.shape[0] != 3 or upper.shape[0] != 3:
        raise ValueError("legacy prism inverse must contain three vertices on each z level")
    return _triangle_area(lower[:, :2]) * float(z_values[1] - z_values[0])


def legacyPrismVertexSets(tetCells):
    """Return one six-vertex set for every legacy prism encoded as three Tet4 cells."""
    tetCells = np.asarray(tetCells, dtype=np.uint32)
    if tetCells.ndim != 2 or tetCells.shape[1] != 4:
        raise ValueError("Tet4 connectivity must have shape (N, 4)")
    if tetCells.shape[0] % 3 != 0:
        raise ValueError("converted legacy wedge fixtures must have three Tet4 cells per prism")
    prism_vertices = []
    for prism_tets in tetCells.reshape((-1, 3, 4)):
        vertices = np.unique(prism_tets.reshape(-1))
        if vertices.size != 6:
            raise ValueError("three Tet4 cells did not reconstruct one six-vertex prism")
        prism_vertices.append(vertices.astype(np.uint32))
    return prism_vertices


def _legacy_triangles_from_tet4(topology):
    points_per_level, levels, thickness = _structured_metadata(topology)
    points = np.asarray(topology.points, dtype=np.float64)
    cells = np.asarray(topology.cellPointIndices, dtype=np.uint32)
    z_min = float(points[:, 2].min())
    z_max = float(points[:, 2].max())

    triangles_per_layer = topology.numberOfCells // (3 * (levels - 1))
    triangles = []
    prism_vertices = []
    for prism_tets in cells.reshape((-1, 3, 4)):
        flattened = prism_tets.reshape(-1)
        lower = [idx for idx in _unique_in_order(flattened) if np.isclose(points[idx, 2], z_min)]
        if len(lower) != 3:
            # For upper layers, use the lower z of that prism, not the global z_min.
            local_z_min = float(points[np.unique(flattened), 2].min())
            lower = [idx for idx in _unique_in_order(flattened) if np.isclose(points[idx, 2], local_z_min)]
        if len(lower) != 3:
            raise ValueError("could not reconstruct one lower triangle from three Tet4 cells")
        triangles.append(np.asarray([idx % points_per_level for idx in lower], dtype=np.uint32))
        prism_vertices.append(np.unique(flattened).astype(np.uint32))

    base_triangles = np.asarray(triangles[:triangles_per_layer], dtype=np.uint32)
    for level in range(1, levels - 1):
        start = level * triangles_per_layer
        stop = start + triangles_per_layer
        if not np.array_equal(base_triangles, np.asarray(triangles[start:stop], dtype=np.uint32)):
            raise ValueError("legacy triangle ordering is not stable across z levels")

    base_points = points[:points_per_level, :2]
    if not np.allclose(points[:points_per_level, 2], z_min):
        raise ValueError("first structured point block is not the first z level")
    if not np.isclose(points[-points_per_level:, 2].max(), z_max):
        raise ValueError("last structured point block is not the final z level")

    return base_points, base_triangles, levels, thickness, prism_vertices


def projectionFromTet4Medium(medium: GainMedium) -> WedgeProjection:
    topology = medium.topology
    base_points, base_triangles, levels, thickness, prism_vertices = _legacy_triangles_from_tet4(topology)
    wedge_topology = MeshTopology(
        points=base_points,
        trianglePointIndices=base_triangles,
        levels=levels,
        thickness=thickness,
        metadata={"source": topology.metadata.get("source", "tet4"), "format": "tet4-inverted-wedge"},
    )

    tet_volumes = np.asarray(topology.cellVolumes, dtype=np.float64).reshape((-1, 3))
    prism_volumes = np.asarray(
        [_extruded_prism_volume(topology.points, vertices) for vertices in prism_vertices],
        dtype=np.float64,
    )
    projected_beta = projectTet4CellValuesToLegacyPrisms(medium.get("betaVolume").value, topology.cellVolumes)
    original_beta = np.asarray(medium.get("betaVolume").value, dtype=np.float64).reshape((-1, 3))[:, 0]
    return WedgeProjection(
        topology=wedge_topology,
        prismVolumes=prism_volumes,
        originalBetaVolume=original_beta,
        projectedBetaVolume=projected_beta,
        betaCells=np.asarray(medium.get("betaCells").value, dtype=np.float64).reshape(-1),
    )


def writeWedgeComparisonArtifacts(
    medium: GainMedium,
    outputDir: Path,
    name: str,
    *,
    tetCellFields=None,
    serialReference=None,
):
    """Write original and Tet4-projected wedge VTKs for one converted fixture."""
    outputDir = Path(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)
    projection = projectionFromTet4Medium(medium)

    original_fields = {
        "betaCells": projection.betaCells,
        "betaVolume": projection.originalBetaVolume,
    }
    roundtrip_fields = {
        "betaCells": projection.betaCells,
        "betaVolume": projection.projectedBetaVolume,
    }
    if serialReference is not None:
        original_fields.update(serialPhiAsePointFields(serialReference, projection.topology))
    for field_name, values in (tetCellFields or {}).items():
        roundtrip_fields[field_name] = projectTet4CellValuesToLegacyPrisms(values, medium.topology.cellVolumes)

    original_path = vtkWedge(
        outputDir / f"{name}_original_wedge.vtk",
        data={},
        geometry=projection.topology,
        fields=original_fields,
    )
    roundtrip_path = vtkWedge(
        outputDir / f"{name}_tet4_roundtrip_wedge.vtk",
        data={},
        geometry=projection.topology,
        fields=roundtrip_fields,
    )
    return {"original": original_path, "roundtrip": roundtrip_path, "projection": projection}


def generateArtifacts(
    repoRoot: Path = REPO_ROOT,
    outputDir: Path | None = None,
    datasets=DATASETS,
    referencePath: Path | None = REFERENCE_PATH,
):
    outputDir = Path(os.environ.get("HASE_COMPARE_SERIAL_ARTIFACT_DIR", "compareSerial_wedge_artifacts")) if outputDir is None else Path(outputDir)
    references = loadCompareSerialReference(referencePath) if referencePath is not None else {}
    written = {}
    for name in datasets:
        if name not in DATASETS:
            raise ValueError(f"unknown compareSerial dataset '{name}'")
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / f"{name}.vtk")
        written[name] = writeWedgeComparisonArtifacts(
            medium,
            outputDir,
            name,
            serialReference=references.get(name),
        )
    return written


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--dataset", action="append", choices=DATASETS, help="Dataset to write; defaults to both.")
    parser.add_argument("--reference-path", type=Path, default=REFERENCE_PATH)
    args = parser.parse_args(argv)

    artifacts = generateArtifacts(
        outputDir=args.output_dir,
        datasets=tuple(args.dataset or DATASETS),
        referencePath=args.reference_path,
    )
    for dataset, paths in artifacts.items():
        print(f"{dataset}: {paths['original']}")
        print(f"{dataset}: {paths['roundtrip']}")


if __name__ == "__main__":
    main()
