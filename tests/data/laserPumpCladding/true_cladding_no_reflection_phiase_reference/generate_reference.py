#!/usr/bin/env python3
# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate the static true-cladding reference with the fixed legacy solver."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np


LEGACY_COMMIT = "effd8077edccef93a68d818e8a5eb2f0ebdc03b4"
CLADDING_NUMBER = 1
PHYSICAL_CLADDING_ABSORPTION = 5.5
INSTRUMENTATION_CLADDING_ABSORPTION = 55.0
CLADDING_ABSORPTIONS = (
    0.0,
    PHYSICAL_CLADDING_ABSORPTION,
    INSTRUMENTATION_CLADDING_ABSORPTION,
)
BETA_GAIN = 0.1
WAVELENGTH = 1030.0e-9
SIGMA_ABSORPTION = 0.5e-20
SIGMA_EMISSION = 1.0e-20
RNG_SEED = 5489


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _legacy_commit(legacy_repo):
    return subprocess.run(
        ["git", "-C", str(legacy_repo), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _load_legacy_example(legacy_repo, legacy_build):
    sys.path[:0] = [str(legacy_build / "python"), str(legacy_repo)]
    import HASEonGPU  # noqa: PLC0415
    import pyInclude  # noqa: PLC0415

    bindings = importlib.import_module("HASEonGPU_Bindings.HASEonGPU")
    for name, module in (
        ("HASEonGPU", HASEonGPU),
        ("pyInclude", pyInclude),
        ("HASEonGPU_Bindings", bindings),
    ):
        module_path = Path(module.__file__).resolve()
        if legacy_repo.resolve() not in module_path.parents:
            raise RuntimeError(f"{name} resolved outside the legacy worktree: {module_path}")

    example_path = legacy_repo / "example" / "python_example" / "laserPumpCladding.py"
    spec = importlib.util.spec_from_file_location("true_cladding_legacy_example", example_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load legacy example from {example_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, HASEonGPU, bindings


def _outer_triangle_mask(triangles):
    triangles = np.asarray(triangles, dtype=np.uint32)
    edges = np.sort(
        np.concatenate(
            [triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]],
            axis=0,
        ),
        axis=1,
    )
    unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
    exterior_edges = {
        tuple(int(value) for value in edge) for edge in unique_edges[counts == 1]
    }
    mask = np.asarray(
        [
            any(
                tuple(sorted((int(triangle[index]), int(triangle[(index + 1) % 3]))))
                in exterior_edges
                for index in range(3)
            )
            for triangle in triangles
        ],
        dtype=bool,
    )
    return mask, len(exterior_edges)


def _extruded_geometry(topology):
    xy = np.asarray(topology.points, dtype=np.float64)
    triangles = np.asarray(topology.trianglePointIndices, dtype=np.uint32)
    number_of_points = topology.numberOfPoints
    levels = int(topology.levels)
    points = np.vstack(
        [
            np.column_stack(
                [xy, np.full(number_of_points, level * topology.thickness, dtype=np.float64)]
            )
            for level in range(levels)
        ]
    )
    cells = np.vstack(
        [
            np.column_stack(
                [
                    triangles + level * number_of_points,
                    triangles + (level + 1) * number_of_points,
                ]
            )
            for level in range(levels - 1)
        ]
    ).astype(np.uint32)
    return points, cells


def _wedge_volumes(points, cells):
    lower = np.asarray(points, dtype=np.float64)[np.asarray(cells, dtype=np.uint32)[:, :3]]
    ab = lower[:, 1, :2] - lower[:, 0, :2]
    ac = lower[:, 2, :2] - lower[:, 0, :2]
    areas = 0.5 * np.abs(ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0])
    vertices = np.asarray(points, dtype=np.float64)[np.asarray(cells, dtype=np.uint32)]
    heights = vertices[:, :, 2].max(axis=1) - vertices[:, :, 2].min(axis=1)
    return areas * heights


def _integrals(points, cells, phi, wedge_cladding_types):
    volumes = _wedge_volumes(points, cells)
    cell_values = np.asarray(phi, dtype=np.float64)[np.asarray(cells, dtype=np.uint32)].mean(axis=1)
    cladding = np.asarray(wedge_cladding_types, dtype=np.uint32) == CLADDING_NUMBER
    weighted = cell_values * volumes
    return {
        "total": float(weighted.sum()),
        "gain": float(weighted[~cladding].sum()),
        "cladding": float(weighted[cladding].sum()),
    }


def generate(args):
    legacy_repo = args.legacy_repo.resolve()
    legacy_build = args.legacy_build.resolve()
    commit = _legacy_commit(legacy_repo)
    if commit != LEGACY_COMMIT:
        raise RuntimeError(f"legacy worktree is at {commit}, expected {LEGACY_COMMIT}")

    example, hase, bindings = _load_legacy_example(legacy_repo, legacy_build)
    medium = example.laserPumpCladdingMedium(cladAbsorption=CLADDING_ABSORPTIONS[0])
    topology = medium.topology
    base_cladding, exterior_edge_count = _outer_triangle_mask(topology.trianglePointIndices)
    if exterior_edge_count != 24 or np.count_nonzero(base_cladding) != 24:
        raise RuntimeError(
            "the audited pt.mat shell must contain 24 exterior edges "
            "and 24 shell triangles"
        )

    wedge_cladding = np.tile(base_cladding, int(topology.levels) - 1).astype(np.uint32)
    beta_volume = np.full(medium.get("betaVolume").expectedShape, BETA_GAIN, dtype=np.float64)
    beta_volume[base_cladding, :] = 0.0
    if np.any(beta_volume[base_cladding, :] != 0.0):
        raise RuntimeError("cladding prisms must carry exactly zero source beta")

    # betaCells is intentionally uniform. Direct PhiASE transport and source
    # selection use betaVolume; cladding segments override gain with the
    # cladding absorption coefficient. This avoids zeroing gain prisms that
    # share interface vertices with cladding prisms.
    medium.get("betaCells").value = np.full(
        medium.get("betaCells").expectedShape,
        BETA_GAIN,
        dtype=np.float64,
    )
    medium.get("betaVolume").value = beta_volume
    medium.get("claddingCellTypes").value = base_cladding.astype(np.uint32)
    medium.get("claddingNumber").value = CLADDING_NUMBER

    cross_sections = hase.CrossSectionData.monochromatic(
        wavelength=WAVELENGTH,
        crossSectionAbsorption=SIGMA_ABSORPTION,
        crossSectionEmission=SIGMA_EMISSION,
    )
    points, cells = _extruded_geometry(topology)
    volumes = _wedge_volumes(points, cells)
    cladding_volume = float(volumes[wedge_cladding == CLADDING_NUMBER].sum())

    phi_by_absorption = []
    mse_by_absorption = []
    rays_by_absorption = []
    diagnostics = []
    for cladding_absorption in CLADDING_ABSORPTIONS:
        medium.get("claddingAbsorption").value = cladding_absorption
        stored_beta_volume = np.asarray(medium.get("betaVolume").value, dtype=np.float64).reshape(
            medium.get("betaVolume").expectedShape,
            order="F",
        )
        if np.any(stored_beta_volume[base_cladding, :] != 0.0):
            raise RuntimeError("cladding source beta changed before the PhiASE call")
        phi_ase = hase.PhiASE(
            crossSections=cross_sections,
            minRaysPerSample=args.rays_per_sample,
            maxRaysPerSample=args.rays_per_sample,
            mseThreshold=1.0,
            repetitions=1,
            adaptiveSteps=1,
            useReflections=False,
            monochromatic=True,
            backend=args.backend,
            rngSeed=RNG_SEED,
        )
        phi_ase.run(gainMedium=medium, crossSections=cross_sections)
        result = phi_ase.getResults()
        phi = np.asarray(result.phiAse, dtype=np.float64)
        mse = np.asarray(result.mse, dtype=np.float64)
        total_rays = np.asarray(result.totalRays, dtype=np.uint32)
        phi_by_absorption.append(phi)
        mse_by_absorption.append(mse)
        rays_by_absorption.append(total_rays)
        diagnostics.append(
            {
                "claddingAbsorption": cladding_absorption,
                "integrals": _integrals(points, cells, phi, wedge_cladding),
                "meanPhi": float(phi.mean()),
                "meanMse": float(mse.mean()),
                "maxMse": float(mse.max()),
                "totalRays": int(total_rays.astype(np.uint64).sum()),
                "droppedRays": int(np.asarray(result.droppedRays, dtype=np.uint64).sum()),
            }
        )

    material_path = legacy_repo / "example" / "python_example" / "legacy" / "pt.mat"
    metadata = {
        "formatVersion": 1,
        "generator": {
            "commit": commit,
            "source": "example/python_example/laserPumpCladding.py",
            "script": Path(__file__).name,
            "legacyFrontend": str(Path(hase.__file__).resolve()),
            "legacyBindings": str(Path(bindings.__file__).resolve()),
            "ptMatSha256": _sha256(material_path),
        },
        "geometry": {
            "definition": (
                "base triangles possessing an exterior boundary edge, "
                "extruded through every wedge layer"
            ),
            "numberOfBasePoints": int(topology.numberOfPoints),
            "numberOfBaseTriangles": int(topology.numberOfTriangles),
            "numberOfLevels": int(topology.levels),
            "numberOfWedges": int(cells.shape[0]),
            "exteriorEdgeCount": exterior_edge_count,
            "claddingBaseTriangles": int(np.count_nonzero(base_cladding)),
            "claddingWedges": int(np.count_nonzero(wedge_cladding)),
            "totalVolume": float(volumes.sum()),
            "gainVolume": float(volumes.sum() - cladding_volume),
            "claddingVolume": cladding_volume,
            "claddingVolumeFraction": cladding_volume / float(volumes.sum()),
        },
        "material": {
            "claddingNumber": CLADDING_NUMBER,
            "claddingAbsorptions": list(CLADDING_ABSORPTIONS),
            "physicalCladdingAbsorption": PHYSICAL_CLADDING_ABSORPTION,
            "instrumentationCladdingAbsorption": INSTRUMENTATION_CLADDING_ABSORPTION,
            "gainBetaVolume": BETA_GAIN,
            "claddingBetaVolume": 0.0,
            "betaCells": BETA_GAIN,
            "sharedInterfaceVertexPolicy": (
                "betaCells remain uniform; authoritative per-volume betaVolume "
                "is zero only in cladding"
            ),
            "wavelength": WAVELENGTH,
            "crossSectionAbsorption": SIGMA_ABSORPTION,
            "crossSectionEmission": SIGMA_EMISSION,
        },
        "parameters": {
            "backend": args.backend,
            "raysPerSample": args.rays_per_sample,
            "mseThreshold": 1.0,
            "rngSeed": RNG_SEED,
            "useReflections": False,
            "monochromatic": True,
        },
        "observable": {
            "field": "phiASE",
            "location": "legacy wedge POINT_DATA",
            "shape": [len(CLADDING_ABSORPTIONS), int(points.shape[0])],
            "projection": "arithmetic mean of six wedge vertices times wedge volume",
            "partitions": ["total", "gain", "cladding"],
        },
        "diagnostics": diagnostics,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        metadata=np.asarray(json.dumps(metadata, sort_keys=True)),
        phiASE=np.asarray(phi_by_absorption, dtype=np.float64),
        mse=np.asarray(mse_by_absorption, dtype=np.float64),
        totalRays=np.asarray(rays_by_absorption, dtype=np.uint32),
        points=points,
        cells=cells,
        cellTypes=np.full(cells.shape[0], 13, dtype=np.uint32),
        baseCladdingCellTypes=base_cladding.astype(np.uint32),
        wedgeCladdingCellTypes=wedge_cladding,
        wedgeBetaVolume=beta_volume.reshape(-1, order="F"),
    )
    print(json.dumps(metadata, indent=2, sort_keys=True))
    print(f"wrote {args.output}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-repo", type=Path, required=True)
    parser.add_argument("--legacy-build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rays-per-sample", type=int, default=2000)
    parser.add_argument("--backend", default="Host_Cpu_CpuOmpBlocks")
    args = parser.parse_args(argv)
    if args.rays_per_sample <= 0:
        parser.error("--rays-per-sample must be positive")
    generate(args)


if __name__ == "__main__":
    main()
