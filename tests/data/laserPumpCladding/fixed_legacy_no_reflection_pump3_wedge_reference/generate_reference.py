# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Pack deterministic legacy wedge VTK snapshots into the regression archive."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np


LEGACY_COMMIT = "effd8077edccef93a68d818e8a5eb2f0ebdc03b4"
IMAGE_ID = "sha256:8cd6743456e7ca8844f28c84e9aedae35429322eb21ac855bb01a2fdf15407a6"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class _TokenReader:
    def __init__(self, path: Path):
        self.tokens = path.read_text(encoding="utf-8").split()
        self.index = 0

    def done(self) -> bool:
        return self.index >= len(self.tokens)

    def next(self) -> str:
        value = self.tokens[self.index]
        self.index += 1
        return value

    def take(self, count: int, dtype=float) -> np.ndarray:
        result = np.asarray(self.tokens[self.index : self.index + count], dtype=dtype)
        self.index += count
        return result


def _parse_vtk(path: Path):
    reader = _TokenReader(path)
    points = None
    cells = []
    cell_types = None
    point_data = {}
    active_data = None
    active_count = 0
    while not reader.done():
        token = reader.next().upper()
        if token == "POINTS":
            count = int(reader.next())
            reader.next()
            points = reader.take(count * 3, np.float64).reshape(count, 3)
        elif token == "CELLS":
            count = int(reader.next())
            reader.next()
            cells = [reader.take(int(reader.next()), np.uint32) for _ in range(count)]
        elif token == "CELL_TYPES":
            cell_types = reader.take(int(reader.next()), np.uint32)
        elif token == "POINT_DATA":
            active_data = point_data
            active_count = int(reader.next())
        elif token == "CELL_DATA":
            active_data = {}
            active_count = int(reader.next())
        elif token == "SCALARS":
            if active_data is None:
                raise ValueError(f"SCALARS without POINT_DATA/CELL_DATA in {path}")
            name = reader.next()
            reader.next()
            components = 1
            if reader.tokens[reader.index].upper() != "LOOKUP_TABLE":
                components = int(reader.next())
            if reader.tokens[reader.index].upper() == "LOOKUP_TABLE":
                reader.next()
                reader.next()
            values = reader.take(active_count * components, np.float64)
            active_data[name] = values.reshape(active_count, components).squeeze()
    if points is None or not cells or cell_types is None:
        raise ValueError(f"incomplete VTK geometry in {path}")
    return points, np.asarray(cells, dtype=np.uint32), cell_types, point_data


def _wedge_volumes(points: np.ndarray, cells: np.ndarray) -> np.ndarray:
    volumes = []
    for cell in cells:
        vertices = points[cell]
        z_values = np.unique(vertices[:, 2])
        lower = vertices[np.isclose(vertices[:, 2], z_values[0])]
        a, b, c = lower
        ab = b[:2] - a[:2]
        ac = c[:2] - a[:2]
        area = 0.5 * abs(float(ab[0] * ac[1] - ab[1] * ac[0]))
        volumes.append(area * float(z_values[1] - z_values[0]))
    return np.asarray(volumes, dtype=np.float64)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vtk-dir", type=Path, required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("phiase_reference.npz"))
    args = parser.parse_args()

    fixture_dir = Path(__file__).resolve().parent
    config_path = fixture_dir / "phiASE-no-reflections.yaml"
    legacy_script = args.legacy_root / "example/python_example/laserPumpCladding.py"
    geometry_source = args.legacy_root / "example/python_example/legacy/pt.mat"
    vtk_paths = [args.vtk_dir / f"laserPumpCladding_{step:03d}.vtk" for step in range(1, 7)]

    actual_commit = subprocess.check_output(
        ["git", "-C", str(args.legacy_root), "rev-parse", "HEAD"],
        text=True,
    ).strip()
    if actual_commit != LEGACY_COMMIT:
        raise RuntimeError(
            f"legacy source is at {actual_commit}, expected fixture commit {LEGACY_COMMIT}"
        )

    snapshots = [_parse_vtk(path) for path in vtk_paths]
    points, cells, cell_types, _ = snapshots[0]
    for step_points, step_cells, step_types, _ in snapshots[1:]:
        np.testing.assert_array_equal(step_points, points)
        np.testing.assert_array_equal(step_cells, cells)
        np.testing.assert_array_equal(step_types, cell_types)
    phi_ase = np.asarray([snapshot[3]["phiASE"] for snapshot in snapshots], dtype=np.float64)
    volumes = _wedge_volumes(points, cells)
    integrals = np.asarray(
        [np.sum(phi[cells].mean(axis=1) * volumes) for phi in phi_ase],
        dtype=np.float64,
    )
    config_text = config_path.read_text(encoding="utf-8")

    metadata = {
        "schema": "hase.laserPumpCladding.fixedLegacyNoReflectionPump3WedgePhiASE.v1",
        "generator": {
            "repository": "detached legacy worktree",
            "commit": LEGACY_COMMIT,
            "commitSubject": "Fix legacy reflection plane height",
            "script": "example/python_example/laserPumpCladding.py",
            "scriptSha256": _sha256(legacy_script),
            "geometrySource": "example/python_example/legacy/pt.mat",
            "geometrySourceSha256": _sha256(geometry_source),
            "config": "phiASE-no-reflections.yaml",
            "configSha256": _sha256(config_path),
            "configText": config_text,
        },
        "docker": {
            "baseImage": "ubuntu:24.04",
            "baseImageId": "sha256:786a8b558f7be160c6c8c4a54f9a57274f3b4fb1491cf65146521ae77ff1dc54",
            "image": "hase-no-reflection-toolchain:ubuntu24.04",
            "imageId": IMAGE_ID,
            "dockerfileSha256": _sha256(fixture_dir / "Dockerfile"),
            "compiler": "gcc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0",
            "cmake": "3.28.3",
            "python": "3.12.3",
        },
        "command": {
            "type": "Docker-mounted Python runExample import",
            "pythonpath": ["/build/python", "/src"],
            "haseModule": "/src/HASEonGPU.py",
            "bindingsPackage": "/build/python/HASEonGPU_Bindings/__init__.py",
            "nativeBindings": "/build/python/HASEonGPU_Bindings/HASEonGPU.cpython-312-x86_64-linux-gnu.so",
            "call": "laserPumpCladding.runExample(phiAseConfigPath=Path('/generation/phiASE-no-reflections.yaml'), backend='Host_Cpu_CpuOmpBlocks', timeSlices=6, pumpSteps=3, vtkOutputDir=Path('/output'), enableAse=True, prePump=True, rngSeed=5489)",
        },
        "parameters": {
            "backend": "Host_Cpu_CpuOmpBlocks",
            "timeSlices": 6,
            "pumpSteps": 3,
            "prePump": True,
            "enableAse": True,
            "rngSeed": 5489,
            "minRaysPerSample": 2000,
            "maxRaysPerSample": 1000000,
            "mseThreshold": 0.01087,
            "adaptiveSteps": 8,
            "useReflections": False,
            "numberOfLevels": 10,
            "thickness": 0.7 / 9.0,
            "timeStepSeconds": 2.0e-5,
            "claddingAbsorption": 5.5,
        },
        "random": {
            "rngSeed": 5489,
            "note": "The legacy CLI had no seed option; runExample received rngSeed explicitly.",
        },
        "observable": {
            "field": "phiASE",
            "location": "legacy wedge POINT_DATA",
            "coverage": "full point buffer",
            "shape": list(phi_ase.shape),
            "pointCount": int(points.shape[0]),
            "wedgeCellCount": int(cells.shape[0]),
            "cellType": 13,
            "stepNumbers": list(range(1, 7)),
        },
        "meanPhi": phi_ase.mean(axis=1).tolist(),
        "minPhi": phi_ase.min(axis=1).tolist(),
        "maxPhi": phi_ase.max(axis=1).tolist(),
        "integrals": integrals.tolist(),
        "totalWedgeVolume": float(volumes.sum()),
        "sourceFiles": {path.name: _sha256(path) for path in vtk_paths},
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        metadata=json.dumps(metadata, sort_keys=True),
        phiASE=phi_ase,
        points=points,
        cells=cells,
        cellTypes=cell_types,
    )
    print(f"wrote {args.output}")
    print(f"integrals={integrals.tolist()}")


if __name__ == "__main__":
    main()
