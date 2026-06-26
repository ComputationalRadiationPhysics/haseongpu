#!/usr/bin/env python3
# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Convert bundled legacy material inputs to frontend-readable VTK fixtures."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from HASEonGPU import GainMedium, MeshTopology, backendFlat  # noqa: E402


def _scalar(path, dtype=float):
    return dtype(np.loadtxt(path))


def _array(path, dtype=float):
    return np.asarray(np.loadtxt(path), dtype=dtype)


def _text_medium(input_dir: Path) -> GainMedium:
    number_of_points = _scalar(input_dir / "numberOfPoints.txt", int)
    number_of_triangles = _scalar(input_dir / "numberOfTriangles.txt", int)
    number_of_levels = _scalar(input_dir / "numberOfLevels.txt", int)
    thickness = _scalar(input_dir / "thickness.txt", float)

    points = _array(input_dir / "points.txt", np.float64)
    if points.size != 2 * number_of_points:
        raise ValueError(f"{input_dir}/points.txt has {points.size} values, expected {2 * number_of_points}")
    points = points.reshape((2, number_of_points)).T

    triangles = _array(input_dir / "trianglePointIndices.txt", np.uint32)
    if triangles.size != 3 * number_of_triangles:
        raise ValueError(
            f"{input_dir}/trianglePointIndices.txt has {triangles.size} values, "
            f"expected {3 * number_of_triangles}"
        )
    triangles = triangles.reshape((3, number_of_triangles)).T

    topology = MeshTopology(
        points=points,
        trianglePointIndices=triangles,
        levels=number_of_levels,
        thickness=thickness,
        metadata={"source": str(input_dir), "format": "legacy-text"},
    )
    return GainMedium(topology=topology).withPhysicalProperties(
        betaCells=backendFlat(_array(input_dir / "betaCells.txt", np.float64)),
        betaVolume=backendFlat(_array(input_dir / "betaVolume.txt", np.float64)),
        claddingCellTypes=_array(input_dir / "claddingCellTypes.txt", np.uint32),
        refractiveIndices=_array(input_dir / "refractiveIndices.txt", np.float32),
        reflectivities=backendFlat(_array(input_dir / "reflectivities.txt", np.float32)),
        nTot=_scalar(input_dir / "nTot.txt", float),
        crystalTFluo=_scalar(input_dir / "crystalTFluo.txt", float),
        claddingNumber=_scalar(input_dir / "claddingNumber.txt", int),
        claddingAbsorption=_scalar(input_dir / "claddingAbsorption.txt", float),
    )


def _pt_mat_medium(material_path: Path, *, number_of_levels: int = 10, length: float = 0.7) -> GainMedium:
    try:
        from scipy.io import loadmat
    except ImportError as exc:
        raise RuntimeError("scipy is required to convert pt.mat") from exc

    material = loadmat(material_path)
    points = np.asarray(material["p"], dtype=np.float64)
    triangles = np.asarray(material["t"], dtype=np.int64) - 1
    if np.any(triangles < 0):
        raise ValueError(f"{material_path} does not contain MATLAB-style 1-based triangle indices")

    topology = MeshTopology(
        points=points[:, :2],
        trianglePointIndices=triangles.astype(np.uint32),
        levels=number_of_levels,
        thickness=length / (number_of_levels - 1),
        metadata={"source": str(material_path), "format": "matlab-pt"},
    )
    return GainMedium(topology=topology).withPhysicalProperties(
        betaCells=np.zeros((topology.numberOfPoints, topology.levels), dtype=np.float64),
        betaVolume=np.zeros((topology.numberOfTriangles, topology.levels - 1), dtype=np.float64),
        claddingCellTypes=np.zeros(topology.numberOfTriangles, dtype=np.uint32),
        refractiveIndices=np.asarray([1.83, 1.0, 1.83, 1.0], dtype=np.float32),
        reflectivities=np.zeros((topology.numberOfTriangles, 2), dtype=np.float32),
        nTot=2.76e20,
        crystalTFluo=9.5e-4,
        claddingNumber=1,
        claddingAbsorption=5.5,
    )


def convert(repo_root: Path = REPO_ROOT, output_dir: Path | None = None) -> list[Path]:
    output_dir = repo_root / "example" / "data" if output_dir is None else output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    inputs = repo_root / "example" / "c_example" / "input"
    outputs = [
        (_pt_mat_medium(repo_root / "example" / "python_example" / "legacy" / "pt.mat"), output_dir / "pt.vtk"),
        (_text_medium(inputs / "cuboid"), output_dir / "cuboid.vtk"),
        (_text_medium(inputs / "cylindrical"), output_dir / "cylindrical.vtk"),
    ]
    return [medium.toVtk(path) for medium, path in outputs]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "example" / "data",
        help="Directory for generated VTK fixtures.",
    )
    args = parser.parse_args(argv)

    for path in convert(REPO_ROOT, args.output_dir):
        print(path)


if __name__ == "__main__":
    main()
