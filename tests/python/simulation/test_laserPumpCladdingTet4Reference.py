# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest


repoRoot = Path(__file__).resolve().parents[3]
exampleDir = repoRoot / "example"
sys.path.insert(0, str(exampleDir))
import laserPumpCladding  # noqa: E402


referenceRoot = repoRoot / "tests" / "data" / "laserPumpCladding"
referencePath = referenceRoot / "upstream_master_reference.json"


def _tokens(path):
    return Path(path).read_text(encoding="utf-8").split()


def _vtkCells(path):
    tokens = _tokens(path)
    index = tokens.index("CELLS")
    count = int(tokens[index + 1])
    cursor = index + 3
    cells = []
    for _ in range(count):
        width = int(tokens[cursor])
        cursor += 1
        cells.append([int(value) for value in tokens[cursor:cursor + width]])
        cursor += width
    return cells


def _vtkCellTypes(path):
    tokens = _tokens(path)
    index = tokens.index("CELL_TYPES")
    count = int(tokens[index + 1])
    start = index + 2
    return np.asarray(tokens[start:start + count], dtype=np.uint32)


def _vtkScalar(path, name):
    tokens = _tokens(path)
    active_count = None
    index = 0
    while index < len(tokens):
        token = tokens[index].upper()
        if token in {"POINT_DATA", "CELL_DATA"}:
            active_count = int(tokens[index + 1])
            index += 2
            continue
        if token == "SCALARS" and tokens[index + 1] == name:
            components = 1
            cursor = index + 4
            if tokens[cursor].upper() != "LOOKUP_TABLE":
                components = int(tokens[cursor])
                cursor += 1
            if tokens[cursor].upper() == "LOOKUP_TABLE":
                cursor += 2
            values = np.asarray(tokens[cursor:cursor + active_count * components], dtype=np.float64)
            return values.reshape(active_count, components).squeeze()
        index += 1
    raise KeyError(name)


@pytest.fixture(scope="module")
def laserPumpCladdingReference():
    if not referencePath.is_file():
        pytest.skip(f"missing generated laserPumpCladding reference data: {referencePath}")
    return json.loads(referencePath.read_text(encoding="utf-8"))


def testLaserPumpCladdingReferenceFilesContainFiveMeanPhiValues(laserPumpCladdingReference):
    steps = laserPumpCladdingReference["steps"]
    assert [step["step"] for step in steps] == [1, 2, 3, 4, 5]
    for step in steps:
        wedge_path = repoRoot / step["wedgeVtk"]
        tet_path = repoRoot / step["tet4Vtk"]
        assert wedge_path.is_file()
        assert tet_path.is_file()
        assert np.isfinite(step["meanPhi"])
        np.testing.assert_allclose(_vtkScalar(wedge_path, "phiASE").mean(), step["meanPhi"], rtol=0.0, atol=1e-14)
        assert len(_vtkCells(tet_path)) == 3 * len(_vtkCells(wedge_path))
        assert np.all(_vtkCellTypes(tet_path) == 10)


def testLaserPumpCladdingTet4InputLoadsLegacyPointSamples(laserPumpCladdingReference):
    first = laserPumpCladdingReference["steps"][0]
    medium = laserPumpCladding.loadLaserPumpCladdingTet4Medium(repoRoot / first["tet4Vtk"])

    assert medium.topology.numberOfPoints == first["points"]
    assert medium.topology.numberOfCells == first["tet4Cells"]
    assert medium.topology.numberOfSamplePoints == first["points"]
    assert medium.get("betaCells").value.size == first["points"]
    assert medium.get("betaVolume").value.size == first["tet4Cells"]
    assert np.isfinite(medium.get("betaCells").value).all()
    assert np.isfinite(medium.get("betaVolume").value).all()
    assert np.count_nonzero(medium.get("betaVolume").value) == first["tet4Cells"]


@pytest.mark.integration
def testLaserPumpCladdingTet4ForwardMeanPhiMatchesReferenceData(laserPumpCladdingReference):
    if os.environ.get("HASE_RUN_LASERPUMP_REFERENCE") != "1":
        pytest.skip("set HASE_RUN_LASERPUMP_REFERENCE=1 to run the full PhiASE reference comparison")

    backend = os.environ.get("HASE_LASERPUMP_REFERENCE_BACKEND", "Host_Cpu_CpuOmpBlocks")
    tolerance = float(os.environ.get("HASE_LASERPUMP_REFERENCE_RTOL", "0.35"))
    forwardRayLength = float(os.environ.get("HASE_LASERPUMP_REFERENCE_FORWARD_RAY_LENGTH", "1.0"))
    observed = []
    expected = []
    for step in laserPumpCladdingReference["steps"]:
        result = laserPumpCladding.runTet4PhiAseInput(
            repoRoot / step["tet4Vtk"],
            backend=backend,
            propagationMode="forward",
            forwardRayLength=forwardRayLength,
            useReflections=False,
            rngSeed=laserPumpCladdingReference.get("rngSeed", 12345),
        )
        phi = np.asarray(result.phiAse, dtype=np.float64).reshape(-1)
        observed.append(float(phi.mean()))
        expected.append(float(step["meanPhi"]))

    np.testing.assert_allclose(observed, expected, rtol=tolerance, atol=0.0)
