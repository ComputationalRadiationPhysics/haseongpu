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
if str(repoRoot) not in sys.path:
    sys.path.insert(0, str(repoRoot))

from pyInclude.geometry.vtk import _parseVtk
from utils.convert_vtk_topology import convertVtk


exampleDir = repoRoot / "example"
sys.path.insert(0, str(exampleDir))
import laserPumpCladding  # noqa: E402


REFERENCE_PATH = (
    repoRoot
    / "tests"
    / "data"
    / "laserPumpCladding"
    / "upstream_master_pump3_wedge_reference"
    / "phiase_reference.npz"
)
POINTWISE_RTOL = 0.35


@pytest.fixture(scope="module")
def laserPumpCladdingReference():
    if not REFERENCE_PATH.is_file():
        pytest.skip(f"missing generated laserPumpCladding reference data: {REFERENCE_PATH}")
    with np.load(REFERENCE_PATH, allow_pickle=False) as data:
        return {
            "metadata": json.loads(str(data["metadata"].item())),
            "phiASE": np.asarray(data["phiASE"], dtype=np.float64),
            "points": np.asarray(data["points"], dtype=np.float64),
            "cells": np.asarray(data["cells"], dtype=np.uint32),
            "cellTypes": np.asarray(data["cellTypes"], dtype=np.uint32),
        }


def testLaserPumpCladdingReferenceDocumentsGeneratorAndSeed(laserPumpCladdingReference):
    metadata = laserPumpCladdingReference["metadata"]

    assert metadata["generator"]["commit"] == "469c87770ed13796f2e82385bcf83528e8aeaf1b"
    assert metadata["parameters"]["backend"] == "Host_Cpu_CpuOmpBlocks"
    assert metadata["parameters"]["timeSlices"] == 6
    assert metadata["parameters"]["pumpSteps"] == 3
    assert metadata["random"]["rngSeed"] == 5489
    assert metadata["observable"]["field"] == "phiASE"
    assert metadata["observable"]["location"] == "legacy wedge POINT_DATA"
    assert metadata["observable"]["coverage"] == "full point buffer"


def testLaserPumpCladdingReferenceStoresFullWedgePointBuffers(laserPumpCladdingReference):
    metadata = laserPumpCladdingReference["metadata"]
    phi = laserPumpCladdingReference["phiASE"]

    assert phi.shape == tuple(metadata["observable"]["shape"])
    assert phi.shape == (6, 4210)
    assert laserPumpCladdingReference["points"].shape == (4210, 3)
    assert laserPumpCladdingReference["cells"].shape == (7308, 6)
    assert np.all(laserPumpCladdingReference["cellTypes"] == 13)
    assert np.isfinite(phi).all()
    np.testing.assert_allclose(phi.mean(axis=1), metadata["meanPhi"], rtol=0.0, atol=0.0)


def testCurrentTet4LaserPumpGeometryConvertsBackToLegacyWedgeOrder(
    laserPumpCladdingReference,
    tmp_path,
):
    wedge_path = convertVtk(
        repoRoot / "example" / "data" / "pt.vtk",
        tmp_path / "pt_wedge.vtk",
        direction="tet4-to-wedge",
    )
    points, cells, cell_types, point_data, cell_data, _fields = _parseVtk(wedge_path)

    np.testing.assert_allclose(points, laserPumpCladdingReference["points"], rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(np.asarray(cells, dtype=np.uint32), laserPumpCladdingReference["cells"])
    np.testing.assert_array_equal(np.asarray(cell_types, dtype=np.uint32), laserPumpCladdingReference["cellTypes"])
    assert "betaCells" in point_data
    assert "betaVolume" in cell_data


@pytest.mark.integration
@pytest.mark.xfail(
    reason=(
        "Current forward Tet4 PhiASE is spatially sparse relative to the "
        "legacy upstream/master wedge reference; keep this as the executable "
        "pointwise comparison until the physics discrepancy is resolved."
    ),
    strict=False,
)
def testCurrentTet4ForwardPhiAseMatchesLegacyWedgeReferencePointwise(
    laserPumpCladdingReference,
    openPmdRuntimeBackend,
    tmp_path,
):
    metadata = laserPumpCladdingReference["metadata"]
    backend = os.environ.get("HASE_LASERPUMP_REFERENCE_BACKEND", metadata["parameters"]["backend"])
    rtol = float(os.environ.get("HASE_LASERPUMP_POINTWISE_RTOL", str(POINTWISE_RTOL)))
    tet4_dir = tmp_path / "current_tet4"
    wedge_dir = tmp_path / "current_wedge"
    tet4_dir.mkdir()
    wedge_dir.mkdir()

    laserPumpCladding.runExample(
        backend=backend,
        openpmdBackend=openPmdRuntimeBackend,
        timeSlices=metadata["parameters"]["timeSlices"],
        pumpSteps=metadata["parameters"]["pumpSteps"],
        vtkOutputDir=tet4_dir,
        enableASE=True,
        rngSeed=metadata["random"]["rngSeed"],
        useReflections=False,
        minRaysPerSample=metadata["parameters"]["minRaysPerSample"],
        maxRaysPerSample=metadata["parameters"]["maxRaysPerSample"],
        mseThreshold=metadata["parameters"]["mseThreshold"],
        adaptiveSteps=metadata["parameters"]["adaptiveSteps"],
    )

    observed = []
    for step in metadata["observable"]["stepNumbers"]:
        tet4_path = tet4_dir / f"laserPumpCladding_{step:03d}.vtk"
        wedge_path = convertVtk(
            tet4_path,
            wedge_dir / tet4_path.name,
            direction="tet4-to-wedge",
        )
        points, cells, cell_types, point_data, _cell_data, _fields = _parseVtk(wedge_path)
        np.testing.assert_allclose(points, laserPumpCladdingReference["points"], rtol=0.0, atol=0.0)
        np.testing.assert_array_equal(np.asarray(cells, dtype=np.uint32), laserPumpCladdingReference["cells"])
        np.testing.assert_array_equal(np.asarray(cell_types, dtype=np.uint32), laserPumpCladdingReference["cellTypes"])
        observed.append(np.asarray(point_data["phiASE"], dtype=np.float64).reshape(-1))

    observed = np.vstack(observed)
    np.testing.assert_allclose(observed, laserPumpCladdingReference["phiASE"], rtol=rtol, atol=0.0)
