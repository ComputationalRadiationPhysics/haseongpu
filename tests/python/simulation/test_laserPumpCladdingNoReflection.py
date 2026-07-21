# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from pyInclude.geometry.vtk import _parseVtk


repoRoot = Path(__file__).resolve().parents[3]
exampleDir = repoRoot / "example"
sys.path.insert(0, str(exampleDir))
import laserPumpCladding  # noqa: E402


REFERENCE_PATH = (
    repoRoot
    / "tests"
    / "data"
    / "laserPumpCladding"
    / "fixed_legacy_no_reflection_pump3_wedge_reference"
    / "phiase_reference.npz"
)
INTEGRAL_RTOL = 0.05


def _triangle_area(points):
    a, b, c = np.asarray(points, dtype=np.float64)
    ab = b[:2] - a[:2]
    ac = c[:2] - a[:2]
    return 0.5 * abs(float(ab[0] * ac[1] - ab[1] * ac[0]))


def _wedge_volume(points, cell):
    vertices = np.asarray(points, dtype=np.float64)[np.asarray(cell, dtype=np.uint32)]
    z_values = np.unique(vertices[:, 2])
    if z_values.size != 2:
        raise ValueError("wedge cell must span exactly two z levels")
    lower = vertices[np.isclose(vertices[:, 2], z_values[0])]
    if lower.shape[0] != 3:
        raise ValueError("wedge cell must have three lower vertices")
    return _triangle_area(lower) * float(z_values[1] - z_values[0])


def _tet_volume(points, cell):
    a, b, c, d = np.asarray(points, dtype=np.float64)[np.asarray(cell, dtype=np.uint32)]
    return abs(float(np.dot(b - a, np.cross(c - a, d - a)))) / 6.0


def _tet_cell_integral(points, cells, values):
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    volumes = np.asarray([_tet_volume(points, cell) for cell in cells], dtype=np.float64)
    if values.shape != volumes.shape:
        raise ValueError(f"Tet4 field has {values.size} values for {volumes.size} cells")
    return float(np.sum(values * volumes))


def _wedge_point_integrals(points, cells, phi_by_step):
    points = np.asarray(points, dtype=np.float64)
    cells = np.asarray(cells, dtype=np.uint32)
    volumes = np.asarray([_wedge_volume(points, cell) for cell in cells], dtype=np.float64)
    return np.asarray(
        [float(np.sum(phi[cells].mean(axis=1) * volumes)) for phi in phi_by_step],
        dtype=np.float64,
    )


@pytest.fixture(scope="module")
def legacyNoReflectionReference():
    if not REFERENCE_PATH.is_file():
        pytest.skip(f"missing generated no-reflection reference data: {REFERENCE_PATH}")
    with np.load(REFERENCE_PATH, allow_pickle=False) as data:
        return {
            "metadata": json.loads(str(data["metadata"].item())),
            "phiASE": np.asarray(data["phiASE"], dtype=np.float64),
            "points": np.asarray(data["points"], dtype=np.float64),
            "cells": np.asarray(data["cells"], dtype=np.uint32),
            "cellTypes": np.asarray(data["cellTypes"], dtype=np.uint32),
        }


def testNoReflectionReferenceDocumentsGeneratorGeometryAndSeed(legacyNoReflectionReference):
    metadata = legacyNoReflectionReference["metadata"]

    assert metadata["generator"]["commit"] == "effd8077edccef93a68d818e8a5eb2f0ebdc03b4"
    assert metadata["parameters"]["timeSlices"] == 6
    assert metadata["parameters"]["pumpSteps"] == 3
    assert metadata["parameters"]["useReflections"] is False
    assert metadata["random"]["rngSeed"] == 5489
    assert "useReflections: false" in metadata["generator"]["configText"]
    assert metadata["generator"]["geometrySource"] == "example/python_example/legacy/pt.mat"
    assert (
        metadata["generator"]["geometrySourceSha256"]
        == "afab3241bb89045a2234006f4c2eca26194bef761e1265e78c5115e6898ed74a"
    )


def testNoReflectionReferenceStoresFullWedgePointBuffers(legacyNoReflectionReference):
    metadata = legacyNoReflectionReference["metadata"]
    phi = legacyNoReflectionReference["phiASE"]

    assert phi.shape == (6, 4210)
    assert phi.shape == tuple(metadata["observable"]["shape"])
    assert legacyNoReflectionReference["points"].shape == (4210, 3)
    assert legacyNoReflectionReference["cells"].shape == (7308, 6)
    assert np.all(legacyNoReflectionReference["cellTypes"] == 13)
    assert np.isfinite(phi).all()
    np.testing.assert_allclose(phi.mean(axis=1), metadata["meanPhi"], rtol=0.0, atol=0.0)


@pytest.mark.integration
def testCurrentTet4WithoutReflectionsMatchesLegacyWedgeIntegral(
    legacyNoReflectionReference,
    alpakaRuntimeBackend,
    openPmdRuntimeBackend,
    tmp_path,
):
    metadata = legacyNoReflectionReference["metadata"]
    assert metadata["parameters"]["useReflections"] is False

    output_dir = tmp_path / "no_reflection_tet4"
    state = laserPumpCladding.runExample(
        backend=alpakaRuntimeBackend,
        openpmdBackend=openPmdRuntimeBackend,
        timeSlices=metadata["parameters"]["timeSlices"],
        pumpSteps=metadata["parameters"]["pumpSteps"],
        vtkOutputDir=output_dir,
        enableASE=True,
        prePump=metadata["parameters"]["prePump"],
        rngSeed=metadata["random"]["rngSeed"],
        useReflections=False,
        minRays=metadata["parameters"]["minRaysPerSample"],
        maxRays=metadata["parameters"]["maxRaysPerSample"],
        relativeStandardErrorThreshold=0.05,
        adaptiveSteps=metadata["parameters"]["adaptiveSteps"],
    )

    assert state.aseResult.srmStatus == "disabled"
    assert state.aseResult.srmPasses == 0
    assert state.aseResult.srmMaxIterations == 0
    # A stale all-zero beta-volume CDF sends every source history to the final
    # Tet. Keep the source visits dispersed as beta evolves from its zero state.
    assert np.max(state.volume_total_rays) < metadata["parameters"]["maxRaysPerSample"] // 20

    observed_integrals = []
    for step in metadata["observable"]["stepNumbers"]:
        points, cells, _cell_types, _point_data, cell_data, _fields = _parseVtk(
            output_dir / f"laserPumpCladding_{step:03d}.vtk"
        )
        observed_integrals.append(_tet_cell_integral(points, cells, cell_data["volumePhiASE"]))

    reference_integrals = _wedge_point_integrals(
        legacyNoReflectionReference["points"],
        legacyNoReflectionReference["cells"],
        legacyNoReflectionReference["phiASE"],
    )
    np.testing.assert_allclose(
        observed_integrals,
        reference_integrals,
        rtol=INTEGRAL_RTOL,
        atol=0.0,
    )
