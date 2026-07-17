# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Static no-reflection regression for actual laserPumpCladding volume cells."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from HASEonGPU import AlpakaBackends, CrossSectionData, PhiASE


repoRoot = Path(__file__).resolve().parents[3]
exampleDir = repoRoot / "example"
sys.path.insert(0, str(exampleDir))
import laserPumpCladding  # noqa: E402


REFERENCE_PATH = (
    repoRoot
    / "tests"
    / "data"
    / "laserPumpCladding"
    / "true_cladding_no_reflection_phiase_reference"
    / "phiase_reference.npz"
)
LEGACY_COMMIT = "effd8077edccef93a68d818e8a5eb2f0ebdc03b4"
CLADDING_NUMBER = 1
PHYSICAL_CLADDING_ABSORPTION = 5.5
INSTRUMENTATION_CLADDING_ABSORPTION = 55.0
CURRENT_FORWARD_RAYS = 4_000_000


def _tet_volume(points, cell):
    a, b, c, d = np.asarray(points, dtype=np.float64)[np.asarray(cell, dtype=np.uint32)]
    return abs(float(np.dot(b - a, np.cross(c - a, d - a)))) / 6.0


def _wedge_volumes(points, cells):
    points = np.asarray(points, dtype=np.float64)
    cells = np.asarray(cells, dtype=np.uint32)
    lower = points[cells[:, :3]]
    ab = lower[:, 1, :2] - lower[:, 0, :2]
    ac = lower[:, 2, :2] - lower[:, 0, :2]
    areas = 0.5 * np.abs(ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0])
    vertices = points[cells]
    return areas * (vertices[:, :, 2].max(axis=1) - vertices[:, :, 2].min(axis=1))


def _partitioned_wedge_integrals(reference, phi):
    cells = reference["cells"]
    volumes = _wedge_volumes(reference["points"], cells)
    values = np.asarray(phi, dtype=np.float64)[cells].mean(axis=1)
    cladding = reference["wedgeCladdingCellTypes"] == CLADDING_NUMBER
    weighted = values * volumes
    return {
        "total": float(weighted.sum()),
        "gain": float(weighted[~cladding].sum()),
        "cladding": float(weighted[cladding].sum()),
    }


def _partitioned_tet_integrals(topology, tet_types, phi):
    values = np.asarray(phi, dtype=np.float64)
    tet_types = np.asarray(tet_types, dtype=np.uint32)
    volumes = np.asarray(topology.cellVolumes, dtype=np.float64)
    if values.shape != (topology.numberOfCells,):
        raise ValueError(
            f"current PhiASE must contain one value per Tet4 cell, got {values.shape}"
        )
    cladding = tet_types == CLADDING_NUMBER
    weighted = values * volumes
    return {
        "total": float(weighted.sum()),
        "gain": float(weighted[~cladding].sum()),
        "cladding": float(weighted[cladding].sum()),
    }


def _serial_backend():
    backends = AlpakaBackends.all()
    for backend in backends:
        if "Host_Cpu_CpuSerial" in backend:
            return backend
    raise RuntimeError(
        "true-cladding regression requires the CPU serial backend, "
        f"available: {backends}"
    )


@pytest.fixture(scope="module")
def trueCladdingReference():
    if not REFERENCE_PATH.is_file():
        pytest.skip(f"missing generated true-cladding reference: {REFERENCE_PATH}")
    with np.load(REFERENCE_PATH, allow_pickle=False) as data:
        return {
            "metadata": json.loads(str(data["metadata"].item())),
            "phiASE": np.asarray(data["phiASE"], dtype=np.float64),
            "mse": np.asarray(data["mse"], dtype=np.float64),
            "totalRays": np.asarray(data["totalRays"], dtype=np.uint32),
            "points": np.asarray(data["points"], dtype=np.float64),
            "cells": np.asarray(data["cells"], dtype=np.uint32),
            "cellTypes": np.asarray(data["cellTypes"], dtype=np.uint32),
            "baseCladdingCellTypes": np.asarray(
                data["baseCladdingCellTypes"], dtype=np.uint32
            ),
            "wedgeCladdingCellTypes": np.asarray(
                data["wedgeCladdingCellTypes"], dtype=np.uint32
            ),
            "wedgeBetaVolume": np.asarray(data["wedgeBetaVolume"], dtype=np.float64),
        }


def _tet_types(reference):
    return np.repeat(reference["wedgeCladdingCellTypes"], 3).astype(np.uint32)


def _make_current_medium(reference, cladding_absorption):
    medium = laserPumpCladding.laserPumpCladdingMedium(cladAbsorption=cladding_absorption)
    topology = medium.topology
    tet_types = _tet_types(reference)
    cladding_cells = np.flatnonzero(tet_types == CLADDING_NUMBER)
    topology = topology.withCellDomains(
        [
            {"domain": 1, "name": "gain_medium", "where": "all"},
            {"domain": 2, "name": "cladding_volume", "cellIndices": cladding_cells},
        ]
    )
    medium.topology = topology

    gain_beta = float(reference["metadata"]["material"]["gainBetaVolume"])
    beta_volume = np.full(topology.numberOfCells, gain_beta, dtype=np.float64)
    beta_volume[cladding_cells] = 0.0

    # PhiASE source selection and propagation use per-cell betaVolume. Keep
    # shared interface betaCells uniform so a shared vertex cannot erase a
    # neighboring gain cell; cladding cells still have exactly zero source.
    medium.get("betaCells").value = np.full(
        medium.get("betaCells").expectedShape,
        gain_beta,
        dtype=np.float64,
    )
    medium.get("betaVolume").value = beta_volume
    medium.get("claddingCellTypes").value = tet_types
    medium.get("claddingNumber").value = CLADDING_NUMBER
    medium.get("claddingAbsorption").value = cladding_absorption
    return medium


def testTrueCladdingReferenceDocumentsStaticLegacyContract(trueCladdingReference):
    reference = trueCladdingReference
    metadata = reference["metadata"]

    assert metadata["generator"]["commit"] == LEGACY_COMMIT
    assert metadata["parameters"]["useReflections"] is False
    assert metadata["parameters"]["monochromatic"] is True
    assert metadata["material"]["physicalCladdingAbsorption"] == PHYSICAL_CLADDING_ABSORPTION
    assert (
        metadata["material"]["instrumentationCladdingAbsorption"]
        == INSTRUMENTATION_CLADDING_ABSORPTION
    )
    assert metadata["material"]["claddingBetaVolume"] == 0.0
    assert reference["phiASE"].shape == (3, 4210)
    assert reference["points"].shape == (4210, 3)
    assert reference["cells"].shape == (7308, 6)
    assert reference["totalRays"].shape == (3, 4210)
    assert np.all(reference["totalRays"] == 2000)
    assert np.all(reference["cellTypes"] == 13)
    assert np.isfinite(reference["phiASE"]).all()
    assert np.max(reference["mse"]) < 0.13
    assert np.count_nonzero(reference["baseCladdingCellTypes"] == CLADDING_NUMBER) == 24
    assert np.count_nonzero(reference["wedgeCladdingCellTypes"] == CLADDING_NUMBER) == 216
    cladding = reference["wedgeCladdingCellTypes"] == CLADDING_NUMBER
    assert np.count_nonzero(reference["wedgeBetaVolume"][cladding]) == 0
    np.testing.assert_array_equal(
        reference["wedgeBetaVolume"][~cladding],
        np.full(np.count_nonzero(~cladding), metadata["material"]["gainBetaVolume"]),
    )

    absorptions = metadata["material"]["claddingAbsorptions"]
    for index, absorption in enumerate(absorptions):
        assert metadata["diagnostics"][index]["claddingAbsorption"] == absorption
        observed = _partitioned_wedge_integrals(reference, reference["phiASE"][index])
        for partition, value in observed.items():
            np.testing.assert_allclose(
                value,
                metadata["diagnostics"][index]["integrals"][partition],
                rtol=1.0e-14,
                atol=0.0,
            )

    no_absorption = metadata["diagnostics"][0]["integrals"]
    physical = metadata["diagnostics"][1]["integrals"]
    instrumentation = metadata["diagnostics"][2]["integrals"]
    assert physical["cladding"] < no_absorption["cladding"] * 0.90
    assert instrumentation["cladding"] < no_absorption["cladding"] * 0.85


def testTrueCladdingShellMapsEveryLegacyWedgeToThreeTetChildren(trueCladdingReference):
    reference = trueCladdingReference
    medium = _make_current_medium(reference, PHYSICAL_CLADDING_ABSORPTION)
    topology = medium.topology
    tet_types = np.asarray(medium.get("claddingCellTypes").value, dtype=np.uint32)
    tet_cells = np.asarray(topology.cellPointIndices, dtype=np.uint32)

    assert tet_cells.shape == (21924, 4)
    assert tet_types.shape == (21924,)
    np.testing.assert_array_equal(
        tet_types.reshape((-1, 3)),
        np.broadcast_to(reference["wedgeCladdingCellTypes"][:, None], (7308, 3)),
    )
    np.testing.assert_array_equal(
        reference["wedgeCladdingCellTypes"].reshape((9, 812)),
        np.broadcast_to(reference["baseCladdingCellTypes"], (9, 812)),
    )
    assert np.count_nonzero(tet_types == CLADDING_NUMBER) == 648
    assert topology.cellDomainNames == {1: "gain_medium", 2: "cladding_volume"}
    np.testing.assert_array_equal(topology.cellDomains == 2, tet_types == CLADDING_NUMBER)

    for wedge, children in zip(reference["cells"], tet_cells.reshape((-1, 3, 4)), strict=True):
        np.testing.assert_array_equal(np.unique(children), np.sort(wedge))

    tet_volumes = np.asarray(
        [_tet_volume(topology.points, cell) for cell in tet_cells],
        dtype=np.float64,
    )
    geometry = reference["metadata"]["geometry"]
    np.testing.assert_allclose(
        tet_volumes.sum(), geometry["totalVolume"], rtol=2.0e-9, atol=0.0
    )
    np.testing.assert_allclose(
        tet_volumes[tet_types == CLADDING_NUMBER].sum(),
        geometry["claddingVolume"],
        rtol=2.0e-9,
        atol=0.0,
    )
    beta_volume = np.asarray(medium.get("betaVolume").value, dtype=np.float64)
    assert np.count_nonzero(beta_volume[tet_types == CLADDING_NUMBER]) == 0
    assert np.all(beta_volume[tet_types != CLADDING_NUMBER] > 0.0)


@pytest.fixture(scope="module")
def currentTrueCladdingResults(trueCladdingReference, openPmdFileBackend):
    reference = trueCladdingReference
    material = reference["metadata"]["material"]
    cross_sections = CrossSectionData.monochromatic(
        wavelength=material["wavelength"],
        crossSectionAbsorption=material["crossSectionAbsorption"],
        crossSectionEmission=material["crossSectionEmission"],
    )
    results = []
    for cladding_absorption in material["claddingAbsorptions"]:
        medium = _make_current_medium(reference, cladding_absorption)
        phi_ase = PhiASE(
            crossSections=cross_sections,
            minRays=CURRENT_FORWARD_RAYS,
            maxRays=CURRENT_FORWARD_RAYS,
            forwardRayCount=CURRENT_FORWARD_RAYS,
            relativeStandardErrorThreshold=1.0,
            repetitions=1,
            adaptiveSteps=1,
            useReflections=False,
            monochromatic=True,
            backend=_serial_backend(),
            openpmdBackend=openPmdFileBackend,
            rngSeed=reference["metadata"]["parameters"]["rngSeed"],
        )
        phi_ase.run(gainMedium=medium, crossSections=cross_sections)
        result = phi_ase.getResults()
        phi = np.asarray(result.phiAse, dtype=np.float64)
        relative_standard_error = np.asarray(result.relativeStandardError, dtype=np.float64)
        tet_types = np.asarray(medium.get("claddingCellTypes").value, dtype=np.uint32)
        assert phi.shape == (medium.topology.numberOfCells,)
        assert phi.shape == (21924,)
        assert relative_standard_error.shape == phi.shape
        results.append(
            {
                "claddingAbsorption": cladding_absorption,
                "phiASE": phi,
                "integrals": _partitioned_tet_integrals(medium.topology, tet_types, phi),
                "relativeStandardError": relative_standard_error,
            }
        )
    return results


@pytest.mark.integration
def testCurrentTrueCladdingPhiAseMatchesLegacyTotalIntegral(
    trueCladdingReference,
    currentTrueCladdingResults,
):
    reference = trueCladdingReference
    for index, current in enumerate(currentTrueCladdingResults[:2]):
        legacy = reference["metadata"]["diagnostics"][index]["integrals"]
        assert np.isfinite(current["relativeStandardError"]).all()
        assert np.max(current["relativeStandardError"]) < 0.05
        np.testing.assert_allclose(
            current["integrals"]["total"],
            legacy["total"],
            rtol=0.05,
            atol=0.0,
            err_msg=f"claddingAbsorption={current['claddingAbsorption']}",
        )

    legacy_no_absorption = reference["metadata"]["diagnostics"][0]["integrals"]["total"]
    legacy_physical = reference["metadata"]["diagnostics"][1]["integrals"]["total"]
    current_no_absorption = currentTrueCladdingResults[0]["integrals"]["total"]
    current_physical = currentTrueCladdingResults[1]["integrals"]["total"]
    legacy_attenuation = legacy_physical / legacy_no_absorption
    current_attenuation = current_physical / current_no_absorption
    assert abs(current_attenuation - legacy_attenuation) < 0.05


@pytest.mark.integration
def testCurrentTrueCladdingHighContrastAbsorptionIsExercised(currentTrueCladdingResults):
    no_absorption, physical, instrumentation = currentTrueCladdingResults
    assert physical["integrals"]["cladding"] < no_absorption["integrals"]["cladding"] * 0.50
    assert instrumentation["integrals"]["cladding"] < no_absorption["integrals"]["cladding"] * 0.10
