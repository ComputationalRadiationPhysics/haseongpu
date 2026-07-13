# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pyInclude import AlpakaBackends
from pyInclude.geometry import GainMedium, VolumeTopology
from pyInclude.laser import CrossSectionData
from pyInclude.openpmd import backendFlat
from pyInclude.openpmd import transport as openpmd_transport
from pyInclude.simulation import PhiASE


FIXTURE_DIR = Path(__file__).resolve().parents[2] / "data" / "juliaASE" / "reflection_surface_reference"
REFERENCE_PATH = FIXTURE_DIR / "reference.json"


def _cpu_backend_or_skip():
    try:
        backends = list(AlpakaBackends.all())
    except Exception as exc:  # pragma: no cover - depends on native build availability
        pytest.skip(f"no Alpaka backend metadata is available: {exc}")
    for preferred in ("Host_Cpu_CpuSerial", "Host_Cpu_CpuOmpBlocks"):
        if preferred in backends:
            return preferred
    for backend in backends:
        if "Cpu" in backend:
            return backend
    pytest.skip(f"no CPU Alpaka backend is available; reported backends: {backends}")


def _native_executable_or_skip():
    try:
        return openpmd_transport.findCalcPhiAse()
    except (FileNotFoundError, RuntimeError) as exc:
        pytest.skip(f"no usable native hase-cpp executable is available: {exc}")


def _medium_from_reference(reference):
    topology = VolumeTopology.fromTetrahedra(
        np.asarray(reference["points"], dtype=np.float64),
        np.asarray(reference["cells"], dtype=np.uint32),
        faceBoundaries=np.asarray(reference["faceBoundaries"], dtype=np.int32),
        metadata={"source": str(FIXTURE_DIR / reference["meshFile"]), "format": "gmsh"},
    )
    material = reference["material"]
    legacy = reference["legacyOpticsFallback"]
    optics = reference["surfaceOptics"]
    return GainMedium(topology).withPhysicalProperties(
        betaVolume=backendFlat(np.asarray(reference["initialBetaVolume"], dtype=np.float64)),
        betaCells=backendFlat(np.asarray(reference["betaCells"], dtype=np.float64)),
        claddingCellTypes=np.asarray(reference["claddingCellTypes"], dtype=np.uint32),
        refractiveIndices=np.asarray(legacy["refractiveIndices"], dtype=np.float32),
        reflectivities=backendFlat(np.asarray(legacy["reflectivities"], dtype=np.float32)),
        surfaceReflectivity=np.asarray(optics["surfaceReflectivity"], dtype=np.float32),
        surfaceRefractiveIndexInside=np.asarray(optics["surfaceRefractiveIndexInside"], dtype=np.float32),
        surfaceRefractiveIndexOutside=np.asarray(optics["surfaceRefractiveIndexOutside"], dtype=np.float32),
        nTot=float(material["nTot"]),
        crystalTFluo=float(material["crystalTFluo"]),
        claddingNumber=int(material["claddingNumber"]),
        claddingAbsorption=float(material["claddingAbsorption"]),
    )


def _cross_sections_from_reference(reference):
    cross_sections = reference["crossSections"]
    return CrossSectionData(
        wavelengthsAbsorption=np.asarray(cross_sections["wavelengthsAbsorption"], dtype=np.float64),
        crossSectionAbsorption=np.asarray(cross_sections["crossSectionAbsorption"], dtype=np.float64),
        wavelengthsEmission=np.asarray(cross_sections["wavelengthsEmission"], dtype=np.float64),
        crossSectionEmission=np.asarray(cross_sections["crossSectionEmission"], dtype=np.float64),
        resolution=int(cross_sections["resolution"]),
    )


def test_hase_forward_reflection_matches_committed_juliaase_surface_fixture():
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    _native_executable_or_skip()
    backend = _cpu_backend_or_skip()

    phi_ase = PhiASE(
        crossSections=_cross_sections_from_reference(reference),
        propagationMode="forward",
        minRaysPerSample=int(reference["rayCount"]),
        maxRaysPerSample=int(reference["rayCount"]),
        forwardRayCount=int(reference["rayCount"]),
        mseThreshold=1.0,
        repetitions=1,
        adaptiveSteps=1,
        useReflections=True,
        reflectionMaxIterations=int(reference["reflectionMaxIterations"]),
        reflectionTolerance=float(reference["reflectionTolerance"]),
        surfaceReservoirSize=int(reference["surfaceReservoirSize"]),
        monochromatic=True,
        backend=backend,
        openpmdBackend="adios",
        parallelMode="single",
        numDevices=1,
        rngSeed=int(reference["seed"]),
    )

    phi_ase.run(gainMedium=_medium_from_reference(reference))

    result = phi_ase.getResults()
    actual_phi = np.asarray(result.phiAse, dtype=np.float64)
    actual_dndt = np.asarray(result.dndtAse, dtype=np.float64)
    actual_final_beta = np.asarray(reference["initialBetaVolume"], dtype=np.float64) - float(reference["timeStep"]) * actual_dndt
    tolerances = reference["tolerances"]

    assert result.srmStatus == "converged"
    assert result.srmPasses == 1
    assert result.srmMaxIterations == int(reference["reflectionMaxIterations"])
    assert result.srmDivergenceStreak == 3
    assert result.srmRemainingFraction == pytest.approx(0.0)

    np.testing.assert_allclose(
        actual_phi,
        np.asarray(reference["phiAse"], dtype=np.float64),
        rtol=float(tolerances["phiAseRtol"]),
        atol=float(tolerances["phiAseAtol"]),
    )
    np.testing.assert_allclose(
        actual_final_beta,
        np.asarray(reference["finalBetaVolume"], dtype=np.float64),
        rtol=float(tolerances["betaVolumeRtol"]),
        atol=float(tolerances["betaVolumeAtol"]),
    )
