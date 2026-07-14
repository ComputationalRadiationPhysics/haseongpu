# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import importlib.util
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

repoRoot = Path(__file__).resolve().parents[3]


from pyInclude.geometry.vtk import _parseVtk
from pyInclude.openpmd import transport


_convert_vtk_topology_path = repoRoot / "utils" / "convert_vtk_topology.py"
_convert_vtk_topology_spec = importlib.util.spec_from_file_location(
    "_hase_convert_vtk_topology",
    _convert_vtk_topology_path,
)
if _convert_vtk_topology_spec is None or _convert_vtk_topology_spec.loader is None:
    raise ImportError(f"cannot load VTK topology converter from {_convert_vtk_topology_path}")
_convert_vtk_topology = importlib.util.module_from_spec(_convert_vtk_topology_spec)
sys.modules[_convert_vtk_topology_spec.name] = _convert_vtk_topology
_convert_vtk_topology_spec.loader.exec_module(_convert_vtk_topology)
convertVtk = _convert_vtk_topology.convertVtk


exampleDir = repoRoot / "example"
sys.path.insert(0, str(exampleDir))
import laserPumpCladding  # noqa: E402


REFERENCE_PATH = (
    repoRoot
    / "tests"
    / "data"
    / "laserPumpCladding"
    / "fixed_legacy_reflection_pump3_wedge_reference"
    / "phiase_reference.npz"
)
# The fixed legacy wedge and Tet4 solvers both reflect at the physical top
# plane, z = (numberOfLevels - 1) * thickness.
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
    cells = np.asarray(cells, dtype=np.uint32)
    volumes = np.asarray([_tet_volume(points, cell) for cell in cells], dtype=np.float64)
    if values.shape != volumes.shape:
        raise ValueError(f"Tet4 field has {values.size} values for {volumes.size} cells")
    return float(np.sum(values * volumes))


def _wedge_point_integrals(points, cells, phi_by_step):
    points = np.asarray(points, dtype=np.float64)
    cells = np.asarray(cells, dtype=np.uint32)
    volumes = np.asarray([_wedge_volume(points, cell) for cell in cells], dtype=np.float64)
    result = []
    for phi in np.asarray(phi_by_step, dtype=np.float64):
        cell_values = phi[cells].mean(axis=1)
        result.append(float(np.sum(cell_values * volumes)))
    return np.asarray(result, dtype=np.float64)


def _vtkScalarNames(path):
    tokens = path.read_text(encoding="utf-8").split()
    return {tokens[index + 1] for index, token in enumerate(tokens) if token.upper() == "SCALARS"}


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

    assert metadata["generator"]["commit"] == "effd8077edccef93a68d818e8a5eb2f0ebdc03b4"
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


def testPtTet4GeometryConvertsBackToLegacyWedgeOrder(
    laserPumpCladdingReference,
    tmp_path,
):
    wedge_path = convertVtk(
        repoRoot / "example" / "data" / "ptTet4.vtk",
        tmp_path / "pt_wedge.vtk",
        direction="tet4-to-wedge",
    )
    points, cells, cell_types, point_data, cell_data, _fields = _parseVtk(wedge_path)

    np.testing.assert_allclose(points, laserPumpCladdingReference["points"], rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(
        np.sort(np.asarray(cells, dtype=np.uint32), axis=1),
        np.sort(laserPumpCladdingReference["cells"], axis=1),
    )
    np.testing.assert_array_equal(np.asarray(cell_types, dtype=np.uint32), laserPumpCladdingReference["cellTypes"])
    assert "betaCells" in point_data
    assert "betaVolume" in cell_data


def testLaserPumpCladdingUsesLegacyEffectiveAseSpectrum():
    spectra = laserPumpCladding.laserPumpCladdingSpectralProperties()
    material_dir = repoRoot / "example" / "input"
    raw_wavelengths = np.loadtxt(material_dir / "lambda_a.txt")
    raw_absorption = np.loadtxt(material_dir / "sigma_a.txt")
    raw_emission = np.loadtxt(material_dir / "sigma_e.txt")
    effective_wavelengths = np.linspace(
        raw_wavelengths[0],
        raw_wavelengths[-1],
        laserPumpCladding.LEGACY_ASE_SPECTRAL_RESOLUTION,
    )

    assert spectra.resolution == 1000
    np.testing.assert_array_equal(spectra.wavelengthsAbsorption, effective_wavelengths)
    np.testing.assert_array_equal(spectra.wavelengthsEmission, effective_wavelengths)
    np.testing.assert_allclose(
        spectra.crossSectionAbsorption,
        np.interp(effective_wavelengths, raw_wavelengths, raw_absorption),
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        spectra.crossSectionEmission,
        np.interp(effective_wavelengths, raw_wavelengths, raw_emission),
        rtol=0.0,
        atol=0.0,
    )


def testLaserPumpCladdingTet4MediumAssignsCylinderSurfaceOptics():
    medium = laserPumpCladding.laserPumpCladdingMedium()
    topology = medium.topology
    boundaries = np.asarray(topology.faceBoundaries, dtype=np.int32)
    points = np.asarray(topology.points, dtype=np.float64)
    face_nodes = np.asarray(topology.facePointIndices, dtype=np.uint32)
    face_z = points[:, 2][face_nodes]

    bottom_id = laserPumpCladding.BOTTOM_ASE_SURFACE_ID
    top_id = laserPumpCladding.TOP_ASE_SURFACE_ID

    cladding_id = laserPumpCladding.CLADDING_SURFACE_ID

    assert topology.surfaceDomainNames == {
        bottom_id: "ase_bottom",
        top_id: "ase_top",
        cladding_id: "cladding",
    }
    assert np.count_nonzero(boundaries == bottom_id) == 812
    assert np.count_nonzero(boundaries == top_id) == 812
    assert np.count_nonzero(boundaries == cladding_id) == 432
    assert np.count_nonzero(topology.neighborCells < 0) == 2056
    np.testing.assert_allclose(face_z[boundaries == bottom_id], np.min(points[:, 2]), rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(face_z[boundaries == top_id], np.max(points[:, 2]), rtol=0.0, atol=1.0e-12)

    surface_reflectivity = np.asarray(medium.get("surfaceReflectivity").value, dtype=np.float32)
    surface_inside = np.asarray(medium.get("surfaceRefractiveIndexInside").value, dtype=np.float32)
    surface_outside = np.asarray(medium.get("surfaceRefractiveIndexOutside").value, dtype=np.float32)

    assert surface_reflectivity.shape == (cladding_id + 1,)
    np.testing.assert_array_equal(surface_reflectivity[[bottom_id, top_id, cladding_id]], np.asarray([0.0, 0.0, 0.0], dtype=np.float32))
    np.testing.assert_array_equal(surface_inside[[bottom_id, top_id, cladding_id]], np.asarray([1.83, 1.83, 1.0], dtype=np.float32))
    np.testing.assert_array_equal(surface_outside[[bottom_id, top_id, cladding_id]], np.asarray([1.0, 1.0, 1.0], dtype=np.float32))


def testLaserPumpCladdingTet4MediumPreservesLegacyTenLayerPumpLayout():
    medium = laserPumpCladding.laserPumpCladdingMedium()
    topology = medium.topology
    points = np.asarray(topology.points, dtype=np.float64)
    z_planes = np.unique(points[:, 2])

    assert topology.cellDomainNames == {1: "gain_medium"}
    np.testing.assert_array_equal(topology.cellDomains, np.ones(topology.numberOfCells, dtype=np.int32))
    assert z_planes.size == laserPumpCladding.NUMBER_OF_Z_LAYERS
    np.testing.assert_allclose(
        z_planes,
        np.linspace(0.0, 0.7, laserPumpCladding.NUMBER_OF_Z_LAYERS),
        rtol=0.0,
        atol=1.0e-12,
    )

    assert topology.structuredNumberOfLevels == laserPumpCladding.NUMBER_OF_Z_LAYERS
    assert topology.numberOfSamplePoints == topology.numberOfPoints
    assert topology.numberOfSamplePoints == topology.structuredNumberOfPoints * topology.structuredNumberOfLevels
    assert medium.get("betaCells").expectedShape == (topology.numberOfSamplePoints,)

    points_by_level = points.reshape(
        (topology.structuredNumberOfLevels, topology.structuredNumberOfPoints, 3)
    )
    np.testing.assert_allclose(
        points_by_level[:, :, 2],
        np.broadcast_to(z_planes[:, None], points_by_level[:, :, 2].shape),
        rtol=0.0,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        points_by_level[:, :, :2],
        np.broadcast_to(points_by_level[0, :, :2], points_by_level[:, :, :2].shape),
        rtol=0.0,
        atol=1.0e-12,
    )


@pytest.mark.integration
def testLaserPumpCladdingRunExampleReflectionToggleChangesPhiAse(
    tmp_path,
    alpakaRuntimeBackend,
    openPmdRuntimeBackend,
):
    def run(use_reflections):
        output_dir = tmp_path / f"reflections_{int(use_reflections)}"
        laserPumpCladding.runExample(
            backend=alpakaRuntimeBackend,
            openpmdBackend=openPmdRuntimeBackend,
            timeSlices=2,
            pumpSteps=1,
            vtkOutputDir=output_dir,
            enableASE=True,
            prePump=True,
            rngSeed=1234,
            useReflections=use_reflections,
        minRays=20_000,
        maxRays=20_000,
            adaptiveSteps=1,
            relativeStandardErrorThreshold=0.1,
            reflectionMaxIterations=17,
            reflectionTolerance=0.0,
            surfaceReservoirSize=32,
        )
        points, cells, _cell_types, _point_data, cell_data, _fields = _parseVtk(
            output_dir / "laserPumpCladding_002.vtk"
        )
        return _tet_cell_integral(points, cells, cell_data["volumePhiASE"])

    without_reflections = run(False)
    with_reflections = run(True)

    assert without_reflections > 0.0
    assert with_reflections > without_reflections * 1.05


@pytest.mark.integration
def testCurrentTet4ForwardPhiAseMatchesLegacyWedgeReferenceIntegral(
    laserPumpCladdingReference,
    openPmdRuntimeBackend,
    tmp_path,
):
    metadata = laserPumpCladdingReference["metadata"]
    backend = os.environ.get("HASE_LASERPUMP_REFERENCE_BACKEND", metadata["parameters"]["backend"])
    rtol = 0.05
    tet4_dir = tmp_path / "current_tet4"
    tet4_dir.mkdir()

    state = laserPumpCladding.runExample(
        backend=backend,
        openpmdBackend=openPmdRuntimeBackend,
        timeSlices=metadata["parameters"]["timeSlices"],
        pumpSteps=metadata["parameters"]["pumpSteps"],
        vtkOutputDir=tet4_dir,
        enableASE=True,
        prePump=metadata["parameters"]["prePump"],
        rngSeed=metadata["random"]["rngSeed"],
        useReflections=metadata["parameters"]["useReflections"],
        minRays=metadata["parameters"]["minRaysPerSample"],
        maxRays=metadata["parameters"]["maxRaysPerSample"],
        relativeStandardErrorThreshold=0.05,
        adaptiveSteps=metadata["parameters"]["adaptiveSteps"],
    )

    relative_standard_error = np.asarray(state.volumeRelativeStandardError, dtype=np.float64)
    defined_relative_standard_error = relative_standard_error[np.isfinite(relative_standard_error)]
    assert defined_relative_standard_error.size > 0
    assert np.max(defined_relative_standard_error) < 0.15

    observed_integrals = []
    for step in metadata["observable"]["stepNumbers"]:
        tet4_path = tet4_dir / f"laserPumpCladding_{step:03d}.vtk"
        tet4_points, tet4_cells, _tet4_cell_types, _tet4_point_data, tet4_cell_data, _tet4_fields = _parseVtk(tet4_path)
        observed_integrals.append(_tet_cell_integral(tet4_points, tet4_cells, tet4_cell_data["volumePhiASE"]))

    reference_integrals = _wedge_point_integrals(
        laserPumpCladdingReference["points"],
        laserPumpCladdingReference["cells"],
        laserPumpCladdingReference["phiASE"],
    )
    np.testing.assert_allclose(np.asarray(observed_integrals), reference_integrals, rtol=rtol, atol=0.0)


def testLaserPumpCladdingMediumUsesPrimitiveReflectivityShape():
    medium = laserPumpCladding.laserPumpCladdingMedium()
    cell_count = medium.topology.numberOfCells

    assert medium.get("reflectivities").expectedShape == (cell_count, 2)
    reflectivities = medium.get("reflectivities").value.reshape((cell_count, 2), order="F")
    assert reflectivities.shape == (cell_count, 2)


@pytest.fixture
def fakeCompiledSnapshots(monkeypatch):
    calls = []

    def fake_run_simulation(simulation, *, steps, pumpSteps=None, transport=None):
        calls.append(
            {
                "simulation": simulation,
                "steps": steps,
                "pumpSteps": pumpSteps,
                "transport": transport,
                "phiASE": simulation.phiASE,
            }
        )
        point_shape = simulation.gainMedium.get("betaCells").expectedShape
        volume_shape = simulation.gainMedium.get("betaVolume").expectedShape
        states = []
        for step in range(1, steps + 1):
            pump_active = pumpSteps is None or step <= pumpSteps
            states.append(
                SimpleNamespace(
                    step=step,
                    time=step * simulation.timeStep,
                    betaCells=np.full(point_shape, 0.05 * step, dtype=np.float64),
                    betaVolume=np.full(volume_shape, 0.025 * step, dtype=np.float64),
                    phiAse=np.ones(point_shape, dtype=np.float64),
                    dndtAse=np.zeros(point_shape, dtype=np.float64),
                    dndtPump=(
                        np.ones(point_shape, dtype=np.float64)
                        if pump_active
                        else np.zeros(point_shape, dtype=np.float64)
                    ),
                    aseResult=object(),
                )
            )
        return states

    monkeypatch.setattr(transport, "runSimulation", fake_run_simulation)
    return calls


def testLaserPumpCladdingExampleWritesVtkFromCompiledSnapshots(
    monkeypatch,
    tmp_path,
    smallGainMedium,
    fakeCompiledSnapshots,
):
    monkeypatch.setattr(laserPumpCladding, "laserPumpCladdingMedium", lambda **kwargs: smallGainMedium)

    state = laserPumpCladding.runExample(timeSlices=2, pumpSteps=1, vtkOutputDir=tmp_path)

    first = tmp_path / "laserPumpCladding_001.vtk"
    second = tmp_path / "laserPumpCladding_002.vtk"
    assert first.is_file()
    assert second.is_file()
    assert state.step == 2
    scalars = _vtkScalarNames(second)
    assert {"betaCells", "phiASE", "dndtAse", "dndtPump", "cladAbs"}.issubset(scalars)
    assert fakeCompiledSnapshots[-1]["pumpSteps"] == 1


def testLaserPumpCladdingExampleWiresOpenPmdBackend(
    monkeypatch,
    tmp_path,
    smallGainMedium,
    fakeCompiledSnapshots,
):
    monkeypatch.setattr(laserPumpCladding, "laserPumpCladdingMedium", lambda **kwargs: smallGainMedium)

    state = laserPumpCladding.runExample(
        timeSlices=2,
        pumpSteps=1,
        vtkOutputDir=tmp_path,
        openpmdBackend="hdf5",
    )

    assert state.step == 2
    assert fakeCompiledSnapshots[-1]["transport"] == "hdf5"
    assert fakeCompiledSnapshots[-1]["phiASE"].openpmdBackend == "hdf5"
    assert np.allclose(state.dndtPump, 0.0)


def testLaserPumpCladdingExampleCanDisableAse(
    monkeypatch,
    tmp_path,
    smallGainMedium,
    fakeCompiledSnapshots,
):
    monkeypatch.setattr(laserPumpCladding, "laserPumpCladdingMedium", lambda **kwargs: smallGainMedium)

    state = laserPumpCladding.runExample(
        timeSlices=1,
        pumpSteps=1,
        vtkOutputDir=tmp_path,
        enableASE=False,
    )

    assert state.step == 1
    assert fakeCompiledSnapshots[-1]["simulation"].enableASE is False


def testLaserPumpCladdingCliAcceptsDisableAse(monkeypatch, tmp_path):
    calls = []

    def fake_run_example(*args, **kwargs):
        calls.append({"args": args, "kwargs": kwargs})
        return SimpleNamespace(phiAse=np.zeros((2, 3)), betaCells=np.zeros((2, 3)))

    monkeypatch.setattr(laserPumpCladding, "runExample", fake_run_example)

    laserPumpCladding.main(
        [
            "--disable-ase",
            "--timeSteps",
            "1",
            "--pumpSteps",
            "1",
            "--vtk-output-dir",
            str(tmp_path),
        ]
    )

    assert calls[-1]["kwargs"]["enableASE"] is False
