# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import os
from pathlib import Path

import numpy as np
import pytest

from HASEonGPU import AlpakaBackends, GainMedium, PhiASE, SpectralDecomposition
from pyInclude.geometry.vtk import _parseVtk
from utils.compare_serial_wedge_projection import (
    projectionFromTet4Medium,
    serialPhiAsePointFields,
    tet4ResultToLegacyPointValues,
    writeWedgeComparisonArtifacts,
)


repoRoot = Path(__file__).resolve().parents[3]
REFERENCE_PATH = repoRoot / "tests" / "data" / "compareSerial" / "phiase_reference.npz"


@pytest.fixture(scope="session")
def compareSerialReference():
    if not REFERENCE_PATH.is_file():
        pytest.skip(f"missing compareSerial reference data: {REFERENCE_PATH}")
    with np.load(REFERENCE_PATH, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata"].item()))
        return {
            "metadata": metadata,
            "datasets": {
                "cuboid": {
                    "phiASE": np.asarray(data["cuboid_phiASE"], dtype=np.float64),
                    "dndtAse": np.asarray(data["cuboid_dndtAse"], dtype=np.float64),
                },
                "cylindrical": {
                    "phiASE": np.asarray(data["cylindrical_phiASE"], dtype=np.float64),
                    "dndtAse": np.asarray(data["cylindrical_dndtAse"], dtype=np.float64),
                },
            },
        }


def testCompareSerialReferenceDataDocumentsGenerator(compareSerialReference):
    metadata = compareSerialReference["metadata"]
    assert metadata["generator"]["commit"] == "469c87770ed13796f2e82385bcf83528e8aeaf1b"
    assert metadata["random"]["serialMt19937Seed"] == 5489
    assert metadata["observable"]["entity"] == "legacy sample point"
    assert metadata["observable"]["coverage"] == "full point buffer"
    assert metadata["datasets"] == ["cuboid", "cylindrical"]


def testCompareSerialReferenceStoresFullPointBuffers(compareSerialReference):
    for name, dataset in compareSerialReference["datasets"].items():
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / f"{name}.vtk")
        expected_size = medium.numberOfPoints * medium.numberOfLevels
        assert dataset["phiASE"].shape == (expected_size,)
        assert dataset["dndtAse"].shape == (expected_size,)
        assert np.isfinite(dataset["phiASE"]).all()
        assert np.isfinite(dataset["dndtAse"]).all()


def testConvertedTet4FixturesRetainInvertibleLegacyPrismGeometry(compareSerialReference):
    for name in compareSerialReference["datasets"]:
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / f"{name}.vtk")
        projection = projectionFromTet4Medium(medium)
        tet_volumes_by_prism = medium.topology.cellVolumes.reshape((-1, 3))

        assert projection.topology.numberOfPrisms * 3 == medium.topology.numberOfCells
        assert projection.topology.numberOfPoints == medium.numberOfPoints
        assert projection.topology.levels == medium.numberOfLevels
        np.testing.assert_allclose(tet_volumes_by_prism.sum(axis=1), projection.prismVolumes, rtol=2.0e-5, atol=0.0)


def testConvertedTet4FixturesCanProjectBetaVolumeBackToLegacyPrisms(compareSerialReference):
    for name in compareSerialReference["datasets"]:
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / f"{name}.vtk")
        projection = projectionFromTet4Medium(medium)

        assert projection.projectedBetaVolume.size * 3 == medium.topology.numberOfCells
        np.testing.assert_allclose(projection.projectedBetaVolume, projection.originalBetaVolume, rtol=0.0, atol=0.0)


def testCompareSerialReferenceMapsToFullWedgePointField(compareSerialReference):
    for name, dataset in compareSerialReference["datasets"].items():
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / f"{name}.vtk")
        projection = projectionFromTet4Medium(medium)
        fields = serialPhiAsePointFields(dataset, projection.topology)

        assert fields["serialPhiASE"].size == projection.topology.numberOfPoints * projection.topology.levels
        np.testing.assert_allclose(fields["serialPhiASE"], dataset["phiASE"], rtol=0.0, atol=0.0)


def testCompareSerialWedgeArtifactWriterCreatesOriginalAndRoundtripVtks(compareSerialReference, tmp_path):
    for name, dataset in compareSerialReference["datasets"].items():
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / f"{name}.vtk")
        artifacts = writeWedgeComparisonArtifacts(medium, tmp_path, name, serialReference=dataset)

        for key in ("original", "roundtrip"):
            path = artifacts[key]
            points, cells, cell_types, point_data, cell_data, _fields = _parseVtk(path)
            assert path.is_file()
            assert np.all(np.asarray(cell_types, dtype=np.uint32) == 13)
            assert len(cells) == artifacts["projection"].topology.numberOfPrisms
            assert points.shape[0] == artifacts["projection"].topology.numberOfPoints * artifacts["projection"].topology.levels
            assert "betaCells" in point_data
            assert "betaVolume" in cell_data

        _points, _cells, _types, original_point_data, original_cell_data, _fields = _parseVtk(artifacts["original"])
        _points, _cells, _types, roundtrip_point_data, roundtrip_cell_data, _fields = _parseVtk(artifacts["roundtrip"])
        np.testing.assert_allclose(original_point_data["betaCells"], roundtrip_point_data["betaCells"], rtol=0.0, atol=0.0)
        np.testing.assert_allclose(original_cell_data["betaVolume"], roundtrip_cell_data["betaVolume"], rtol=0.0, atol=0.0)
        assert "serialPhiASE" in original_point_data
        assert "serialPhiASE" not in roundtrip_point_data
        np.testing.assert_allclose(
            original_point_data["serialPhiASE"],
            dataset["phiASE"],
            rtol=0.0,
            atol=0.0,
        )


def testTet4CellResultCanMapToLegacyPointBuffer(compareSerialReference):
    for name in compareSerialReference["datasets"]:
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / f"{name}.vtk")
        projection = projectionFromTet4Medium(medium)
        direct = np.arange(projection.topology.numberOfPoints * projection.topology.levels, dtype=np.float64)
        np.testing.assert_array_equal(tet4ResultToLegacyPointValues(medium, direct), direct)

        synthetic = np.arange(medium.topology.numberOfCells, dtype=np.float64)
        point_values = tet4ResultToLegacyPointValues(medium, synthetic)

        assert point_values.shape == (projection.topology.numberOfPoints * projection.topology.levels,)
        if name == "cylindrical":
            assert np.count_nonzero(~np.isfinite(point_values)) == 20
        else:
            assert np.isfinite(point_values).all()


def _default_backend():
    backends = AlpakaBackends.all()
    for backend in backends:
        if "CpuSerial" in backend:
            return backend
    if backends:
        return backends[0]
    pytest.skip("no Alpaka backend is available")


def _legacy_spectral_properties(metadata, name):
    spectra = metadata["spectra"][name]
    return SpectralDecomposition(
        wavelengthsAbsorption=spectra["lambdaA"],
        crossSectionAbsorption=spectra["sigmaA"],
        wavelengthsEmission=spectra["lambdaE"],
        crossSectionEmission=spectra["sigmaE"],
        resolution=int(metadata["parameters"]["spectralResolution"]),
    )


@pytest.mark.integration
def testCurrentTet4ForwardPhiAseCanBeComparedPointwiseWithCompareSerialReference(
    compareSerialReference,
    openPmdRuntimeBackend,
):
    if os.environ.get("HASE_COMPARE_SERIAL_RUN_FORWARD") != "1":
        pytest.skip("set HASE_COMPARE_SERIAL_RUN_FORWARD=1 to run the expensive pointwise physics comparison")
    if "HASE_COMPARE_SERIAL_FORWARD_RAY_LENGTH" not in os.environ:
        pytest.skip("set HASE_COMPARE_SERIAL_FORWARD_RAY_LENGTH for the forward Tet4 comparison")

    metadata = compareSerialReference["metadata"]
    ray_count = int(os.environ.get("HASE_COMPARE_SERIAL_FORWARD_RAY_COUNT", metadata["parameters"]["experiment"]["maxRays"]))
    forward_ray_length = float(os.environ["HASE_COMPARE_SERIAL_FORWARD_RAY_LENGTH"])
    rtol = float(os.environ.get("HASE_COMPARE_SERIAL_POINTWISE_RTOL", "0.35"))
    backend = os.environ.get("HASE_COMPARE_SERIAL_BACKEND", _default_backend())

    for name, dataset in compareSerialReference["datasets"].items():
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / f"{name}.vtk")
        phi_ase = PhiASE(
            spectralProperties=_legacy_spectral_properties(metadata, name),
            minRaysPerSample=ray_count,
            maxRaysPerSample=ray_count,
            forwardRayCount=ray_count,
            forwardRayLength=forward_ray_length,
            repetitions=1,
            adaptiveSteps=1,
            mseThreshold=float(metadata["parameters"]["experiment"]["mseThreshold"]),
            useReflections=False,
            backend=backend,
            openpmdBackend=openPmdRuntimeBackend,
            parallelMode="single",
            numDevices=1,
            monochromatic=True,
            rngSeed=int(metadata["random"]["serialMt19937Seed"]),
        )
        phi_ase.run(gainMedium=medium)
        current_point_phi = tet4ResultToLegacyPointValues(medium, phi_ase.getResults().phiAse)
        legacy_point_phi = dataset["phiASE"]

        assert current_point_phi.shape == legacy_point_phi.shape
        np.testing.assert_allclose(current_point_phi, legacy_point_phi, rtol=rtol, atol=0.0)
