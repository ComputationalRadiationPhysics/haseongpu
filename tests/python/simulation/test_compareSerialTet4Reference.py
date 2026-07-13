# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

repoRoot = Path(__file__).resolve().parents[3]


from HASEonGPU import AlpakaBackends, GainMedium, PhiASE, SpectralDecomposition
from pyInclude.geometry.vtk import _parseVtk

_compare_serial_helper_path = repoRoot / "utils" / "compare_serial_wedge_projection.py"
_compare_serial_helper_spec = importlib.util.spec_from_file_location(
    "_hase_compare_serial_wedge_projection",
    _compare_serial_helper_path,
)
if _compare_serial_helper_spec is None or _compare_serial_helper_spec.loader is None:
    raise ImportError(f"cannot load compareSerial helper from {_compare_serial_helper_path}")
_compare_serial_helper = importlib.util.module_from_spec(_compare_serial_helper_spec)
sys.modules[_compare_serial_helper_spec.name] = _compare_serial_helper
_compare_serial_helper_spec.loader.exec_module(_compare_serial_helper)
projectionFromTet4Medium = _compare_serial_helper.projectionFromTet4Medium
serialPhiAsePointFields = _compare_serial_helper.serialPhiAsePointFields
tet4ResultToLegacyPointValues = _compare_serial_helper.tet4ResultToLegacyPointValues
writeWedgeComparisonArtifacts = _compare_serial_helper.writeWedgeComparisonArtifacts
legacyWedgePointIntegral = getattr(_compare_serial_helper, "legacyWedgePointIntegral", None)
tet4CellVolumeIntegral = getattr(_compare_serial_helper, "tet4CellVolumeIntegral", None)

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
    assert metadata["parameters"]["experiment"] == {
        "minRays": 1_000_000,
        "maxRays": 1_000_000,
        "mseThreshold": 0.1,
        "reflection": False,
        "monochromatic": True,
    }
    assert metadata["parameters"]["compute"] == {
        "parallelMode": "single",
        "repetitions": 1,
        "adaptiveSteps": 1,
        "numDevices": 1,
        "minSampleIndex": 0,
        "maxSampleIndex": "full point buffer",
    }
    assert metadata["parameters"]["spectralResolution"] == 1


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
    spectral_length = len(spectra["lambdaE"])
    assert spectral_length == len(spectra["lambdaA"])
    assert spectral_length == len(spectra["sigmaA"])
    assert spectral_length == len(spectra["sigmaE"])
    assert metadata["parameters"]["experiment"]["monochromatic"] is True
    assert metadata["parameters"]["spectralResolution"] == 1

    # The legacy parser reduced monochromatic input tables to sigmaA.front()
    # and sigmaE.front() before constructing ExperimentParameters.  Sending
    # the complete tables would make the current forward kernel sample every
    # stored wavelength even when the run is marked monochromatic.
    return SpectralDecomposition(
        wavelengthsAbsorption=spectra["lambdaA"][:1],
        crossSectionAbsorption=spectra["sigmaA"][:1],
        wavelengthsEmission=spectra["lambdaE"][:1],
        crossSectionEmission=spectra["sigmaE"][:1],
        resolution=1,
    )


def testLegacyMonochromaticSpectralPropertiesUseFirstStoredPair(compareSerialReference):
    metadata = compareSerialReference["metadata"]
    for name in compareSerialReference["datasets"]:
        spectra = metadata["spectra"][name]
        properties = _legacy_spectral_properties(metadata, name)

        assert properties.resolution == 1
        np.testing.assert_array_equal(properties.wavelengthsAbsorption, spectra["lambdaA"][:1])
        np.testing.assert_array_equal(properties.crossSectionAbsorption, spectra["sigmaA"][:1])
        np.testing.assert_array_equal(properties.wavelengthsEmission, spectra["lambdaE"][:1])
        np.testing.assert_array_equal(properties.crossSectionEmission, spectra["sigmaE"][:1])


@pytest.mark.integration
def testCurrentTet4ForwardPhiAseVolumeIntegralMatchesCompareSerialReference(
    compareSerialReference,
    openPmdRuntimeBackend,
):
    if legacyWedgePointIntegral is None or tet4CellVolumeIntegral is None:
        pytest.skip("volume-centered Tet4 comparison helpers are not enabled in this trimmed runtime")

    metadata = compareSerialReference["metadata"]
    experiment = metadata["parameters"]["experiment"]
    compute = metadata["parameters"]["compute"]
    ray_count = int(os.environ.get("HASE_COMPARE_SERIAL_FORWARD_RAY_COUNT", experiment["maxRays"]))
    rtol = float(os.environ.get("HASE_COMPARE_SERIAL_INTEGRAL_RTOL", "0.05"))
    backend = os.environ.get("HASE_COMPARE_SERIAL_BACKEND", _default_backend())

    for name, dataset in compareSerialReference["datasets"].items():
        if name == "cuboid":
            relative_standard_error_threshold = 0.14
        elif name == "cylindrical":
            relative_standard_error_threshold = 0.10
        else:
            raise ValueError(f"no calibrated RSE threshold for compareSerial dataset '{name}'")
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / f"{name}.vtk")
        projection = projectionFromTet4Medium(medium)
        phi_ase = PhiASE(
            spectralProperties=_legacy_spectral_properties(metadata, name),
            minRaysPerSample=int(experiment["minRays"]),
            maxRaysPerSample=int(experiment["maxRays"]),
            forwardRayCount=ray_count,
            repetitions=int(compute["repetitions"]),
            adaptiveSteps=int(compute["adaptiveSteps"]),
            relativeStandardErrorThreshold=relative_standard_error_threshold,
            useReflections=bool(experiment["reflection"]),
            backend=backend,
            openpmdBackend=openPmdRuntimeBackend,
            parallelMode=str(compute["parallelMode"]),
            numDevices=int(compute["numDevices"]),
            monochromatic=bool(experiment["monochromatic"]),
            rngSeed=int(metadata["random"]["serialMt19937Seed"]),
        )
        phi_ase.run(gainMedium=medium)
        np.testing.assert_array_less(
            np.asarray(phi_ase.getResults().relativeStandardError, dtype=np.float64),
            relative_standard_error_threshold,
        )
        current_integral = tet4CellVolumeIntegral(medium, phi_ase.getResults().phiAse)
        legacy_integral = legacyWedgePointIntegral(dataset["phiASE"], projection)

        np.testing.assert_allclose(current_integral, legacy_integral, rtol=rtol, atol=0.0)
