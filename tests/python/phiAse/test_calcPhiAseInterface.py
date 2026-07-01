# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


from HASEonGPU import AlpakaBackends
from pyInclude import calcPhiASE
import numpy as np
import copy
import importlib
import pytest
from types import SimpleNamespace

calcPhiASEModule = importlib.import_module("pyInclude.calcPhiASE")


alpakaBackends = AlpakaBackends.all()
REGRESSION_RNG_SEED = 5489


def testFindCalcPhiASENearModuleFallsBackToInstalledBindingsPath(monkeypatch, tmp_path):
    localBindings = tmp_path / "repo" / "HASEonGPU_Bindings"
    installedBindings = tmp_path / "site-packages" / "HASEonGPU_Bindings"
    localBindings.mkdir(parents=True)
    installedBindings.mkdir(parents=True)

    expected = installedBindings / "calcPhiASE"
    expected.write_text("", encoding="utf-8")

    fakeBindings = SimpleNamespace(
        __file__=str(localBindings / "__init__.py"),
        __path__=[str(localBindings), str(installedBindings)],
    )
    monkeypatch.setattr(calcPhiASEModule, "HASEonGPU_Bindings", fakeBindings)

    assert calcPhiASEModule.findCalcPhiASENearModule().resolve() == expected.resolve()


def toFlat(arr, width=None, dtype=None):
    arr = np.asarray(arr, dtype=dtype)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2:
        if width is not None and arr.shape[1] != width:
            raise ValueError(f"Expected shape (N, {width}), got {arr.shape}")
        return arr.reshape(-1, order="F")
    raise ValueError(f"Unsupported shape {arr.shape}")


def deepToList(x):
    if isinstance(x, dict):
        return {k: deepToList(v) for k, v in x.items()}
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x


def createDummyData(phiAseConfig):
    meshZ = 3  # number of point levels
    experiment = phiAseConfig["experiment"]
    compute = phiAseConfig["compute"]

    p = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
    ], dtype=np.float64)

    trianglePointIndices = np.array([
        [0, 1, 2],
        [0, 2, 3],
    ], dtype=np.uint32)

    betaCell = np.array([
        [1.0, 1.1, 1.2],
        [1.3, 1.4, 1.5],
        [1.6, 1.7, 1.8],
        [1.9, 2.0, 2.1],
    ], dtype=np.float64)

    betaVolume = np.array([
        [0.10, 0.20],
        [0.30, 0.40],
    ], dtype=np.float64)

    claddingCellTypes = np.array([0, 1], dtype=np.uint32)

    refractiveIndices = np.array([1.8, 1.0, 1.8, 1.0], dtype=np.float32)
    reflectivities = np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float32)

    normalsX = np.array([
        [0.0, 1.0, -1.0],
        [0.0, 1.0, -1.0],
    ], dtype=np.float64)

    normalsY = np.array([
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ], dtype=np.float64)

    orderedInt = np.array([
        [-1, -1, 1],
        [0, -1, -1],
    ], dtype=np.int32)

    surface = np.array([0.5, 0.5], dtype=np.float32)

    xCenter = np.array([2.0 / 3.0, 1.0 / 3.0], dtype=np.float64)
    yCenter = np.array([1.0 / 3.0, 2.0 / 3.0], dtype=np.float64)

    normalPoints = np.array([
        [0, 1, 2],
        [0, 2, 3],
    ], dtype=np.uint32)

    forbidden = np.array([
        [-1, -1, 0],
        [2, -1, -1],
    ], dtype=np.int32)

    laser = {
        "l_abs": np.array([900.0, 910.0], dtype=np.float64),
        "l_ems": np.array([1000.0, 1010.0], dtype=np.float64),
        "s_abs": np.array([0.01, 0.02], dtype=np.float64),
        "s_ems": np.array([0.03, 0.04], dtype=np.float64),
        "l_res": 2,
    }

    crystal = {"tfluo": 1.23}

    scalars = {
        "clad_number": 2,
        "clad_abs": 0.05,
        "useReflections": bool(experiment["useReflections"]),
        "minRaysPerSample": int(experiment["minRaysPerSample"]),
        "maxRaysPerSample": int(experiment["maxRaysPerSample"]),
        "mseThreshold": float(experiment["mseThreshold"]),
        "repetitions": int(compute["repetitions"]),
        "adaptiveSteps": int(compute["adaptiveSteps"]),
        "N_tot": 5.0,
        "z_mesh": 0.25,
        "mesh_z": meshZ,
        "backend": str(compute["backend"]),
        "parallelMode": str(compute["parallelMode"]),
        "numDevices": int(compute["numDevices"]),
        "nPerNode": 1
    }

    arrays = {
        "p": p,
        "t_int": trianglePointIndices,
        "beta_cell": betaCell,
        "beta_vol": betaVolume,
        "clad_int": claddingCellTypes,
        "refractiveIndices": refractiveIndices,
        "reflectivities": reflectivities,
        "normals_x": normalsX,
        "normals_y": normalsY,
        "ordered_int": orderedInt,
        "surface": surface,
        "x_center": xCenter,
        "y_center": yCenter,
        "normals_p": normalPoints,
        "forbidden": forbidden,
        "laser": laser,
        "crystal": crystal,
    }

    return {**arrays, **scalars}


def constructCartesianProduct(layoutOptions, containerOptions, parallelModes):
    return [(layout, container, parallel)
            for layout in layoutOptions
            for container in containerOptions
            for parallel in parallelModes]


def convertCase(data, layout, container):
    case = copy.deepcopy(data)
    meshZ = int(case["mesh_z"])

    if layout == "flat":
        case["p"] = toFlat(case["p"], width=2, dtype=np.float64)
        case["t_int"] = toFlat(case["t_int"], width=3, dtype=np.uint32)
        case["beta_cell"] = toFlat(case["beta_cell"], dtype=np.float64)
        case["beta_vol"] = toFlat(case["beta_vol"], width=meshZ - 1, dtype=np.float64)
        case["clad_int"] = toFlat(case["clad_int"], dtype=np.uint32)
        case["refractiveIndices"] = toFlat(case["refractiveIndices"], dtype=np.float32)
        case["reflectivities"] = toFlat(case["reflectivities"], dtype=np.float32)
        case["normals_x"] = toFlat(case["normals_x"], width=3, dtype=np.float64)
        case["normals_y"] = toFlat(case["normals_y"], width=3, dtype=np.float64)
        case["ordered_int"] = toFlat(case["ordered_int"], width=3, dtype=np.int32)
        case["surface"] = toFlat(case["surface"], dtype=np.float32)
        case["x_center"] = toFlat(case["x_center"], dtype=np.float64)
        case["y_center"] = toFlat(case["y_center"], dtype=np.float64)
        case["normals_p"] = toFlat(case["normals_p"], width=3, dtype=np.uint32)
        case["forbidden"] = toFlat(case["forbidden"], width=3, dtype=np.int32)
        case["laser"]["l_abs"] = toFlat(case["laser"]["l_abs"], dtype=np.float64)
        case["laser"]["l_ems"] = toFlat(case["laser"]["l_ems"], dtype=np.float64)
        case["laser"]["s_abs"] = toFlat(case["laser"]["s_abs"], dtype=np.float64)
        case["laser"]["s_ems"] = toFlat(case["laser"]["s_ems"], dtype=np.float64)
    elif layout != "matrix":
        raise ValueError(f"Unsupported layout: {layout}")

    if container == "list":
        case = deepToList(case)
    elif container != "ndarray":
        raise ValueError(f"Unsupported container: {container}")

    return case


def runWithOptions(layout, container, parallel, backend, data, rngSeed=None):
    case = convertCase(data, layout, container)
    case["backend"] = backend

    phiAseValues, mseValues, rayCounts = calcPhiASE(
        case["p"],
        case["t_int"],
        case["beta_cell"],
        case["beta_vol"],
        case["clad_int"],
        case["clad_number"],
        case["clad_abs"],
        case["useReflections"],
        case["refractiveIndices"],
        case["reflectivities"],
        case["normals_x"],
        case["normals_y"],
        case["ordered_int"],
        case["surface"],
        case["x_center"],
        case["y_center"],
        case["normals_p"],
        case["forbidden"],
        case["minRaysPerSample"],
        case["maxRaysPerSample"],
        case["mseThreshold"],
        case["repetitions"],
        case["N_tot"],
        case["z_mesh"],
        case["laser"],
        case["crystal"],
        case["mesh_z"],
        case["backend"],
        parallel,
        case["numDevices"],
        case["adaptiveSteps"],
        case["nPerNode"],
        rngSeed=rngSeed,
    )

    return phiAseValues, mseValues, rayCounts


def normalizeOutput(values, numberOfPoints, numberOfLevels):
    arr = np.asarray(values)
    if arr.ndim == 1:
        return arr.reshape((numberOfPoints, numberOfLevels), order="F")
    return arr


@pytest.fixture(scope="module")
def dummyData(phiAseTestConfig):
    return createDummyData(phiAseTestConfig)


@pytest.fixture
def referenceOutputs(dummyData, backend):
    phiAseValues, mseValues, rayCounts = runWithOptions(
        layout="matrix",
        container="ndarray",
        parallel="single",
        backend=backend,
        data=dummyData,
        rngSeed=REGRESSION_RNG_SEED,
    )

    nPoints = np.asarray(dummyData["p"]).shape[0]
    nLevels = int(dummyData["mesh_z"])

    return (
        normalizeOutput(phiAseValues, nPoints, nLevels),
        normalizeOutput(mseValues, nPoints, nLevels),
        normalizeOutput(rayCounts, nPoints, nLevels),
    )


@pytest.mark.parametrize("backend", alpakaBackends)
@pytest.mark.parametrize("layout", ["matrix", "flat"])
@pytest.mark.parametrize("container", ["ndarray", "list"])
@pytest.mark.parametrize("parallel", ["single", "debugFileIOPath"])
def testCalcPhiAseCases(backend, layout, container, parallel, dummyData, referenceOutputs):
    phiAseValues, mseValues, rayCounts = runWithOptions(
        layout,
        container,
        parallel,
        backend,
        dummyData,
        rngSeed=REGRESSION_RNG_SEED,
    )

    nPoints = np.asarray(dummyData["p"]).shape[0]
    nLevels = int(dummyData["mesh_z"])

    phiNorm = normalizeOutput(phiAseValues, nPoints, nLevels)
    mseNorm = normalizeOutput(mseValues, nPoints, nLevels)
    raysNorm = normalizeOutput(rayCounts, nPoints, nLevels)

    referencePhi, referenceMse, referenceRays = referenceOutputs

    assert np.allclose(phiNorm, referencePhi)
    assert np.allclose(mseNorm, referenceMse)
    assert np.array_equal(raysNorm, referenceRays)

    if container == "ndarray":
        assert isinstance(phiAseValues, np.ndarray)
        assert isinstance(mseValues, np.ndarray)
        assert isinstance(rayCounts, np.ndarray)
    else:
        assert isinstance(phiAseValues, list)
        assert isinstance(mseValues, list)
        assert isinstance(rayCounts, list)

    if layout == "matrix":
        if container == "ndarray":
            assert phiAseValues.ndim == 2
            assert mseValues.ndim == 2
            assert rayCounts.ndim == 2
        else:
            assert isinstance(phiAseValues[0], list)
            assert isinstance(mseValues[0], list)
            assert isinstance(rayCounts[0], list)
    else:
        if container == "ndarray":
            assert phiAseValues.ndim == 1
            assert mseValues.ndim == 1
            assert rayCounts.ndim == 1
        else:
            assert not isinstance(phiAseValues[0], list)
            assert not isinstance(mseValues[0], list)
            assert not isinstance(rayCounts[0], list)


@pytest.mark.parametrize("backend", alpakaBackends)
def testDirectRepeat(backend, dummyData):
    firstPhiAseValues, firstMseValues, firstRayCounts = runWithOptions(
        layout="matrix",
        container="ndarray",
        parallel="single",
        backend=backend,
        data=dummyData,
        rngSeed=REGRESSION_RNG_SEED,
    )
    secondPhiAseValues, secondMseValues, secondRayCounts = runWithOptions(
        layout="matrix",
        container="ndarray",
        parallel="single",
        backend=backend,
        data=dummyData,
        rngSeed=REGRESSION_RNG_SEED,
    )

    assert np.allclose(firstPhiAseValues, secondPhiAseValues)
    assert np.allclose(firstMseValues, secondMseValues)
    assert np.array_equal(firstRayCounts, secondRayCounts)
