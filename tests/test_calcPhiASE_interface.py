from HASEonGPU import calcPhiASE
import numpy as np
import copy
import pytest
from pathlib import Path
def to_flat(arr, width=None, dtype=None):
    arr = np.asarray(arr, dtype=dtype)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2:
        if width is not None and arr.shape[1] != width:
            raise ValueError(f"Expected shape (N, {width}), got {arr.shape}")
        return arr.reshape(-1, order="F")
    raise ValueError(f"Unsupported shape {arr.shape}")


def deep_to_list(x):
    if isinstance(x, dict):
        return {k: deep_to_list(v) for k, v in x.items()}
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x
def createDummyData():
    mesh_z = 3  # number of point levels

    p = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
    ], dtype=np.float64)

    t_int = np.array([
        [0, 1, 2],
        [0, 2, 3],
    ], dtype=np.uint32)

    beta_cell = np.array([
        [1.0, 1.1, 1.2],
        [1.3, 1.4, 1.5],
        [1.6, 1.7, 1.8],
        [1.9, 2.0, 2.1],
    ], dtype=np.float64)

    beta_vol = np.array([
        [0.10, 0.20],
        [0.30, 0.40],
    ], dtype=np.float64)

    clad_int = np.array([0, 1], dtype=np.uint32)

    refractiveIndices = np.array([1.8, 1.0, 1.8, 1.0], dtype=np.float32)
    reflectivities = np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float32)

    normals_x = np.array([
        [0.0, 1.0, -1.0],
        [0.0, 1.0, -1.0],
    ], dtype=np.float64)

    normals_y = np.array([
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ], dtype=np.float64)

    ordered_int = np.array([
        [-1, -1, 1],
        [0, -1, -1],
    ], dtype=np.int32)

    surface = np.array([0.5, 0.5], dtype=np.float32)

    x_center = np.array([2.0 / 3.0, 1.0 / 3.0], dtype=np.float64)
    y_center = np.array([1.0 / 3.0, 2.0 / 3.0], dtype=np.float64)

    normals_p = np.array([
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
        "useReflections": True,
        "minRaysPerSample": 10,
        "maxRaysPerSample": 20,
        "mseThreshold": 1e-3,
        "repetitions": 3,
        "N_tot": 5.0,
        "z_mesh": 0.25,
        "mesh_z": mesh_z,
        "deviceMode": "cpu",
        "parallelMode": "threaded",
        "maxGPUs": 1,
        "nPerNode": 1
    }

    arrays = {
        "p": p,
        "t_int": t_int,
        "beta_cell": beta_cell,
        "beta_vol": beta_vol,
        "clad_int": clad_int,
        "refractiveIndices": refractiveIndices,
        "reflectivities": reflectivities,
        "normals_x": normals_x,
        "normals_y": normals_y,
        "ordered_int": ordered_int,
        "surface": surface,
        "x_center": x_center,
        "y_center": y_center,
        "normals_p": normals_p,
        "forbidden": forbidden,
        "laser": laser,
        "crystal": crystal,
    }

    return {**arrays, **scalars}

def constructCartesianProcuct(layoutOpts,containerOpts,parallelModes):
    return [(layout, container, parallel)
            for layout in layoutOpts
            for container in containerOpts
            for parallel in parallelModes]

def convert_case(data, layout, container):
    case = copy.deepcopy(data)
    mesh_z = int(case["mesh_z"])

    if layout == "flat":
        case["p"] = to_flat(case["p"], width=2, dtype=np.float64)
        case["t_int"] = to_flat(case["t_int"], width=3, dtype=np.uint32)
        case["beta_cell"] = to_flat(case["beta_cell"], dtype=np.float64)
        case["beta_vol"] = to_flat(case["beta_vol"], width=mesh_z - 1, dtype=np.float64)
        case["clad_int"] = to_flat(case["clad_int"], dtype=np.uint32)
        case["refractiveIndices"] = to_flat(case["refractiveIndices"], dtype=np.float32)
        case["reflectivities"] = to_flat(case["reflectivities"], dtype=np.float32)
        case["normals_x"] = to_flat(case["normals_x"], width=3, dtype=np.float64)
        case["normals_y"] = to_flat(case["normals_y"], width=3, dtype=np.float64)
        case["ordered_int"] = to_flat(case["ordered_int"], width=3, dtype=np.int32)
        case["surface"] = to_flat(case["surface"], dtype=np.float32)
        case["x_center"] = to_flat(case["x_center"], dtype=np.float64)
        case["y_center"] = to_flat(case["y_center"], dtype=np.float64)
        case["normals_p"] = to_flat(case["normals_p"], width=3, dtype=np.uint32)
        case["forbidden"] = to_flat(case["forbidden"], width=3, dtype=np.int32)
        case["laser"]["l_abs"] = to_flat(case["laser"]["l_abs"], dtype=np.float64)
        case["laser"]["l_ems"] = to_flat(case["laser"]["l_ems"], dtype=np.float64)
        case["laser"]["s_abs"] = to_flat(case["laser"]["s_abs"], dtype=np.float64)
        case["laser"]["s_ems"] = to_flat(case["laser"]["s_ems"], dtype=np.float64)
    elif layout != "matrix":
        raise ValueError(f"Unsupported layout: {layout}")

    if container == "list":
        case = deep_to_list(case)
    elif container != "ndarray":
        raise ValueError(f"Unsupported container: {container}")

    return case
def run_withOptions(layout, container, parallel, data):
    case = convert_case(data, layout, container)

    phi_ASE, mse_values, N_rays = calcPhiASE(
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
        case["deviceMode"],
        parallel,
        case["maxGPUs"],
        case["nPerNode"]
    )

    return phi_ASE, mse_values, N_rays
def normalize_output(values, numberOfPoints, numberOfLevels):
    arr = np.asarray(values)
    if arr.ndim == 1:
        return arr.reshape((numberOfPoints, numberOfLevels), order="F")
    return arr

@pytest.fixture(scope="module")
def dummy_data():
    return createDummyData()


@pytest.fixture(scope="module")
def reference_outputs(dummy_data):
    phi_ASE, mse_values, N_rays = run_withOptions(
        layout="matrix",
        container="ndarray",
        parallel="threaded",
        data=dummy_data,
    )

    n_points = np.asarray(dummy_data["p"]).shape[0]
    n_levels = int(dummy_data["mesh_z"])

    return (
        normalize_output(phi_ASE, n_points, n_levels),
        normalize_output(mse_values, n_points, n_levels),
        normalize_output(N_rays, n_points, n_levels),
    )


@pytest.mark.parametrize("layout", ["matrix", "flat"])
@pytest.mark.parametrize("container", ["ndarray", "list"])
@pytest.mark.parametrize("parallel", ["threaded", "debugFileIOPath"])
def test_calcPhiASE_cases(layout, container, parallel, dummy_data, reference_outputs):
    phi_ASE, mse_values, N_rays = run_withOptions(layout, container, parallel, dummy_data)

    n_points = np.asarray(dummy_data["p"]).shape[0]
    n_levels = int(dummy_data["mesh_z"])

    phi_norm = normalize_output(phi_ASE, n_points, n_levels)
    mse_norm = normalize_output(mse_values, n_points, n_levels)
    rays_norm = normalize_output(N_rays, n_points, n_levels)

    reference_phi, reference_mse, reference_rays = reference_outputs

    assert np.allclose(phi_norm, reference_phi)
    assert np.allclose(mse_norm, reference_mse)
    assert np.array_equal(rays_norm, reference_rays)

    if container == "ndarray":
        assert isinstance(phi_ASE, np.ndarray)
        assert isinstance(mse_values, np.ndarray)
        assert isinstance(N_rays, np.ndarray)
    else:
        assert isinstance(phi_ASE, list)
        assert isinstance(mse_values, list)
        assert isinstance(N_rays, list)

    if layout == "matrix":
        if container == "ndarray":
            assert phi_ASE.ndim == 2
            assert mse_values.ndim == 2
            assert N_rays.ndim == 2
        else:
            assert isinstance(phi_ASE[0], list)
            assert isinstance(mse_values[0], list)
            assert isinstance(N_rays[0], list)
    else:
        if container == "ndarray":
            assert phi_ASE.ndim == 1
            assert mse_values.ndim == 1
            assert N_rays.ndim == 1
        else:
            assert not isinstance(phi_ASE[0], list)
            assert not isinstance(mse_values[0], list)
            assert not isinstance(N_rays[0], list)
