import shutil
import warnings

import HASEonGPU_Bindings
import numpy as np
import os
from pathlib import Path
Mesh=HASEonGPU_Bindings.HostMesh

############################## clean_IO_files #################################
def clean_IO_files(TMP_FOLDER):
    if os.path.exists(TMP_FOLDER) and os.path.isdir(TMP_FOLDER):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shutil.rmtree(TMP_FOLDER)
########################### create_calcPhiASE_input ###########################
def create_calcPhiASE_input(
        packed,
        thickness,
        numberOfLevels,
        nTot,
        crystal,
        claddingNumber,
        claddingAbsorption,
        FOLDER):

    CURRENT_DIR = os.getcwd()
    os.makedirs(FOLDER, exist_ok=True)
    os.chdir(FOLDER)

    try:
        np.savetxt('points.txt', packed["p_flat"], delimiter='\n', fmt='%.50f')
        np.savetxt('triangleNormalsX.txt', packed["normals_x_flat"], delimiter='\n', fmt='%.50f')
        np.savetxt('triangleNormalsY.txt', packed["normals_y_flat"], delimiter='\n', fmt='%.50f')

        # signed integer arrays
        np.savetxt('forbiddenEdge.txt', packed["forbidden_flat"], delimiter='\n', fmt='%d')
        np.savetxt('triangleNeighbors.txt', packed["ordered_int_flat"], delimiter='\n', fmt='%d')

        # unsigned index arrays
        np.savetxt('triangleNormalPoint.txt', packed["normals_p_flat"], delimiter='\n', fmt='%d')
        np.savetxt('trianglePointIndices.txt', packed["t_int_flat"], delimiter='\n', fmt='%d')

        with open('thickness.txt', 'w') as f:
            f.write(str(thickness) + '\n')
        with open('numberOfLevels.txt', 'w') as f:
            f.write(str(numberOfLevels) + '\n')
        with open('numberOfTriangles.txt', 'w') as f:
            f.write(str(packed["numberOfTriangles"]) + '\n')
        with open('numberOfPoints.txt', 'w') as f:
            f.write(str(packed["numberOfPoints"]) + '\n')

        with open('nTot.txt', 'w') as f:
            f.write(str(float(nTot)) + '\n')

        np.savetxt('betaVolume.txt', packed["beta_vol_flat"], delimiter='\n', fmt='%.50f')
        np.savetxt('sigmaA.txt', packed["laser"]["s_abs"], delimiter='\n', fmt='%.50f')
        np.savetxt('sigmaE.txt', packed["laser"]["s_ems"], delimiter='\n', fmt='%.50f')
        np.savetxt('lambdaA.txt', packed["laser"]["l_abs"], delimiter='\n', fmt='%.50f')
        np.savetxt('lambdaE.txt', packed["laser"]["l_ems"], delimiter='\n', fmt='%.50f')

        with open('crystalTFluo.txt', 'w') as f:
            f.write(str(crystal['tfluo']) + '\n')

        np.savetxt('betaCells.txt', packed["beta_cell_flat"], delimiter='\n', fmt='%.50f')
        np.savetxt('triangleSurfaces.txt', packed["surface_flat"], delimiter='\n', fmt='%.50f')
        np.savetxt('triangleCenterX.txt', packed["x_center_flat"], delimiter='\n', fmt='%.50f')
        np.savetxt('triangleCenterY.txt', packed["y_center_flat"], delimiter='\n', fmt='%.50f')

        np.savetxt('claddingCellTypes.txt', packed["clad_int_flat"], delimiter='\n', fmt='%d')
        with open('claddingNumber.txt', 'w') as f:
            f.write(str(claddingNumber) + '\n')
        with open('claddingAbsorption.txt', 'w') as f:
            f.write(str(claddingAbsorption) + '\n')

        np.savetxt('refractiveIndices.txt', packed["refractiveIndices_flat"], delimiter='\n', fmt='%3.5f')
        np.savetxt('reflectivities.txt', packed["reflectivities_flat"], delimiter='\n', fmt='%.50f')

    finally:
        os.chdir(CURRENT_DIR)



######################### parse_calcPhiASE_output #############################
def parse_calcPhiASE_output(FOLDER, layout):
    CURRENT_DIR = os.getcwd()
    os.chdir(FOLDER)

    try:
        def _read_array(fname, dtype=float):
            with open(fname, "r") as fid:
                shape_tokens = next(fid).strip().replace(",", "").split()
                arraySize = [int(x) for x in shape_tokens]

                line = next(fid).strip().replace(",", "")
                vals = np.asarray([dtype(x) for x in line.split()], dtype=dtype)

                if layout == "matrix":
                    return np.reshape(vals, arraySize, order="F").squeeze()

                return vals

        phiASE = _read_array("phi_ASE.txt", float)
        mseValues = _read_array("mse_values.txt", float)
        raysUsedPerSample = _read_array("N_rays.txt", float)

    finally:
        os.chdir(CURRENT_DIR)

    return phiASE, mseValues, raysUsedPerSample

def find_calcPhiASE_near_module():
    so_dir = Path(HASEonGPU_Bindings.__file__).resolve().parent

    # direct sibling
    direct = so_dir / "calcPhiASE"
    if direct.is_file():
        return direct

    # one level below
    for sub in so_dir.iterdir():
        if sub.is_dir():
            candidate = sub / "calcPhiASE"
            if candidate.is_file():
                return candidate

    raise FileNotFoundError(
        f"Could not find calcPhiASE near module location: {so_dir}"
    )


################################ calcPhiASE ########################################
def calcPhiASE_mpi(
        packed,
        claddingNumber,
        claddingAbsorption,
        useReflections,
        minRaysPerSample,
        maxRaysPerSample,
        mseThreshold,
        repetitions,
        nTot,
        thickness,
        crystal,
        numberOfLevels,
        deviceMode,
        parallelMode,
        maxGPUs,
        nPerNode
):

    minSample = 0
    nP = int(packed["numberOfPoints"])
    maxSample = (int(numberOfLevels) * nP) - 1

    REFLECT = ' --reflection=1' if useReflections else ' --reflection=0'

    Prefix = ''
    if parallelMode == 'mpi' or parallelMode == 'graybat':
        Prefix = f"mpiexec -npernode {nPerNode} "
        maxGPUs = 1

    CALCPHIASE_DIR = Path(__file__).resolve().parent
    TMP_FOLDER = os.path.join(CALCPHIASE_DIR, 'input_tmp')

    create_calcPhiASE_input(
        packed,
        thickness,
        numberOfLevels,
        nTot,
        crystal,
        claddingNumber,
        claddingAbsorption,
        TMP_FOLDER
    )
    exec_path = find_calcPhiASE_near_module().resolve()
    cmd = (
            Prefix + str(exec_path)
            + f' --parallel-mode={parallelMode}'
            + f' --device-mode={deviceMode}'
            + f' --min-rays={int(minRaysPerSample)}'
            + f' --max-rays={int(maxRaysPerSample)}'
            + REFLECT
            + f' --input-path={TMP_FOLDER}'
            + f' --output-path={TMP_FOLDER}'
            + f' --min-sample-i={minSample}'
            + f' --max-sample-i={maxSample}'
            + f' --ngpus={maxGPUs}'
            + f' --repetitions={repetitions}'
            + f' --mse-threshold={mseThreshold}'
            + f' --spectral-resolution={packed["laser"]["l_res"]}'
    )
    print(cmd)
    status = os.system(cmd)

    if status != 0:
        print('This step of the raytracing computation did NOT finish successfully. Aborting.')
        exit()

    phiASE, mseValues, raysUsedPerSample = parse_calcPhiASE_output(
        TMP_FOLDER,
        packed["layout"]
    )
    print("bef delete")
    clean_IO_files(TMP_FOLDER)

    if packed["container"] == "list" and packed["layout"] != "matrix":
        return phiASE.tolist(), mseValues.tolist(), raysUsedPerSample.tolist()

    if packed["container"] == "list":
        return phiASE.tolist(), mseValues.tolist(), raysUsedPerSample.tolist()

    return phiASE, mseValues, raysUsedPerSample


def calcPhiASE(
        p,
        t_int,
        beta_cell,
        beta_vol,
        clad_int,
        clad_number,
        clad_abs,
        useReflections,
        refractiveIndices,
        reflectivities,
        normals_x,
        normals_y,
        ordered_int,
        surface,
        x_center,
        y_center,
        normals_p,
        forbidden,
        minRaysPerSample,
        maxRaysPerSample,
        mseThreshold,
        repetitions,
        N_tot,
        z_mesh,
        laser,
        crystal,
        mesh_z,
        deviceMode,
        parallelMode,
        maxGPUs,
        nPerNode
):
    def transform_inputs(
            p,
            t_int,
            beta_cell,
            beta_vol,
            clad_int,
            refractiveIndices,
            reflectivities,
            normals_x,
            normals_y,
            ordered_int,
            surface,
            x_center,
            y_center,
            normals_p,
            forbidden,
            laser,
            mesh_z,
    ):
        def is_numpy(x):
            return isinstance(x, np.ndarray)

        def is_list(x):
            return isinstance(x, list)

        def as_array(x, name):
            if not is_numpy(x) and not is_list(x):
                raise ValueError(
                    f"{name} must be either a numpy.ndarray or a list, got {type(x)}"
                )
            return np.asarray(x)

        def pack(x, name, width=None, dtype=None):
            """
            Accept either:
            - flat README/MATLAB layout
            - matrix layout (N, width)

            Returns backend-flat NumPy array.
            """
            arr = as_array(x, name)

            if arr.ndim == 1:
                flat = arr
            elif arr.ndim == 2:
                if width is not None:
                    if arr.shape[1] != width:
                        raise ValueError(
                            f"{name} must have shape (N, {width}) or be flat, got {arr.shape}"
                        )
                    flat = arr.reshape(-1, order="F")
                else:
                    flat = arr.reshape(-1)
            else:
                raise ValueError(f"{name} has unsupported shape {arr.shape}")

            if dtype is not None:
                flat = flat.astype(dtype, copy=False)

            return flat

        p_arr = as_array(p, "p")
        if p_arr.ndim == 2:
            if p_arr.shape[1] != 2:
                raise ValueError(f"p must have shape (numberOfPoints, 2), got {p_arr.shape}")
            layout = "matrix"
            numberOfPoints = int(p_arr.shape[0])
        elif p_arr.ndim == 1:
            if p_arr.size % 2 != 0:
                raise ValueError(f"flat p must have even length, got {p_arr.size}")
            layout = "flattened"
            numberOfPoints = int(p_arr.size // 2)
        else:
            raise ValueError(f"Unsupported p shape {p_arr.shape}")

        t_int_arr = as_array(t_int, "t_int")
        if t_int_arr.ndim == 2:
            if t_int_arr.shape[1] != 3:
                raise ValueError(
                    f"t_int must have shape (numberOfTriangles, 3), got {t_int_arr.shape}"
                )
            numberOfTriangles = int(t_int_arr.shape[0])
        elif t_int_arr.ndim == 1:
            if t_int_arr.size % 3 != 0:
                raise ValueError(
                    f"flat t_int length must be divisible by 3, got {t_int_arr.size}"
                )
            numberOfTriangles = int(t_int_arr.size // 3)
        else:
            raise ValueError(f"Unsupported t_int shape {t_int_arr.shape}")

        container = "ndarray" if is_numpy(p) else "list"

        return {
            "layout": layout,
            "container": container,
            "numberOfPoints": numberOfPoints,
            "numberOfTriangles": numberOfTriangles,
            "numberOfLevels": int(mesh_z),

            "p_flat": pack(p, "p", width=2, dtype=np.float64),
            "t_int_flat": pack(t_int, "t_int", width=3, dtype=np.uint32),
            "beta_cell_flat": pack(beta_cell, "beta_cell", dtype=np.float64),
            "beta_vol_flat": pack(beta_vol, "beta_vol", width=int(mesh_z) - 1, dtype=np.float64),
            "clad_int_flat": pack(clad_int, "clad_int", dtype=np.uint32),
            "refractiveIndices_flat": pack(refractiveIndices, "refractiveIndices", dtype=np.float32),
            "reflectivities_flat": pack(reflectivities, "reflectivities", dtype=np.float32),
            "normals_x_flat": pack(normals_x, "normals_x", width=3, dtype=np.float64),
            "normals_y_flat": pack(normals_y, "normals_y", width=3, dtype=np.float64),
            "ordered_int_flat": pack(ordered_int, "ordered_int", width=3, dtype=np.int32),
            "surface_flat": pack(surface, "surface", dtype=np.float32),
            "x_center_flat": pack(x_center, "x_center", dtype=np.float64),
            "y_center_flat": pack(y_center, "y_center", dtype=np.float64),
            "normals_p_flat": pack(normals_p, "normals_p", width=3, dtype=np.uint32),
            "forbidden_flat": pack(forbidden, "forbidden", width=3, dtype=np.int32),

            "laser": {
                "l_abs": pack(laser["l_abs"], 'laser["l_abs"]', dtype=np.float64),
                "l_ems": pack(laser["l_ems"], 'laser["l_ems"]', dtype=np.float64),
                "s_abs": pack(laser["s_abs"], 'laser["s_abs"]', dtype=np.float64),
                "s_ems": pack(laser["s_ems"], 'laser["s_ems"]', dtype=np.float64),
                "l_res": laser["l_res"],
            }
        }
    packed = transform_inputs(
        p,
        t_int,
        beta_cell,
        beta_vol,
        clad_int,
        refractiveIndices,
        reflectivities,
        normals_x,
        normals_y,
        ordered_int,
        surface,
        x_center,
        y_center,
        normals_p,
        forbidden,
        laser,
        mesh_z,
    )
    if parallelMode=="mpi" or parallelMode=="graybat" or parallelMode=="debugFileIOPath":
        if parallelMode=="debugFileIOPath":
            parallelMode="threaded"
        return calcPhiASE_mpi(
            packed,
            clad_number,
            clad_abs,
            useReflections,
            minRaysPerSample,
            maxRaysPerSample,
            mseThreshold,
            repetitions,
            N_tot,
            z_mesh,
            crystal,
            mesh_z,
            deviceMode,
            parallelMode,
            maxGPUs,
            nPerNode
        )

    numberOfPoints = int(packed["numberOfPoints"])
    numberOfTriangles = int(packed["numberOfTriangles"])
    numberOfLevels = int(packed["numberOfLevels"])
    thicknessOfPrism = float(z_mesh)

    experiment = HASEonGPU_Bindings.ExperimentParameters(
        minRaysPerSample=int(minRaysPerSample),
        maxRaysPerSample=int(maxRaysPerSample),
        lambdaA=packed["laser"]["l_abs"],
        lambdaE=packed["laser"]["l_ems"],
        sigmaA=packed["laser"]["s_abs"],
        sigmaE=packed["laser"]["s_ems"],
        maxSigmaA=float(np.max(packed["laser"]["s_abs"])),
        maxSigmaE=float(np.max(packed["laser"]["s_ems"])),
        mseThreshold=float(mseThreshold),
        useReflections=bool(useReflections),
        spectral=int(packed["laser"]["l_res"]),
    )

    compute = HASEonGPU_Bindings.ComputeParameters(
        maxRepetitions=int(repetitions),
        adaptiveSteps=5,
        maxGpus=int(maxGPUs),
        gpu_i=0,
        deviceMode=str(deviceMode),
        parallelMode=str(parallelMode),
        writeVtk=False,
        devices=[],
        minSampleRange=0,
        maxSampleRange=int((numberOfPoints * numberOfLevels) - 1),
    )

    host_mesh = Mesh(
        triangleIndices=packed["t_int_flat"],
        numberOfTriangles=numberOfTriangles,
        numberOfLevels=numberOfLevels,
        numberOfPoints=numberOfPoints,
        thicknessOfPrism=thicknessOfPrism,
        pointsVector=packed["p_flat"],
        xOfTriangleCenter=packed["x_center_flat"],
        yOfTriangleCenter=packed["y_center_flat"],
        positionsOfNormalVectors=packed["normals_p_flat"],
        xOfNormals=packed["normals_x_flat"],
        yOfNormals=packed["normals_y_flat"],
        forbiddenVector=packed["forbidden_flat"],
        neighborsVector=packed["ordered_int_flat"],
        surfacesVector=packed["surface_flat"],
        betaValuesVector=packed["beta_vol_flat"],
        betaCells=packed["beta_cell_flat"],
        cellTypes=packed["clad_int_flat"],
        refractiveIndices=packed["refractiveIndices_flat"],
        reflectivities=packed["reflectivities_flat"],
        nTot=float(N_tot),
        crystalTFluo=float(crystal["tfluo"]),
        claddingNumber=int(clad_number),
        claddingAbsorption=float(clad_abs),
    )

    result = HASEonGPU_Bindings.calcPhiASE(experiment, compute, host_mesh)

    if packed["container"] == "list" and packed["layout"] != "matrix":
        return list(result.phiAse), list(result.mse), list(result.totalRays)

    phi_ASE = np.asarray(result.phiAse, dtype=np.float32)
    mse_values = np.asarray(result.mse, dtype=np.float64)
    N_rays = np.asarray(result.totalRays, dtype=np.uint32)

    if packed["layout"] == "matrix":
        phi_ASE = phi_ASE.reshape((numberOfPoints, numberOfLevels), order="F")
        mse_values = mse_values.reshape((numberOfPoints, numberOfLevels), order="F")
        N_rays = N_rays.reshape((numberOfPoints, numberOfLevels), order="F")

    if packed["container"] == "list":
        return phi_ASE.tolist(), mse_values.tolist(), N_rays.tolist()

    return phi_ASE, mse_values, N_rays
