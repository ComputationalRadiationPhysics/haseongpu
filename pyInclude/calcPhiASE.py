import shutil
import warnings

import HASEonGPU_Bindings
import numpy as np
import os
Mesh=HASEonGPU_Bindings.HostMesh

############################## clean_IO_files #################################
def clean_IO_files(TMP_FOLDER):
    if os.path.exists(TMP_FOLDER) and os.path.isdir(TMP_FOLDER):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            shutil.rmtree(TMP_FOLDER)
########################### create_calcPhiASE_input ###########################
def create_calcPhiASE_input(points,
                            triangleNormalsX,
                            triangleNormalsY,
                            forbiddenEdge,
                            triangleNormalPoint,
                            triangleNeighbors,
                            trianglePointIndices,
                            thickness,
                            numberOfLevels,
                            nTot,
                            betaVolume,
                            laserParameter,
                            crystal,
                            betaCells,
                            triangleSurfaces,
                            triangleCenterX,
                            triangleCenterY,
                            claddingCellTypes,
                            claddingNumber,
                            claddingAbsorption,
                            refractiveIndices,
                            reflectivities,
                            FOLDER):

    CURRENT_DIR = os.getcwd()
    os.makedirs(FOLDER, exist_ok=True)
    os.chdir(FOLDER)

    # ensure arrays are transposed properly for saving (MATLAB/Fortran style flattening)
    points = np.transpose(points)  # (2,N)
    triangleNormalsX = np.transpose(triangleNormalsX)
    triangleNormalsY = np.transpose(triangleNormalsY)
    forbiddenEdge = np.transpose(forbiddenEdge)
    triangleNormalPoint = np.transpose(triangleNormalPoint)
    triangleNeighbors = np.transpose(triangleNeighbors)
    trianglePointIndices = np.transpose(trianglePointIndices)
    betaVolume = np.transpose(betaVolume)

    laserParameter['s_abs'] = np.transpose(laserParameter['s_abs'])
    laserParameter['s_ems'] = np.transpose(laserParameter['s_ems'])
    laserParameter['l_abs'] = np.transpose(laserParameter['l_abs'])
    laserParameter['l_ems'] = np.transpose(laserParameter['l_ems'])

    betaCells = np.transpose(betaCells)
    triangleSurfaces = np.transpose(triangleSurfaces)
    triangleCenterX = np.transpose(triangleCenterX)
    triangleCenterY = np.transpose(triangleCenterY)
    claddingCellTypes = np.transpose(claddingCellTypes)
    refractiveIndices = np.transpose(refractiveIndices)
    reflectivities = np.transpose(reflectivities)

    # save arrays as text files
    np.savetxt('points.txt', points, delimiter='\n', fmt='%.50f')
    np.savetxt('triangleNormalsX.txt', triangleNormalsX, delimiter='\n', fmt='%.50f')
    np.savetxt('triangleNormalsY.txt', triangleNormalsY, delimiter='\n', fmt='%.50f')

    # IMPORTANT TYPES:
    # forbiddenEdge and triangleNeighbors are allowed to contain -1 (parser expects that), so keep them signed.
    np.savetxt('forbiddenEdge.txt', forbiddenEdge, delimiter='\n', fmt='%d')
    np.savetxt('triangleNeighbors.txt', triangleNeighbors, delimiter='\n', fmt='%d')

    # These are unsigned indices, but saving as %d is fine as long as values are non-negative.
    np.savetxt('triangleNormalPoint.txt', triangleNormalPoint, delimiter='\n', fmt='%d')
    np.savetxt('trianglePointIndices.txt', trianglePointIndices, delimiter='\n', fmt='%d')

    with open('thickness.txt', 'w') as f:
        f.write(str(thickness) + '\n')
    with open('numberOfLevels.txt', 'w') as f:
        f.write(str(numberOfLevels) + '\n')
    with open('numberOfTriangles.txt', 'w') as f:
        f.write(str(trianglePointIndices.shape[1]) + '\n')
    with open('numberOfPoints.txt', 'w') as f:
        f.write(str(points.shape[1]) + '\n')

    with open('nTot.txt', 'w') as f:
        f.write(str(float(nTot)) + '\n')

    np.savetxt('betaVolume.txt', betaVolume, delimiter='\n', fmt='%.50f')
    np.savetxt('sigmaA.txt', laserParameter['s_abs'], delimiter='\n', fmt='%.50f')
    np.savetxt('sigmaE.txt', laserParameter['s_ems'], delimiter='\n', fmt='%.50f')
    np.savetxt('lambdaA.txt', laserParameter['l_abs'], delimiter='\n', fmt='%.50f')
    np.savetxt('lambdaE.txt', laserParameter['l_ems'], delimiter='\n', fmt='%.50f')

    with open('crystalTFluo.txt', 'w') as f:
        f.write(str(crystal['tfluo']) + '\n')

    np.savetxt('betaCells.txt', betaCells, delimiter='\n', fmt='%.50f')
    np.savetxt('triangleSurfaces.txt', triangleSurfaces, delimiter='\n', fmt='%.50f')
    np.savetxt('triangleCenterX.txt', triangleCenterX, delimiter='\n', fmt='%.50f')
    np.savetxt('triangleCenterY.txt', triangleCenterY, delimiter='\n', fmt='%.50f')

    np.savetxt('claddingCellTypes.txt', claddingCellTypes, delimiter='\n', fmt='%d')
    with open('claddingNumber.txt', 'w') as f:
        f.write(str(claddingNumber) + '\n')
    with open('claddingAbsorption.txt', 'w') as f:
        f.write(str(claddingAbsorption) + '\n')

    np.savetxt('refractiveIndices.txt', refractiveIndices, delimiter='\n', fmt='%3.5f')
    np.savetxt('reflectivities.txt', reflectivities, delimiter='\n', fmt='%.50f')

    os.chdir(CURRENT_DIR)



######################### parse_calcPhiASE_output #############################
def parse_calcPhiASE_output(FOLDER):
    import re
    CURRENT_DIR = os.getcwd()
    os.chdir(FOLDER)

    def _read_array(fname, dtype=float):
        with open(fname, "r") as fid:
            # First line: shape (may contain commas too)
            shape_tokens = next(fid).strip().replace(",", "").split()
            arraySize = [int(x) for x in shape_tokens]
            # Second line: values (may contain commas as thousands separators)
            line = next(fid).strip()
            # Remove commas inside numbers: "1,234" -> "1234"
            line = line.replace(",", "")
            # Split on whitespace
            vals = [dtype(x) for x in line.split()]
            arr = np.reshape(vals, arraySize, order="F")
            return arr

    phiASE = _read_array("phi_ASE.txt", float)
    mseValues = _read_array("mse_values.txt", float)
    raysUsedPerSample = _read_array("N_rays.txt", float)

    os.chdir(CURRENT_DIR)
    return phiASE, mseValues, raysUsedPerSample



################################ calcPhiASE ########################################
def calcPhiASE_mpi(
        points,
        trianglePointIndices,
        betaCells,
        betaVolume,
        claddingCellTypes,
        claddingNumber,
        claddingAbsorption,
        useReflections,
        refractiveIndices,
        reflectivities,
        triangleNormalsX,
        triangleNormalsY,
        triangleNeighbors,
        triangleSurfaces,
        triangleCenterX,
        triangleCenterY,
        triangleNormalPoint,
        forbiddenEdge,
        minRaysPerSample,
        maxRaysPerSample,
        mseThreshold,
        repetitions,
        nTot,
        thickness,
        laserParameter,
        crystal,
        numberOfLevels,
        deviceMode,
        parallelMode,
        maxGPUs,
        nPerNode
):


    # # -------------------------
    # # Ensure parser-safe mesh I/O by compacting/remapping points
    # # -------------------------
    # points2, tri2, tnp2, betaCells2 = _compact_mesh_for_parser(
    #     points=points,
    #     trianglePointIndices=trianglePointIndices,
    #     triangleNormalPoint=triangleNormalPoint,
    #     betaCells=betaCells
    # )
    #
    # # Replace with compacted arrays for file writing and for maxSample calculation
    # points = points2
    # trianglePointIndices = tri2
    # triangleNormalPoint = tnp2
    # betaCells = betaCells2

    minSample = 0
    nP = points.shape[0]  # points is (N,2)
    maxSample = (numberOfLevels * nP) - 1

    REFLECT = ' --reflection=1' if useReflections else ' --reflection=0'

    Prefix = ''
    if parallelMode == 'mpi' or parallelMode=='graybat':
        Prefix = f"mpiexec -npernode {nPerNode} "
        maxGPUs = 1

    CALCPHIASE_DIR = os.getcwd()
    TMP_FOLDER = os.path.join(CALCPHIASE_DIR, 'input_tmp')

    create_calcPhiASE_input(
        points,
        triangleNormalsX,
        triangleNormalsY,
        forbiddenEdge,
        triangleNormalPoint,
        triangleNeighbors,
        trianglePointIndices,
        thickness,
        numberOfLevels,
        nTot,
        betaVolume,
        laserParameter,
        crystal,
        betaCells,
        triangleSurfaces,
        triangleCenterX,
        triangleCenterY,
        claddingCellTypes,
        claddingNumber,
        claddingAbsorption,
        refractiveIndices,
        reflectivities,
        TMP_FOLDER
    )
    cmd = (
            Prefix + CALCPHIASE_DIR + '/../../build/calcPhiASE'
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
            + f' --spectral-resolution={laserParameter["l_res"]}'
    )
    print(cmd)
    status = os.system(cmd)

    if status != 0:
        print('This step of the raytracing computation did NOT finish successfully. Aborting.')
        exit()

    phiASE, mseValues, raysUsedPerSample = parse_calcPhiASE_output(TMP_FOLDER)
    print("bef delete")
    clean_IO_files(TMP_FOLDER)

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

    """
    Compute amplified spontaneous emission (ASE) flux on an extruded mesh.

    This function is a low-level Python interface to the HASEonGPU backend.
    It performs Monte Carlo ray tracing to estimate ASE photon flux based on
    mesh geometry, material properties, and gain distribution.

    ------------------------------------------------------------------------
    Input Data Requirements
    ------------------------------------------------------------------------

    All array inputs must be provided as **NumPy arrays (np.ndarray)** with
    consistent shapes.

    General rule:
        Arrays must be shaped as (N, components)

        Example:
            Points must be shaped (N_points, 2):
                [[x0, y0],
                 [x1, y1],
                 ...]

    ------------------------------------------------------------------------
    Mesh Geometry
    ------------------------------------------------------------------------

    p : np.ndarray, shape (N_points, 2)
        Coordinates of mesh points (x, y).

    t_int : np.ndarray, shape (N_triangles, 3)
        Triangle connectivity (indices into `p`).

    mesh_z : int
        Number of layers in z-direction (extrusion of the 2D mesh).

    z_mesh : float
        Physical thickness of the extruded prism.

    x_center, y_center : np.ndarray, shape (N_triangles,)
        Coordinates of triangle centers.

    surface : np.ndarray, shape (N_triangles,)
        Area of each triangle.

    ------------------------------------------------------------------------
    Mesh Topology / Connectivity
    ------------------------------------------------------------------------

    ordered_int : np.ndarray
        Neighbor indices for each triangle (topological connectivity).

    normals_p : np.ndarray
        Indices mapping normal vectors to mesh elements.

    forbidden : np.ndarray
        Mask indicating forbidden transitions (e.g. invalid neighbors).

    ------------------------------------------------------------------------
    Normals
    ------------------------------------------------------------------------

    normals_x : np.ndarray
        x-components of normal vectors.

    normals_y : np.ndarray
        y-components of normal vectors.

    ------------------------------------------------------------------------
    Material & Optical Properties
    ------------------------------------------------------------------------

    beta_cell : np.ndarray
        Gain / absorption values per cell (typically from pumping step).

    beta_vol : np.ndarray
        Volumetric gain values.

    clad_int : np.ndarray
        Cell type identifiers (e.g. core vs cladding).

    clad_number : int
        Number of cladding regions.

    clad_abs : float
        Absorption coefficient of cladding.

    refractiveIndices : np.ndarray
        Refractive indices of materials.

    reflectivities : np.ndarray
        Reflectivity values for boundaries.

    useReflections : bool
        Enable / disable reflection handling.

    ------------------------------------------------------------------------
    Physical Parameters
    ------------------------------------------------------------------------

    N_tot : float
        Total dopant / particle density.

    crystal : dict
        Crystal properties with required key:
            "tfluo" : fluorescence lifetime

    laser : dict
        Spectral and cross-section data:

            "l_abs" : absorption wavelengths
            "l_ems" : emission wavelengths
            "s_abs" : absorption cross sections
            "s_ems" : emission cross sections
            "l_res" : spectral resolution (int)

    ------------------------------------------------------------------------
    Monte Carlo / Sampling Parameters
    ------------------------------------------------------------------------

    minRaysPerSample : int
        Minimum number of rays per sampling step.

    maxRaysPerSample : int
        Maximum number of rays per sampling step.

    mseThreshold : float
        Convergence criterion for adaptive sampling.

    repetitions : int
        Maximum number of adaptive iterations.

    ------------------------------------------------------------------------
    Compute / Parallelization Parameters
    ------------------------------------------------------------------------

    deviceMode : str
        Execution mode, e.g. "CPU" or "GPU".

    parallelMode : str
        Parallel backend configuration.

    maxGPUs : int
        Maximum number of GPUs used.

    nPerNode : int
        Number of Devices per Node (is simply ignored when using parallel-mode: "threaded")

    ------------------------------------------------------------------------
    Returns
    ------------------------------------------------------------------------

    phi_ASE : list
        Computed ASE photon flux.

    mse_values : list
        Mean squared error values (convergence behavior).

    N_rays : list
        Total number of simulated rays.

    ------------------------------------------------------------------------
    Notes
    ------------------------------------------------------------------------

    - All inputs must have consistent sizes and indexing.
    - Incorrect shapes will typically not raise immediate errors but will
      result in incorrect physical results.
    - It is recommended to look into the provided example script 'laserPumpCladdingExample.py' to generate inputs.

    MPI Execution
    -------------

    If ``parallelMode`` is set to "mpi" or "graybat", the computation is executed
    via the standalone ``calcPhiASE`` executable using ``mpiexec``.

    In this mode:
    - Input data is written to a temporary folder
    - The external executable is launched
    - Results are read back into Python

    This is required because distributed MPI execution cannot be controlled
    directly from the Python binding (TODO : think about moving mpi layer (broadcast and gather) into python).

    For single-node execution, the Python binding is used directly without
    disk I/O.
    """

    if parallelMode=="mpi" or parallelMode=="graybat":
        return calcPhiASE_mpi(p,
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
                              nPerNode)
    p = np.asarray(p)
    t_int = np.asarray(t_int)
    beta_cell = np.asarray(beta_cell)
    beta_vol = np.asarray(beta_vol)
    clad_int = np.asarray(clad_int)
    refractiveIndices = np.asarray(refractiveIndices)
    reflectivities = np.asarray(reflectivities)
    normals_x = np.asarray(normals_x)
    normals_y = np.asarray(normals_y)
    ordered_int = np.asarray(ordered_int)
    surface = np.asarray(surface)
    x_center = np.asarray(x_center)
    y_center = np.asarray(y_center)
    normals_p = np.asarray(normals_p)
    forbidden = np.asarray(forbidden)

    numberOfPoints = int(p.shape[0])
    numberOfTriangles = int(t_int.shape[0])
    numberOfLevels = int(mesh_z)
    thicknessOfPrism = float(z_mesh)

    experiment = HASEonGPU_Bindings.ExperimentParameters(
        minRaysPerSample=int(minRaysPerSample),
        maxRaysPerSample=int(maxRaysPerSample),
        lambdaA=np.asarray(laser["l_abs"], dtype=float).reshape(-1).tolist(),
        lambdaE=np.asarray(laser["l_ems"], dtype=float).reshape(-1).tolist(),
        sigmaA=np.asarray(laser["s_abs"], dtype=float).reshape(-1).tolist(),
        sigmaE=np.asarray(laser["s_ems"], dtype=float).reshape(-1).tolist(),
        maxSigmaA=float(np.max(laser["s_abs"])),
        maxSigmaE=float(np.max(laser["s_ems"])),
        mseThreshold=float(mseThreshold),
        useReflections=bool(useReflections),
        spectral=int(laser["l_res"]),
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
        maxSampleRange=int((numberOfPoints*numberOfLevels)-1),
    )

    host_mesh = Mesh(
        triangleIndices=np.asarray(t_int, dtype=np.uint32).reshape(-1, order="F").tolist(),
        numberOfTriangles=numberOfTriangles,
        numberOfLevels=numberOfLevels,
        numberOfPoints=numberOfPoints,
        thicknessOfPrism=thicknessOfPrism,
        pointsVector=np.asarray(p, dtype=float).reshape(-1, order="F").tolist(),
        xOfTriangleCenter=np.asarray(x_center, dtype=float).reshape(-1).tolist(),
        yOfTriangleCenter=np.asarray(y_center, dtype=float).reshape(-1).tolist(),
        positionsOfNormalVectors=np.asarray(normals_p, dtype=np.uint32).reshape(-1, order="F").tolist(),
        xOfNormals=np.asarray(normals_x, dtype=float).reshape(-1, order="F").tolist(),
        yOfNormals=np.asarray(normals_y, dtype=float).reshape(-1, order="F").tolist(),
        forbiddenVector=np.asarray(forbidden, dtype=int).reshape(-1, order="F").tolist(),
        neighborsVector=np.asarray(ordered_int, dtype=int).reshape(-1, order="F").tolist(),
        surfacesVector=np.asarray(surface, dtype=np.float32).reshape(-1).tolist(),
        betaValuesVector=np.asarray(beta_vol, dtype=float).reshape(-1, order="F").tolist(),
        betaCells=np.asarray(beta_cell, dtype=float).reshape(-1, order="F").tolist(),
        cellTypes=np.asarray(clad_int, dtype=np.uint32).reshape(-1, order="F").tolist(),
        refractiveIndices=np.asarray(refractiveIndices, dtype=np.float32).reshape(-1, order="F").tolist(),
        reflectivities=np.asarray(reflectivities, dtype=np.float32).reshape(-1, order="F").tolist(),
        nTot=float(N_tot),
        crystalTFluo=float(crystal["tfluo"]),
        claddingNumber=int(clad_number),
        claddingAbsorption=float(clad_abs),
    )

    result = HASEonGPU_Bindings.calcPhiASE(experiment, compute, host_mesh)
    phi_ASE = np.asarray(result.phiAse, dtype=np.float32).reshape(
        (numberOfPoints, numberOfLevels), order="F"
    )
    mse_values = np.asarray(result.mse, dtype=np.float64).reshape(
        (numberOfPoints, numberOfLevels), order="F"
    )
    N_rays = np.asarray(result.totalRays, dtype=np.uint32).reshape(
        (numberOfPoints, numberOfLevels), order="F"
    )

    return phi_ASE, mse_values, N_rays
