from HASEonGPU_Bindings import *
from pyInclude import *
del HASEonGPU_Bindings.HostMesh
import numpy as np
def calcPhiASE_legacy(
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
        nPerNode,
):
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

    experiment = ExperimentParameters(
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

    compute = ComputeParameters(
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
    print(f'spectral: {int(laser["l_res"])}')
    print(f'number of points: {int(numberOfPoints - 1)}')

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

    result = calcPhiASE(experiment, compute, host_mesh)

    phi_ASE = np.asarray(result.phiAse, dtype=np.float32).reshape(-1, order="F").tolist()
    mse_values = np.asarray(result.mse, dtype=np.float64).reshape(-1, order="F").tolist()
    N_rays = np.asarray(result.totalRays, dtype=np.uint32).reshape(-1, order="F").tolist()

    return phi_ASE, mse_values, N_rays
