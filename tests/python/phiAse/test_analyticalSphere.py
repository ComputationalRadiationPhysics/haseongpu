# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import os
import tempfile

import numpy as np
import pytest

from HASEonGPU import AlpakaBackends, GainMedium, PhiASE, SpectralDecomposition, VolumeTopology
from openpmd_backend_matrix import openpmd_runtime_test_backends


def analyticalPhiAseSphereCenter(gain, radius, beta, nTot, tauRad):
    gain = float(gain)
    radius = float(radius)
    beta = float(beta)

    if abs(gain) < 1.0e-12:
        return nTot * (beta / tauRad) * radius

    return nTot * (beta / tauRad) * np.expm1(gain * radius) / gain


def calcBetaFromGain(gain, nTot, sigmaA, sigmaE):
    return (gain / nTot + sigmaA) / (sigmaA + sigmaE)


def testAnalyticalPhiAseSphereCenterHasCorrectZeroGainLimit():
    radius = 2.0
    beta = 0.25
    nTot = 3.0
    tauRad = 0.5
    expected = nTot * beta * radius / tauRad

    assert analyticalPhiAseSphereCenter(0.0, radius, beta, nTot, tauRad) == expected
    assert np.isclose(analyticalPhiAseSphereCenter(1.0e-11, radius, beta, nTot, tauRad), expected)


def _requireGmsh():
    try:
        import gmsh as gmshApi
    except (ImportError, OSError) as exc:
        if "libGLU.so.1" in str(exc):
            pytest.xfail(f"gmsh runtime dependency is unavailable: {exc}")
        pytest.fail(f"gmsh is required for analytical sphere topology generation: {exc}")
    return gmshApi


def constructExplicitSphereTopology(radius, *, samplePoints=None, meshSizeDivisor=8.0):
    gmshApi = _requireGmsh()
    center = np.zeros(3, dtype=np.float64)
    with tempfile.TemporaryDirectory() as tmpdir:
        msh = f"{tmpdir}/sphere_tet4.msh"
        gmshApi.initialize()
        try:
            gmshApi.option.setNumber("General.Terminal", 0)
            gmshApi.clear()
            gmshApi.model.add("sphere_tet4")
            sphere = gmshApi.model.occ.addSphere(float(center[0]), float(center[1]), float(center[2]), float(radius))

            # Preserve one regular central tetrahedron as a distinct CAD volume. Its
            # signed-coordinate vertices have an exact centroid at the sphere center.
            centralScale = float(radius) / (3.0 * float(meshSizeDivisor))
            centralCoordinates = (
                (+centralScale, +centralScale, +centralScale),
                (+centralScale, -centralScale, -centralScale),
                (-centralScale, +centralScale, -centralScale),
                (-centralScale, -centralScale, +centralScale),
            )
            points = [gmshApi.model.occ.addPoint(*coordinate) for coordinate in centralCoordinates]
            line01 = gmshApi.model.occ.addLine(points[0], points[1])
            line02 = gmshApi.model.occ.addLine(points[0], points[2])
            line03 = gmshApi.model.occ.addLine(points[0], points[3])
            line12 = gmshApi.model.occ.addLine(points[1], points[2])
            line13 = gmshApi.model.occ.addLine(points[1], points[3])
            line23 = gmshApi.model.occ.addLine(points[2], points[3])
            faces = [
                gmshApi.model.occ.addPlaneSurface(
                    [gmshApi.model.occ.addCurveLoop([line02, -line12, -line01])]
                ),
                gmshApi.model.occ.addPlaneSurface(
                    [gmshApi.model.occ.addCurveLoop([line01, line13, -line03])]
                ),
                gmshApi.model.occ.addPlaneSurface(
                    [gmshApi.model.occ.addCurveLoop([line03, -line23, -line02])]
                ),
                gmshApi.model.occ.addPlaneSurface(
                    [gmshApi.model.occ.addCurveLoop([line12, line23, -line13])]
                ),
            ]
            centralSurfaceLoop = gmshApi.model.occ.addSurfaceLoop(faces)
            centralTetrahedron = gmshApi.model.occ.addVolume([centralSurfaceLoop])
            fragments, _ = gmshApi.model.occ.fragment(
                [(3, sphere)],
                [(3, centralTetrahedron)],
                removeObject=True,
                removeTool=True,
            )
            gmshApi.model.occ.synchronize()
            volumeTags = [tag for dim, tag in fragments if dim == 3]
            gmshApi.model.addPhysicalGroup(3, volumeTags, 1)
            gmshApi.model.setPhysicalName(3, 1, "gain")
            combinedBoundary = gmshApi.model.getBoundary(
                [(3, tag) for tag in volumeTags],
                combined=True,
                oriented=False,
                recursive=False,
            )
            surfaces = [tag for dim, tag in combinedBoundary if dim == 2]
            if surfaces:
                gmshApi.model.addPhysicalGroup(2, surfaces, 2)
                gmshApi.model.setPhysicalName(2, 2, "outer")
            meshSize = max(float(radius) / float(meshSizeDivisor), 1.0e-3)
            gmshApi.option.setNumber("Mesh.CharacteristicLengthMin", meshSize)
            gmshApi.option.setNumber("Mesh.CharacteristicLengthMax", meshSize)
            gmshApi.model.mesh.generate(3)
            gmshApi.write(msh)
        finally:
            gmshApi.finalize()
        topology = VolumeTopology.fromFile(msh)
    if samplePoints is not None:
        topology.samplePoints = np.asarray(samplePoints, dtype=np.float64).reshape((-1, 3))
    return topology


def nearestVolumeIndex(topology, point):
    point = np.asarray(point, dtype=np.float64)
    distances = np.linalg.norm(np.asarray(topology.cellCenters, dtype=np.float64) - point, axis=1)
    return int(np.argmin(distances))


def centeredVolumeIndex(topology, radius):
    centerVolume = nearestVolumeIndex(topology, np.zeros(3, dtype=np.float64))
    center = np.asarray(topology.cellCenters[centerVolume], dtype=np.float64)
    np.testing.assert_allclose(
        center,
        np.zeros(3, dtype=np.float64),
        rtol=0.0,
        atol=64.0 * np.finfo(np.float64).eps * max(1.0, float(radius)),
    )
    return centerVolume


nTot = np.float64(1.38e20 * 1.0)
sigmaA = np.float64(0.11e-20)
sigmaE = np.float64(2.1e-20)
sphereCases = [
    (np.float64(radiusValue), np.float64(g0Value / 100))
    for radiusValue in np.geomspace(0.1, 100, num=8)
    for g0Value in np.geomspace(5, 400, num=8)
    if 5.0 >= np.float64(radiusValue) * np.float64(g0Value / 100) >= 1.0 >= calcBetaFromGain(g0Value / 100, nTot, sigmaA, sigmaE) >= 0.0
]


sphereCaseIds = [f"R{float(radius):g}_g0_{float(g0):.2f}" for radius, g0 in sphereCases]
alpakaBackends = AlpakaBackends.all()
_NO_ANALYTICAL_SPHERE_BACKEND = "__no_analytical_sphere_backend__"


def analyticalSphereBackends():
    for preferred in ("Host_Cpu_CpuOmpBlocks", "Host_Cpu_CpuSerial"):
        if preferred in alpakaBackends:
            return [preferred]
    cpuBackends = [backend for backend in alpakaBackends if "Cpu" in backend]
    if cpuBackends:
        return [cpuBackends[0]]
    if alpakaBackends:
        return [alpakaBackends[0]]
    return [_NO_ANALYTICAL_SPHERE_BACKEND]


def analyticalSphereRayCount():
    return int(os.environ.get("HASE_ANALYTICAL_SPHERE_RAYS", "2000000"))


def analyticalSphereMeshSizeDivisor():
    divisor = float(os.environ.get("HASE_ANALYTICAL_SPHERE_MESH_SIZE_DIVISOR", "10.0"))
    if not np.isfinite(divisor) or divisor <= 0.0:
        raise ValueError("HASE_ANALYTICAL_SPHERE_MESH_SIZE_DIVISOR must be finite and positive")
    return divisor


@pytest.mark.parametrize("backend", analyticalSphereBackends())
@pytest.mark.parametrize("openpmdBackend", openpmd_runtime_test_backends())
@pytest.mark.parametrize(("radius", "gain"), sphereCases, ids=sphereCaseIds)
def testForwardSphereCenterVolumeMatchesAnalyticalSolution(radius, gain, openpmdBackend, backend):
    if backend == _NO_ANALYTICAL_SPHERE_BACKEND:
        pytest.fail("analytical sphere test requires at least one Alpaka backend")

    nTot = np.float64(1.38e20)
    sigmaA = np.float64(0.11e-20)
    sigmaE = np.float64(2.1e-20)
    beta = calcBetaFromGain(gain, nTot, sigmaA=sigmaA, sigmaE=sigmaE)
    flourescenceLifetime = np.float64(9.41e-4)
    topology = constructExplicitSphereTopology(radius, meshSizeDivisor=analyticalSphereMeshSizeDivisor())
    assert topology.numberOfCells >= 1_000

    centerVolume = centeredVolumeIndex(topology, radius)
    medium = GainMedium(topology=topology)
    medium.withPhysicalProperties(
        betaVolume=np.full(topology.numberOfCells, beta, dtype=np.float64),
        claddingCellTypes=np.zeros(topology.numberOfCells, dtype=np.uint32),
        nTot=nTot,
        crystalTFluo=flourescenceLifetime,
    )
    crossSections = SpectralDecomposition.monochromatic(
        wavelength=np.float64(1030e-9),
        crossSectionAbsorption=sigmaA,
        crossSectionEmission=sigmaE,
    )
    rayCount = analyticalSphereRayCount()
    phiAse = PhiASE(
        spectralProperties=crossSections,
        maxRays=rayCount,
        forwardRayCount=rayCount,
        repetitions=1,
        adaptiveSteps=1,
        relativeStandardErrorThreshold=0.05,
        useReflections=False,
        backend=backend,
        openpmdBackend=openpmdBackend,
        parallelMode="single",
        numDevices=1,
        monochromatic=True,
        rngSeed=1234,
    )

    phiAse.run(gainMedium=medium)

    result = phiAse.getResults()
    phiAseValues = np.asarray(result.phiAse, dtype=np.float64).reshape(-1)
    totalRays = np.asarray(result.totalRays, dtype=np.uint32).reshape(-1)
    assert phiAseValues.size == topology.numberOfCells
    assert totalRays[centerVolume] > 0

    numerical = phiAseValues[centerVolume]
    expected = analyticalPhiAseSphereCenter(
        gain=gain,
        radius=radius,
        beta=beta,
        nTot=nTot,
        tauRad=flourescenceLifetime,
    )
    print(
        f"forward center volume: tets={topology.numberOfCells}, "
        f"centerVolume={centerVolume}, visits={int(totalRays[centerVolume])}, "
        f"expected={expected}, numerical={numerical}"
    )
    assert np.isfinite(numerical)
    assert numerical > 0.0
    assert np.isclose(numerical, expected, rtol=0.05)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
