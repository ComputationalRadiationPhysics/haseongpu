# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import tempfile

import numpy as np
import pytest

from HASEonGPU import AlpakaBackends, GainMedium, PhiASE, SpectralDecomposition, VolumeTopology


def analyticalPhiAseSphereCenter(gain, radius, beta, nTot, tauRad):
    gain = float(gain)
    radius = float(radius)
    beta = float(beta)

    if abs(gain) < 1.0e-12:
        return beta * radius

    return nTot * (beta / tauRad) * np.expm1(gain * radius) / gain


def calcBetaFromGain(gain, nTot, sigmaA, sigmaE):
    return (gain / nTot + sigmaA) / (sigmaA + sigmaE)


def constructExplicitSphereTopology(radius, *, samplePoints=None, meshSizeDivisor=8.0):
    gmshApi = pytest.importorskip("gmsh")
    center = np.asarray((radius, radius, radius), dtype=np.float64)
    with tempfile.TemporaryDirectory() as tmpdir:
        msh = f"{tmpdir}/sphere_tet4.msh"
        gmshApi.initialize()
        try:
            gmshApi.option.setNumber("General.Terminal", 0)
            gmshApi.clear()
            gmshApi.model.add("sphere_tet4")
            sphere = gmshApi.model.occ.addSphere(float(center[0]), float(center[1]), float(center[2]), float(radius))
            gmshApi.model.occ.synchronize()
            gmshApi.model.addPhysicalGroup(3, [sphere], 1)
            gmshApi.model.setPhysicalName(3, 1, "gain")
            surfaces = [tag for dim, tag in gmshApi.model.getEntities(2)]
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
if not alpakaBackends:
    raise RuntimeError("analytical sphere tests require at least one Alpaka backend")


def _requireAvailableAlpakaBackend(backend):
    if backend not in alpakaBackends:
        available = ", ".join(alpakaBackends)
        raise AssertionError(
            f"requested Alpaka backend {backend!r} is not reported by "
            f"libHaseAlpakaBackendNames; available backends: {available}"
        )


@pytest.mark.parametrize("backend", alpakaBackends)
@pytest.mark.parametrize(("radius", "g0"), sphereCases, ids=sphereCaseIds)
def test_centerPointIntegralMatchesAnalyticalSolution(radius, g0, backend, phiAseTestConfigPath, openPmdRuntimeBackend):
    _requireAvailableAlpakaBackend(backend)
    xDim = radius * 2.0
    nTot = np.float64(1.38e20 * 1.0)
    sigmaA = np.float64(0.11e-20)
    sigmaE = np.float64(2.1e-20)
    gain = g0
    beta = calcBetaFromGain(gain, nTot, sigmaA=sigmaA, sigmaE=sigmaE)
    print(f' running with: g0: {g0} and radius: {radius} and beta: {beta}')
    flourescenceLifetime = np.float64(9.41e-4)

    crossSections = SpectralDecomposition.monochromatic(
        wavelength=np.float64(1030e-9),
        crossSectionAbsorption=sigmaA,
        crossSectionEmission=sigmaE,
    )
    center = (radius, radius, radius)

    topology = constructExplicitSphereTopology(radius, samplePoints=[center])
    medium = GainMedium(topology=topology)
    betaCells = np.full(topology.numberOfSamplePoints, beta, dtype=np.float64)
    betaVolume = np.full(topology.numberOfCells, beta, dtype=np.float64)
    flatBetaVolume = betaVolume.reshape(-1)
    assert np.any(flatBetaVolume > 0.0)
    assert np.any(betaCells > 0.0)
    cells = medium.get("betaCells").expectedShape
    volume = medium.get("betaVolume").expectedShape
    print(f'betaCells: {cells}, betaVolume: {volume}')
    medium.withPhysicalProperties(
        betaCells=betaCells,
        betaVolume=betaVolume,
        nTot=nTot,
        crystalTFluo=flourescenceLifetime,
    )

    centerSample = 0
    phiAse = PhiASE.fromYaml(
        phiAseTestConfigPath,
        spectralProperties=crossSections,
        propagationMode="backward",
        minRaysPerSample=10000,
        maxRaysPerSample=100000,
        repetitions=2,
        adaptiveSteps=3,
        mseThreshold=0.05,
        useReflections=False,
        backend=backend,
        parallelMode="single",
        numDevices=1,
        minSampleRange=centerSample,
        maxSampleRange=centerSample,
        monochromatic=True,
        openpmdBackend=openPmdRuntimeBackend,
    )
    phiAse.run(gainMedium=medium)
    result = phiAse.getResults()
    numerical = np.array(result.phiAse, dtype=np.float64).reshape(-1)[centerSample]

    expected = analyticalPhiAseSphereCenter(
        gain=gain,
        radius=radius,
        beta=beta,
        nTot=nTot,
        tauRad=flourescenceLifetime,
    )
    print(f'expected: {expected}, numerical {numerical}')
    assert np.isfinite(numerical)
    assert numerical > 0.0
    assert np.isclose(numerical, expected, rtol=0.05)


@pytest.mark.parametrize("backend", ["Host_Cpu_CpuOmpBlocks"])
def testForwardSphereCenterVolumeMatchesAnalyticalSolutionOnOmp(backend):
    radius = np.float64(0.7196856730011522)
    gain = np.float64(2.14)
    nTot = np.float64(1.38e20)
    sigmaA = np.float64(0.11e-20)
    sigmaE = np.float64(2.1e-20)
    beta = calcBetaFromGain(gain, nTot, sigmaA=sigmaA, sigmaE=sigmaE)
    flourescenceLifetime = np.float64(9.41e-4)
    center = np.asarray((radius, radius, radius), dtype=np.float64)

    topology = constructExplicitSphereTopology(radius, meshSizeDivisor=17.0)
    assert topology.numberOfCells >= 80_000

    centerVolume = nearestVolumeIndex(topology, center)
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
    phiAse = PhiASE(
        spectralProperties=crossSections,
        maxRaysPerSample=1_000_000,
        forwardRayCount=1_000_000,
        forwardRayLength=2.0 * float(radius),
        repetitions=1,
        adaptiveSteps=1,
        mseThreshold=0.05,
        useReflections=False,
        backend=backend,
        parallelMode="single",
        numDevices=1,
        monochromatic=True,
        rngSeed=1234,
    )

    try:
        phiAse.run(gainMedium=medium)
    except RuntimeError as exc:
        if "return code 1" in str(exc):
            pytest.skip(f"backend {backend} is not available in this build")
        raise

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
    assert np.isclose(numerical, expected, rtol=0.2)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
