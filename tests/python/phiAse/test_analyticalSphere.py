# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np
import pytest

from HASEonGPU import AlpakaBackends, GainMedium, Grid, MeshTopology, PhiASE, SpectralDecomposition, vtkWedge


def analyticalPhiAseSphereCenter(gain, radius, beta, nTot, tauRad):
    gain = float(gain)
    radius = float(radius)
    beta = float(beta)

    if abs(gain) < 1.0e-12:
        return beta * radius

    return nTot * (beta / tauRad) * np.expm1(gain * radius) / gain



def _pointInTriangle(point, trianglePoints, tol=1e-12):
    p = np.asarray(point, dtype=np.float64)
    a, b, c = np.asarray(trianglePoints, dtype=np.float64)
    v0 = c - a
    v1 = b - a
    v2 = p - a
    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)
    denom = dot00 * dot11 - dot01 * dot01
    if abs(denom) <= tol:
        return False
    inv = 1.0 / denom
    u = (dot11 * dot02 - dot01 * dot12) * inv
    v = (dot00 * dot12 - dot01 * dot02) * inv
    return u >= -tol and v >= -tol and (u + v) <= 1.0 + tol


def _triangleIndexAt(topology, x, y):
    point = np.asarray([x, y], dtype=np.float64)
    for triangleIndex, triangle in enumerate(topology.trianglePointIndices):
        if _pointInTriangle(point, topology.points[np.asarray(triangle, dtype=np.uint32)]):
            return int(triangleIndex)
    raise ValueError(f"no triangle contains ({x}, {y})")


def _betaVolumeFlatIndex(topology, x, y, z):
    topology._require_levels()
    topology._require_thickness()
    layerIndex = int(np.floor(float(z) / float(topology.thickness)))
    layerIndex = min(max(layerIndex, 0), topology.levels - 2)
    triangleIndex = _triangleIndexAt(topology, x, y)
    return triangleIndex + layerIndex * topology.numberOfTriangles


def _triangleCenters(topology):
    centers = []
    for triangle in topology.trianglePointIndices:
        triangle = np.asarray(triangle, dtype=np.uint32)
        centers.append(np.mean(topology.points[triangle], axis=0))
    return np.asarray(centers, dtype=np.float64)


def _prismCenters(topology):
    topology._require_levels()
    topology._require_thickness()

    triangleCenters = _triangleCenters(topology)
    centers = np.empty((topology.numberOfPrisms, 3), dtype=np.float64)
    levelCoordinates = topology.levelCoordinates()

    for layerIndex in range(topology.levels - 1):
        zCenter = 0.5 * (levelCoordinates[layerIndex] + levelCoordinates[layerIndex + 1])
        for triangleIndex, xyCenter in enumerate(triangleCenters):
            flatIndex = triangleIndex + layerIndex * topology.numberOfTriangles
            centers[flatIndex] = [xyCenter[0], xyCenter[1], zCenter]

    return centers


def _flatPrismMaskSphere(topology, *, center=(5.0, 5.0, 5.0), radius=5.0):
    centers = _prismCenters(topology)
    center = np.asarray(center, dtype=np.float64)
    return np.linalg.norm(centers - center, axis=1) < radius


def constructBetaVolumeSphere(topology, *, center=(5.0, 5.0, 5.0), radius=5.0, beta):
    """
    Set prism gain volumes inside a sphere.

    betaVolume uses the documented flat layout internally:

        flatIndex = triangleIndex + layerIndex * numberOfTriangles

    The returned matrix layout is accepted by GainMedium.withPhysicalProperties.
    """
    flat = np.zeros(topology.numberOfPrisms, dtype=np.float64)
    flat[_flatPrismMaskSphere(topology, center=center, radius=radius)] = beta
    return flat.reshape((topology.numberOfTriangles, topology.levels - 1), order="F")


def calcBetaFromGain(gain, nTot, sigmaA, sigmaE):
    return (gain / nTot + sigmaA) / (sigmaA + sigmaE)


def constructBetaCellsSphere(topology, *, center=(5.0, 5.0, 5.0), radius=5.0, beta):
    center = np.asarray(center, dtype=np.float64)
    betaCells = np.zeros((topology.numberOfPoints, topology.levels), dtype=np.float64)
    for pointIndex, (x, y) in enumerate(topology.points):
        for levelIndex, z in enumerate(topology.levelCoordinates()):
            r = np.linalg.norm(np.asarray([x, y, z]) - center)
            if r < radius:
                betaCells[pointIndex, levelIndex] = float(beta)
    return betaCells
nTot = np.float64(1.38e20 * 1.0)
sigmaA = np.float64(0.11e-20)
sigmaE = np.float64(2.1e-20)
sphereCases = [
    (np.float64(radiusValue), np.float64(g0Value / 100))
    for radiusValue in np.geomspace(0.1, 100, num=8)
    for g0Value in np.geomspace(5, 400, num=8)
    if 5.0 >= np.float64(radiusValue) * np.float64(g0Value / 100) >= 1.0 >= calcBetaFromGain(g0Value/ 100, nTot, sigmaA, sigmaE) >= 0.0
]


sphereCaseIds = [f"R{float(radius):g}_g0_{float(g0):.2f}" for radius, g0 in sphereCases]
alpakaBackends = AlpakaBackends.all()
@pytest.mark.parametrize("backend", alpakaBackends)
@pytest.mark.parametrize(("radius", "g0"), sphereCases, ids=sphereCaseIds)
def testCenterPointIntegralMatchesAnalyticalSolution(radius, g0, backend, phiAseTestConfigPath):
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

    grid = Grid(xExtent=xDim, yExtent=xDim, zExtent=xDim, tileSizeX=xDim / 100)
    topology = MeshTopology.fromGrid(grid)
    medium = GainMedium(topology=topology)
    betaCells = constructBetaCellsSphere(topology, center=center, radius=radius, beta=beta)
    betaVolume = constructBetaVolumeSphere(topology, center=center, radius=radius, beta=beta)
    flatBetaVolume = betaVolume.reshape(-1, order="F")
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

    centerSample = medium.betaCellIndexAt(*center, flat=True)
    phiAse = PhiASE.fromYaml(
        phiAseTestConfigPath,
        spectralProperties=crossSections,
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
    )
    try:
        phiAse.run(gainMedium=medium)
    except RuntimeError as exc:
        if "return code 1" in str(exc):
            pytest.skip(f"backend {backend} is not available in this build")
        raise
    result = phiAse.getResults()
    shape = medium.get("betaCells").expectedShape
    phiAseValues = np.array(result.phiAse, dtype=np.float64).reshape(shape, order="F")
    numerical = phiAseValues.reshape(-1, order="F")[centerSample]

    expected = analyticalPhiAseSphereCenter(
        gain=gain,
        radius=radius,
        beta=beta,
        nTot=nTot,
        tauRad=flourescenceLifetime
    )
    print(f'expected: {expected}, numerical {numerical}')
    # assert np.isfinite(numerical)
    # assert numerical > 0.0
    #
    # vtkWedge("phi0.vtk", data=phiAse, geometry=medium.topology,field="phiAse")
    assert np.isclose(numerical, expected, rtol=0.05)
def make_grid_within_radius(radius: float):
    return [
        (x, y)
        for x in [i * radius / 10 for i in range(-10, 11)]
        for y in [j * radius / 10 for j in range(-10, 11)]
        if x * x + y * y <= radius * radius
    ]
diskCases = [
    (np.float64(radiusValue), np.float64(g0Value / 100))
    for radiusValue in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 50.0, 70.0, 100.0]
    for g0Value in range(10, 405, 10)
    if 1.0 >= calcBetaFromGain(g0Value / 100, nTot, sigmaA, sigmaE) >= 0.0 and (g0Value / 100) * radiusValue == 2.0
]

diskCaseIds = [f"R{float(radius):g}_g0_{float(g0):.2f}" for radius, g0 in diskCases]
@pytest.mark.parametrize(("radius", "g0"), diskCases, ids=diskCaseIds)
def testDiskPointsMatchAnalyticalSolutionForG02(radius, g0, phiAseTestConfigPath):
    global expected
    backend=alpakaBackends[-1]
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

    grid = Grid(xExtent=xDim, yExtent=xDim, zExtent=xDim, tileSizeX=xDim / 100)
    topology = MeshTopology.fromGrid(grid)
    medium = GainMedium(topology=topology)
    betaCells = constructBetaCellsSphere(topology, center=center, radius=radius, beta=beta)
    betaVolume = constructBetaVolumeSphere(topology, center=center, radius=radius, beta=beta)
    flatBetaVolume = betaVolume.reshape(-1, order="F")
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
    expectedValues=[]
    for offsetX,offsetY in make_grid_within_radius(radius):
        center = (radius+offsetX, radius+offsetY, radius)
        centerSample = medium.betaCellIndexAt(*center, flat=True)
        phiAse = PhiASE.fromYaml(
            phiAseTestConfigPath,
            spectralProperties=crossSections,
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
        )
        try:
            phiAse.run(gainMedium=medium)
        except RuntimeError as exc:
            if "return code 1" in str(exc):
                pytest.skip(f"backend {backend} is not available in this build")
            raise
        result = phiAse.getResults()
        shape = medium.get("betaCells").expectedShape
        phiAseValues = np.array(result.phiAse, dtype=np.float64).reshape(shape, order="F")
        numerical = phiAseValues.reshape(-1, order="F")[centerSample]

        expected = analyticalPhiAseSphereCenter(
            gain=gain,
            radius=radius,
            beta=beta,
            nTot=nTot,
            tauRad=flourescenceLifetime
        )
        print(f'expected: {expected}, numerical {numerical}')
        # assert np.isfinite(numerical)
        # assert numerical > 0.0
        #
        # vtkWedge("phi0.vtk", data=phiAse, geometry=medium.topology,field="phiAse")
        assert np.isclose(numerical, expected, rtol=0.05)
        expectedValues.append(expected)
    assert all(np.allclose(a, expected[0], rtol=0.2) for a in expected[1:])

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
