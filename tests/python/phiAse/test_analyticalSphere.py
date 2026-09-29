# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import os
import tempfile
from pathlib import Path

import numpy as np
import pytest

from HASEonGPU import Domain, GainMedium, OpticalComponent, PhiASE, VolumeTopology
from hase_units import units
from material_library import CrossSectionTable, Material
from openpmd_backend_matrix import openpmd_runtime_test_backends
from alpaka_backend_matrix import alpaka_runtime_backend
from pyInclude.geometry.vtk import _parseVtk


repoRoot = Path(__file__).resolve().parents[3]
TERRA_DIAMOND_SPHERE_PATH = (
    repoRoot / "tests" / "data" / "analyticalSphere" / "terraDiamondSphere.vtk"
)
TERRA_DIAMOND_SPHERE_RADIUS = np.float64(0.1)
TERRA_DIAMOND_SPHERE_GAIN = np.float64(10.0)


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
            meshSize = max(float(radius) / float(meshSizeDivisor), 1.0e-5)
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


nTot = np.float64(1.38e26)
sigmaA = np.float64(0.11e-24)
sigmaE = np.float64(2.1e-24)
sphereCases = [
    (np.float64(radiusValue), np.float64(gainValue))
    for radiusValue in np.geomspace(0.001, 1.0, num=8)
    for gainValue in np.geomspace(5, 400, num=8)
    if 5.0 >= np.float64(radiusValue) * np.float64(gainValue) >= 1.0 >= calcBetaFromGain(gainValue, nTot, sigmaA, sigmaE) >= 0.0
]


sphereCaseIds = [f"R{float(radius):g}_g0_{float(g0):.2f}" for radius, g0 in sphereCases]
_NO_ANALYTICAL_SPHERE_BACKEND = "__no_analytical_sphere_backend__"


def analyticalSphereBackends():
    try:
        return [alpaka_runtime_backend()]
    except RuntimeError:
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

    nTot = np.float64(1.38e26)
    sigmaA = np.float64(0.11e-24)
    sigmaE = np.float64(2.1e-24)
    beta = calcBetaFromGain(gain, nTot, sigmaA=sigmaA, sigmaE=sigmaE)
    flourescenceLifetime = np.float64(9.41e-4)
    topology = constructExplicitSphereTopology(radius, meshSizeDivisor=analyticalSphereMeshSizeDivisor())
    assert topology.numberOfCells >= 1_000

    centerVolume = centeredVolumeIndex(topology, radius)
    material = Material(
        materialName="analytical sphere material",
        temperature=293.15 * units.K,
        refractiveIndex=1.0,
        fluorescenceLifetime=flourescenceLifetime * units.s,
        crossSections=CrossSectionTable.monochromatic(
            wavelength=np.float64(1030e-9) * units.m,
            absorption=sigmaA * units.m**2,
            emission=sigmaE * units.m**2,
        ),
        active=True,
        activeIonDensity=nTot / units.m**3,
    )
    component = OpticalComponent(domain=Domain.fromTopology(topology), material=material)
    medium = GainMedium([component])
    rayCount = analyticalSphereRayCount()
    phiAse = PhiASE(
        maxRays=rayCount,
        forwardRayCount=rayCount,
        repetitions=1,
        adaptiveSteps=1,
        relativeStandardErrorThreshold=0.05,
        enableDiagnostics=True,
        useReflections=False,
        backend=backend,
        openpmdBackend=openpmdBackend,
        parallelMode="single",
        numDevices=1,
        monochromatic=True,
        rngSeed=1234,
    )

    phiAse.run(gainMedium=medium, initialExcitation=beta)

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


@pytest.mark.parametrize("backend", analyticalSphereBackends())
@pytest.mark.parametrize("openpmdBackend", openpmd_runtime_test_backends())
def testForwardTerraDiamondSphereCenterMatchesMonolithicAndAnalyticalSolutions(
    openpmdBackend,
    backend,
):
    if backend == _NO_ANALYTICAL_SPHERE_BACKEND:
        pytest.fail("analytical sphere test requires at least one Alpaka backend")

    topology = VolumeTopology.fromFile(TERRA_DIAMOND_SPHERE_PATH)
    _points, _cells, _types, _point_data, cell_data, _fields = _parseVtk(
        TERRA_DIAMOND_SPHERE_PATH
    )
    diamond_ids = np.asarray(cell_data["diamondId"], dtype=np.uint32)
    cell_domains = np.asarray(topology.cellDomains, dtype=np.uint32)
    assert topology.numberOfCells == 40_745
    np.testing.assert_array_equal(diamond_ids, cell_domains)
    np.testing.assert_array_equal(np.unique(diamond_ids), np.arange(10, dtype=np.uint32))
    assert np.all(np.bincount(diamond_ids, minlength=10) > 0)

    beta = calcBetaFromGain(
        TERRA_DIAMOND_SPHERE_GAIN,
        nTot,
        sigmaA=sigmaA,
        sigmaE=sigmaE,
    )
    fluorescence_lifetime = np.float64(9.41e-4)
    material = Material(
        materialName="analytical TERRA diamond sphere material",
        temperature=293.15 * units.K,
        refractiveIndex=1.0,
        fluorescenceLifetime=fluorescence_lifetime * units.s,
        crossSections=CrossSectionTable.monochromatic(
            wavelength=np.float64(1030e-9) * units.m,
            absorption=sigmaA * units.m**2,
            emission=sigmaE * units.m**2,
        ),
        active=True,
        activeIonDensity=nTot / units.m**3,
    )
    diamond_domains = [
        Domain.fromGmsh(topology, diamond_id, entityKind="volume")
        for diamond_id in range(10)
    ]
    coverage = np.zeros(topology.numberOfCells, dtype=np.uint32)
    for domain in diamond_domains:
        coverage += domain.maskFor(topology)
    np.testing.assert_array_equal(coverage, np.ones_like(coverage))

    ray_count = analyticalSphereRayCount()

    def run(components):
        phi_ase = PhiASE(
            maxRays=ray_count,
            forwardRayCount=ray_count,
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
        phi_ase.run(
            gainMedium=GainMedium(components),
            initialExcitation=beta,
        )
        return phi_ase.getResults()

    monolithic = run(
        [OpticalComponent(domain=Domain.fromTopology(topology), material=material)]
    )
    decomposed = run(
        [
            OpticalComponent(domain=domain, material=material)
            for domain in diamond_domains
        ]
    )
    center_volume = centeredVolumeIndex(topology, TERRA_DIAMOND_SPHERE_RADIUS)
    monolithic_value = np.asarray(monolithic.phiAse, dtype=np.float64)[center_volume]
    decomposed_value = np.asarray(decomposed.phiAse, dtype=np.float64)[center_volume]
    decomposed_rse = np.asarray(decomposed.relativeStandardError, dtype=np.float64)[
        center_volume
    ]
    expected = analyticalPhiAseSphereCenter(
        gain=TERRA_DIAMOND_SPHERE_GAIN,
        radius=TERRA_DIAMOND_SPHERE_RADIUS,
        beta=beta,
        nTot=nTot,
        tauRad=fluorescence_lifetime,
    )
    print(
        f"TERRA diamond sphere center: tets={topology.numberOfCells}, "
        f"expected={expected}, monolithic={monolithic_value}, "
        f"decomposed={decomposed_value}, decomposedRse={decomposed_rse}"
    )
    assert np.isfinite(decomposed_rse)
    assert np.isclose(monolithic_value, expected, rtol=0.05)
    assert np.isclose(decomposed_value, expected, rtol=0.05)
    assert np.isclose(decomposed_value, monolithic_value, rtol=0.05)


@pytest.mark.parametrize("backend", analyticalSphereBackends())
@pytest.mark.parametrize("openpmdBackend", openpmd_runtime_test_backends())
def testForwardSphereCenterRefinementStudy(openpmdBackend, backend):
    if backend == _NO_ANALYTICAL_SPHERE_BACKEND:
        pytest.fail("analytical sphere test requires at least one Alpaka backend")

    radius = np.float64(0.01)
    gain = np.float64(100.0)
    beta = calcBetaFromGain(gain, nTot, sigmaA, sigmaE)
    fluorescenceLifetime = np.float64(9.41e-4)
    expected = analyticalPhiAseSphereCenter(gain, radius, beta, nTot, fluorescenceLifetime)
    material = Material(
        materialName="analytical sphere material",
        temperature=293.15 * units.K,
        refractiveIndex=1.0,
        fluorescenceLifetime=fluorescenceLifetime * units.s,
        crossSections=CrossSectionTable.monochromatic(
            wavelength=np.float64(1030e-9) * units.m,
            absorption=sigmaA * units.m**2,
            emission=sigmaE * units.m**2,
        ),
        active=True,
        activeIonDensity=nTot / units.m**3,
    )

    for meshSizeDivisor in (8.0, 12.0):
        topology = constructExplicitSphereTopology(radius, meshSizeDivisor=meshSizeDivisor)
        centerVolume = centeredVolumeIndex(topology, radius)
        medium = GainMedium([OpticalComponent(domain=Domain.fromTopology(topology), material=material)])
        for rayCount in (2_000_000, 4_000_000):
            phiAse = PhiASE(
                maxRays=rayCount,
                forwardRayCount=rayCount,
                repetitions=1,
                adaptiveSteps=1,
                relativeStandardErrorThreshold=0.05,
                enableDiagnostics=True,
                useReflections=False,
                backend=backend,
                openpmdBackend=openpmdBackend,
                parallelMode="single",
                numDevices=1,
                monochromatic=True,
                rngSeed=1234,
            )
            phiAse.run(gainMedium=medium, initialExcitation=beta)
            result = phiAse.getResults()
            numerical = float(np.asarray(result.phiAse).reshape(-1)[centerVolume])
            rse = float(np.asarray(result.relativeStandardError).reshape(-1)[centerVolume])
            visits = int(np.asarray(result.totalRays).reshape(-1)[centerVolume])
            relativeError = abs(numerical / expected - 1.0)
            print(
                f"sphere refinement: backend={backend}, meshDivisor={meshSizeDivisor}, "
                f"tets={topology.numberOfCells}, rays={rayCount}, visits={visits}, "
                f"relativeError={relativeError:.6g}, rse={rse:.6g}"
            )
            assert visits > 0
            assert np.isfinite(rse)
            assert np.isclose(numerical, expected, rtol=0.05)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
