# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


from pathlib import Path

import numpy as np
import pytest

from HASEonGPU import GainMedium, Grid, MeshTopology, backendFlat, unitDimension


def _lineCount(path):
    return sum(1 for _ in Path(path).open("r", encoding="utf-8"))


def _assertTopologyLaunchContract(topology):
    topology._require_levels()
    topology._require_thickness()

    numberOfPoints = topology.numberOfPoints
    numberOfTriangles = topology.numberOfTriangles
    numberOfLevels = int(topology.levels)

    assert topology.points.shape == (numberOfPoints, 2)
    assert topology.points.reshape(-1, order="F").size == 2 * numberOfPoints

    assert topology.trianglePointIndices.shape == (numberOfTriangles, 3)
    assert topology.trianglePointIndices.reshape(-1, order="F").size == 3 * numberOfTriangles
    assert np.min(topology.trianglePointIndices) >= 0
    assert np.max(topology.trianglePointIndices) < numberOfPoints

    assert numberOfLevels >= 2
    assert topology.thickness > 0.0
    assert topology.levelCoordinates().shape == (numberOfLevels,)
    assert topology.numberOfPrisms == numberOfTriangles * (numberOfLevels - 1)

    derived = topology._topology()
    assert derived["triangleCenterX"].shape == (numberOfTriangles,)
    assert derived["triangleCenterY"].shape == (numberOfTriangles,)
    assert derived["triangleSurfaces"].shape == (numberOfTriangles,)
    assert np.all(derived["triangleSurfaces"] > 0.0)

    for name in ("triangleNormalsX", "triangleNormalsY"):
        values = np.asarray(derived[name])
        assert values.shape == (3 * numberOfTriangles,)
        assert np.isfinite(values).all()

    for name in ("triangleNeighbors", "forbiddenEdge", "triangleNormalPoint"):
        assert np.asarray(derived[name]).shape == (3 * numberOfTriangles,)

    neighbors = np.asarray(derived["triangleNeighbors"])
    assert np.all((neighbors == -1) | ((0 <= neighbors) & (neighbors < numberOfTriangles)))

    forbiddenEdges = np.asarray(derived["forbiddenEdge"])
    assert np.all((forbiddenEdges == -1) | ((0 <= forbiddenEdges) & (forbiddenEdges <= 2)))

    normalPoints = np.asarray(derived["triangleNormalPoint"])
    assert np.all((0 <= normalPoints) & (normalPoints < numberOfPoints))


def _assertGainMediumLaunchContract(medium):
    topology = medium.topology
    numberOfPoints = topology.numberOfPoints
    numberOfTriangles = topology.numberOfTriangles
    numberOfLevels = int(topology.levels)

    assert medium.get("betaCells").expectedShape == (numberOfPoints, numberOfLevels)
    assert medium.get("betaVolume").expectedShape == (numberOfTriangles, numberOfLevels - 1)
    assert medium.get("claddingCellTypes").expectedShape == (numberOfTriangles,)
    assert medium.get("reflectivities").expectedShape == (numberOfTriangles, 2)
    assert medium.get("refractiveIndices").expectedShape == (4,)

    betaCells = np.asarray(medium.get("betaCells").value)
    betaVolume = np.asarray(medium.get("betaVolume").value)
    cladding = np.asarray(medium.get("claddingCellTypes").value)
    reflectivities = np.asarray(medium.get("reflectivities").value)
    refractive = np.asarray(medium.get("refractiveIndices").value)

    assert betaCells.size == numberOfPoints * numberOfLevels
    assert betaVolume.size == numberOfTriangles * (numberOfLevels - 1)
    assert betaVolume.size == topology.numberOfPrisms
    assert cladding.size == numberOfTriangles
    assert reflectivities.size == 2 * numberOfTriangles
    assert refractive.size == 4

    assert betaCells.reshape((numberOfPoints, numberOfLevels), order="F").shape == (
        numberOfPoints,
        numberOfLevels,
    )
    assert betaVolume.reshape((numberOfTriangles, numberOfLevels - 1), order="F").shape == (
        numberOfTriangles,
        numberOfLevels - 1,
    )


def _mediumForGrid(grid, *, flat):
    topology = MeshTopology.fromGrid(grid)
    medium = GainMedium(topology=topology)

    betaCells = np.linspace(0.0, 1.0, topology.numberOfPoints * topology.levels, dtype=np.float64)
    betaVolume = np.linspace(0.0, 1.0, topology.numberOfPrisms, dtype=np.float64)
    cladding = np.zeros(topology.numberOfTriangles, dtype=np.uint32)
    reflectivities = np.zeros(2 * topology.numberOfTriangles, dtype=np.float32)

    if not flat:
        betaCells = betaCells.reshape((topology.numberOfPoints, topology.levels), order="F")
        betaVolume = betaVolume.reshape((topology.numberOfTriangles, topology.levels - 1), order="F")
        reflectivities = reflectivities.reshape((topology.numberOfTriangles, 2), order="F")

    medium.withPhysicalProperties(
        betaCells=backendFlat(betaCells) if flat else betaCells,
        betaVolume=backendFlat(betaVolume) if flat else betaVolume,
        claddingCellTypes=cladding,
        refractiveIndices=[1.8, 1.0, 1.8, 1.0],
        reflectivities=backendFlat(reflectivities) if flat else reflectivities,
        nTot=1.0,
        crystalTFluo=1.0,
        claddingNumber=0,
        claddingAbsorption=0.0,
    )
    return medium


@pytest.mark.parametrize(
    ("grid", "flat"),
    [
        (Grid(xExtent=1.0, yExtent=1.0, zExtent=0.5, tileSizeX=1.0, tileSizeZ=0.25), False),
        (Grid(xExtent=2.0, yExtent=1.0, zExtent=1.0, tileSizeX=0.5, tileSizeY=1.0, tileSizeZ=0.5), True),
        (Grid(xExtent=3.0, yExtent=2.0, zExtent=1.5, tileSizeX=1.0, tileSizeY=0.5, tileSizeZ=0.75), False),
    ],
)
def testGridGeometrySatisfiesPhiAseGainMediumShapeContracts(grid, flat):
    medium = _mediumForGrid(grid, flat=flat)

    _assertTopologyLaunchContract(medium.topology)
    _assertGainMediumLaunchContract(medium)


def testPointIndexAtFallsBackToNearestTopologyPoint():
    topology = MeshTopology.fromGrid(Grid(xExtent=1.0, yExtent=1.0, zExtent=1.0, tileSizeX=0.5))

    exact = topology.pointIndexAt(0.5, 0.5)
    nearest = topology.pointIndexAt(0.49, 0.49)

    assert nearest == exact



def testGainMediumPrimitiveViewsShareCanonicalFlatStorage():
    topology = MeshTopology.fromGrid(Grid(xExtent=1.0, yExtent=1.0, zExtent=0.5, tileSizeZ=0.25))
    medium = GainMedium(topology=topology).withPhysicalProperties(
        betaCells=backendFlat(np.arange(topology.numberOfPoints * topology.levels, dtype=np.float64)),
        betaVolume=backendFlat(np.array([0.1, 0.3, 0.2, 0.4], dtype=np.float64)),
        claddingCellTypes=np.array([0, 1], dtype=np.uint32),
        refractiveIndices=[1.8, 1.0, 1.8, 1.0],
        reflectivities=backendFlat(np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float32)),
    )

    prismView = medium.getPrisms()["betaVolume"]
    pointView = medium.getPoints()["betaCells"]
    triangleView = medium.getTriangles()["reflectivities"]

    assert prismView.shape == (topology.numberOfTriangles, topology.levels - 1)
    assert pointView.shape == (topology.numberOfPoints, topology.levels)
    assert triangleView.shape == (topology.numberOfTriangles, 2)
    assert np.shares_memory(prismView, medium.get("betaVolume").value)
    assert np.shares_memory(pointView, medium.get("betaCells").value)
    assert np.shares_memory(triangleView, medium.get("reflectivities").value)

    prismView[1, 1] = 9.0
    pointView[2, 1] = 8.0
    triangleView[0, 1] = 0.75

    assert medium.get("betaVolume").value[1 + topology.numberOfTriangles] == 9.0
    assert medium.get("betaCells").value[2 + topology.numberOfPoints] == 8.0
    assert medium.get("reflectivities").value[topology.numberOfTriangles] == np.float32(0.75)


def testGainMediumRejectsAmbiguousFlatArraysUnlessMarkedBackendFlat():
    topology = MeshTopology.fromGrid(Grid(xExtent=1.0, yExtent=1.0, zExtent=0.5, tileSizeZ=0.25))
    medium = GainMedium(topology=topology)
    flat = np.arange(topology.numberOfPrisms, dtype=np.float64)

    with pytest.raises(ValueError, match="ambiguous flat array"):
        medium.withPhysicalProperties(betaVolume=flat)

    medium.withPhysicalProperties(betaVolume=backendFlat(flat))
    np.testing.assert_array_equal(medium.get("betaVolume").value, flat)


def testGainMediumRejectsArraysThatBreakGeometryDependencies():
    topology = MeshTopology.fromGrid(Grid(xExtent=1.0, yExtent=1.0, zExtent=0.5, tileSizeZ=0.25))
    medium = GainMedium(topology=topology)

    with pytest.raises(ValueError, match="betaCells expects"):
        medium.withPhysicalProperties(betaCells=np.zeros(topology.numberOfPoints * topology.levels - 1))

    with pytest.raises(ValueError, match="betaVolume expects"):
        medium.withPhysicalProperties(betaVolume=np.zeros(topology.numberOfPrisms + 1))

    with pytest.raises(ValueError, match="claddingCellTypes expects"):
        medium.withPhysicalProperties(claddingCellTypes=np.zeros(topology.numberOfTriangles + 1, dtype=np.uint32))

    with pytest.raises(ValueError, match="reflectivities expects"):
        medium.withPhysicalProperties(reflectivities=np.zeros((topology.numberOfTriangles + 1, 2)))



def testGainMediumDefineFieldAttachesCustomFieldsToPrimitiveViews():
    topology = MeshTopology.fromGrid(Grid(xExtent=1.0, yExtent=1.0, zExtent=0.5, tileSizeZ=0.25))
    medium = GainMedium(topology=topology)
    temperature = np.array(
        [[300.0, 301.0], [310.0, 311.0]],
        dtype=np.float64,
    )
    marker = np.array([4, 5], dtype=np.int32)

    medium.defineField(
        "temperature",
        entity=("cell", "layer"),
        values=temperature,
        unit="K",
        unitDimension=unitDimension.tFluo,
    )
    medium.defineField("triangleMarker", entity="triangle", values=marker, dtype=np.int32, unit="1")

    prismView = medium.getPrisms()["temperature"]
    triangleView = medium.getTriangles()["triangleMarker"]

    assert prismView.shape == (topology.numberOfTriangles, topology.levels - 1)
    assert triangleView.shape == (topology.numberOfTriangles,)
    np.testing.assert_array_equal(prismView, temperature)
    np.testing.assert_array_equal(triangleView, marker)
    assert np.shares_memory(prismView, medium.customFields["temperature"].value)
    assert medium.customFields["temperature"].spec.unitDimension == unitDimension.tFluo

    prismView[1, 1] = 333.0
    assert medium.customFields["temperature"].value[1 + topology.numberOfTriangles] == 333.0
    assert medium.customFields["temperature"].dynamic is False
    assert medium.customFields["temperature"].spec.backendRequired is False


def testGainMediumDefineFieldValidatesNameEntityDtypeUnitAndFlatOrder():
    topology = MeshTopology.fromGrid(Grid(xExtent=1.0, yExtent=1.0, zExtent=0.5, tileSizeZ=0.25))
    medium = GainMedium(topology=topology)

    with pytest.raises(ValueError, match="field name"):
        medium.defineField("1bad", entity="point", values=np.zeros(topology.numberOfPoints))

    with pytest.raises(ValueError, match="unknown field entity"):
        medium.defineField("foo", entity="unknown", values=np.zeros(topology.numberOfPoints))

    with pytest.raises(ValueError, match="non-empty string"):
        medium.defineField("foo", entity="point", values=np.zeros(topology.numberOfPoints), unit="")

    with pytest.raises(ValueError, match="unitDimension"):
        medium.defineField("foo", entity="point", values=np.zeros(topology.numberOfPoints), unitDimension=(0.0, 0.0))

    with pytest.raises(ValueError, match="not supported"):
        medium.defineField("labels", entity="point", values=np.array(["a"] * topology.numberOfPoints))

    flat = np.arange(topology.numberOfTriangles * (topology.levels - 1), dtype=np.float64)
    with pytest.raises(ValueError, match="ambiguous flat array"):
        medium.defineField("temperature", entity="prism", values=flat, unit="K")

    medium.defineField("temperature", entity="prism", values=backendFlat(flat), unit="K", dynamic=True)
    assert medium.customFields["temperature"].dynamic is True
    np.testing.assert_array_equal(medium.getPrisms()["temperature"].reshape(-1, order="F"), flat)
