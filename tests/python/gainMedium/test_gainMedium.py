# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np
from HASEonGPU import GainMedium, Grid, MeshTopology


def testMeshTopologyFromPointsConstructsTriangles():
    points = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ],
        dtype=np.float64,
    )

    topology = MeshTopology.fromPoints(points).numberOfLevels(3)

    assert topology.numberOfPoints == 4
    assert topology.numberOfTriangles == 2
    assert topology.numberOfPrisms == 4
    assert topology.trianglePointIndices.shape == (2, 3)


def testGridDefersPointConstructionUntilUsed():
    grid = Grid(xExtent=2, yExtent=1, zExtent=1, tileSizeX=0.5)

    points = grid.constructPoints()
    topology = MeshTopology.fromGrid(grid)

    assert points.shape == (15, 2)
    assert np.allclose(points[:, 0].max(), 2.0)
    assert np.allclose(points[:, 1].max(), 1.0)
    assert topology.levels == 3
    assert topology.thickness == 0.5
    assert topology.numberOfPoints == 15
    assert topology.numberOfTriangles == 16


def testGridTopologyTriangulatesXyPlaneOnly():
    shallow = MeshTopology.fromGrid(Grid(xExtent=2, yExtent=1, zExtent=1, tileSizeX=1.0, tileSizeZ=0.5))
    deep = MeshTopology.fromGrid(Grid(xExtent=2, yExtent=1, zExtent=10, tileSizeX=1.0, tileSizeZ=0.25))

    expectedTriangles = np.array(
        [
            [0, 1, 3],
            [1, 2, 4],
            [3, 1, 4],
            [4, 2, 5],
        ],
        dtype=np.uint32,
    )

    assert np.array_equal(shallow.points, deep.points)
    assert np.array_equal(shallow.trianglePointIndices, expectedTriangles)
    assert np.array_equal(deep.trianglePointIndices, expectedTriangles)
    assert shallow.levels != deep.levels
    assert shallow.thickness != deep.thickness


def testMeshTopologyFromAsciiStl(tmp_path):
    stl = tmp_path / "surface.stl"
    stl.write_text(
        """
solid planar
  facet normal 0 0 1
    outer loop
      vertex 0 0 0
      vertex 1 0 0
      vertex 0 1 0
    endloop
  endfacet
endsolid planar
""",
        encoding="utf-8",
    )

    topology = MeshTopology.fromFile(stl, format="stl", numberOfLevels=2)

    assert topology.numberOfPoints == 3
    assert topology.numberOfTriangles == 1
    assert topology.metadata["format"] == "stl"


def testGainMediumOwnsPhysicalProperties():
    topology = MeshTopology.fromGrid(Grid(xExtent=1, yExtent=1, zExtent=2, tileSizeX=1.0, tileSizeZ=0.25))
    gainMedium = (
        GainMedium(topology=topology)
        .withPhysicalProperties(
            betaVolume=np.ones((topology.numberOfTriangles, topology.levels - 1), dtype=np.float64),
            betaCells=np.ones((topology.numberOfPoints, topology.levels), dtype=np.float64),
            claddingCellTypes=np.zeros(2, dtype=np.uint32),
            refractiveIndices=np.array([1.8, 1.0, 1.8, 1.0], dtype=np.float32),
            reflectivities=np.ones(4, dtype=np.float32),
            nTot=5.0,
            crystalTFluo=1.23,
            claddingNumber=2,
            claddingAbsorption=0.05,
        )
    )

    assert gainMedium.numberOfPrisms == 16
    assert not hasattr(gainMedium, "toHostMesh")


def testGainMediumPhysicalPropertiesAreDiscoverableAndAssignable():
    topology = MeshTopology.fromGrid(Grid(xExtent=1, yExtent=1, zExtent=2, tileSizeZ=0.25))
    gainMedium = GainMedium(topology=topology)

    betaVolume = gainMedium.get("betaVolume")
    refractiveIndices = gainMedium.get("refractiveIndicies")

    assert betaVolume.expectedShape == (2, 8)
    assert betaVolume.dtype == np.dtype(np.float64)
    assert refractiveIndices.name == "refractiveIndices"
    assert refractiveIndices.expectedShape == (4,)

    betaVolume.value = np.ones(betaVolume.expectedShape)
    refractiveIndices.value = [1.8, 1.0, 1.8, 1.0]

    assert gainMedium.physical["betaVolume"].shape == (16,)
    assert np.allclose(gainMedium.get("refractiveIndices").value, [1.8, 1.0, 1.8, 1.0])
    assert any(item["name"] == "reflectivities" for item in gainMedium.listProperties())
