# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path

import numpy as np

from HASEonGPU import GainMedium, MeshTopology


repoRoot = Path(__file__).resolve().parents[3]


def testGainMediumRoundTripsThroughLegacyVtk(tmp_path, smallGainMedium):
    vtk_path = tmp_path / "gain_medium.vtk"

    smallGainMedium.toVtk(vtk_path)
    loaded = GainMedium.fromVtk(vtk_path)

    assert loaded.topology.numberOfPoints == smallGainMedium.topology.numberOfPoints
    assert loaded.topology.numberOfTriangles == smallGainMedium.topology.numberOfTriangles
    assert loaded.topology.levels == smallGainMedium.topology.levels
    assert loaded.topology.thickness == smallGainMedium.topology.thickness
    assert np.allclose(loaded.topology.points, smallGainMedium.topology.points)
    assert np.array_equal(loaded.topology.trianglePointIndices, smallGainMedium.topology.trianglePointIndices)
    assert np.allclose(loaded.get("betaCells").value, smallGainMedium.get("betaCells").value.reshape(-1, order="F"))
    assert np.allclose(loaded.get("betaVolume").value, smallGainMedium.get("betaVolume").value.reshape(-1, order="F"))
    assert np.array_equal(loaded.get("claddingCellTypes").value, smallGainMedium.get("claddingCellTypes").value)
    assert np.allclose(loaded.get("refractiveIndices").value, smallGainMedium.get("refractiveIndices").value)
    assert np.allclose(loaded.get("reflectivities").value, smallGainMedium.get("reflectivities").value.reshape(-1, order="F"))
    assert loaded.get("nTot").value == smallGainMedium.get("nTot").value
    assert loaded.get("crystalTFluo").value == smallGainMedium.get("crystalTFluo").value
    assert loaded.get("claddingNumber").value == smallGainMedium.get("claddingNumber").value
    assert loaded.get("claddingAbsorption").value == smallGainMedium.get("claddingAbsorption").value


def testMeshTopologyCanBeLoadedFromVtk(tmp_path, smallGainMedium):
    vtk_path = tmp_path / "topology.vtk"

    smallGainMedium.toVtk(vtk_path)
    topology = MeshTopology.fromFile(vtk_path)

    assert topology.numberOfPoints == smallGainMedium.topology.numberOfPoints
    assert topology.numberOfTriangles == smallGainMedium.topology.numberOfTriangles
    assert topology.levels == smallGainMedium.topology.levels
    assert topology.thickness == smallGainMedium.topology.thickness


def testBundledExampleVtkFixturesExposeFrontendFields():
    fixtures = {
        "pt.vtk": (421, 812, 10),
        "cuboid.vtk": (321, 600, 10),
        "cylindrical.vtk": (421, 812, 10),
    }

    for filename, (points, triangles, levels) in fixtures.items():
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / filename)

        assert medium.numberOfPoints == points
        assert medium.numberOfTriangles == triangles
        assert medium.numberOfLevels == levels
        assert np.asarray(medium.get("betaCells").value).size == points * levels
        assert np.asarray(medium.get("betaVolume").value).size == triangles * (levels - 1)
        assert np.asarray(medium.get("claddingCellTypes").value).shape == (triangles,)
        assert np.asarray(medium.get("refractiveIndices").value).shape == (4,)
        assert np.asarray(medium.get("reflectivities").value).size == triangles * 2
        assert np.isfinite(medium.get("nTot").value)
        assert np.isfinite(medium.get("crystalTFluo").value)
        assert medium.get("claddingNumber").value >= 1
        assert np.isfinite(medium.get("claddingAbsorption").value)
