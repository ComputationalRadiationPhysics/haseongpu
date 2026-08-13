# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path

import numpy as np
import pytest

from HASEonGPU import GainMedium, MeshTopology, SurfaceOptics, VolumeTopology


repoRoot = Path(__file__).resolve().parents[3]


def _smallVolumeGainMedium():
    topology = VolumeTopology.fromTetrahedra(
        np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        ),
        np.asarray([[0, 1, 2, 3]], dtype=np.uint32),
        faceBoundaries=np.asarray([[7, -1, -1, -1]], dtype=np.int32),
        metadata={"structured": {"numberOfPoints": 1, "numberOfLevels": 1, "thickness": 0.0}},
    )
    return GainMedium(topology=topology).withPhysicalProperties(
        betaVolume=np.asarray([0.5], dtype=np.float64),
        claddingCellTypes=np.asarray([2], dtype=np.uint32),
        nTot=1.0,
        crystalTFluo=2.0,
        claddingNumber=3,
        claddingAbsorption=4.0,
    ).with_surface_optics(
        {7: SurfaceOptics(reflectivity=0.7, n_inside=1.2, n_outside=1.0)}
    )


def test_gainMediumRoundTripsThroughTet4Vtk(tmp_path):
    medium = _smallVolumeGainMedium()
    vtk_path = tmp_path / "gain_medium.vtk"

    medium.toVtk(vtk_path)
    loaded = GainMedium.fromVtk(vtk_path)

    assert loaded.topology.numberOfPoints == medium.topology.numberOfPoints
    assert loaded.topology.numberOfCells == medium.topology.numberOfCells
    assert np.allclose(loaded.topology.points, medium.topology.points)
    assert np.array_equal(loaded.topology.cellPointIndices, medium.topology.cellPointIndices)
    assert np.array_equal(loaded.topology.faceBoundaries, medium.topology.faceBoundaries)
    assert np.allclose(loaded.get("betaVolume").value, medium.get("betaVolume").value)
    assert np.array_equal(loaded.get("claddingCellTypes").value, medium.get("claddingCellTypes").value)
    assert set(loaded.surface_optics) == {7}
    assert loaded.surface_optics[7].reflectivity == pytest.approx(0.7)
    assert loaded.surface_optics[7].n_inside == pytest.approx(1.2)
    assert loaded.surface_optics[7].n_outside == pytest.approx(1.0)
    assert loaded.get("nTot").value == medium.get("nTot").value
    assert loaded.get("crystalTFluo").value == medium.get("crystalTFluo").value
    assert loaded.get("claddingNumber").value == medium.get("claddingNumber").value
    assert loaded.get("claddingAbsorption").value == medium.get("claddingAbsorption").value


def test_meshTopologyRejectsVtkInput():
    with pytest.raises(NotImplementedError, match="Tet4 VolumeTopology"):
        MeshTopology.fromFile(repoRoot / "example" / "data" / "ptTet4.vtk")


def test_gainMediumRejectsWedgeVtkInput(tmp_path):
    vtk_path = tmp_path / "wedge.vtk"
    vtk_path.write_text(
        "\n".join(
            [
                "# vtk DataFile Version 2.0",
                "legacy wedge",
                "ASCII",
                "DATASET UNSTRUCTURED_GRID",
                "POINTS 6 double",
                "0 0 0  1 0 0  0 1 0  0 0 1  1 0 1  0 1 1",
                "CELLS 1 7",
                "6 0 1 2 3 4 5",
                "CELL_TYPES 1",
                "13",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Tet4 cells"):
        GainMedium.fromVtk(vtk_path)


def test_gainMediumRejectsNonVolumeTopologyVtkWrite(tmp_path, smallGainMedium):
    with pytest.raises(TypeError, match="Tet4 VolumeTopology"):
        smallGainMedium.toVtk(tmp_path / "legacy.vtk")


def test_bundledExampleVtkFixturesExposeFrontendFields():
    fixtures = {
        "ptTet4.vtk": (421, 21924, 10),
        "cuboid.vtk": (321, 16200, 10),
        "cylindrical.vtk": (421, 21924, 10),
    }

    for filename, (points, triangles, levels) in fixtures.items():
        medium = GainMedium.fromVtk(repoRoot / "example" / "data" / filename)

        assert medium.numberOfPoints == points
        assert medium.numberOfTriangles == triangles
        assert medium.numberOfLevels == levels
        assert np.asarray(medium.get("betaVolume").value).size == triangles
        assert np.asarray(medium.get("claddingCellTypes").value).shape == (triangles,)
        with pytest.raises(KeyError, match="unknown gain medium property"):
            medium.get("refractiveIndices")
        with pytest.raises(KeyError, match="unknown gain medium property"):
            medium.get("reflectivities")
        assert np.isfinite(medium.get("nTot").value)
        assert np.isfinite(medium.get("crystalTFluo").value)
        assert medium.get("claddingNumber").value >= 1
        assert np.isfinite(medium.get("claddingAbsorption").value)
