# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import pytest

from HASEonGPU import Gmsh, VolumeTopology
from pyInclude.geometry import GmshElement
import pyInclude.openpmd.transport as transport
from pyInclude.geometry.volume import BOUND_STOP, GMSH_TET4, GMSH_TRI3, VTK_TETRA

try:
    import gmsh as gmshApi
except (ImportError, OSError) as exc:
    gmshApi = None
    GMSH_SKIP_REASON = f"gmsh is not importable: {exc}"
else:
    GMSH_SKIP_REASON = ""


def _twoTetPoints():
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )


def _oneTetTopology():
    return VolumeTopology.fromTetrahedra(_twoTetPoints()[:4], np.array([[0, 1, 2, 3]], dtype=np.uint32))


def testVolumeTopologyBuildsExplicitTetNeighborsWithoutLayoutOrdering():
    topology = VolumeTopology.fromTetrahedra(
        _twoTetPoints(),
        np.array(
            [
                [0, 1, 2, 3],
                [1, 2, 3, 4],
            ],
            dtype=np.uint32,
        ),
        cellDomains=np.array([7, 3], dtype=np.int32),
    )

    assert topology.numberOfCells == 2
    assert topology.numberOfFacesPerCell == 4
    assert np.all(topology.cellTypes == VTK_TETRA)
    assert np.array_equal(topology.cellDomains, np.array([7, 3], dtype=np.int32))

    sharedA = np.argwhere(topology.neighborCells[0] == 1)
    sharedB = np.argwhere(topology.neighborCells[1] == 0)
    assert sharedA.shape == (1, 1)
    assert sharedB.shape == (1, 1)
    faceA = int(sharedA[0, 0])
    faceB = int(sharedB[0, 0])
    assert topology.neighborLocalFaces[0, faceA] == faceB
    assert topology.neighborLocalFaces[1, faceB] == faceA
    assert topology.faceBoundaries[0, faceA] == 0
    assert topology.faceBoundaries[1, faceB] == 0
    assert np.count_nonzero(topology.faceBoundaries == BOUND_STOP) == 6
    assert np.all(topology.cellVolumes > 0.0)


def testVolumeTopologyImportsSyntheticGmshTetPhysicalGroups():
    nodes = {
        10: (0.0, 0.0, 0.0),
        11: (1.0, 0.0, 0.0),
        12: (0.0, 1.0, 0.0),
        13: (0.0, 0.0, 1.0),
    }
    gmsh = Gmsh(
        nodes=nodes,
        elements=[
            GmshElement(1, GMSH_TET4, (10, 11, 12, 13), physical_tag=4),
            GmshElement(2, GMSH_TRI3, (10, 11, 12), physical_tag=30),
        ],
        physical_names={3: {4: "gain"}, 2: {30: "side_wall"}},
        source="synthetic.msh",
    )

    topology = VolumeTopology.fromGmsh(gmsh)

    assert topology.metadata["dimension"] == 3
    assert topology.cellDomains.tolist() == [4]
    assert np.count_nonzero(topology.faceBoundaries == 30) == 1
    assert np.count_nonzero(topology.faceBoundaries == BOUND_STOP) == 3


def testVolumeTopologyImportsGenerated3dGmshTets(tmp_path):
    if gmshApi is None:
        pytest.skip(GMSH_SKIP_REASON)
    msh = tmp_path / "tet4_3d.msh"

    gmshApi.initialize()
    try:
        gmshApi.option.setNumber("General.Terminal", 0)
        gmshApi.model.add("tet4_3d")
        box = gmshApi.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        gmshApi.model.occ.synchronize()
        gmshApi.model.addPhysicalGroup(3, [box], 4)
        gmshApi.model.setPhysicalName(3, 4, "gain")
        surfaces = [tag for dim, tag in gmshApi.model.getEntities(2)]
        gmshApi.model.addPhysicalGroup(2, surfaces, 30)
        gmshApi.model.setPhysicalName(2, 30, "outer")
        gmshApi.model.mesh.generate(3)
        gmshApi.write(str(msh))
    finally:
        gmshApi.finalize()

    topology = VolumeTopology.fromFile(msh)

    assert topology.metadata["dimension"] == 3
    assert topology.numberOfCells >= 1
    assert np.all(topology.cellTypes == VTK_TETRA)
    assert set(topology.cellDomains.tolist()) == {4}
    assert np.count_nonzero(topology.faceBoundaries == 30) > 0


def testVolumeTopologyImportsTet4VtkLookupTables(tmp_path):
    vtkPath = tmp_path / "tet4.vtk"
    vtkPath.write_text(
        "\n".join(
            [
                "# vtk DataFile Version 2.0",
                "HASEonGPU Tet4 test",
                "ASCII",
                "DATASET UNSTRUCTURED_GRID",
                "POINTS 4 double",
                "0 0 0  1 0 0  0 1 0  0 0 1",
                "CELLS 1 5",
                "4 0 1 2 3",
                "CELL_TYPES 1",
                "10",
                "FIELD HASEonGPU 2",
                "cellDomains 1 1 int",
                "7",
                "faceBoundaries 4 1 int",
                "-1 -1 -1 30",
            ]
        ),
        encoding="utf-8",
    )

    topology = VolumeTopology.fromFile(vtkPath)

    assert topology.numberOfCells == 1
    np.testing.assert_array_equal(topology.cellDomains, np.array([7], dtype=np.int32))
    np.testing.assert_array_equal(topology.faceBoundaries, np.array([[-1, -1, -1, 30]], dtype=np.int32))


def testVolumeTopologyRejectsUnsupportedGmshVolumes():
    gmsh = Gmsh(
        nodes={
            1: (0.0, 0.0, 0.0),
            2: (1.0, 0.0, 0.0),
            3: (0.0, 1.0, 0.0),
            4: (0.0, 0.0, 1.0),
            5: (1.0, 0.0, 1.0),
            6: (0.0, 1.0, 1.0),
        },
        elements=[GmshElement(1, 6, (1, 2, 3, 4, 5, 6), physical_tag=1)],
    )

    with pytest.raises(NotImplementedError, match="Tet4 only"):
        VolumeTopology.fromGmsh(gmsh)


def testExplicitOpenPmdTopologySpecsUseTet4Shapes():
    topology = _oneTetTopology()
    context = transport._explicit_topology_context(topology)

    assert transport.CANONICAL_CONNECTIVITY_SPEC.expectedShape(context) == (1, 4)
    assert transport.EXPLICIT_CELL_FACES_SPEC.expectedShape(context) == (1, 4, 3)
    assert transport.EXPLICIT_CELL_NEIGHBORS_SPEC.expectedShape(context) == (1, 4)
    assert transport.EXPLICIT_SAMPLE_POINTS_SPEC.expectedShape(context) == (3, 1)
    np.testing.assert_array_equal(topology.cellsOffsets(), np.array([0, 4], dtype=np.uint32))
    assert topology.faceConnectivityFlat().shape == (12,)


def _readOpenPmdScalar(series, iteration, name):
    io = transport._io()
    chunk = iteration.meshes[name][io.Mesh_Record_Component.SCALAR].load_chunk()
    series.flush()
    return np.array(chunk, copy=True).reshape(-1)


def testExplicitOpenPmdStaticTopologyWriterStoresFaceLookupTables(tmp_path):
    topology = _oneTetTopology()
    path = tmp_path / ("explicit_volume" + transport._backend_spec("adios").suffix)

    series = transport._open_input_series(path, backend="adios")
    iteration = series.snapshots()[0]
    try:
        transport._write_explicit_static_topology(iteration, topology)
        iteration.close()
    finally:
        series.close()

    io = transport._io()
    series = io.Series(str(path), io.Access.read_only)
    iteration = series.iterations[0]
    try:
        assert "core_points" in iteration.meshes
        assert "core_sample_points" in iteration.meshes
        assert "core_cell_faces" in iteration.meshes
        assert "core_cell_neighbor_cells" in iteration.meshes
        assert "core_cell_neighbor_local_faces" in iteration.meshes
        assert "core_cell_face_boundaries" in iteration.meshes
        assert "core_cell_domains" in iteration.meshes
        assert iteration.meshes["core_points"].get_attribute("geometryParameters") == "topology=explicit_tet4_volume"
        assert iteration.meshes["core_sample_points"].get_attribute("geometryParameters") == "topology=explicit_tet4_volume"
        np.testing.assert_array_equal(_readOpenPmdScalar(series, iteration, "core_cell_faces"), topology.facePointIndices.reshape(-1))
        np.testing.assert_array_equal(_readOpenPmdScalar(series, iteration, "core_cell_neighbor_cells"), topology.neighborCells.reshape(-1))
        np.testing.assert_array_equal(_readOpenPmdScalar(series, iteration, "core_cell_face_boundaries"), topology.faceBoundaries.reshape(-1))
    finally:
        series.close()
