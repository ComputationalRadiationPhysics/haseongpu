# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import pytest

from HASEonGPU import Gmsh, MeshTopology, VolumeTopology
from pyInclude.geometry import GmshElement
import pyInclude.openpmd.transport as transport
from pyInclude.geometry.volume import BOUND_STOP, GMSH_PRISM6, GMSH_QUAD4, GMSH_TET4, VTK_WEDGE

try:
    import gmsh as gmsh_api
except (ImportError, OSError) as exc:
    gmsh_api = None
    GMSH_SKIP_REASON = f"gmsh is not importable: {exc}"
else:
    GMSH_SKIP_REASON = ""


def _two_prism_points():
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )


def testVolumeTopologyBuildsExplicitPrismNeighborsWithoutLayoutOrdering():
    # Two prisms share the quad face {1,2,4,5}.  The neighbor relation is found
    # by explicit face keys, not by any triangle/layer index arithmetic.
    topology = VolumeTopology.fromPrisms(
        _two_prism_points(),
        np.array(
            [
                [1, 6, 2, 4, 7, 5],
                [0, 1, 2, 3, 4, 5],
            ],
            dtype=np.uint32,
        ),
        cellDomains=np.array([7, 3], dtype=np.int32),
    )

    assert topology.numberOfCells == 2
    assert topology.numberOfFacesPerCell == 5
    assert np.all(topology.cellTypes == VTK_WEDGE)
    assert np.array_equal(topology.cellDomains, np.array([7, 3], dtype=np.int32))

    shared_a = np.argwhere(topology.neighborCells[0] == 1)
    shared_b = np.argwhere(topology.neighborCells[1] == 0)
    assert shared_a.shape == (1, 1)
    assert shared_b.shape == (1, 1)
    face_a = int(shared_a[0, 0])
    face_b = int(shared_b[0, 0])
    assert topology.neighborLocalFaces[0, face_a] == face_b
    assert topology.neighborLocalFaces[1, face_b] == face_a
    assert topology.faceBoundaries[0, face_a] == 0
    assert topology.faceBoundaries[1, face_b] == 0
    assert np.count_nonzero(topology.faceBoundaries == BOUND_STOP) == 8
    assert np.all(topology.cellVolumes > 0.0)


def testVolumeTopologyFromLegacyExtrudedMeshCreatesExplicitWedges():
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    triangles = np.array([[0, 1, 2]], dtype=np.uint32)
    legacy = MeshTopology(points, triangles, levels=3, thickness=0.5)

    topology = VolumeTopology.fromExtrudedTopology(legacy)

    assert topology.numberOfPoints == 9
    assert topology.numberOfCells == 2
    np.testing.assert_array_equal(
        topology.cellPointIndices,
        np.array([[0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]], dtype=np.uint32),
    )
    assert np.count_nonzero(topology.neighborCells[0] == 1) == 1
    assert np.count_nonzero(topology.neighborCells[1] == 0) == 1


def testVolumeTopologyImportsSyntheticGmshPrismPhysicalGroups():
    nodes = {
        10: (0.0, 0.0, 0.0),
        11: (1.0, 0.0, 0.0),
        12: (0.0, 1.0, 0.0),
        13: (0.0, 0.0, 1.0),
        14: (1.0, 0.0, 1.0),
        15: (0.0, 1.0, 1.0),
    }
    gmsh = Gmsh(
        nodes=nodes,
        elements=[
            GmshElement(1, GMSH_PRISM6, (10, 11, 12, 13, 14, 15), physical_tag=4),
            GmshElement(2, GMSH_QUAD4, (10, 11, 14, 13), physical_tag=30),
        ],
        physical_names={3: {4: "gain"}, 2: {30: "side_wall"}},
        source="synthetic.msh",
    )

    topology = VolumeTopology.fromGmsh(gmsh)

    assert topology.metadata["dimension"] == 3
    assert topology.cellDomains.tolist() == [4]
    assert np.count_nonzero(topology.faceBoundaries == 30) == 1
    assert np.count_nonzero(topology.faceBoundaries == BOUND_STOP) == 4


def testVolumeTopologyImportsGenerated3dGmshPrisms(tmp_path):
    if gmsh_api is None:
        pytest.skip(GMSH_SKIP_REASON)
    msh = tmp_path / "prism3d.msh"

    gmsh_api.initialize()
    try:
        gmsh_api.option.setNumber("General.Terminal", 0)
        gmsh_api.model.add("prism3d")
        geo = gmsh_api.model.geo
        p1 = geo.addPoint(0.0, 0.0, 0.0, 1.0)
        p2 = geo.addPoint(1.0, 0.0, 0.0, 1.0)
        p3 = geo.addPoint(0.0, 1.0, 0.0, 1.0)
        l1 = geo.addLine(p1, p2)
        l2 = geo.addLine(p2, p3)
        l3 = geo.addLine(p3, p1)
        surface = geo.addPlaneSurface([geo.addCurveLoop([l1, l2, l3])])
        extrusion = geo.extrude([(2, surface)], 0.0, 0.0, 1.0, numElements=[1], recombine=True)
        geo.synchronize()
        volumes = [tag for dim, tag in extrusion if dim == 3]
        surfaces = [tag for dim, tag in gmsh_api.model.getEntities(2)]
        gmsh_api.model.addPhysicalGroup(3, volumes, 4)
        gmsh_api.model.setPhysicalName(3, 4, "gain")
        gmsh_api.model.addPhysicalGroup(2, surfaces, 30)
        gmsh_api.model.setPhysicalName(2, 30, "outer")
        gmsh_api.model.mesh.generate(3)
        gmsh_api.write(str(msh))
    finally:
        gmsh_api.finalize()

    topology = VolumeTopology.fromFile(msh)

    assert topology.metadata["dimension"] == 3
    assert topology.numberOfCells >= 1
    assert np.all(topology.cellTypes == VTK_WEDGE)
    assert set(topology.cellDomains.tolist()) == {4}
    assert np.count_nonzero(topology.faceBoundaries == 30) > 0


def testVolumeTopologyRejectsTetGmshVolumesForLaterMilestone():
    gmsh = Gmsh(
        nodes={1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (0.0, 1.0, 0.0), 4: (0.0, 0.0, 1.0)},
        elements=[GmshElement(1, GMSH_TET4, (1, 2, 3, 4), physical_tag=1)],
    )

    with pytest.raises(NotImplementedError, match="Prism6 only"):
        VolumeTopology.fromGmsh(gmsh)


def testExplicitOpenPmdTopologySpecsUseVolumeCellShapes():
    topology = VolumeTopology.fromPrisms(_two_prism_points(), np.array([[0, 1, 2, 3, 4, 5]], dtype=np.uint32))
    context = transport._explicit_topology_context(topology)

    assert transport.CANONICAL_CONNECTIVITY_SPEC.expectedShape(context) == (1, 6)
    assert transport.EXPLICIT_CELL_FACES_SPEC.expectedShape(context) == (1, 5, 4)
    assert transport.EXPLICIT_CELL_NEIGHBORS_SPEC.expectedShape(context) == (1, 5)
    assert transport.EXPLICIT_SAMPLE_POINTS_SPEC.expectedShape(context) == (3, 1)
    np.testing.assert_array_equal(topology.cellsOffsets(), np.array([0, 6], dtype=np.uint32))
    assert topology.faceConnectivityFlat().shape == (20,)


def _read_openpmd_scalar(series, iteration, name):
    io = transport._io()
    chunk = iteration.meshes[name][io.Mesh_Record_Component.SCALAR].load_chunk()
    series.flush()
    return np.array(chunk, copy=True).reshape(-1)


def testExplicitOpenPmdStaticTopologyWriterStoresFaceLookupTables(tmp_path):
    topology = VolumeTopology.fromPrisms(_two_prism_points(), np.array([[0, 1, 2, 3, 4, 5]], dtype=np.uint32))
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
        assert iteration.meshes["core_points"].get_attribute("geometryParameters") == "topology=explicit_volume_cells"
        assert iteration.meshes["core_sample_points"].get_attribute("geometryParameters") == "topology=explicit_volume_cells"
        np.testing.assert_array_equal(
            _read_openpmd_scalar(series, iteration, "core_cell_faces"),
            topology.facePointIndices.reshape(-1, order="F"),
        )
        np.testing.assert_array_equal(
            _read_openpmd_scalar(series, iteration, "core_cell_neighbor_cells"),
            topology.neighborCells.reshape(-1, order="F"),
        )
        np.testing.assert_array_equal(
            _read_openpmd_scalar(series, iteration, "core_cell_face_boundaries"),
            topology.faceBoundaries.reshape(-1, order="F"),
        )
    finally:
        series.close()
