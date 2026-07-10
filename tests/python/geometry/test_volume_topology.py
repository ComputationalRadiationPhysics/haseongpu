# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
from pathlib import Path

import numpy as np
import pytest

repoRoot = Path(__file__).resolve().parents[3]


def _module_from_repo(module):
    filename = getattr(module, "__file__", None)
    if filename is None:
        return True
    try:
        Path(filename).resolve().relative_to(repoRoot)
        return True
    except ValueError:
        return False


for module_name, module in list(sys.modules.items()):
    if (
        module_name == "HASEonGPU"
        or module_name.startswith("HASEonGPU.")
        or module_name == "pyInclude"
        or module_name.startswith("pyInclude.")
    ) and not _module_from_repo(module):
        del sys.modules[module_name]

if str(repoRoot) not in sys.path:
    sys.path.insert(0, str(repoRoot))

from HASEonGPU import AlpakaBackends, PhiASE, SpectralDecomposition
from pyInclude.geometry import GainMedium, Gmsh, GmshElement, SurfaceOptics, VolumeTopology
import pyInclude.openpmd.transport as transport
from pyInclude.geometry.volume import BOUND_STOP, GMSH_TET4, GMSH_TRI3, VTK_TETRA

try:
    import gmsh as gmshApi
except (ImportError, OSError) as exc:
    gmshApi = None
    GMSH_IMPORT_ERROR = f"gmsh is required for gmsh-backed volume topology tests: {exc}"
else:
    GMSH_IMPORT_ERROR = ""


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


def _slabTopology():
    points = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )
    cells = np.asarray(
        [
            [0, 1, 2, 3],
            [1, 2, 4, 3],
            [2, 4, 5, 3],
        ],
        dtype=np.uint32,
    )
    return VolumeTopology.fromTetrahedra(points, cells)


def _writeClosedCubeStl(path):
    vertices = {
        "000": (0.0, 0.0, 0.0),
        "100": (1.0, 0.0, 0.0),
        "110": (1.0, 1.0, 0.0),
        "010": (0.0, 1.0, 0.0),
        "001": (0.0, 0.0, 1.0),
        "101": (1.0, 0.0, 1.0),
        "111": (1.0, 1.0, 1.0),
        "011": (0.0, 1.0, 1.0),
    }
    triangles = [
        ("000", "010", "110"),
        ("000", "110", "100"),
        ("001", "101", "111"),
        ("001", "111", "011"),
        ("000", "100", "101"),
        ("000", "101", "001"),
        ("010", "011", "111"),
        ("010", "111", "110"),
        ("000", "001", "011"),
        ("000", "011", "010"),
        ("100", "110", "111"),
        ("100", "111", "101"),
    ]
    lines = ["solid cube"]
    for triangle in triangles:
        lines.extend(["  facet normal 0 0 0", "    outer loop"])
        for vertex in triangle:
            x, y, z = vertices[vertex]
            lines.append(f"      vertex {x} {y} {z}")
        lines.extend(["    endloop", "  endfacet"])
    lines.append("endsolid cube")
    path.write_text("\n".join(lines), encoding="utf-8")


def _runtimeAlpakaBackend():
    backends = AlpakaBackends.all()
    for preferred in ("Host_Cpu_CpuOmpBlocks", "Host_Cpu_CpuSerial"):
        if preferred in backends:
            return preferred
    for backend in backends:
        if "Cpu" in backend:
            return backend
    if backends:
        return backends[0]
    pytest.skip("no Alpaka backend is available in this build")


def _runtimeOpenPmdBackend():
    from openpmd_backend_matrix import openpmd_runtime_backend

    try:
        return openpmd_runtime_backend()
    except RuntimeError as exc:
        pytest.skip(str(exc))


def _requireRuntimeBackendExecutable(monkeypatch):
    try:
        executable = transport.findCalcPhiAse()
    except FileNotFoundError as exc:
        pytest.skip(str(exc))
    monkeypatch.setenv("HASE_CPP_EXECUTABLE", str(executable))


def _requireOpenPmdTransportBackend(backend):
    try:
        transport._ensure_backend_available(backend)
    except (FileNotFoundError, ImportError) as exc:
        pytest.skip(str(exc))


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
    assert topology.cellDomainNames == {4: "gain"}
    assert topology.surfaceDomainNames == {30: "side_wall"}


def testVolumeTopologyAssignsCellAndSurfaceDomainsBySelector():
    topology = _slabTopology()
    assigned = topology.withDomains(
        cellDomains={"where": "all", "domain": 6},
        surfaceDomains=[
            {"where": "z_min", "domain": 1},
            {"where": "z_max", "domain": 2},
        ],
    )

    np.testing.assert_array_equal(topology.cellDomains, np.ones(topology.numberOfCells, dtype=np.int32))
    np.testing.assert_array_equal(assigned.cellDomains, np.full(topology.numberOfCells, 6, dtype=np.int32))
    assert np.count_nonzero(assigned.faceBoundaries == 1) > 0
    assert np.count_nonzero(assigned.faceBoundaries == 2) > 0
    assert np.count_nonzero(topology.faceBoundaries > 0) == 0


def testVolumeTopologyRemapsGmshDomainsByName():
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
        physical_names={3: {4: "gain"}, 2: {30: "outer"}},
        source="synthetic.msh",
    )
    topology = VolumeTopology.fromGmsh(gmsh).withDomains(
        cellDomains={"gmshName": "gain", "domain": 9},
        surfaceDomains={"gmshName": "outer", "domain": 11},
    )

    np.testing.assert_array_equal(topology.cellDomains, np.asarray([9], dtype=np.int32))
    assert np.count_nonzero(topology.faceBoundaries == 11) == 1
    assert np.count_nonzero(topology.faceBoundaries == BOUND_STOP) == 3
    assert topology.cellDomainNames[9] == "gain"
    assert topology.surfaceDomainNames[11] == "outer"

    medium = GainMedium(topology).withSurfaceOptics(
        {"outer": SurfaceOptics(reflectivity=0.5, n_inside=1.4, n_outside=1.0)}
    )
    assert medium.get("surfaceReflectivity").expectedShape == (12,)
    assert medium.get("surfaceReflectivity").value[11] == np.float32(0.5)


def testGainMediumSurfaceOpticsUsesAssignedSurfaceDomains():
    topology = _slabTopology().withSurfaceDomains(
        [
            {"where": "z_min", "domain": 1},
            {"where": "z_max", "domain": 2},
        ]
    )
    medium = GainMedium(topology).withSurfaceOptics(
        {
            1: SurfaceOptics(reflectivity=0.0, n_inside=1.83, n_outside=1.0),
            2: SurfaceOptics(reflectivity=0.25, n_inside=1.5, n_outside=1.0),
        }
    )

    np.testing.assert_allclose(medium.get("surfaceReflectivity").value, np.asarray([0.0, 0.0, 0.25], dtype=np.float32))
    np.testing.assert_allclose(
        medium.get("surfaceRefractiveIndexInside").value,
        np.asarray([1.0, 1.83, 1.5], dtype=np.float32),
    )
    np.testing.assert_allclose(
        medium.get("surfaceRefractiveIndexOutside").value,
        np.asarray([1.0, 1.0, 1.0], dtype=np.float32),
    )


@pytest.mark.integration
def testVolumeTopologyImportsClosed3dStlAndRunsBackendOnce(tmp_path, monkeypatch):
    if gmshApi is None:
        pytest.fail(GMSH_IMPORT_ERROR)
    stl = tmp_path / "closed_cube.stl"
    _writeClosedCubeStl(stl)

    with pytest.warns(RuntimeWarning, match="STL volume import assumes a closed 3D surface"):
        topology = VolumeTopology.fromFile(stl, format="stl", meshSize=0.35)

    assert topology.numberOfPrisms >= 10
    assert topology.numberOfCells == topology.numberOfPrisms
    assert np.all(topology.cellTypes == VTK_TETRA)
    _requireRuntimeBackendExecutable(monkeypatch)

    medium = GainMedium(topology=topology).withPhysicalProperties(
        betaVolume=np.zeros(topology.numberOfCells, dtype=np.float64),
        betaCells=np.zeros(topology.numberOfSamplePoints, dtype=np.float64),
        claddingCellTypes=np.zeros(topology.numberOfCells, dtype=np.uint32),
        refractiveIndices=np.array([1.5, 1.0, 1.5, 1.0], dtype=np.float32),
        reflectivities=np.zeros((topology.numberOfCells, 2), dtype=np.float32),
        nTot=1.0,
        crystalTFluo=1.0,
        claddingNumber=99,
        claddingAbsorption=0.0,
    )
    crossSections = SpectralDecomposition.monochromatic(
        wavelength=1000e-9,
        crossSectionAbsorption=0.0,
        crossSectionEmission=0.0,
    )
    phiAse = PhiASE(
        spectralProperties=crossSections,
        minRaysPerSample=1,
        maxRaysPerSample=1,
        forwardRayCount=1,
        forwardRayLength=1.0,
        repetitions=1,
        adaptiveSteps=1,
        mseThreshold=0.25,
        useReflections=False,
        backend=_runtimeAlpakaBackend(),
        openpmdBackend=_runtimeOpenPmdBackend(),
        parallelMode="single",
        numDevices=1,
        minSampleRange=0,
        maxSampleRange=0,
        monochromatic=True,
        rngSeed=1234,
    )
    _requireOpenPmdTransportBackend(phiAse.openpmdBackend)

    phiAse.run(gainMedium=medium)

    result = phiAse.getResults()
    assert result.phiAse.shape == (topology.numberOfSamplePoints,)
    assert result.mse.shape == (topology.numberOfSamplePoints,)
    assert result.totalRays.shape == (topology.numberOfSamplePoints,)
    assert result.dndtAse.shape == (topology.numberOfSamplePoints,)
    assert np.all(np.isfinite(result.phiAse))
    assert np.all(result.totalRays >= 0)


def testVolumeTopologyImportsGenerated3dGmshTets(tmp_path):
    if gmshApi is None:
        pytest.fail(GMSH_IMPORT_ERROR)
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
        assert "core_sample_points" not in iteration.meshes
        assert "core_cell_faces" in iteration.meshes
        assert "core_cell_neighbor_cells" in iteration.meshes
        assert "core_cell_neighbor_local_faces" in iteration.meshes
        assert "core_cell_face_boundaries" in iteration.meshes
        assert "core_cell_domains" in iteration.meshes
        assert iteration.meshes["core_points"].get_attribute("geometryParameters") == "topology=explicit_tet4_volume"
        np.testing.assert_array_equal(_readOpenPmdScalar(series, iteration, "core_cell_faces"), topology.facePointIndices.reshape(-1))
        np.testing.assert_array_equal(_readOpenPmdScalar(series, iteration, "core_cell_neighbor_cells"), topology.neighborCells.reshape(-1))
        np.testing.assert_array_equal(_readOpenPmdScalar(series, iteration, "core_cell_face_boundaries"), topology.faceBoundaries.reshape(-1))
    finally:
        series.close()


def testForwardOpenPmdInputWritesVolumeRecords(tmp_path):
    from HASEonGPU import GainMedium, PhiASE, SpectralDecomposition

    topology = _oneTetTopology()
    medium = GainMedium(topology=topology)
    medium.withPhysicalProperties(betaVolume=np.ones(topology.numberOfCells, dtype=np.float64))
    crossSections = SpectralDecomposition.monochromatic(
        wavelength=1.0,
        crossSectionAbsorption=0.0,
        crossSectionEmission=0.0,
    )
    phiAse = PhiASE(spectralProperties=crossSections, forwardRayLength=1.0)
    path = tmp_path / ("forward_volume" + transport._backend_spec("adios").suffix)

    with transport.OpenPmdInputSeries(path, backend="adios") as series:
        series.write(phiAse, medium, crossSections)

    io = transport._io()
    series = io.Series(str(path), io.Access.read_only)
    iteration = series.iterations[0]
    try:
        assert iteration.get_attribute("propagation_mode") == "forward"
        assert "core_points" in iteration.meshes
        assert "core_beta_volume" in iteration.meshes
    finally:
        series.close()
