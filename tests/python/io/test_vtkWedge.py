# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np

from HASEonGPU import Grid, MeshTopology, TimeStepState, vtkWedge


def _scalarNamesByDataKind(path):
    tokens = path.read_text(encoding="utf-8").split()
    names = {"POINT_DATA": set(), "CELL_DATA": set()}
    active = None
    active_count = 0
    index = 0
    while index < len(tokens):
        token = tokens[index].upper()
        if token in names:
            active = token
            active_count = int(tokens[index + 1])
            index += 2
        elif token == "SCALARS" and active is not None:
            names[active].add(tokens[index + 1])
            components = 1
            index += 3
            if tokens[index].upper() != "LOOKUP_TABLE":
                components = int(tokens[index])
                index += 1
            if tokens[index].upper() == "LOOKUP_TABLE":
                index += 2
            index += active_count * components
        else:
            index += 1
    return names


def testVtkWedgeWritesPointDataFromState(tmp_path):
    topology = MeshTopology.fromGrid(Grid(xExtent=1, yExtent=1, zExtent=0.5, tileSizeZ=0.25))
    values = np.arange(topology.numberOfPoints * topology.levels, dtype=np.float64).reshape(
        (topology.numberOfPoints, topology.levels),
        order="F",
    )
    state = TimeStepState(
        step=2,
        time=2e-5,
        betaCells=values,
        betaVolume=np.zeros((topology.numberOfTriangles, topology.levels - 1)),
        phiAse=values + 1.0,
        dndtAse=values,
        dndtPump=values,
        aseResult=None,
    )

    callback = vtkWedge(tmp_path / "phi_{step:03d}", topology)
    path = callback(state)

    text = path.read_text(encoding="utf-8")
    assert path.name == "phi_002.vtk"
    assert "DATASET UNSTRUCTURED_GRID" in text
    assert f"POINTS {topology.numberOfPoints * topology.levels} float" in text
    assert f"CELLS {topology.numberOfPrisms} {topology.numberOfPrisms * 7}" in text
    assert "CELL_TYPES 4" in text
    assert "POINT_DATA 12" in text
    assert "phiAse" in _scalarNamesByDataKind(path)["POINT_DATA"]


def testVtkWedgeWritesMultiplePointFieldsFromState(tmp_path):
    topology = MeshTopology.fromGrid(Grid(xExtent=1, yExtent=1, zExtent=0.5, tileSizeZ=0.25))
    values = np.arange(topology.numberOfPoints * topology.levels, dtype=np.float64).reshape(
        (topology.numberOfPoints, topology.levels),
        order="F",
    )
    state = TimeStepState(
        step=3,
        time=3e-5,
        betaCells=values,
        betaVolume=np.zeros((topology.numberOfTriangles, topology.levels - 1)),
        phiAse=values + 1.0,
        dndtAse=values + 2.0,
        dndtPump=values + 3.0,
        aseResult=None,
    )

    callback = vtkWedge(tmp_path / "out_{field}_{step:03d}", topology, field=["phiAse", "dndtAse"])
    path = callback(state)

    text = path.read_text(encoding="utf-8")
    assert path.name == "out_phiAse_dndtAse_003.vtk"
    assert text.count(f"POINT_DATA {topology.numberOfPoints * topology.levels}") == 1
    assert "phiAse" in _scalarNamesByDataKind(path)["POINT_DATA"]
    assert "dndtAse" in _scalarNamesByDataKind(path)["POINT_DATA"]


def testVtkWedgeWritesAliasedFieldMappingFromState(tmp_path):
    topology = MeshTopology.fromGrid(Grid(xExtent=1, yExtent=1, zExtent=0.5, tileSizeZ=0.25))
    values = np.arange(topology.numberOfPoints * topology.levels, dtype=np.float64).reshape(
        (topology.numberOfPoints, topology.levels),
        order="F",
    )
    state = TimeStepState(
        step=4,
        time=4e-5,
        betaCells=values,
        betaVolume=np.zeros((topology.numberOfTriangles, topology.levels - 1)),
        phiAse=values + 1.0,
        dndtAse=values + 2.0,
        dndtPump=values + 3.0,
        aseResult=None,
    )

    callback = vtkWedge(tmp_path / "named_{step:03d}", topology, field={"phi": "phiAse", "dn": "dndtAse"})
    path = callback(state)

    text = path.read_text(encoding="utf-8")
    assert path.name == "named_004.vtk"
    assert "phi" in _scalarNamesByDataKind(path)["POINT_DATA"]
    assert "dn" in _scalarNamesByDataKind(path)["POINT_DATA"]


def testVtkWedgeWritesDirectlyFromStateTopology(tmp_path):
    topology = MeshTopology.fromGrid(Grid(xExtent=1, yExtent=1, zExtent=0.5, tileSizeZ=0.25))
    values = np.arange(topology.numberOfPoints * topology.levels, dtype=np.float64).reshape(
        (topology.numberOfPoints, topology.levels),
        order="F",
    )
    state = TimeStepState(
        step=5,
        time=5e-5,
        betaCells=values,
        betaVolume=np.zeros((topology.numberOfTriangles, topology.levels - 1)),
        phiAse=values + 1.0,
        dndtAse=values + 2.0,
        dndtPump=values + 3.0,
        aseResult=None,
        topology=topology,
    )

    path = vtkWedge(
        tmp_path / "state_{step:03d}",
        state,
        fields={
            "betaCells": state.betaCells,
            "phiASE": state.phiAse,
            "dndtAse": state.dndtAse,
            "cladAbs": state.phiAse * 5.5,
        },
    )

    text = path.read_text(encoding="utf-8")
    assert path.name == "state_005.vtk"
    assert "betaCells" in _scalarNamesByDataKind(path)["POINT_DATA"]
    assert "phiASE" in _scalarNamesByDataKind(path)["POINT_DATA"]
    assert "dndtAse" in _scalarNamesByDataKind(path)["POINT_DATA"]
    assert "cladAbs" in _scalarNamesByDataKind(path)["POINT_DATA"]


def testVtkWedgeWritesMixedPointAndCellFields(tmp_path):
    topology = MeshTopology.fromGrid(Grid(xExtent=1, yExtent=1, zExtent=0.5, tileSizeZ=0.25))
    point_values = np.arange(topology.numberOfPoints * topology.levels, dtype=np.float64).reshape(
        (topology.numberOfPoints, topology.levels),
        order="F",
    )
    cell_values = np.arange(topology.numberOfPrisms, dtype=np.float64).reshape(
        (topology.numberOfTriangles, topology.levels - 1),
        order="F",
    )

    path = vtkWedge(
        tmp_path / "mixed",
        geometry=topology,
        fields={"phiAse": point_values, "betaVolume": cell_values},
    )

    text = path.read_text(encoding="utf-8")
    assert f"POINT_DATA {topology.numberOfPoints * topology.levels}" in text
    assert "phiAse" in _scalarNamesByDataKind(path)["POINT_DATA"]
    assert f"CELL_DATA {topology.numberOfPrisms}" in text
    assert "betaVolume" in _scalarNamesByDataKind(path)["CELL_DATA"]
    assert text.index("POINT_DATA") < text.index("CELL_DATA")
