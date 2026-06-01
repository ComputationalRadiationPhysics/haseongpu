# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np

from HASEonGPU import Grid, MeshTopology, TimeStepState, vtkWedge


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
    assert "SCALARS phiAse float 1" in text
