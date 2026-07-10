# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

import subprocess
import sys
from pathlib import Path

import numpy as np

from pyInclude.geometry.vtk import _parseVtk


repoRoot = Path(__file__).resolve().parents[3]
converter = repoRoot / "utils" / "convert_vtk_topology.py"


def _writeTwoWedgeVtk(path):
    path.write_text(
        "\n".join(
            [
                "# vtk DataFile Version 2.0",
                "two wedge prisms",
                "ASCII",
                "DATASET UNSTRUCTURED_GRID",
                "POINTS 8 double",
                "0 0 0",
                "1 0 0",
                "0 1 0",
                "1 1 0",
                "0 0 1",
                "1 0 1",
                "0 1 1",
                "1 1 1",
                "CELLS 2 14",
                "6 0 1 2 4 5 6",
                "6 1 3 2 5 7 6",
                "CELL_TYPES 2",
                "13",
                "13",
                "FIELD HASEonGPU 2",
                "cellMaterial 1 2 unsigned_int",
                "7 8",
                "globalValue 1 1 double",
                "42",
                "POINT_DATA 8",
                "SCALARS betaCells double 1",
                "LOOKUP_TABLE default",
                "0.0",
                "0.1",
                "0.2",
                "0.3",
                "0.4",
                "0.5",
                "0.6",
                "0.7",
                "CELL_DATA 2",
                "SCALARS betaVolume double 1",
                "LOOKUP_TABLE default",
                "10.0",
                "20.0",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _runConverter(direction, inputPath, outputPath):
    completed = subprocess.run(
        [
            sys.executable,
            str(converter),
            str(inputPath),
            str(outputPath),
            "--direction",
            direction,
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout


def testConvertWedgeVtkToTet4AndBackTransfersCellAdjacentFields(tmp_path):
    wedge_path = tmp_path / "wedge.vtk"
    tet4_path = tmp_path / "tet4.vtk"
    roundtrip_path = tmp_path / "roundtrip.vtk"
    _writeTwoWedgeVtk(wedge_path)

    _runConverter("wedge-to-tet4", wedge_path, tet4_path)

    points, cells, cell_types, point_data, cell_data, fields = _parseVtk(tet4_path)
    assert points.shape == (8, 3)
    assert len(cells) == 6
    assert np.all(np.asarray(cell_types) == 10)
    np.testing.assert_allclose(point_data["betaCells"], np.arange(8, dtype=np.float64) / 10.0)
    np.testing.assert_allclose(cell_data["betaVolume"], [10.0, 10.0, 10.0, 20.0, 20.0, 20.0])
    np.testing.assert_array_equal(fields["cellMaterial"], [7, 7, 7, 8, 8, 8])
    np.testing.assert_allclose(fields["globalValue"], 42.0)
    np.testing.assert_array_equal(fields["structuredNumberOfPoints"], [4])
    np.testing.assert_array_equal(fields["structuredNumberOfLevels"], [2])
    np.testing.assert_allclose(fields["structuredThickness"], [1.0])

    _runConverter("tet4-to-wedge", tet4_path, roundtrip_path)

    round_points, round_cells, round_types, round_point_data, round_cell_data, round_fields = _parseVtk(roundtrip_path)
    assert round_points.shape == (8, 3)
    assert len(round_cells) == 2
    assert np.all(np.asarray(round_types) == 13)
    np.testing.assert_array_equal(np.asarray(round_cells, dtype=np.uint32), np.asarray([[0, 1, 2, 4, 5, 6], [1, 3, 2, 5, 7, 6]], dtype=np.uint32))
    np.testing.assert_allclose(round_point_data["betaCells"], point_data["betaCells"])
    np.testing.assert_allclose(round_cell_data["betaVolume"], [10.0, 20.0])
    np.testing.assert_array_equal(round_fields["cellMaterial"], [7, 8])
    np.testing.assert_allclose(round_fields["globalValue"], 42.0)
