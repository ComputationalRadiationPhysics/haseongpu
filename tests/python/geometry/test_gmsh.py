# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later


import numpy as np
import pytest

from HASEonGPU import GainMedium, Gmsh, MeshTopology

try:
    import gmsh
except (ImportError, OSError) as exc:
    pytest.skip(f"gmsh is not importable: {exc}", allow_module_level=True)


def _withGmshModel(name, path, dim, build):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(name)
        build()
        gmsh.model.mesh.generate(dim)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


def _write2dTriangleMesh(path):
    def build():
        core, cladding = _cylindricalCoreCladdingSurfaces(0.5, 0.8, meshSize=0.4)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [core], 20)
        gmsh.model.setPhysicalName(2, 20, "Core")
        gmsh.model.addPhysicalGroup(2, [cladding], 21)
        gmsh.model.setPhysicalName(2, 21, "CladdingShell")

    _withGmshModel("triangle_2d", path, 2, build)


def _cylindricalCoreCladdingSurfaces(coreRadius, claddingRadius, *, meshSize):
    geo = gmsh.model.geo
    center = geo.addPoint(0.0, 0.0, 0.0, meshSize)
    rings = []
    for radius in (coreRadius, claddingRadius):
        rings.append(
            [
                geo.addPoint(radius, 0.0, 0.0, meshSize),
                geo.addPoint(0.0, radius, 0.0, meshSize),
                geo.addPoint(-radius, 0.0, 0.0, meshSize),
                geo.addPoint(0.0, -radius, 0.0, meshSize),
            ]
        )
    corePoints, claddingPoints = rings
    coreArcs = [geo.addCircleArc(corePoints[i], center, corePoints[(i + 1) % 4]) for i in range(4)]
    claddingArcs = [geo.addCircleArc(claddingPoints[i], center, claddingPoints[(i + 1) % 4]) for i in range(4)]
    core = geo.addPlaneSurface([geo.addCurveLoop(coreArcs)])
    cladding = geo.addPlaneSurface([geo.addCurveLoop(claddingArcs), geo.addCurveLoop([-arc for arc in coreArcs])])
    return core, cladding


def test_gmshBuildsTopology(tmp_path):
    msh = tmp_path / "surface.msh"
    _write2dTriangleMesh(msh)

    with pytest.raises(ValueError, match="numberOfLevels and thickness"):
        MeshTopology.fromFile(msh, format="gmsh")

    topology = MeshTopology.fromFile(msh, format="gmsh", numberOfLevels=3, thickness=0.25)

    assert topology.numberOfPoints >= 3
    assert topology.numberOfTriangles >= 1
    assert topology.numberOfPrisms == topology.numberOfTriangles * 2
    assert topology.thickness == 0.25
    assert topology.metadata["format"] == "gmsh"
    assert topology.metadata["dimension"] == 2
    assert isinstance(topology.metadata["gmsh"], Gmsh)


def test_gmshMapsCladdingSurface(tmp_path):
    msh = tmp_path / "core_cladding.msh"
    _write2dTriangleMesh(msh)

    topology = MeshTopology.fromFile(msh, format="gmsh", numberOfLevels=6, thickness=0.25)
    medium = GainMedium(topology=topology)

    assert topology.numberOfPoints > 4
    assert topology.numberOfTriangles > 2
    assert topology.numberOfPrisms == topology.numberOfTriangles * 5
    assert topology.levels == 6
    assert topology.thickness == pytest.approx(0.25)
    cladding = medium.get("claddingCellTypes").value
    assert cladding.dtype == np.uint32
    assert np.count_nonzero(cladding == 21) > 0
    assert np.count_nonzero(cladding == 0) > 0
