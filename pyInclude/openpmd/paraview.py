from __future__ import annotations

from pathlib import Path

import numpy as np


def _access(io, name):
    if hasattr(io, "Access_Type"):
        return getattr(io.Access_Type, name)
    return getattr(io.Access, name)


def _unit_dimensionless(io):
    return {}


def _reset_scalar_record(iteration, name, values, primitive_shape, axis_labels):
    import openpmd_api as io

    data = np.ascontiguousarray(np.asarray(values).reshape(-1))
    record = iteration.meshes[name]
    record.set_attribute("geometry", "other")
    record.set_attribute("geometryParameters", "topology=extruded_triangular_prism")
    record.set_attribute("dataOrder", "C")
    record.set_attribute("hasePrimitiveShape", list(primitive_shape))
    record.axis_labels = ["flatIndex"]
    record.grid_spacing = [1.0]
    record.grid_global_offset = [0.0]
    record.grid_unit_SI = 1.0
    record.unit_dimension = _unit_dimensionless(io)

    component = record[io.Mesh_Record_Component.SCALAR]
    component.unit_SI = 1.0
    component.position = [0.0]
    component.reset_dataset(io.Dataset(data.dtype, data.shape))
    component.store_chunk(data)
    record.set_attribute("haseAxes", list(axis_labels))


def writeParaviewState(
    state,
    outputDir,
    claddingAbsorption=1.0,
    pattern="laserPumpCladding_%06T.bp",
    handleName="laserPumpCladding.pmd",
):
    """Append one simulation state to an openPMD series and write a ParaView handle."""
    import openpmd_api as io

    output_dir = Path(outputDir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / handleName).write_text(pattern + "\n", encoding="utf-8")

    series_path = output_dir / pattern
    has_existing_series = any(output_dir.glob(pattern.replace("%06T", "*")))
    access = _access(io, "append") if has_existing_series else _access(io, "create")
    series = io.Series(str(series_path), access)
    series.set_software("HASEonGPU")

    iteration = series.iterations[int(state.step)]
    iteration.time = float(state.time)
    iteration.dt = 1.0
    iteration.time_unit_SI = 1.0

    topology = state.topology
    beta_cells = np.asarray(state.betaCells)
    primitive_shape = beta_cells.shape
    _reset_scalar_record(iteration, "beta_cells", beta_cells, primitive_shape, ["point", "level"])

    if state.phiAse is not None:
        phi_ase = np.asarray(state.phiAse)
        _reset_scalar_record(iteration, "phi_ase", phi_ase, phi_ase.shape, ["point", "level"])
        _reset_scalar_record(
            iteration,
            "cladding_absorption",
            phi_ase * np.float64(claddingAbsorption),
            phi_ase.shape,
            ["point", "level"],
        )
    if state.dndtAse is not None:
        dndt_ase = np.asarray(state.dndtAse)
        _reset_scalar_record(iteration, "dndt_ase", dndt_ase, dndt_ase.shape, ["point", "level"])
    if state.dndtPump is not None:
        dndt_pump = np.asarray(state.dndtPump)
        _reset_scalar_record(iteration, "dndt_pump", dndt_pump, dndt_pump.shape, ["point", "level"])
    if state.betaVolume is not None:
        beta_volume = np.asarray(state.betaVolume)
        _reset_scalar_record(iteration, "beta_volume", beta_volume, beta_volume.shape, ["cell", "layer"])

    if topology is not None:
        _reset_scalar_record(iteration, "points", np.asarray(topology.points), topology.points.shape, ["point", "coordinate"])
        _reset_scalar_record(
            iteration,
            "triangle_point_indices",
            np.asarray(topology.trianglePointIndices, dtype=np.uint32),
            topology.trianglePointIndices.shape,
            ["cell", "local_vertex"],
        )

    iteration.close()
    series.close()
    return output_dir / handleName
