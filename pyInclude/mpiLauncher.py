# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""MPI launcher support for ``PhiASE(parallelMode="mpi")``."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import signal
import subprocess
import tempfile
import sys

import numpy as np

from .laser import CrossSectionData
from .openpmd import transport
from .geometry import _flat


@dataclass
class PhiAseMpiResult:
    """Result object shaped like the direct pybind ASE result."""

    phiAse: np.ndarray
    r"""Flattened ASE flux :math:`\Phi_i` for each beta sample."""
    mse: np.ndarray
    """Monte Carlo mean-squared-error estimate per sample."""
    totalRays: np.ndarray
    """Number of rays used per sample."""
    droppedRays: np.ndarray
    """Number of non-finite ray contributions dropped per accepted sample."""
    dndtAse: np.ndarray
    """ASE contribution to ``d beta / dt`` reconstructed in Python."""


_PROBED_MPI_IO_ROOTS: set[tuple[Path, int]] = set()


def _runInterruptible(cmd: list[str]) -> subprocess.CompletedProcess:
    process = subprocess.Popen(cmd, start_new_session=True)
    try:
        return subprocess.CompletedProcess(cmd, process.wait())
    except KeyboardInterrupt:
        try:
            os.killpg(process.pid, signal.SIGTERM)
            process.wait(timeout=5)
        except ProcessLookupError:
            pass
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        raise


def _repoRoot() -> Path:
    return Path(__file__).resolve().parents[1]


def _mpiIoRoot() -> Path:
    root = _repoRoot() / "IO" / "phiase_mpi"
    root.mkdir(parents=True, exist_ok=True)
    return root


def _ensureMpiIoVisible(ioRoot: Path, nPerNode: int) -> None:
    key = (ioRoot.resolve(), int(nPerNode))
    if key in _PROBED_MPI_IO_ROOTS:
        return

    with tempfile.TemporaryDirectory(prefix="probe_", dir=ioRoot) as probe_dir:
        probe_path = Path(probe_dir)
        marker = probe_path / "rank_visibility_marker.txt"
        marker.write_text("hase mpi io visibility probe\n", encoding="utf-8")

        code = """
from pathlib import Path
import os
import sys

probe_dir = Path(sys.argv[1])
marker = probe_dir / "rank_visibility_marker.txt"
if not marker.is_file():
    raise SystemExit(f"MPI rank cannot read shared IO marker: {marker}")
rank = (
    os.environ.get("OMPI_COMM_WORLD_RANK")
    or os.environ.get("PMI_RANK")
    or os.environ.get("PMIX_RANK")
    or os.environ.get("SLURM_PROCID")
    or "0"
)
(probe_dir / f"rank_{rank}.txt").write_text("ok\\n", encoding="utf-8")
"""
        cmd = [
            "mpiexec",
            "-npernode",
            str(int(nPerNode)),
            sys.executable,
            "-c",
            code,
            str(probe_path),
        ]
        status = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if status.returncode != 0:
            details = "\n".join(part for part in (status.stdout, status.stderr) if part).strip()
            raise RuntimeError(
                "PhiASE MPI IO preflight failed. The repository IO directory must be visible "
                f"to every MPI rank: {ioRoot}. "
                f"Command: {' '.join(cmd)}"
                + (f"\n{details}" if details else "")
            )

    _PROBED_MPI_IO_ROOTS.add(key)


def _packedFromModernInputs(gainMedium, laser):
    topology = gainMedium.topology
    topology._require_levels()
    if topology.thickness is None:
        raise ValueError("topology thickness is required before running a simulation")

    derived = topology._topology()
    return {
        "layout": "flattened",
        "container": "ndarray",
        "numberOfPoints": int(topology.numberOfPoints),
        "numberOfTriangles": int(topology.numberOfTriangles),
        "numberOfLevels": int(topology.levels),
        "points_flat": _flat(topology.points, 2, np.float64, "points"),
        "trianglePointIndices_flat": _flat(topology.trianglePointIndices, 3, np.uint32, "trianglePointIndices"),
        "betaCells_flat": _flat(gainMedium.get("betaCells").value, topology.levels, np.float64, "betaCells"),
        "betaVolume_flat": _flat(gainMedium.get("betaVolume").value, topology.levels - 1, np.float64, "betaVolume"),
        "claddingCellTypes_flat": _flat(gainMedium.get("claddingCellTypes").value, None, np.uint32, "claddingCellTypes"),
        "refractiveIndices_flat": _flat(gainMedium.get("refractiveIndices").value, None, np.float32, "refractiveIndices"),
        "reflectivities_flat": _flat(gainMedium.get("reflectivities").value, None, np.float32, "reflectivities"),
        "triangleNormalsX_flat": derived["triangleNormalsX"],
        "triangleNormalsY_flat": derived["triangleNormalsY"],
        "triangleNeighbors_flat": derived["triangleNeighbors"],
        "triangleSurfaces_flat": derived["triangleSurfaces"],
        "triangleCenterX_flat": derived["triangleCenterX"],
        "triangleCenterY_flat": derived["triangleCenterY"],
        "triangleNormalPoint_flat": derived["triangleNormalPoint"],
        "forbiddenEdge_flat": derived["forbiddenEdge"],
        "laserParameter": {
            "l_abs": np.asarray(laser["l_abs"], dtype=np.float64),
            "l_ems": np.asarray(laser["l_ems"], dtype=np.float64),
            "s_abs": np.asarray(laser["s_abs"], dtype=np.float64),
            "s_ems": np.asarray(laser["s_ems"], dtype=np.float64),
            "l_res": int(laser["l_res"]),
        },
    }


def runPhiaseMPI(phiAse, gainMedium, laser, laserProperties):
    """Run ``calcPhiASE`` through ``mpiexec`` and return arrays to Python.

    The launcher writes an openPMD input series below ``IO/phiase_mpi`` so every
    MPI rank can see the same directory, then reads the openPMD result series
    back into a ``PhiAseMpiResult``.
    """
    io_root = _mpiIoRoot()
    _ensureMpiIoVisible(io_root, int(phiAse.nPerNode))

    cross_sections = CrossSectionData(
        wavelengthsAbsorption=laser["l_abs"],
        crossSectionAbsorption=laser["s_abs"],
        wavelengthsEmission=laser["l_ems"],
        crossSectionEmission=laser["s_ems"],
        resolution=laser["l_res"],
    )
    result = transport.runPhiASE(
        phiAse,
        gainMedium,
        cross_sections,
        command_prefix=["mpiexec", "-npernode", str(int(phiAse.nPerNode))],
        workspace_dir=io_root,
    )

    return PhiAseMpiResult(
        phiAse=np.asarray(result.phiAse, dtype=np.float64).reshape(-1),
        mse=np.asarray(result.mse, dtype=np.float64).reshape(-1),
        totalRays=np.asarray(result.totalRays, dtype=np.uint32).reshape(-1),
        droppedRays=np.asarray(getattr(result, "droppedRays", np.zeros_like(result.totalRays)), dtype=np.uint32).reshape(-1),
        dndtAse=np.asarray(result.dndtAse, dtype=np.float64).reshape(-1),
    )
