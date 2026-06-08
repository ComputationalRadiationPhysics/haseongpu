# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""High-level Python simulation wrapper around pump, ASE, and time stepping."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import HASEonGPU_Bindings.HASEonGPU as HASEonGPU_Bindings
import numpy as np

from .geometry import GainMedium, _flat
from .laser import CrossSectionData, LaserProperties, PumpProperties, SpectralDecomposition
from . import mpiLauncher
from .pumping import BetaIntegrationGaussianSolver, Constants
from .timeIntegration import TimeDerivative, TimeIntegrationSolver


def _constructHostMesh(gainMedium):
    """Build the pybind host mesh from a ``GainMedium``.

    Arrays are flattened in Fortran order because the lower-level C++/legacy
    interfaces index beta first by transverse sample and then by z level.
    """
    topology = gainMedium.topology
    topology._require_levels()
    if topology.thickness is None:
        raise ValueError("topology thickness is required before running a simulation")

    derived = topology._topology()
    beta_volume = gainMedium.get("betaVolume").value
    beta_cells = gainMedium.get("betaCells").value
    cladding = gainMedium.get("claddingCellTypes").value
    refractive = gainMedium.get("refractiveIndices").value
    reflective = gainMedium.get("reflectivities").value

    return HASEonGPU_Bindings.HostMesh(
        _flat(topology.trianglePointIndices, 3, np.uint32, "trianglePointIndices"),
        topology.numberOfTriangles,
        int(topology.levels),
        topology.numberOfPoints,
        float(topology.thickness),
        _flat(topology.points, 2, np.float64, "points"),
        derived["triangleCenterX"],
        derived["triangleCenterY"],
        derived["triangleNormalPoint"],
        derived["triangleNormalsX"],
        derived["triangleNormalsY"],
        derived["forbiddenEdge"],
        derived["triangleNeighbors"],
        derived["triangleSurfaces"],
        _flat(beta_volume, topology.levels - 1, np.float64, "betaVolume"),
        _flat(beta_cells, topology.levels, np.float64, "betaCells"),
        _flat(cladding, None, np.uint32, "claddingCellTypes"),
        _flat(refractive, None, np.float32, "refractiveIndices"),
        _flat(reflective, None, np.float32, "reflectivities"),
        float(gainMedium.get("nTot").value),
        float(gainMedium.get("crystalTFluo").value),
        int(gainMedium.get("claddingNumber").value),
        float(gainMedium.get("claddingAbsorption").value),
    )


@dataclass
class PhiASE:
    """Configure and run the ASE flux calculation for one gain-medium state.

    ``Simulation`` normally owns this object and calls ``run(...)`` during each
    time-step derivative evaluation. Advanced users can also call ``run``
    directly with a ``GainMedium`` and ``CrossSectionData`` object.
    """

    config: object | None = None
    """Optional YAML filename or mapping with PhiASE run-control settings."""
    crossSections: CrossSectionData | None = None
    """Absorption/emission spectra used by the ASE calculation."""
    spectralProperties: SpectralDecomposition | None = None
    """Alias for ``crossSections`` kept for the public spectral API."""
    laserProperties: LaserProperties | None = None
    """Lower-level laser property store accepted by legacy workflows."""
    gainMedium: GainMedium | None = None
    """Optional medium stored for direct ``run()`` calls."""

    minRaysPerSample: int = 100000
    """Minimum Monte Carlo rays launched for each beta sample."""
    maxRaysPerSample: int = 100000
    """Maximum rays per sample allowed during adaptive refinement."""
    mseThreshold: float = 0.1
    """Target mean-squared-error threshold for adaptive ASE sampling."""
    repetitions: int = 4
    """Maximum repeated ASE estimates at a fixed ray count."""
    adaptiveSteps: int = 4
    """Number of ray-count increases between min and max rays."""
    useReflections: bool = True
    """Whether top/bottom surface reflectivities affect ray propagation."""
    monochromatic: bool = False
    """Use only the first spectral samples instead of wavelength integration."""

    backend: str = None
    """Alpaka backend name; inspect valid strings with ``AlpakaBackends.all()``."""
    parallelMode: str = "single"
    """Execution mode: direct binding ``single`` or MPI launcher ``mpi``."""
    numDevices: int = 1
    """Maximum compute devices made available to the lower-level run."""
    nPerNode: int = 1
    """MPI launcher ranks per node when ``parallelMode`` is ``mpi``."""
    writeVtk: bool = False
    """Request VTK output from lower-level compute paths when supported."""
    devices: list[int] = field(default_factory=list)
    """Optional explicit device ids passed to the lower-level compute path."""
    minSampleRange: int | None = None
    """Inclusive first flattened beta sample processed by ASE."""
    maxSampleRange: int | None = None
    """Inclusive last flattened beta sample processed by ASE."""
    rngSeed: int | None = None
    """Optional RNG seed for reproducible Monte Carlo sampling."""

    _experiment: object | None = field(default=None, init=False, repr=False)
    _compute: object | None = field(default=None, init=False, repr=False)
    _hostMesh: object | None = field(default=None, init=False, repr=False)
    _result: object | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        if isinstance(self.config, (str, Path)):
            self._applyConfig(self._loadConfig(self.config))
        elif isinstance(self.config, dict):
            self._applyConfig(self.config)

        self._syncCrossSections()

    def _syncCrossSections(self):
        if self.crossSections is None and self.spectralProperties is not None:
            self.crossSections = self.spectralProperties
        if self.spectralProperties is None and self.crossSections is not None:
            self.spectralProperties = self.crossSections
        if self.crossSections is None and self.laserProperties is not None:
            laser = self.laserProperties.toDict()
            self.crossSections = CrossSectionData(
                wavelengthsAbsorption=laser["l_abs"],
                crossSectionAbsorption=laser["s_abs"],
                wavelengthsEmission=laser["l_ems"],
                crossSectionEmission=laser["s_ems"],
                resolution=laser["l_res"],
            )
            self.spectralProperties = self.crossSections
        return self

    @classmethod
    def fromYaml(cls, filename, **overrides):
        """Create a ``PhiASE`` configuration from YAML plus Python overrides."""
        obj = cls(filename)
        for name, value in overrides.items():
            setattr(obj, name, value)
        return obj._syncCrossSections()

    @staticmethod
    def addArguments(parser):
        """Add command-line arguments that map to ``PhiASE`` settings."""
        parser.add_argument("--phi-ase-config", default=None, help="YAML file with PhiASE compute/experiment settings")
        parser.add_argument("--min-rays-per-sample", type=int, default=None)
        parser.add_argument("--max-rays-per-sample", type=int, default=None)
        parser.add_argument("--mse-threshold", type=float, default=None)
        parser.add_argument("--repetitions", type=int, default=None)
        parser.add_argument("--adaptive-steps", type=int, default=None)
        parser.add_argument("--backend", default=None)
        parser.add_argument("--parallel-mode", default=None)
        parser.add_argument("--max-gpus", type=int, default=None)
        parser.add_argument("--n-per-node", type=int, default=None)
        parser.add_argument("--rng-seed", type=int, default=None)
        return parser

    @classmethod
    def fromArgs(cls, args, **overrides):
        """Create a ``PhiASE`` configuration from parsed argparse results."""
        config = getattr(args, "phi_ase_config", None)
        obj = cls(config) if config else cls()
        mapping = {
            "min_rays_per_sample": "minRaysPerSample",
            "max_rays_per_sample": "maxRaysPerSample",
            "mse_threshold": "mseThreshold",
            "repetitions": "repetitions",
            "adaptive_steps": "adaptiveSteps",
            "backend": "backend",
            "parallel_mode": "parallelMode",
            "max_gpus": "numDevices",
            "n_per_node": "nPerNode",
            "rng_seed": "rngSeed",
        }
        for arg_name, attr_name in mapping.items():
            value = getattr(args, arg_name, None)
            if value is not None:
                setattr(obj, attr_name, value)
        for name, value in overrides.items():
            setattr(obj, name, value)
        return obj

    @staticmethod
    def _loadConfig(filename):
        path = Path(filename)
        try:
            import yaml
        except ImportError as exc:
            raise ImportError("PhiASE YAML configuration requires PyYAML") from exc
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
        if not isinstance(data, dict):
            raise ValueError(f"PhiASE config '{filename}' must contain a mapping")
        return data

    def _applyConfig(self, config):
        sections = []
        for key in ("phiASE", "phi_ase", "experiment", "compute"):
            value = config.get(key)
            if isinstance(value, dict):
                sections.append(value)
        sections.append(config)
        aliases = {
            "minRays": "minRaysPerSample",
            "maxRays": "maxRaysPerSample",
            "min_rays_per_sample": "minRaysPerSample",
            "max_rays_per_sample": "maxRaysPerSample",
            "mse_threshold": "mseThreshold",
            "adaptive_steps": "adaptiveSteps",
            "use_reflections": "useReflections",
            "parallel_mode": "parallelMode",
            "max_gpus": "numDevices",
            "n_per_node": "nPerNode",
            "write_vtk": "writeVtk",
            "min_sample_range": "minSampleRange",
            "max_sample_range": "maxSampleRange",
            "rng_seed": "rngSeed",
        }
        allowed = {
            "minRaysPerSample", "maxRaysPerSample", "mseThreshold", "repetitions",
            "adaptiveSteps", "useReflections", "monochromatic", "backend", "parallelMode",
            "numDevices", "nPerNode", "writeVtk", "devices", "minSampleRange", "maxSampleRange", "rngSeed",
        }
        for section in sections:
            for name, value in section.items():
                attr = aliases.get(name, name)
                if attr in allowed:
                    setattr(self, attr, value)
        return self

    def run(self, gainMedium=None, crossSections=None):
        """Run ASE for the supplied or configured ``GainMedium``.

        Returns ``self``. Use ``getResults()`` afterwards to access the raw
        lower-level result, including ``phiAse``.
        """
        medium = gainMedium if gainMedium is not None else self.gainMedium
        if medium is None:
            raise ValueError("PhiASE.run requires gainMedium; pass it through Simulation or run(gainMedium=...)")
        cross_sections = crossSections if crossSections is not None else self.crossSections
        if cross_sections is None and self.laserProperties is not None:
            cross_sections = self.crossSections
        if cross_sections is None:
            raise ValueError("PhiASE.run requires crossSections")

        self._hostMesh = _constructHostMesh(medium)
        laser_properties = LaserProperties(crossSections=cross_sections)
        laser = laser_properties.toDict()

        self._experiment = HASEonGPU_Bindings.ExperimentParameters(
            minRaysPerSample=int(self.minRaysPerSample),
            maxRaysPerSample=int(self.maxRaysPerSample),
            lambdaA=laser["l_abs"],
            lambdaE=laser["l_ems"],
            sigmaA=laser["s_abs"],
            sigmaE=laser["s_ems"],
            maxSigmaA=laser_properties.maxSigmaA,
            maxSigmaE=laser_properties.maxSigmaE,
            mseThreshold=float(self.mseThreshold),
            useReflections=bool(self.useReflections),
            spectral=int(laser["l_res"]),
            monochromatic=bool(self.monochromatic),
        )

        topology = medium.topology
        min_sample = 0 if self.minSampleRange is None else int(self.minSampleRange)
        max_sample = (
            topology.numberOfPoints * topology.levels - 1
            if self.maxSampleRange is None
            else int(self.maxSampleRange)
        )

        self._compute = HASEonGPU_Bindings.ComputeParameters(
            maxRepetitions=int(self.repetitions),
            adaptiveSteps=int(self.adaptiveSteps),
            numDevices=int(self.numDevices),
            backend=str(self.backend),
            parallelMode=str(self.parallelMode),
            writeVtk=bool(self.writeVtk),
            devices=list(self.devices),
            minSampleRange=min_sample,
            maxSampleRange=max_sample,
        )
        if self.rngSeed is not None:
            self._compute.rngSeed = int(self.rngSeed)

        if str(self.parallelMode).lower() == "mpi":
            self._result = mpiLauncher.runPhiaseMPI(self, medium, laser, laser_properties)
            return self

        self._result = HASEonGPU_Bindings.calcPhiASE(
            self._experiment, self._compute, self._hostMesh
        )
        return self

    def getResults(self):
        """Return the raw ASE result from the most recent ``run(...)`` call."""
        if self._result is None:
            raise RuntimeError("simulation has not been run yet")
        return self._result

    @property
    def experimentParameters(self):
        """Low-level experiment parameters built for the last ``run(...)``."""
        return self._experiment

    @property
    def computeParameters(self):
        """Low-level compute parameters built for the last ``run(...)``."""
        return self._compute

    @property
    def hostMesh(self):
        """Low-level host mesh built from the last ``GainMedium`` input."""
        return self._hostMesh


class LegacyGridDataBetaVolumeMapper:
    """Map point-centered ``betaCells`` to prism-centered ``betaVolume``."""

    def __init__(self, method="linear"):
        """Create a mapper using a scipy ``griddata`` interpolation method.

        The mapper samples point-level beta values at prism centers so ASE can
        use one representative excited-state fraction per wedge volume.
        """
        self.method = method

    def map(self, medium):
        """Return prism-centered ``betaVolume`` for the supplied medium."""
        try:
            from scipy.interpolate import griddata
        except ImportError as exc:
            raise ImportError("LegacyGridDataBetaVolumeMapper requires scipy") from exc

        topology = medium.topology
        derived = topology._topology()
        beta_cells = np.asarray(medium.get("betaCells").value, dtype=np.float64).reshape(
            (topology.numberOfPoints, topology.levels),
            order="F",
        )

        x_1 = topology.points[:, 0]
        y_1 = topology.points[:, 1]
        x_2 = topology.points[:, 0]
        y_2 = topology.points[:, 1]
        x = np.concatenate((x_1, x_2))
        y = np.concatenate((y_1, y_2))
        z_1 = np.zeros(x_1.shape[0])
        z_2 = np.zeros(x_2.shape[0]) + topology.thickness
        z = np.concatenate((z_1, z_2))

        xi = np.asarray(derived["triangleCenterX"], dtype=np.float64)
        yi = np.asarray(derived["triangleCenterY"], dtype=np.float64)
        zi = np.zeros(xi.shape) + topology.thickness / 2.0

        beta_volume = np.zeros((topology.numberOfTriangles, topology.levels - 1), dtype=np.float64)
        for level in range(topology.levels - 1):
            values = np.concatenate((beta_cells[:, level], beta_cells[:, level + 1]))
            beta_volume[:, level] = griddata((x, y, z), values, (xi, yi, zi), method=self.method)
            z = z + topology.thickness
            zi = zi + topology.thickness
        return beta_volume


@dataclass
class TimeStepState:
    """Snapshot handed to ``onStep`` callbacks after a completed time step.

    The arrays are copies of the simulation outputs at ``step``/``time``.
    ``betaCells`` and ``phiAse`` are point-by-level arrays with shape
    ``(numberOfPoints, numberOfLevels)``. ``betaVolume`` is prism-centered
    with shape ``(numberOfTriangles, numberOfLevels - 1)``.
    """

    step: int
    """Completed one-based step index."""
    time: float
    """Physical simulation time after the step, in seconds."""
    betaCells: np.ndarray
    """Excited-state fraction at mesh points and z-levels."""
    betaVolume: np.ndarray
    """Excited-state fraction interpolated to wedge-prism centers."""
    phiAse: np.ndarray | None
    """ASE flux at mesh points and z-levels, or ``None`` if unavailable."""
    dndtAse: np.ndarray
    """ASE depletion contribution to ``d beta / dt``."""
    dndtPump: np.ndarray
    """Pump contribution to ``d beta / dt``."""
    aseResult: object | None
    """Raw lower-level ASE result object for advanced inspection."""


@dataclass
class Simulation:
    """High-level pump, ASE, fluorescence, and time-integration loop.

    Register ``onInit`` or ``beforeStep`` callbacks when a function needs the
    live ``Simulation`` object and may inspect or mutate it. Register
    ``onStep`` callbacks when a function needs the completed ``TimeStepState``
    snapshot produced by each step.
    """

    gainMedium: GainMedium
    """Geometry, material data, and current beta arrays."""
    pump: PumpProperties
    """Pump model that adds population to ``betaCells``."""
    phiASE: PhiASE
    """ASE configuration and execution handle."""
    timeIntegrationSolver: TimeIntegrationSolver
    """Solver object that advances ``betaCells`` using ``d beta / dt``."""
    timeStep: float
    """Physical time increment per step, in seconds."""
    crossSections: CrossSectionData | None = None
    """Shared spectra used by pump and ASE when not set on either object."""
    endTime: float | None = None
    """Optional target physical time for ``runUntil()``."""
    constants: Constants = field(default_factory=Constants)
    """Physical constants used by the pump integrator."""
    updateTerminalLevel: bool = True
    """Whether the last z-level beta values are advanced by the solver."""

    _time: float = field(default=0.0, init=False, repr=False)
    _step: int = field(default=0, init=False, repr=False)
    _initialized: bool = field(default=False, init=False, repr=False)
    _initCallbacks: list = field(default_factory=list, init=False, repr=False)
    _beforeStepCallbacks: list = field(default_factory=list, init=False, repr=False)
    _callbacks: list = field(default_factory=list, init=False, repr=False)
    _states: list[TimeStepState] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self):
        if self.timeIntegrationSolver is None or not hasattr(self.timeIntegrationSolver, "step"):
            raise ValueError("Simulation requires a timeIntegrationSolver with a step(rhs, betaCells, time, timeStep) method")
        if self.timeStep <= 0.0:
            raise ValueError("timeStep must be positive")
        if self.crossSections is None:
            self.crossSections = self._resolveSpectralProperties()
        if self.phiASE.crossSections is None:
            self.phiASE.crossSections = self.crossSections
        if self.phiASE.spectralProperties is None:
            self.phiASE.spectralProperties = self.crossSections
        self._ensureStateArrays()

    def _resolveSpectralProperties(self):
        if self.phiASE.spectralProperties is not None:
            return self.phiASE.spectralProperties
        if self.phiASE.crossSections is not None:
            return self.phiASE.crossSections
        if self.pump.spectralProperties is not None:
            return self.pump.spectralProperties
        if self.pump.crossSections is not None:
            return self.pump.crossSections
        raise ValueError("Simulation requires spectral properties via Simulation.crossSections, phiASE, or pump")

    def onStep(self, callback):
        """Register ``callback(state)`` to run after every completed step.

        ``state`` is a ``TimeStepState`` containing copied result arrays such as
        ``betaCells``, ``betaVolume``, ``phiAse``, ``dndtPump``, and
        ``dndtAse``. Callback return values are ignored.
        """
        self._callbacks.append(callback)
        return self

    def onInit(self, callback):
        """Register ``callback(simulation)`` to run once before the first step.

        The callback receives this live ``Simulation`` object, so it can read
        or change ``gainMedium``, ``pump``, ``phiASE``, ``timeStep``, and other
        simulation settings before any derivative is evaluated.
        """
        self._initCallbacks.append(callback)
        return self

    def beforeStep(self, callback):
        """Register ``callback(simulation)`` to run before every step.

        The callback receives this live ``Simulation`` object at the current
        ``time`` and ``stepIndex``. Use this for controlled changes to inputs
        before the next time-step update.
        """
        self._beforeStepCallbacks.append(callback)
        return self

    def runUntil(self, endtime=None, endTime=None):
        """Advance steps until the configured or supplied end time is reached."""
        target = self.endTime if endtime is None and endTime is None else (endtime if endtime is not None else endTime)
        if target is None:
            raise ValueError("runUntil requires endtime or an endTime configured on construction")
        while self._time < float(target) - 0.5 * self.timeStep:
            self.step()
        return self

    def runSteps(self, steps):
        """Run exactly ``steps`` calls to ``step()`` and return ``self``."""
        for _ in range(int(steps)):
            self.step()
        return self

    def step(self):
        """Advance one time step and return the completed ``TimeStepState``."""
        self._runInitCallbacks()
        for callback in self._beforeStepCallbacks:
            callback(self)

        beta_cells = np.asarray(self.gainMedium.get("betaCells").value, dtype=np.float64).reshape(
            self.gainMedium.get("betaCells").expectedShape,
            order="F",
        )

        integration_result = self.timeIntegrationSolver.step(
            self._timeDerivative,
            beta_cells.copy(),
            self._time,
            self.timeStep,
        )
        updated_beta = integration_result.betaCells
        derivative = integration_result.evaluation

        beta_next = updated_beta if self.updateTerminalLevel else beta_cells.copy()
        if not self.updateTerminalLevel:
            beta_next[:, :-1] = updated_beta[:, :-1]
        self.gainMedium.get("betaCells").value = np.clip(beta_next, 0.0, 1.0)
        self._updateBetaVolumeFromCells()

        self._step += 1
        self._time += self.timeStep

        state = TimeStepState(
            step=self._step,
            time=self._time,
            betaCells=np.asarray(self.gainMedium.get("betaCells").value, dtype=np.float64).reshape(
                beta_cells.shape,
                order="F",
            ).copy(),
            betaVolume=np.asarray(self.gainMedium.get("betaVolume").value, dtype=np.float64).reshape(
                self.gainMedium.get("betaVolume").expectedShape,
                order="F",
            ).copy(),
            phiAse=None if derivative.phiAse is None else derivative.phiAse.copy(),
            dndtAse=derivative.dndtAse.copy(),
            dndtPump=derivative.dndtPump.copy(),
            aseResult=derivative.aseResult,
        )
        self._states.append(state)
        for callback in self._callbacks:
            callback(state)
        return state

    def getResults(self):
        """Return all stored ``TimeStepState`` snapshots in step order."""
        return self._states

    @property
    def time(self):
        """Current physical simulation time in seconds."""
        return self._time

    @property
    def stepIndex(self):
        """Number of completed time steps."""
        return self._step

    def _ensureStateArrays(self):
        topology = self.gainMedium.topology
        if "betaCells" not in self.gainMedium.physical:
            self.gainMedium.get("betaCells").value = np.zeros((topology.numberOfPoints, topology.levels))
        if "betaVolume" not in self.gainMedium.physical:
            self._updateBetaVolumeFromCells()

    def _runInitCallbacks(self):
        if self._initialized:
            return
        self._initialized = True
        for callback in self._initCallbacks:
            callback(self)

    def _dndtPump(self, beta_cells):
        pump_duration = self.pump.activeDuration(self.timeStep)
        solver = self.pump.solver if self.pump.solver is not None else BetaIntegrationGaussianSolver()
        updated = solver.step(
            {
                "betaCell": beta_cells,
                "_medium": self.gainMedium,
                "_timeStep": self.timeStep,
                "_constants": self.constants,
                "_substeps": None,
            },
            self.pump,
        )
        return (updated - beta_cells) / pump_duration

    def _timeDerivative(self, beta_cells, time):
        beta_cells = np.asarray(beta_cells, dtype=np.float64).reshape(
            self.gainMedium.get("betaCells").expectedShape,
            order="F",
        )
        self.gainMedium.get("betaCells").value = beta_cells
        self._updateBetaVolumeFromCells()
        self.phiASE.run(gainMedium=self.gainMedium, crossSections=self.crossSections)
        ase_result = self.phiASE.getResults()
        phi_ase = np.asarray(ase_result.phiAse, dtype=np.float64).reshape(beta_cells.shape, order="F")
        dndt_ase = self._aseDerivative(phi_ase, betaCells=beta_cells)
        dndt_pump = self._dndtPump(beta_cells)
        tau = max(float(self.gainMedium.get("crystalTFluo").value), np.finfo(float).tiny)
        total_derivative = dndt_pump - dndt_ase - beta_cells / tau
        return TimeDerivative(
            betaCells=beta_cells.copy(),
            dndtPump=dndt_pump.copy(),
            dndtAse=dndt_ase.copy(),
            derivative=total_derivative,
            tau=tau,
            phiAse=phi_ase.copy(),
            aseResult=ase_result,
        )

    def _aseDerivative(self, phi_ase, betaCells=None):
        laser = self.phiASE.laserProperties
        laser = LaserProperties(crossSections=self.crossSections) if laser is None else laser
        sigma_e = laser.maxSigmaE
        sigma_a = laser.absorptionAtEmissionPeak
        if betaCells is None:
            betaCells = np.asarray(self.gainMedium.get("betaCells").value, dtype=np.float64).reshape(phi_ase.shape, order="F")
        gain_per_density = betaCells * (sigma_e + sigma_a) - sigma_a
        active_mask = self._activePointMask()[:, None]
        return np.where(active_mask, gain_per_density * phi_ase, 0.0)

    def _activePointMask(self):
        topology = self.gainMedium.topology
        cladding_types = np.asarray(self.gainMedium.get("claddingCellTypes").value, dtype=np.uint32).reshape(-1)
        cladding_number = int(self.gainMedium.get("claddingNumber").value)
        active = np.zeros(topology.numberOfPoints, dtype=bool)
        for tri_id, tri in enumerate(topology.trianglePointIndices):
            if tri_id >= len(cladding_types) or cladding_types[tri_id] != cladding_number:
                active[np.asarray(tri, dtype=np.uint32)] = True
        return active

    def _updateBetaVolumeFromCells(self):
        beta_volume = LegacyGridDataBetaVolumeMapper().map(self.gainMedium)
        self.gainMedium.get("betaVolume").value = beta_volume


TimeSteppedSimulation = Simulation
