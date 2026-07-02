# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""High-level Python simulation wrapper around pump, ASE, and time stepping."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .alpakaUtils import AlpakaBackends
from .geometry import GainMedium
from .laser import CrossSectionData, LaserProperties, PumpProperties, SpectralDecomposition
from .openpmd import transport
from .timeIntegration import TimeIntegrationSolver


HASE_CONFIGURE_HINT = "Run `hase-configure` to generate a matching backend/openPMD setup."


def _preferredDefaultBackend():
    try:
        from .alpakaUtils import AlpakaBackends

        backends = AlpakaBackends.all()
    except Exception as exc:
        raise RuntimeError(
            "PhiASE.backend is not set and HASEonGPU could not query installed Alpaka "
            f"backends. {HASE_CONFIGURE_HINT}"
        ) from exc
    if not backends:
        raise RuntimeError(f"PhiASE.backend is not set and no Alpaka backend is available. {HASE_CONFIGURE_HINT}")
    for marker in ("Host_Cpu_CpuSerial", "CpuSerial"):
        for backend in backends:
            if marker in backend:
                return backend
    return backends[0]


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
    openpmdBackend: str | None = "auto"
    """openPMD backend; ``auto`` prefers SST, ADIOS, then HDF5 when supported."""
    parallelMode: str = "single"
    """Execution mode: local ``single`` execution or the MPI launcher ``mpi``."""
    numDevices: int = 1
    """Maximum compute devices made available to the lower-level run."""
    nPerNode: int = 1
    """MPI ranks per node launched automatically when ``parallelMode`` is ``mpi``."""
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

    _result: object | None = field(default=None, init=False, repr=False)
    _openpmdSession: object | None = field(default=None, init=False, repr=False)

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
        parser.add_argument("--openpmd-backend", default=None)
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
            "openpmd_backend": "openpmdBackend",
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
            "openpmd_backend": "openpmdBackend",
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
            "adaptiveSteps", "useReflections", "monochromatic", "backend", "openpmdBackend",
            "parallelMode", "numDevices", "nPerNode", "writeVtk", "devices",
            "minSampleRange", "maxSampleRange", "rngSeed",
        }
        for section in sections:
            for name, value in section.items():
                attr = aliases.get(name, name)
                if attr in allowed:
                    setattr(self, attr, value)
        return self

    def openPmdAttributes(self, *, numberOfSamples):
        min_sample = 0 if self.minSampleRange is None else int(self.minSampleRange)
        max_sample = int(numberOfSamples) - 1 if self.maxSampleRange is None else int(self.maxSampleRange)
        attributes = {
            "minRaysPerSample": self.minRaysPerSample,
            "maxRaysPerSample": self.maxRaysPerSample,
            "mseThreshold": self.mseThreshold,
            "repetitions": self.repetitions,
            "adaptiveSteps": self.adaptiveSteps,
            "useReflections": self.useReflections,
            "monochromatic": self.monochromatic,
            "backend": _preferredDefaultBackend() if self.backend is None else self.backend,
            "maxGpus": self.numDevices,
            "parallelMode": self.parallelMode,
            "minSampleRange": min_sample,
            "maxSampleRange": max_sample,
        }
        if self.rngSeed is not None:
            attributes["rngSeed"] = int(self.rngSeed)
        return attributes

    def _transportLaunchOptions(self):
        if str(self.parallelMode).strip().lower() != "mpi":
            return {}
        if isinstance(self.nPerNode, bool) or not isinstance(self.nPerNode, (int, np.integer)):
            raise ValueError("nPerNode must be a positive integer for MPI execution")
        ranks_per_node = int(self.nPerNode)
        if ranks_per_node < 1:
            raise ValueError("nPerNode must be a positive integer for MPI execution")
        return {
            "command_prefix": ["mpiexec", "-npernode", str(ranks_per_node)],
            # A scheduler allocation commonly spans nodes. Keep file-based
            # openPMD artifacts below the launch directory instead of /tmp so
            # they are visible when that directory is on shared storage.
            "workspace_dir": Path.cwd() / "IO" / "phiase_mpi",
        }

    def openStream(self, **kwargs):
        """Open a persistent openPMD transport session owned by this ``PhiASE``."""
        if self._openpmdSession is None:
            for name, value in self._transportLaunchOptions().items():
                kwargs.setdefault(name, value)
            if self.openpmdBackend is not None and "transport" not in kwargs:
                kwargs["transport"] = self.openpmdBackend
            self._openpmdSession = transport.openStream(**kwargs)
        return self._openpmdSession

    def closeStream(self):
        """Close this ``PhiASE`` object's persistent openPMD transport session."""
        session = self._openpmdSession
        self._openpmdSession = None
        return transport.closeStream(session)

    def run(self, gainMedium=None, crossSections=None, *, openpmdSession=None):
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

        if openpmdSession == "persistent":
            openpmdSession = self.openStream()
        elif openpmdSession == "interval":
            openpmdSession = None

        launch_options = {} if openpmdSession is not None else self._transportLaunchOptions()
        self._result = transport.runPhiASE(
            self,
            medium,
            cross_sections,
            transport=self.openpmdBackend,
            openpmdSession=openpmdSession,
            **launch_options,
        )
        return self

    def getResults(self):
        """Return the raw ASE result from the most recent ``run(...)`` call."""
        if self._result is None:
            raise RuntimeError("simulation has not been run yet")
        return self._result



class ConnectivityAverageBetaVolumeMapper:
    """Map point-centered ``betaCells`` to prism beta by vertex averaging.

    This matches the C++/Alpaka prism-beta kernel: each prism value is the
    arithmetic mean of the three triangle vertices on the lower and upper
    z-levels.
    """

    def map(self, medium):
        """Return prism-centered ``betaVolume`` for the supplied medium."""
        topology = medium.topology
        beta_cells = np.asarray(medium.get("betaCells").value, dtype=np.float64).reshape(
            (topology.numberOfPoints, topology.levels),
            order="F",
        )
        triangles = np.asarray(topology.trianglePointIndices, dtype=np.int64)
        if triangles.shape[0] != topology.numberOfTriangles:
            triangles = triangles.reshape((topology.numberOfTriangles, 3), order="F")
        beta_volume = np.empty((topology.numberOfTriangles, topology.levels - 1), dtype=np.float64)
        for level in range(topology.levels - 1):
            lower = beta_cells[triangles, level]
            upper = beta_cells[triangles, level + 1]
            beta_volume[:, level] = (lower.sum(axis=1) + upper.sum(axis=1)) / 6.0
        return beta_volume


LegacyGridDataBetaVolumeMapper = ConnectivityAverageBetaVolumeMapper


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
    topology: object | None = None
    """Static mesh topology used by geometry-aware state callbacks."""


@dataclass
class Simulation:
    """High-level Python wrapper for compiled C++/Alpaka simulation runs.

    Python sends the initial setup to the compiled backend and receives
    ``TimeStepState`` snapshots after completed steps. Register ``onInit`` for
    one-time Python setup before launch and ``onStep`` for snapshot consumers.
    ``beforeStep`` is retained only to report that per-step Python mutation is
    unsupported by compiled runs.
    """

    gainMedium: GainMedium
    """Geometry, material data, and current beta arrays."""
    pump: PumpProperties
    """Pump model that adds population to ``betaCells``."""
    phiASE: PhiASE
    """ASE configuration and execution handle."""
    timeIntegrationSolver: TimeIntegrationSolver | str
    """Compiled integrator name or descriptor with a ``name`` attribute."""
    timeStep: float
    """Physical time increment per step, in seconds."""
    crossSections: CrossSectionData | None = None
    """Shared spectra used by pump and ASE when not set on either object."""
    endTime: float | None = None
    """Optional target physical time for ``runUntil()``."""
    enableASE: bool = True
    """Whether compiled time-stepped runs include ASE depletion."""
    _time: float = field(default=0.0, init=False, repr=False)
    _step: int = field(default=0, init=False, repr=False)
    _initialized: bool = field(default=False, init=False, repr=False)
    _initCallbacks: list = field(default_factory=list, init=False, repr=False)
    _beforeStepCallbacks: list = field(default_factory=list, init=False, repr=False)
    _callbacks: list = field(default_factory=list, init=False, repr=False)
    _lastState: TimeStepState | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        if self.timeIntegrationSolver is None:
            raise ValueError("Simulation requires a timeIntegrationSolver")
        if not isinstance(self.timeIntegrationSolver, str) and not hasattr(self.timeIntegrationSolver, "name"):
            raise ValueError("Simulation requires a compiled time integrator name or descriptor with a .name attribute")
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

    def onStep(self, callback, *args, **kwargs):
        """Register a post-step callback.

        The callback signature is ``callback(state, *args, **kwargs)``.
        ``Simulation`` always supplies the completed ``TimeStepState`` as the
        first argument, then appends the positional and keyword arguments passed
        to ``onStep``. For example,
        ``simulation.onStep(write_vtk, output_dir, scale=5.5)`` calls
        ``write_vtk(state, output_dir, scale=5.5)`` after every completed step.

        Use this hook for logging, writing VTK files, explicit state storage,
        or other work that should consume the immutable step snapshot. Callback
        return values are ignored. The method returns ``self`` so registrations can be
        chained.
        """
        self._callbacks.append((callback, args, kwargs))
        return self

    def onInit(self, callback, *args, **kwargs):
        """Register a one-time initialization callback.

        The callback signature is ``callback(simulation, *args, **kwargs)``.
        ``Simulation`` supplies the live simulation object as the first
        argument, then appends the user arguments passed to ``onInit``. The hook
        runs once, immediately before the first step is evaluated.

        Use this hook to initialize or normalize mutable simulation inputs such
        as ``gainMedium``, ``pump``, ``phiASE``, or ``timeStep``. Callback return
        values are ignored. The method returns ``self`` for chaining.
        """
        self._initCallbacks.append((callback, args, kwargs))
        return self

    def beforeStep(self, callback, *args, **kwargs):
        """Register an unsupported legacy pre-step callback.

        Compiled simulations cannot call Python between C++-owned steps. The
        registration is stored so ``runSteps`` can raise a clear error instead
        of silently ignoring the callback.
        """
        self._beforeStepCallbacks.append((callback, args, kwargs))
        return self

    def runUntil(self, endtime=None, endTime=None, *, openpmdSession=None):
        """Advance steps until the configured or supplied end time is reached."""
        target = self.endTime if endtime is None and endTime is None else (endtime if endtime is not None else endTime)
        if target is None:
            raise ValueError("runUntil requires endtime or an endTime configured on construction")
        steps = 0
        while self._time + steps * self.timeStep < float(target) - 0.5 * self.timeStep:
            steps += 1
        if steps:
            self.runSteps(steps, openpmdSession=openpmdSession)
        return self

    def runSteps(self, steps, pumpSteps=None, *, openpmdSession=None):
        """Run exactly ``steps`` compiled C++/Alpaka time steps and return ``self``.

        ``pumpSteps`` optionally limits the pump contribution to the first
        ``pumpSteps`` calls in this run. When omitted, ``PumpProperties`` may
        provide a ``pumpSteps`` custom property. The complete time loop is
        executed by the C++ backend; Python only sends the initial setup and
        receives streamed snapshots.
        """
        if openpmdSession not in (None, "interval"):
            raise ValueError("compiled Simulation owns its C++ openPMD lifetime; openpmdSession is no longer supported")
        steps = int(steps)
        if steps <= 0:
            raise ValueError("steps must be positive")
        if self._beforeStepCallbacks:
            raise ValueError("compiled Simulation does not support Python beforeStep callbacks during C++-owned steps")
        self._runInitCallbacks()
        if pumpSteps is None and hasattr(self, "pump"):
            pumpSteps = self.pump.getProperty("pumpSteps")
        if pumpSteps is not None and int(pumpSteps) < 0:
            raise ValueError("pumpSteps must be non-negative")

        previous_step = self._step
        previous_time = self._time
        states = transport.runSimulation(
            self,
            steps=steps,
            pumpSteps=pumpSteps,
            transport=self.phiASE.openpmdBackend,
        )
        for raw_state in states:
            state = TimeStepState(
                step=previous_step + int(raw_state.step),
                time=previous_time + float(raw_state.time),
                betaCells=np.asarray(raw_state.betaCells, dtype=np.float64).copy(),
                betaVolume=np.asarray(raw_state.betaVolume, dtype=np.float64).copy(),
                phiAse=np.asarray(raw_state.phiAse, dtype=np.float64).copy(),
                dndtAse=np.asarray(raw_state.dndtAse, dtype=np.float64).copy(),
                dndtPump=np.asarray(raw_state.dndtPump, dtype=np.float64).copy(),
                aseResult=raw_state.aseResult,
                topology=self.gainMedium.topology,
            )
            self.gainMedium.get("betaCells").value = state.betaCells
            self.gainMedium.get("betaVolume").value = state.betaVolume
            self._lastState = state
            self._step = state.step
            self._time = state.time
            for callback, args, kwargs in self._callbacks:
                callback(state, *args, **kwargs)
        return self

    def step(self, *, openpmdSession=None):
        """Advance one compiled C++/Alpaka time step and return the completed state."""
        return self.runSteps(1, openpmdSession=openpmdSession).getLastState()

    def getLastState(self):
        """Return the most recent completed ``TimeStepState`` snapshot."""
        if self._lastState is None:
            raise RuntimeError("simulation has not completed a time step yet")
        return self._lastState

    def getResults(self):
        """Return the most recent completed ``TimeStepState`` snapshot.

        ``Simulation`` does not retain a full time-step history. Register an
        ``onStep`` callback to write or store per-step state explicitly.
        """
        return self.getLastState()

    @property
    def lastState(self):
        """Most recent completed ``TimeStepState`` snapshot."""
        return self.getLastState()

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
        for callback, args, kwargs in self._initCallbacks:
            callback(self, *args, **kwargs)

    def _updateBetaVolumeFromCells(self):
        beta_volume = ConnectivityAverageBetaVolumeMapper().map(self.gainMedium)
        self.gainMedium.get("betaVolume").value = beta_volume


TimeSteppedSimulation = Simulation
