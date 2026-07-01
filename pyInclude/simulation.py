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

from .geometry import GainMedium
from .laser import CrossSectionData, LaserProperties, PumpProperties, SpectralDecomposition
from .openpmd import transport
from .pumping import BetaIntegrationGaussianSolver, Constants
from .timeIntegration import TimeDerivative, TimeIntegrationSolver


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
    topology: object | None = None
    """Static mesh topology used by geometry-aware state callbacks."""


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
    enableAse: bool = True
    """Whether each derivative evaluation runs ASE depletion through ``PhiASE``."""
    prePump: bool = False
    """When true, the first time step advances without ASE so pump can seed beta."""

    _time: float = field(default=0.0, init=False, repr=False)
    _step: int = field(default=0, init=False, repr=False)
    _initialized: bool = field(default=False, init=False, repr=False)
    _initCallbacks: list = field(default_factory=list, init=False, repr=False)
    _beforeStepCallbacks: list = field(default_factory=list, init=False, repr=False)
    _callbacks: list = field(default_factory=list, init=False, repr=False)
    _lastState: TimeStepState | None = field(default=None, init=False, repr=False)
    _pumpEnabled: bool = field(default=True, init=False, repr=False)
    _openpmdSession: object | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        if self.timeIntegrationSolver is None:
            raise ValueError("Simulation requires a timeIntegrationSolver")
        if not isinstance(self.timeIntegrationSolver, str):
            has_step = hasattr(self.timeIntegrationSolver, "step")
            has_simulation_step = hasattr(self.timeIntegrationSolver, "stepSimulation")
            if not has_step and not has_simulation_step:
                raise ValueError("Simulation requires a timeIntegrationSolver with step(...) or stepSimulation(...)")
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
        """Register a pre-step callback.

        The callback signature is ``callback(simulation, *args, **kwargs)``.
        ``Simulation`` supplies the live simulation object as the first
        argument, then appends the user arguments passed to ``beforeStep``. The
        hook runs before every step, after one-time ``onInit`` callbacks.

        Use this hook for controlled changes that must happen before derivative
        evaluation, such as time-dependent pump settings. Callback return values
        are ignored. The method returns ``self`` for chaining.
        """
        self._beforeStepCallbacks.append((callback, args, kwargs))
        return self

    def _withOpenPmdSession(self, openpmdSession):
        if openpmdSession is None:
            selected = "auto" if self.phiASE.openpmdBackend is None else str(self.phiASE.openpmdBackend).strip().lower()
            if selected in {"auto", "adios-sst"}:
                return self.phiASE.openStream(), True
            return None, False
        if openpmdSession == "interval":
            return None, False
        if openpmdSession == "persistent":
            return self.phiASE.openStream(), True
        return openpmdSession, False

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
        """Run exactly ``steps`` calls to ``step()`` and return ``self``.

        ``pumpSteps`` optionally limits the pump contribution to the first
        ``pumpSteps`` calls in this run. When omitted, ``PumpProperties`` may
        provide a ``pumpSteps`` custom property. If neither is set, the pump is
        active for every step, matching the historical ``runSteps(steps)``
        behavior. ``PumpProperties.pumpSubsteps`` is separate: it controls the
        internal pump integration resolution inside one pumped simulation step.
        """
        if not self._usesCompiledBackend():
            return self._runStepsLegacy(steps, pumpSteps=pumpSteps, openpmdSession=openpmdSession)
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
        """Advance one time step and return the completed ``TimeStepState``."""
        if not self._usesCompiledBackend():
            return self._legacyPythonStep(openpmdSession=openpmdSession)
        return self.runSteps(1, openpmdSession=openpmdSession).getLastState()

    def _usesCompiledBackend(self):
        if not isinstance(self.phiASE, PhiASE):
            return False
        return "run" not in vars(self.phiASE) and "getResults" not in vars(self.phiASE)

    def _runStepsLegacy(self, steps, pumpSteps=None, *, openpmdSession=None):
        steps = int(steps)
        session, close_session = self._withOpenPmdSession(openpmdSession)
        previous_session = self._openpmdSession
        self._openpmdSession = session
        if pumpSteps is None and hasattr(self, "pump"):
            pumpSteps = self.pump.getProperty("pumpSteps")
        previousPumpEnabled = self._pumpEnabled
        try:
            if pumpSteps is None:
                for _ in range(steps):
                    self.step(openpmdSession=session)
            else:
                pumpSteps = int(pumpSteps)
                if pumpSteps < 0:
                    raise ValueError("pumpSteps must be non-negative")
                for index in range(steps):
                    self._pumpEnabled = index < pumpSteps
                    self.step(openpmdSession=session)
        finally:
            self._openpmdSession = previous_session
            self._pumpEnabled = previousPumpEnabled
            if close_session:
                self.phiASE.closeStream()
        return self

    def _legacyPythonStep(self, *, openpmdSession=None):
        """Legacy Python-side step implementation retained for focused internal tests."""
        previous_session = self._openpmdSession
        if openpmdSession is not None:
            self._openpmdSession = openpmdSession
        try:
            self._runInitCallbacks()
            for callback, args, kwargs in self._beforeStepCallbacks:
                callback(self, *args, **kwargs)

            beta_cells = np.asarray(self.gainMedium.get("betaCells").value, dtype=np.float64).reshape(
                self.gainMedium.get("betaCells").expectedShape,
                order="F",
            )

            if hasattr(self.timeIntegrationSolver, "stepSimulation"):
                integration_result = self.timeIntegrationSolver.stepSimulation(
                    self,
                    beta_cells.copy(),
                    self._time,
                    self.timeStep,
                )
            else:
                integration_result = self.timeIntegrationSolver.step(
                    self._timeDerivative,
                    beta_cells.copy(),
                    self._time,
                    self.timeStep,
                )
        finally:
            self._openpmdSession = previous_session
        updated_beta = integration_result.betaCells
        derivative = integration_result.evaluation

        self.gainMedium.get("betaCells").value = np.clip(updated_beta, 0.0, 1.0)
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
            topology=self.gainMedium.topology,
        )
        self._lastState = state
        for callback, args, kwargs in self._callbacks:
            callback(state, *args, **kwargs)
        return state

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

    def _dndtPump(self, beta_cells):
        if not self._pumpEnabled:
            return np.zeros_like(beta_cells, dtype=np.float64)
        pump_duration = self.pump.activeDuration(self.timeStep)
        solver = self.pump.solver if self.pump.solver is not None else BetaIntegrationGaussianSolver()
        # Pump solvers return β(t + Δt); this converts the update into ∂β/∂t|pump.
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
        phi_ase, ase_result = self._runPhiAseForBeta(beta_cells)
        return self._timeDerivativeWithFrozenPhiAse(beta_cells, time, phi_ase, ase_result)

    def _aseEnabledForCurrentStep(self):
        return self.enableAse and not (self.prePump and self._step == 0)

    def _runPhiAseForBeta(self, beta_cells):
        beta_cells = np.asarray(beta_cells, dtype=np.float64).reshape(
            self.gainMedium.get("betaCells").expectedShape,
            order="F",
        )
        self.gainMedium.get("betaCells").value = beta_cells
        self._updateBetaVolumeFromCells()
        if not self._aseEnabledForCurrentStep():
            return np.zeros_like(beta_cells, dtype=np.float64), None

        self.phiASE.run(
            gainMedium=self.gainMedium,
            crossSections=self.crossSections,
            openpmdSession=self._openpmdSession,
        )
        ase_result = self.phiASE.getResults()
        phi_ase = np.asarray(ase_result.phiAse, dtype=np.float64).reshape(beta_cells.shape, order="F")
        return phi_ase.copy(), ase_result

    def _timeDerivativeWithFrozenPhiAse(self, beta_cells, time, phi_ase, ase_result):
        beta_cells = np.asarray(beta_cells, dtype=np.float64).reshape(
            self.gainMedium.get("betaCells").expectedShape,
            order="F",
        )
        phi_ase = np.asarray(phi_ase, dtype=np.float64).reshape(beta_cells.shape, order="F")
        self.gainMedium.get("betaCells").value = beta_cells
        self._updateBetaVolumeFromCells()
        dndt_ase = (
            self._aseDerivative(phi_ase, betaCells=beta_cells)
            if self._aseEnabledForCurrentStep()
            else np.zeros_like(beta_cells)
        )
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
