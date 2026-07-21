# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""High-level Python simulation wrapper around pump, ASE, and time stepping."""

from __future__ import annotations

import os
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter

import numpy as np

from .alpakaUtils import AlpakaBackends
from .geometry import GainMedium
from .laser import (
    CrossSectionData,
    LaserProperties,
    MonteCarloPumpSolver,
    PlanarPumpRelay,
    Pump,
    SpectralDecomposition,
    SurfacePumpInjector,
    _PumpProperties,
    _PumpSource,
)
from .openpmd import backendFlat, transport
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

    propagationMode: str = "forward"
    """ASE propagation mode; only ``forward`` is supported."""
    minRays: int = 100000
    """Initial number of globally launched Monte Carlo rays."""
    maxRays: int = 100000
    """Maximum total number of globally launched rays during adaptive refinement."""
    forwardRayCount: int | None = None
    """Explicit fixed forward-ray count; disables adaptive refinement when set."""
    relativeStandardErrorThreshold: float = 0.1
    """Target one-sigma relative sampling uncertainty for ASE flux estimates."""
    repetitions: int = 4
    """Maximum repeated ASE estimates at a fixed ray count."""
    adaptiveSteps: int = 4
    """Maximum geometric ray-count increases from ``minRays`` to ``maxRays``."""
    useReflections: bool = False
    """Whether surface reflectivities affect forward propagation."""
    reflectionMaxIterations: int = 40
    """Maximum reflected-source passes after the direct volume-source pass."""
    reflectionTolerance: float = 1e-4
    """Stop reflected passes when their source-weight fraction is below this value."""
    surfaceReservoirSize: int = 32
    """Maximum reflected source records retained per physical boundary face."""
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
        parser.add_argument("--min-rays", "--min-rays-per-sample", dest="min_rays", type=int, default=None)
        parser.add_argument("--max-rays", "--max-rays-per-sample", dest="max_rays", type=int, default=None)
        parser.add_argument("--propagation-mode", choices=("forward",), default=None)
        parser.add_argument("--forward-ray-count", type=int, default=None)
        parser.add_argument("--relative-standard-error-threshold", type=float, default=None)
        parser.add_argument("--reflection-max-iterations", type=int, default=None)
        parser.add_argument("--reflection-tolerance", type=float, default=None)
        parser.add_argument("--surface-reservoir-size", type=int, default=None)
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
            "min_rays": "minRays",
            "max_rays": "maxRays",
            "min_rays_per_sample": "minRays",
            "max_rays_per_sample": "maxRays",
            "propagation_mode": "propagationMode",
            "forward_ray_count": "forwardRayCount",
            "relative_standard_error_threshold": "relativeStandardErrorThreshold",
            "reflection_max_iterations": "reflectionMaxIterations",
            "reflection_tolerance": "reflectionTolerance",
            "surface_reservoir_size": "surfaceReservoirSize",
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
            "minRaysPerSample": "minRays",
            "maxRaysPerSample": "maxRays",
            "min_rays": "minRays",
            "max_rays": "maxRays",
            "min_rays_per_sample": "minRays",
            "max_rays_per_sample": "maxRays",
            "propagation_mode": "propagationMode",
            "forward_ray_count": "forwardRayCount",
            "relative_standard_error_threshold": "relativeStandardErrorThreshold",
            "reflection_max_iterations": "reflectionMaxIterations",
            "reflection_tolerance": "reflectionTolerance",
            "surface_reservoir_size": "surfaceReservoirSize",
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
            "minRays", "maxRays", "propagationMode", "forwardRayCount",
            "relativeStandardErrorThreshold", "reflectionMaxIterations", "reflectionTolerance",
            "surfaceReservoirSize", "repetitions", "adaptiveSteps", "useReflections", "monochromatic",
            "backend", "openpmdBackend", "parallelMode", "numDevices", "nPerNode", "writeVtk", "devices",
            "minSampleRange", "maxSampleRange", "rngSeed",
        }
        for section in sections:
            for name, value in section.items():
                if name in {"forwardRayLength", "forward_ray_length"}:
                    raise ValueError(
                        "forward_ray_length is retired; forward rays now propagate to their physical boundary"
                    )
                if name in {"mseThreshold", "mse_threshold"}:
                    raise ValueError(
                        "mse_threshold is retired; configure relative_standard_error_threshold instead"
                    )
                attr = aliases.get(name, name)
                if attr in allowed:
                    setattr(self, attr, value)
        return self

    def openPmdAttributes(self, *, numberOfSamples):
        if str(self.propagationMode).strip().lower() != "forward":
            raise ValueError("PhiASE.propagationMode must be 'forward'")
        min_rays = int(self.minRays)
        max_rays = int(self.maxRays)
        adaptive_steps = int(self.adaptiveSteps)
        forward_ray_count = 0 if self.forwardRayCount is None else int(self.forwardRayCount)
        if min_rays == 0:
            raise ValueError("PhiASE.minRays must be greater than zero")
        if max_rays < min_rays:
            raise ValueError("PhiASE.maxRays must be greater than or equal to PhiASE.minRays")
        if adaptive_steps < 0:
            raise ValueError("PhiASE.adaptiveSteps must not be negative")
        if forward_ray_count < 0:
            raise ValueError("PhiASE.forwardRayCount must not be negative")
        min_sample = 0 if self.minSampleRange is None else int(self.minSampleRange)
        max_sample = int(numberOfSamples) - 1 if self.maxSampleRange is None else int(self.maxSampleRange)
        attributes = {
            "minRays": min_rays,
            "maxRays": max_rays,
            "propagationMode": self.propagationMode,
            "forwardRayCount": forward_ray_count,
            "relativeStandardErrorThreshold": self.relativeStandardErrorThreshold,
            "reflectionMaxIterations": self.reflectionMaxIterations,
            "reflectionTolerance": self.reflectionTolerance,
            "surfaceReservoirSize": self.surfaceReservoirSize,
            "repetitions": self.repetitions,
            "adaptiveSteps": adaptive_steps,
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
        mpiexec_extra_args = shlex.split(os.environ.get("HASE_MPIEXEC_EXTRA_ARGS", ""))
        return {
            "command_prefix": [
                "mpiexec",
                *mpiexec_extra_args,
                "-npernode",
                str(ranks_per_node),
            ],
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
    volumePhiAse: np.ndarray | None = None
    """Native volume-centered ASE flux, when provided by the compiled backend."""
    volumeDndtAse: np.ndarray | None = None
    """Native volume-centered ASE depletion contribution, when available."""
    volumeStandardError: np.ndarray | None = None
    """Native volume-centered absolute one-sigma sampling uncertainty, when available."""
    volumeRelativeStandardError: np.ndarray | None = None
    """Native volume-centered relative one-sigma sampling uncertainty, when available."""
    volumeTotalRays: np.ndarray | None = None
    """Native volume-centered ray visit counts, when available."""

    @property
    def beta_cells(self):
        return self.betaCells

    @property
    def beta_volume(self):
        return self.betaVolume

    @property
    def phi_ase(self):
        return self.phiAse

    @property
    def dndt_ase(self):
        return self.dndtAse

    @property
    def dndt_pump(self):
        return self.dndtPump

    @property
    def ase_result(self):
        return self.aseResult

    @property
    def volume_phi_ase(self):
        return self.volumePhiAse

    @property
    def volume_dndt_ase(self):
        return self.volumeDndtAse

    @property
    def volume_standard_error(self):
        return self.volumeStandardError

    @property
    def volume_relative_standard_error(self):
        return self.volumeRelativeStandardError

    @property
    def volume_total_rays(self):
        return self.volumeTotalRays


@dataclass(init=False)
class Simulation:
    """High-level Python wrapper for compiled C++/Alpaka simulation runs.

    Python sends the initial setup to the compiled backend and receives
    ``TimeStepState`` snapshots after completed steps. Register ``on_init`` for
    one-time Python setup before launch and ``on_step`` for snapshot consumers.
    ``beforeStep`` is retained only to report that per-step Python mutation is
    unsupported by compiled runs.
    """

    gainMedium: GainMedium
    pump: _PumpProperties | None
    phiASE: PhiASE
    timeIntegrationSolver: TimeIntegrationSolver | str
    timeStep: float
    crossSections: CrossSectionData | None
    endTime: float | None
    enableASE: bool
    prePump: bool
    pumpSolver: MonteCarloPumpSolver
    maxSteps: int | None
    reportTimings: bool
    _pumpRegistrations: list
    _time: float
    _step: int
    _initialized: bool
    _initCallbacks: list
    _beforeStepCallbacks: list
    _callbacks: list
    _lastState: TimeStepState | None

    def __init__(
        self,
        *,
        gain_medium,
        phi_ase,
        time_integrator,
        time_step_size,
        pump_solver=None,
        cross_sections=None,
        max_steps=None,
        max_time=None,
        enable_ase=True,
        pre_pump=False,
        report_timings=False,
    ):
        self.gainMedium = gain_medium
        self.pump = None
        self.phiASE = phi_ase
        self.timeIntegrationSolver = time_integrator
        self.timeStep = float(time_step_size)
        self.crossSections = cross_sections
        self.endTime = max_time
        self.enableASE = bool(enable_ase)
        self.prePump = bool(pre_pump)
        self.pumpSolver = MonteCarloPumpSolver() if pump_solver is None else pump_solver
        self.maxSteps = None if max_steps is None else int(max_steps)
        self.reportTimings = bool(report_timings)
        self._pumpRegistrations = []
        self._time = 0.0
        self._step = 0
        self._initialized = False
        self._initCallbacks = []
        self._beforeStepCallbacks = []
        self._callbacks = []
        self._lastState = None
        self.__post_init__()

    def __post_init__(self):
        if self.timeIntegrationSolver is None:
            raise ValueError("Simulation requires a time_integrator")
        if not isinstance(self.timeIntegrationSolver, str) and not hasattr(self.timeIntegrationSolver, "name"):
            raise ValueError("Simulation requires a compiled time integrator name or descriptor with a .name attribute")
        if self.timeStep <= 0.0:
            raise ValueError("time_step_size must be positive")
        if not isinstance(self.pumpSolver, MonteCarloPumpSolver):
            raise TypeError("pump_solver must be MonteCarloPumpSolver")
        if self.maxSteps is not None and self.maxSteps <= 0:
            raise ValueError("max_steps must be positive")
        if self.crossSections is None and (
            self.phiASE.spectralProperties is not None or self.phiASE.crossSections is not None
        ):
            self.crossSections = self._resolveSpectralProperties()
        if self.phiASE.crossSections is None and self.crossSections is not None:
            self.phiASE.crossSections = self.crossSections
        if self.phiASE.spectralProperties is None and self.crossSections is not None:
            self.phiASE.spectralProperties = self.crossSections
        self._ensureStateArrays()

    def _resolveSpectralProperties(self):
        if self.phiASE.spectralProperties is not None:
            return self.phiASE.spectralProperties
        if self.phiASE.crossSections is not None:
            return self.phiASE.crossSections
        if self.pump is not None and self.pump.sources:
            return self.pump.sources[0].crossSections
        raise ValueError("Simulation requires spectral properties via Simulation.crossSections, phiASE, or pump")

    def add_pump(self, pump, injection_method, *, relays=()):
        """Register a physical pump and its numerical injection method."""
        if self._initialized:
            raise RuntimeError("pumps must be added before the simulation is initialized")
        if not isinstance(pump, Pump):
            raise TypeError("pump must be a Pump")
        if not isinstance(injection_method, SurfacePumpInjector):
            raise TypeError("injection_method must be SurfacePumpInjector")
        relays = tuple(relays)
        if not all(isinstance(relay, PlanarPumpRelay) for relay in relays):
            raise TypeError("relays must contain PlanarPumpRelay values")
        self._pumpRegistrations.append((pump, injection_method, relays))
        self.pump = _PumpProperties(
            sources=tuple(
                _PumpSource(
                    surfaceDomains=injector.surface_domains,
                    totalPower=physical.total_power,
                    spectrum=physical.spectrum,
                    crossSections=physical.cross_sections,
                    angularDistribution=physical.angular_distribution,
                    profile=physical.profile,
                    relays=registered_relays,
                )
                for physical, injector, registered_relays in self._pumpRegistrations
            ),
            rayCount=self.pumpSolver.ray_count,
            rngSeed=self.pumpSolver.seed,
            pumpSteps=self.pumpSolver.max_steps,
        )
        if self.crossSections is None:
            self.crossSections = self._resolveSpectralProperties()
        if self.phiASE.crossSections is None:
            self.phiASE.crossSections = self.crossSections
        if self.phiASE.spectralProperties is None:
            self.phiASE.spectralProperties = self.crossSections
        return self

    def on_step(self, callback, *args, **kwargs):
        """Register a callback that receives each completed state snapshot."""
        return self.onStep(callback, *args, **kwargs)

    def on_init(self, callback, *args, **kwargs):
        """Register a callback that runs once before compiled execution."""
        return self.onInit(callback, *args, **kwargs)

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

    def run_until(self, max_time=None):
        """Advance to ``max_time`` or the constructor's configured maximum."""
        target = self.endTime if max_time is None else max_time
        if target is None:
            raise ValueError("run_until requires max_time or a configured max_time")
        steps = 0
        while self._time + steps * self.timeStep < float(target) - 0.5 * self.timeStep:
            steps += 1
        if steps:
            self.step(steps)
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

        Internal transport helper. ``pumpSteps`` limits pump contribution to
        the first compiled steps and defaults to the registered pump solver. The complete time loop is
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
        if self.pump is None:
            raise ValueError("Simulation requires at least one pump registered with add_pump")
        self._runInitCallbacks()
        if pumpSteps is None:
            pumpSteps = getattr(self.pump, "pumpSteps", None)
        if pumpSteps is not None and int(pumpSteps) < 0:
            raise ValueError("pumpSteps must be non-negative")

        previous_step = self._step
        previous_time = self._time
        run_started = perf_counter() if self.reportTimings else None
        transport_started = perf_counter() if self.reportTimings else None
        states = transport.runSimulation(
            self,
            steps=steps,
            pumpSteps=pumpSteps,
            transport=self.phiASE.openpmdBackend,
            **self.phiASE._transportLaunchOptions(),
        )
        transport_seconds = perf_counter() - transport_started if self.reportTimings else 0.0
        state_materialization_seconds = 0.0
        callback_seconds = {}
        for raw_state in states:
            state_started = perf_counter() if self.reportTimings else None
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
                volumePhiAse=(
                    None if getattr(raw_state, "volumePhiAse", None) is None
                    else np.asarray(raw_state.volumePhiAse, dtype=np.float64).copy()
                ),
                volumeDndtAse=(
                    None if getattr(raw_state, "volumeDndtAse", None) is None
                    else np.asarray(raw_state.volumeDndtAse, dtype=np.float64).copy()
                ),
                volumeStandardError=(
                    None if getattr(raw_state, "volumeStandardError", None) is None
                    else np.asarray(raw_state.volumeStandardError, dtype=np.float64).copy()
                ),
                volumeRelativeStandardError=(
                    None if getattr(raw_state, "volumeRelativeStandardError", None) is None
                    else np.asarray(raw_state.volumeRelativeStandardError, dtype=np.float64).copy()
                ),
                volumeTotalRays=(
                    None if getattr(raw_state, "volumeTotalRays", None) is None
                    else np.asarray(raw_state.volumeTotalRays, dtype=np.uint32).copy()
                ),
            )
            if hasattr(self.gainMedium.topology, "cellPointIndices"):
                self.gainMedium.get("betaCells").value = backendFlat(state.betaCells.reshape(-1, order="F"))
                self.gainMedium.get("betaVolume").value = backendFlat(state.betaVolume.reshape(-1, order="F"))
            else:
                self.gainMedium.get("betaCells").value = state.betaCells
                self.gainMedium.get("betaVolume").value = state.betaVolume
            self._lastState = state
            self._step = state.step
            self._time = state.time
            if self.reportTimings:
                state_materialization_seconds += perf_counter() - state_started
            for callback, args, kwargs in self._callbacks:
                callback_started = perf_counter() if self.reportTimings else None
                callback(state, *args, **kwargs)
                if self.reportTimings:
                    callback_name = getattr(callback, "__qualname__", getattr(callback, "__name__", type(callback).__name__))
                    callback_seconds[callback_name] = callback_seconds.get(callback_name, 0.0) + (
                        perf_counter() - callback_started
                    )
        if self.reportTimings:
            total_seconds = perf_counter() - run_started
            callbacks_total = sum(callback_seconds.values())
            print(
                "HASE frontend timing: "
                f"steps={len(states)} total={total_seconds:.6f}s "
                f"compiled_transport_and_decode={transport_seconds:.6f}s "
                f"snapshot_materialization={state_materialization_seconds:.6f}s "
                f"callbacks={callbacks_total:.6f}s"
            )
            for callback_name, seconds in sorted(callback_seconds.items(), key=lambda item: item[1], reverse=True):
                print(f"HASE frontend callback timing: {callback_name}={seconds:.6f}s")
        return self

    def step(self, nsteps=1, *, pump_steps=None):
        """Advance ``nsteps`` time steps, following the PICMI call pattern."""
        if int(nsteps) <= 0:
            raise ValueError("nsteps must be positive")
        if pump_steps is not None and int(pump_steps) < 0:
            raise ValueError("pump_steps must be non-negative")
        self.runSteps(int(nsteps), pumpSteps=pump_steps)
        return self

    def get_last_state(self):
        return self.getLastState()

    @property
    def current_step(self):
        return self._step

    @property
    def current_time(self):
        return self._time

    @property
    def gain_medium(self):
        return self.gainMedium

    @property
    def phi_ase(self):
        return self.phiASE

    @property
    def time_integrator(self):
        return self.timeIntegrationSolver

    @property
    def time_step_size(self):
        return self.timeStep

    @property
    def pump_solver(self):
        return self.pumpSolver

    @property
    def cross_sections(self):
        return self.crossSections

    @property
    def enable_ase(self):
        return self.enableASE

    @property
    def pre_pump(self):
        return self.prePump

    @property
    def max_steps(self):
        return self.maxSteps

    @property
    def max_time(self):
        return self.endTime

    @property
    def pumps(self):
        return tuple(physical for physical, _injector, _relays in self._pumpRegistrations)

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
