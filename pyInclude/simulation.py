# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""High-level Python simulation wrapper around pump, ASE, and time stepping."""

from __future__ import annotations

from collections.abc import Mapping
import os
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter

import numpy as np

from hase_transport import PrimitiveDescription, field as transportField, reference
from .alpakaUtils import AlpakaBackends
from .fields import backendFlat
from .physical import Domain, GainMedium, OpticalComponent, validateComponentOverlap
from .frontendState import projectCellMask, projectFrontendState, projectSurfaceDomain
from .laser import (
    PlanarPumpRelay,
    Pump,
    SurfacePumpInjector,
    _PumpProperties,
    _PumpSource,
)
from .openpmd import transport
from .timeIntegration import TimeIntegrationSolver


HASE_CONFIGURE_HINT = "Run `hase-configure` to generate a matching backend/openPMD setup."

SIMULATION_OUTPUT_FIELDS = (
    "beta_volume",
    "phi_ase",
    "standard_error",
    "relative_standard_error",
    "total_rays",
    "dndt_ase",
    "dndt_pump",
)
SIMULATION_CONTROL_FIELDS = ("beta_volume", "cross_sections")


def autonomous_final(number_of_steps):
    """Return the explicit output index list for a final-state-only autonomous run."""
    final_step = int(number_of_steps)
    if final_step <= 0:
        raise ValueError("number_of_steps must be positive")
    return (final_step,)


def _optional_state_array(values, dtype):
    return None if values is None else np.asarray(values, dtype=dtype).copy()


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


def _validate_alpaka_backend(phi_ase):
    """Validate process-local compute availability for a non-MPI launch."""
    if str(getattr(phi_ase, "parallelMode", "single")).strip().lower() == "mpi":
        return
    try:
        available = AlpakaBackends.all()
    except Exception as exc:
        raise RuntimeError(
            "HASEonGPU could not query Alpaka backends available to the Python "
            f"launcher. {HASE_CONFIGURE_HINT}"
        ) from exc
    configured = getattr(phi_ase, "backend", None)
    selected = _preferredDefaultBackend() if configured is None else str(configured)
    if selected in available:
        return
    available_text = "\n".join(f"  {backend}" for backend in available) or "  (none)"
    raise RuntimeError(
        f"Alpaka backend '{selected}' is unavailable in the Python launcher process. "
        f"Available Alpaka backends:\n{available_text}\n{HASE_CONFIGURE_HINT}"
    )


def _validate_launch_backends(phi_ase, *, openpmd_session=None, openpmd_backend=None):
    """Reject unavailable compute and transport selections before backend launch."""
    _validate_alpaka_backend(phi_ase)
    if openpmd_session is None:
        selected_openpmd = phi_ase.openpmdBackend if openpmd_backend is None else openpmd_backend
        transport._ensure_backend_available(selected_openpmd)


@dataclass
class PhiASE:
    """Configure and run the ASE flux calculation for one gain-medium state.

    ``Simulation`` normally owns this object and calls ``run(...)`` during each
    time-step derivative evaluation. Advanced users can also call ``run``
    directly with a ``GainMedium``. Cross sections are owned by each active
    component's ``Material``.
    """

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
    trackRayVisits: bool = False
    """Count per-cell ray visits in ``totalRays`` for diagnostics."""
    repetitions: int = 4
    """Maximum repeated ASE estimates at a fixed ray count."""
    adaptiveSteps: int = 4
    """Maximum geometric ray-count increases from ``minRays`` to ``maxRays``."""
    useReflections: bool = False
    """Whether surface reflectivities affect forward propagation."""
    reflectionMode: str = "direct"
    """Reflection source representation: ``direct`` exact hits or bounded ``srm`` bins."""
    surfaceReservoirSize: int = 64
    """Number of retained statistical ray records per boundary face in SRM mode."""
    srmPositionMode: str = "exact"
    """SRM relaunch position: the retained ``exact`` hit or the face ``centroid``."""
    reflectionMaxIterations: int = 40
    """Maximum reflected-source passes after the direct volume-source pass."""
    reflectionTolerance: float = 1e-4
    """Stop reflected passes when their source-weight fraction is below this value."""
    monochromatic: bool = False
    """Use only the first spectral samples instead of wavelength integration."""

    backend: str = None
    """Alpaka backend name; inspect valid strings with ``AlpakaBackends.all()``."""
    openpmdBackend: str | None = "auto"
    """openPMD backend; ``auto`` prefers ADIOS, SST, then HDF5 when supported."""
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
    ase_steps: int | None = None
    """Outer steps that include ASE; ``None`` or zero disables ASE in ``Simulation``."""

    _result: object | None = field(default=None, init=False, repr=False)
    _openpmdSession: object | None = field(default=None, init=False, repr=False)

    def _transportDescription(self):
        return PrimitiveDescription(
            "phiAse",
            fields=(
                transportField("propagationMode"),
                transportField("minRays"),
                transportField("maxRays"),
                transportField("forwardRayCount", optional=True),
                transportField("relativeStandardErrorThreshold"),
                transportField("trackRayVisits"),
                transportField("repetitions"),
                transportField("adaptiveSteps"),
                transportField("useReflections"),
                transportField("reflectionMode"),
                transportField("surfaceReservoirSize"),
                transportField("srmPositionMode"),
                transportField("reflectionMaxIterations"),
                transportField("reflectionTolerance"),
                transportField("monochromatic"),
                transportField(
                    "backend",
                    lambda owner: _preferredDefaultBackend() if owner.backend is None else owner.backend,
                ),
                transportField("parallelMode"),
                transportField("numDevices"),
                transportField("writeVtk"),
                transportField(
                    "devices",
                    lambda owner: owner.devices if owner.devices else None,
                    axes=("device",),
                    optional=True,
                ),
                transportField("minSampleRange", optional=True),
                transportField("maxSampleRange", optional=True),
                transportField("rngSeed", optional=True),
                transportField("aseSteps", "ase_steps", optional=True),
            ),
        )

    def __post_init__(self):
        if self.reflectionMode not in {"direct", "srm"}:
            raise ValueError("PhiASE.reflectionMode must be 'direct' or 'srm'")
        if isinstance(self.surfaceReservoirSize, bool) or not isinstance(
            self.surfaceReservoirSize, (int, np.integer)
        ):
            raise TypeError("PhiASE.surfaceReservoirSize must be an integer")
        if self.surfaceReservoirSize <= 0:
            raise ValueError("PhiASE.surfaceReservoirSize must be positive")
        if self.srmPositionMode not in {"exact", "centroid"}:
            raise ValueError("PhiASE.srmPositionMode must be 'exact' or 'centroid'")
        if self.ase_steps is not None:
            if isinstance(self.ase_steps, bool) or not isinstance(self.ase_steps, (int, np.integer)):
                raise TypeError("PhiASE.ase_steps must be an integer or None")
            if self.ase_steps < 0:
                raise ValueError("PhiASE.ase_steps must be non-negative")

    @classmethod
    def fromYaml(cls, filename, **overrides):
        """Read ``simulation.phi_ase`` from a schema-v3 configuration."""
        from .configuration import _loadYaml, _phiAse

        _path, data = _loadYaml(filename)
        simulation = data.get("simulation")
        if not isinstance(simulation, dict) or not isinstance(simulation.get("phi_ase"), dict):
            raise ValueError("schema-v3 PhiASE YAML requires simulation.phi_ase")
        obj = _phiAse(simulation["phi_ase"])
        for name, value in overrides.items():
            setattr(obj, name, value)
        obj.__post_init__()
        return obj

    @staticmethod
    def addArguments(parser):
        """Add command-line arguments that map to ``PhiASE`` settings."""
        parser.add_argument("--phi-ase-config", default=None, help="YAML file with PhiASE compute/experiment settings")
        parser.add_argument("--min-rays", "--min-rays-per-sample", dest="min_rays", type=int, default=None)
        parser.add_argument("--max-rays", "--max-rays-per-sample", dest="max_rays", type=int, default=None)
        parser.add_argument("--propagation-mode", choices=("forward",), default=None)
        parser.add_argument("--forward-ray-count", type=int, default=None)
        parser.add_argument("--relative-standard-error-threshold", type=float, default=None)
        ray_visits = parser.add_mutually_exclusive_group()
        ray_visits.add_argument("--track-ray-visits", dest="track_ray_visits", action="store_true", default=None)
        ray_visits.add_argument("--no-track-ray-visits", dest="track_ray_visits", action="store_false")
        parser.add_argument("--reflection-max-iterations", type=int, default=None)
        parser.add_argument("--reflection-tolerance", type=float, default=None)
        reflections = parser.add_mutually_exclusive_group()
        reflections.add_argument("--use-reflections", dest="use_reflections", action="store_true", default=None)
        reflections.add_argument("--no-reflections", dest="use_reflections", action="store_false")
        parser.add_argument("--reflection-mode", choices=("direct", "srm"), default=None)
        parser.add_argument("--surface-reservoir-size", type=int, default=None)
        parser.add_argument("--srm-position-mode", choices=("exact", "centroid"), default=None)
        parser.add_argument("--repetitions", type=int, default=None)
        parser.add_argument("--adaptive-steps", type=int, default=None)
        spectral_mode = parser.add_mutually_exclusive_group()
        spectral_mode.add_argument("--monochromatic", dest="monochromatic", action="store_true", default=None)
        spectral_mode.add_argument("--polychromatic", dest="monochromatic", action="store_false")
        parser.add_argument("--backend", default=None)
        parser.add_argument("--openpmd-backend", default=None)
        parser.add_argument("--parallel-mode", default=None)
        parser.add_argument("--max-gpus", type=int, default=None)
        parser.add_argument("--n-per-node", type=int, default=None)
        vtk_output = parser.add_mutually_exclusive_group()
        vtk_output.add_argument("--write-vtk", dest="write_vtk", action="store_true", default=None)
        vtk_output.add_argument("--no-write-vtk", dest="write_vtk", action="store_false")
        parser.add_argument("--devices", nargs="+", type=int, default=None)
        parser.add_argument("--min-sample-range", type=int, default=None)
        parser.add_argument("--max-sample-range", type=int, default=None)
        parser.add_argument("--rng-seed", type=int, default=None)
        parser.add_argument("--ase-steps", type=int, default=None)
        return parser

    @classmethod
    def fromArgs(cls, args, **overrides):
        """Create a ``PhiASE`` configuration from parsed argparse results."""
        config = getattr(args, "phi_ase_config", None)
        obj = cls.fromYaml(config) if config else cls()
        mapping = {
            "min_rays": "minRays",
            "max_rays": "maxRays",
            "min_rays_per_sample": "minRays",
            "max_rays_per_sample": "maxRays",
            "propagation_mode": "propagationMode",
            "forward_ray_count": "forwardRayCount",
            "relative_standard_error_threshold": "relativeStandardErrorThreshold",
            "track_ray_visits": "trackRayVisits",
            "reflection_max_iterations": "reflectionMaxIterations",
            "reflection_tolerance": "reflectionTolerance",
            "use_reflections": "useReflections",
            "reflection_mode": "reflectionMode",
            "surface_reservoir_size": "surfaceReservoirSize",
            "srm_position_mode": "srmPositionMode",
            "repetitions": "repetitions",
            "adaptive_steps": "adaptiveSteps",
            "monochromatic": "monochromatic",
            "backend": "backend",
            "openpmd_backend": "openpmdBackend",
            "parallel_mode": "parallelMode",
            "max_gpus": "numDevices",
            "n_per_node": "nPerNode",
            "write_vtk": "writeVtk",
            "devices": "devices",
            "min_sample_range": "minSampleRange",
            "max_sample_range": "maxSampleRange",
            "rng_seed": "rngSeed",
            "ase_steps": "ase_steps",
        }
        for arg_name, attr_name in mapping.items():
            value = getattr(args, arg_name, None)
            if value is not None:
                setattr(obj, attr_name, value)
        for name, value in overrides.items():
            setattr(obj, name, value)
        obj.__post_init__()
        return obj

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
            "trackRayVisits": self.trackRayVisits,
            "reflectionMaxIterations": self.reflectionMaxIterations,
            "reflectionTolerance": self.reflectionTolerance,
            "reflectionMode": self.reflectionMode,
            "surfaceReservoirSize": int(self.surfaceReservoirSize),
            "srmPositionMode": self.srmPositionMode,
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
            _validate_launch_backends(self, openpmd_backend=kwargs.get("transport", self.openpmdBackend))
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

    def run(
        self,
        gainMedium=None,
        *,
        opticalComponents=None,
        initialExcitation=None,
        openpmdSession=None,
    ):
        """Run ASE for the supplied or configured ``GainMedium``.

        ``opticalComponents`` may include passive components traversed by the
        rays; it defaults to the active components in ``gainMedium``.
        Returns ``self``. Use ``getResults()`` afterwards to access the raw
        lower-level result, including ``phiAse``.
        """
        medium = gainMedium if gainMedium is not None else self.gainMedium
        if medium is None:
            raise ValueError("PhiASE.run requires gainMedium; pass it through Simulation or run(gainMedium=...)")
        if not isinstance(medium, GainMedium):
            raise TypeError("PhiASE.run requires the physical GainMedium frontend primitive")
        if initialExcitation is None:
            initialExcitation = 0.0

        if openpmdSession == "persistent":
            openpmdSession = self.openStream()
        elif openpmdSession == "interval":
            openpmdSession = None

        _validate_launch_backends(self, openpmd_session=openpmdSession)

        launch_options = {} if openpmdSession is not None else self._transportLaunchOptions()
        self._result = transport.runPhiASE(
            _AseSimulationRequest(self, medium, initialExcitation, opticalComponents),
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


@dataclass(frozen=True)
class ExcitationState:
    """Excited-state values associated with physical volume domains."""

    domains: tuple[Domain, ...]
    values: tuple[object, ...]

    @classmethod
    def fromValue(cls, value, gainMedium):
        if isinstance(value, Mapping):
            domains = tuple(value)
            values = tuple(value[domain] for domain in domains)
        else:
            domains = (gainMedium.domain,)
            values = (value,)
        parsed = []
        for item in values:
            array = np.asarray(item, dtype=np.float64)
            if np.any(~np.isfinite(array)) or np.any((array < 0.0) | (array > 1.0)):
                raise ValueError("initial excitation must be finite and within [0, 1]")
            parsed.append(array.reshape(-1) if array.ndim else np.asarray([float(array)]))
        return cls(domains=domains, values=tuple(parsed))

    def _transportDescription(self):
        return PrimitiveDescription(
            "excitationState",
            fields=(transportField("values", axes=("cell",), dynamic=True, encoding="ragged"),),
            references=(reference("domains", many=True),),
        )


@dataclass(frozen=True)
class PumpRegistration:
    """One physical pump, its injection method, and optional relay chain."""

    pump: Pump
    injectionMethod: SurfacePumpInjector
    relays: tuple[PlanarPumpRelay, ...] = ()

    def _transportDescription(self):
        return PrimitiveDescription(
            "pumpRegistration",
            references=(
                reference("pump"),
                reference("injectionMethod"),
                reference("relays", many=True),
            ),
        )


def _simulationTransportDescription():
    return PrimitiveDescription(
        "simulation",
        fields=(
            transportField("timeStep"),
            transportField(
                "simulationSteps",
                lambda owner: getattr(owner, "_requestedSteps", owner.simulationSteps),
                optional=True,
            ),
            transportField("endTime", optional=True),
            transportField("prePump"),
            transportField("executionMode"),
            transportField("outputSteps", optional=True),
            transportField("outputFields"),
            transportField("controlFields"),
            transportField("currentStep", dynamic=True),
            transportField("currentTime", dynamic=True),
        ),
        references=(
            reference("opticalComponents", many=True),
            reference("gainMedium"),
            reference("exteriorSurface"),
            reference("excitationState"),
            reference("phiAse", "phiASE"),
            reference("timeIntegrationSolver"),
            reference("pumpRegistrations", "_pumpRegistrations", many=True),
        ),
    )


class _AseSimulationRequest:
    """Internal projection of direct ``PhiASE.run`` onto the Simulation graph."""

    def __init__(self, phiAse, gainMedium, initialExcitation, opticalComponents=None):
        if not gainMedium.components:
            raise ValueError("PhiASE.run requires a non-empty GainMedium")
        self.opticalComponents = tuple(
            gainMedium.components if opticalComponents is None else opticalComponents
        )
        if any(component not in self.opticalComponents for component in gainMedium.components):
            raise ValueError("every GainMedium component must belong to the direct ASE optical assembly")
        validateComponentOverlap(self.opticalComponents)
        self.gainMedium = gainMedium
        occupiedDomain = Domain(entityKind="volume")
        for component in self.opticalComponents:
            occupiedDomain += component.domain
        self.exteriorSurface = occupiedDomain.boundary()
        self.excitationState = ExcitationState.fromValue(initialExcitation, gainMedium)
        self.phiASE = phiAse
        self.timeIntegrationSolver = TimeIntegrationSolver("explicit-euler")
        self._pumpRegistrations = ()
        self.timeStep = 1.0
        self.simulationSteps = 1
        self.endTime = None
        self.prePump = False
        self.executionMode = "autonomous"
        self.outputSteps = None
        self.outputFields = (
            "phi_ase",
            "standard_error",
            "relative_standard_error",
            "total_rays",
            "dndt_ase",
        )
        self.controlFields = ()
        self.currentStep = 0
        self.currentTime = 0.0

    def _transportDescription(self):
        return _simulationTransportDescription()



@dataclass
class TimeStepState:
    """Snapshot handed to ``onStep`` callbacks after a completed time step.

    The arrays are copies of the simulation outputs at ``step``/``time`` and
    contain exactly one value per Tet4 cell. The compiled simulation has no
    point-centered state representation.
    """

    step: int
    """Completed one-based step index."""
    time: float
    """Physical simulation time after the step, in seconds."""
    betaVolume: np.ndarray | None
    """Authoritative excited-state fraction for every Tet4 cell, when selected."""
    phiAse: np.ndarray | None
    """ASE flux for every Tet4 cell, when selected."""
    dndtAse: np.ndarray | None
    """Cell-centered ASE depletion contribution to ``d beta / dt``, when selected."""
    dndtPump: np.ndarray | None
    """Cell-centered pump contribution to ``d beta / dt``, when selected."""
    aseResult: object | None
    """Raw lower-level ASE result object for advanced inspection."""
    standardError: np.ndarray | None = None
    """Absolute one-sigma error of the cell ASE estimate, when selected."""
    relativeStandardError: np.ndarray | None = None
    """Relative one-sigma error of the cell ASE estimate, when selected."""
    totalRays: np.ndarray | None = None
    """Ray visit count for each cell, when selected."""
    topology: object | None = None
    """Static mesh topology used by geometry-aware state callbacks."""

@dataclass(init=False)
class Simulation:
    """High-level Python wrapper for compiled C++/Alpaka simulation runs.

    Python sends the initial setup to the compiled backend and receives
    ``TimeStepState`` snapshots after completed steps. ``onInit`` prepares the
    initial state before step 1, ``onStep`` consumes selected completed-step
    snapshots, and ``beforeStep`` updates selected control fields between two
    steps. Autonomous runs never call Python between backend steps.
    Synchronized-debug runs call ``beforeStep`` after each nonfinal
    ``onStep`` callback and before allowing the next backend step to begin.

    Parameters
    ----------
    opticalComponents
        Complete set of disjoint optical volume components traversed by rays.
        Their domains provide all executable geometry; ``Simulation`` owns no
        separate topology collection or aggregate optical domain.
    gainMedium
        Active subset of ``opticalComponents`` participating in population and
        stimulated-emission dynamics.
    exteriorSurface
        Optional user-declared surface :class:`Domain`. When omitted or
        ``None``, it is inferred as the boundary of the union of all component
        domains.
    initialExcitation
        Scalar active-region excitation or an exact mapping from volume domains
        to scalar/array values. Passive cells are initialized to zero.
    """

    gainMedium: GainMedium
    exteriorSurface: Domain
    """Exterior surface, inferred from all component domains unless supplied."""
    pump: _PumpProperties | None
    phiASE: PhiASE
    timeIntegrationSolver: TimeIntegrationSolver | str
    timeStep: float
    endTime: float | None
    simulationSteps: int | None
    prePump: bool
    reportTimings: bool
    executionMode: str
    outputSteps: tuple[int, ...] | None
    outputFields: tuple[str, ...]
    controlFields: tuple[str, ...]
    _pumpRegistrations: list
    _time: float
    _step: int
    _initialized: bool
    _init_callbacks: list
    _before_step_callbacks: list
    _step_callbacks: list
    _lastState: TimeStepState | None

    def __init__(
        self,
        *,
        opticalComponents,
        gainMedium,
        phiASE,
        timeIntegrator,
        timeStepSize,
        initialExcitation=0.0,
        exteriorSurface=None,
        simulationSteps=None,
        maxTime=None,
        prePump=False,
        reportTimings=False,
        executionMode="autonomous",
        outputSteps=None,
        outputFields=None,
        controlFields=(),
    ):
        if not isinstance(gainMedium, GainMedium):
            raise TypeError("gainMedium must be GainMedium")
        self.opticalComponents = tuple(opticalComponents)
        if not all(
            isinstance(component, OpticalComponent)
            for component in self.opticalComponents
        ):
            raise TypeError("opticalComponents must contain OpticalComponent values")
        if any(component not in self.opticalComponents for component in gainMedium.components):
            raise ValueError("every GainMedium component must be owned by Simulation")
        validateComponentOverlap(self.opticalComponents)
        if exteriorSurface is None:
            occupiedDomain = Domain(entityKind="volume")
            for component in self.opticalComponents:
                occupiedDomain += component.domain
            selectedSurface = occupiedDomain.boundary()
        else:
            selectedSurface = (
                exteriorSurface
                if isinstance(exteriorSurface, Domain)
                else Domain(exteriorSurface)
            )
        if selectedSurface.entityKind != "surface" or selectedSurface.isEmpty:
            raise ValueError("exteriorSurface must resolve to a non-empty surface Domain")
        self.exteriorSurface = selectedSurface
        self.gainMedium = gainMedium
        self.initialExcitation = initialExcitation
        self.pump = None
        self.phiASE = phiASE
        self.timeIntegrationSolver = (
            TimeIntegrationSolver(str(timeIntegrator))
            if isinstance(timeIntegrator, str)
            else timeIntegrator
        )
        self.timeStep = float(timeStepSize)
        self.endTime = maxTime
        self.simulationSteps = None if simulationSteps is None else int(simulationSteps)
        self.prePump = bool(prePump)
        self.reportTimings = bool(reportTimings)
        self.executionMode = str(executionMode)
        self.outputSteps = None if outputSteps is None else tuple(int(step) for step in outputSteps)
        self.outputFields = tuple(
            SIMULATION_OUTPUT_FIELDS if outputFields is None else (str(field) for field in outputFields)
        )
        self.controlFields = tuple(str(field) for field in controlFields)
        self._pumpRegistrations = []
        self._projectedPumpRegistrations = []
        self._time = 0.0
        self._step = 0
        self._initialized = False
        self._init_callbacks = []
        self._before_step_callbacks = []
        self._step_callbacks = []
        self._lastState = None
        self.__post_init__()

    @classmethod
    def fromYaml(cls, filename, **objects):
        """Construct a schema-v3 simulation, resolving missing named objects."""
        from .configuration import simulationFromYaml

        return simulationFromYaml(filename, simulationCls=cls, **objects)

    def __post_init__(self):
        if self.timeIntegrationSolver is None:
            raise ValueError("Simulation requires a time_integrator")
        if not isinstance(self.timeIntegrationSolver, str) and not hasattr(self.timeIntegrationSolver, "name"):
            raise ValueError("Simulation requires a compiled time integrator name or descriptor with a .name attribute")
        if self.timeStep <= 0.0:
            raise ValueError("time_step_size must be positive")
        if self.simulationSteps is not None and self.simulationSteps <= 0:
            raise ValueError("simulation_steps must be positive")
        if self.simulationSteps is not None and self.endTime is not None:
            raise ValueError("configure either simulation_steps or max_time, not both")
        if self.executionMode not in {"autonomous", "synchronized-debug"}:
            raise ValueError("execution_mode must be 'autonomous' or 'synchronized-debug'")
        if self.outputSteps is not None and not self.outputSteps:
            raise ValueError("output_steps must be omitted or contain at least one completed-step index")
        if self.outputSteps is not None and any(step <= 0 for step in self.outputSteps):
            raise ValueError("output_steps must contain positive completed-step indices")
        if self.outputSteps is not None and tuple(sorted(set(self.outputSteps))) != self.outputSteps:
            raise ValueError("output_steps must be strictly increasing and unique")
        if not self.outputFields:
            raise ValueError("output_fields must contain at least one field")
        unknown_output_fields = sorted(set(self.outputFields) - set(SIMULATION_OUTPUT_FIELDS))
        if unknown_output_fields:
            raise ValueError(f"unsupported output_fields: {unknown_output_fields}")
        if len(set(self.outputFields)) != len(self.outputFields):
            raise ValueError("output_fields must be unique")
        if self.executionMode == "autonomous" and self.controlFields:
            raise ValueError("control_fields require execution_mode='synchronized-debug'")
        if self.executionMode == "synchronized-debug" and self.outputSteps is not None:
            raise ValueError("synchronized-debug emits every completed step; output_steps must be omitted")
        unknown_control_fields = sorted(set(self.controlFields) - set(SIMULATION_CONTROL_FIELDS))
        if unknown_control_fields:
            raise ValueError(f"unsupported control_fields: {unknown_control_fields}")
        if len(set(self.controlFields)) != len(self.controlFields):
            raise ValueError("control_fields must be unique")
        self._simulationState = projectFrontendState(
            self.gainMedium,
            self.initialExcitation,
            opticalComponents=self.opticalComponents,
        )
        self.phiASE.gainMedium = self.gainMedium
        self._ensureStateArrays()
        self._refreshExcitationState()

    def _transportDescription(self):
        return _simulationTransportDescription()

    def addPump(self, pump, injectionMethod, *, relays=()):
        """Register a physical pump and its numerical injection method."""
        if self._initialized:
            raise RuntimeError("pumps must be added before the simulation is initialized")
        if not isinstance(pump, Pump):
            raise TypeError("pump must be a Pump")
        if not isinstance(injectionMethod, SurfacePumpInjector):
            raise TypeError("injectionMethod must be SurfacePumpInjector")
        relays = tuple(relays)
        if not all(isinstance(relay, PlanarPumpRelay) for relay in relays):
            raise TypeError("relays must contain PlanarPumpRelay values")
        def lowerDomains(values):
            return tuple(
                projectSurfaceDomain(self._simulationState, value)
                if hasattr(value, "entityKind")
                else value
                for value in values
            )

        projectedInjection = SurfacePumpInjector(lowerDomains(injectionMethod.surface_domains))
        projectedRelays = tuple(
            PlanarPumpRelay(
                lowerDomains(relay.exit_domains),
                lowerDomains(relay.entry_domains),
                flip_u=relay.flip_u,
                flip_v=relay.flip_v,
                rotation=relay.rotation,
                offset=relay.offset,
                tilt=relay.tilt,
                magnification=relay.magnification,
                transmission=relay.transmission,
            )
            for relay in relays
        )
        self._pumpRegistrations.append(PumpRegistration(pump, injectionMethod, relays))
        self._projectedPumpRegistrations.append((pump, projectedInjection, projectedRelays))
        self.pump = _PumpProperties(
            sources=tuple(
                _PumpSource(
                    surfaceDomains=injector.surface_domains,
                    totalPower=physical.total_power,
                    spectrum=physical.spectrum,
                    angularDistribution=physical.angular_distribution,
                    profile=physical.profile,
                    relays=registered_relays,
                    rayCount=int(physical.ray_count),
                    pumpSteps=0 if physical.pump_steps is None else int(physical.pump_steps),
                    rngSeed=int(physical.rng_seed),
                )
                for physical, injector, registered_relays in self._projectedPumpRegistrations
            ),
        )
        return self

    def cellMask(self, domain):
        """Map a physical volume domain to the callback state's cell order."""
        return projectCellMask(self._simulationState, domain)

    def onStep(self, callback, *args, **kwargs):
        """Register a post-snapshot callback.

        The callback signature is ``callback(state, *args, **kwargs)``.
        ``Simulation`` supplies the completed ``TimeStepState`` first, followed
        by the arguments passed during registration. With ``output_steps``
        omitted, the callback runs after every completed step. An explicit
        ``output_steps`` schedule restricts it to those snapshots.

        State arrays are copies and changing them does not update the backend.
        With a streaming openPMD backend, callbacks run as snapshots arrive;
        with a file backend, they run after the compiled process finishes. Use
        this hook for logging, output, or explicit state retention. Callback
        return values are ignored. The method returns ``self`` so registrations
        can be chained.
        """
        self._step_callbacks.append((callback, args, kwargs))
        return self

    def onInit(self, callback, *args, **kwargs):
        """Register a one-time initialization callback.

        The callback signature is ``callback(simulation, *args, **kwargs)``.
        ``Simulation`` supplies the live simulation object as the first
        argument, then appends the user arguments passed to ``onInit``. The
        hook runs once, immediately before the first compiled run. Use it to
        modify the initial state consumed by step 1; it does not run again for
        later calls to ``step`` on the same object.

        Callback return values are ignored. The method returns ``self`` for
        chaining.
        """
        self._init_callbacks.append((callback, args, kwargs))
        return self

    def beforeStep(self, callback, *args, **kwargs):
        """Register a synchronized-debug callback that runs between steps.

        The callback signature is ``callback(simulation, *args, **kwargs)``.
        It first runs after step 1: Python has called ``onStep`` for the
        completed snapshot, then calls ``beforeStep`` to prepare selected
        ``control_fields`` before step 2 starts. It repeats after every
        nonfinal step and is not called after the final step. Use ``onInit``
        to modify state before step 1.

        This hook requires ``execution_mode="synchronized-debug"`` and a
        streaming openPMD backend. ``beta_volume`` and material-owned
        ``cross_sections`` are the supported control fields. Synchronized-debug
        emits every step and does
        not accept ``output_steps``. Callback return values are ignored. The
        method returns ``self`` for chaining.
        """
        self._before_step_callbacks.append((callback, args, kwargs))
        return self

    def runUntil(self, maxTime=None):
        """Advance to ``max_time`` or the constructor's configured maximum."""
        target = self.endTime if maxTime is None else maxTime
        if target is None:
            raise ValueError("runUntil requires maxTime or a configured maxTime")
        steps = 0
        while self._time + steps * self.timeStep < float(target) - 0.5 * self.timeStep:
            steps += 1
        if steps:
            self.step(steps)
        return self

    def runSteps(self, steps, *, openpmdSession=None):
        """Run exactly ``steps`` compiled C++/Alpaka time steps and return ``self``.

        Internal transport helper. Each registered pump and ``PhiASE`` owns its
        activity window. The complete time loop is executed by the C++ backend;
        Python only sends the initial setup and receives streamed snapshots.
        """
        if openpmdSession not in (None, "interval"):
            raise ValueError("compiled Simulation owns its C++ openPMD lifetime; openpmdSession is no longer supported")
        steps = int(steps)
        if steps <= 0:
            raise ValueError("steps must be positive")
        if self._before_step_callbacks and self.executionMode != "synchronized-debug":
            raise ValueError(
                "beforeStep callbacks require executionMode='synchronized-debug'; "
                "autonomous runs do not contact Python between steps"
            )
        if self.pump is None:
            raise ValueError("Simulation requires at least one pump registered with addPump")
        self._run_init_callbacks()
        _validate_launch_backends(self.phiASE)

        previous_step = self._step
        previous_time = self._time
        run_started = perf_counter() if self.reportTimings else None
        transport_started = perf_counter() if self.reportTimings else None
        state_materialization_seconds = 0.0
        callback_seconds = {}
        received_states = []

        def consume_raw_state(raw_state):
            nonlocal state_materialization_seconds
            state_started = perf_counter() if self.reportTimings else None
            state = TimeStepState(
                step=previous_step + int(raw_state.step),
                time=previous_time + float(raw_state.time),
                betaVolume=_optional_state_array(raw_state.betaVolume, np.float64),
                phiAse=_optional_state_array(raw_state.phiAse, np.float64),
                standardError=_optional_state_array(getattr(raw_state, "standardError", None), np.float64),
                relativeStandardError=_optional_state_array(
                    getattr(raw_state, "relativeStandardError", None), np.float64
                ),
                totalRays=_optional_state_array(getattr(raw_state, "totalRays", None), np.uint32),
                dndtAse=_optional_state_array(raw_state.dndtAse, np.float64),
                dndtPump=_optional_state_array(raw_state.dndtPump, np.float64),
                aseResult=raw_state.aseResult,
                topology=self._simulationState.topology,
            )
            explicit_topology = hasattr(self._simulationState.topology, "cellPointIndices")
            if state.betaVolume is not None:
                self._simulationState.get("betaVolume").value = (
                    backendFlat(state.betaVolume.reshape(-1, order="F"))
                    if explicit_topology
                    else state.betaVolume
                )
            self._lastState = state
            self._step = state.step
            self._time = state.time
            received_states.append(state)
            if self.reportTimings:
                state_materialization_seconds += perf_counter() - state_started
            for callback, args, kwargs in self._step_callbacks:
                callback_started = perf_counter() if self.reportTimings else None
                callback(state, *args, **kwargs)
                if self.reportTimings:
                    callback_name = getattr(callback, "__qualname__", getattr(callback, "__name__", type(callback).__name__))
                    callback_seconds[callback_name] = callback_seconds.get(callback_name, 0.0) + (
                        perf_counter() - callback_started
                    )
            if self.executionMode == "synchronized-debug" and int(raw_state.step) < steps:
                for callback, args, kwargs in self._before_step_callbacks:
                    callback(self, *args, **kwargs)

        self._requestedSteps = steps
        try:
            states = transport.runSimulation(
                self,
                steps=steps,
                transport=self.phiASE.openpmdBackend,
                on_state=consume_raw_state,
                **self.phiASE._transportLaunchOptions(),
            )
        finally:
            del self._requestedSteps
        transport_seconds = perf_counter() - transport_started if self.reportTimings else 0.0
        if not received_states:
            for raw_state in states:
                consume_raw_state(raw_state)
        self._step = previous_step + steps
        self._time = previous_time + steps * self.timeStep
        if self.reportTimings:
            total_seconds = perf_counter() - run_started
            callbacks_total = sum(callback_seconds.values())
            print(
                "HASE frontend timing: "
                f"snapshots={len(received_states)} total={total_seconds:.6f}s "
                f"compiled_transport_and_decode={transport_seconds:.6f}s "
                f"snapshot_materialization={state_materialization_seconds:.6f}s "
                f"callbacks={callbacks_total:.6f}s"
            )
            for callback_name, seconds in sorted(callback_seconds.items(), key=lambda item: item[1], reverse=True):
                print(f"HASE frontend callback timing: {callback_name}={seconds:.6f}s")
        return self

    def _derived_simulation_steps(self):
        activity_steps = [0 if self.phiASE.ase_steps is None else int(self.phiASE.ase_steps)]
        activity_steps.extend(
            0 if registration.pump.pump_steps is None else int(registration.pump.pump_steps)
            for registration in self._pumpRegistrations
        )
        derived = max(activity_steps, default=0)
        if derived <= 0:
            raise ValueError(
                "simulation_steps is required when all pump_steps and ase_steps are disabled"
            )
        return derived

    def step(self, nsteps=None):
        """Advance ``nsteps`` time steps, following the PICMI call pattern."""
        if nsteps is None:
            if self.simulationSteps is not None:
                nsteps = self.simulationSteps
            elif self.endTime is not None:
                return self.runUntil()
            else:
                nsteps = self._derived_simulation_steps()
        if int(nsteps) <= 0:
            raise ValueError("nsteps must be positive")
        self.runSteps(int(nsteps))
        return self

    @property
    def currentStep(self):
        return self._step

    @property
    def currentTime(self):
        return self._time

    @property
    def pumps(self):
        return tuple(registration.pump for registration in self._pumpRegistrations)

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
        if "betaVolume" not in self._simulationState.physical:
            self._simulationState.get("betaVolume").value = np.zeros(
                self._simulationState.get("betaVolume").expectedShape,
                dtype=np.float64,
            )

    def _refreshExcitationState(self):
        values = np.asarray(
            self._simulationState.get("betaVolume").value,
            dtype=np.float64,
        ).reshape(-1, order="F")
        domains = []
        domainValues = []
        cellMaps = self._simulationState.__dict__["_domainCellMaps"]
        for component in self.gainMedium.components:
            selectedValues = []
            for topology, mask in component.domain._shards:
                owner, mapping = cellMaps[id(topology)]
                if owner is not topology:
                    raise RuntimeError("gain-medium topology identity changed after state projection")
                indices = mapping[mask]
                if np.any(indices < 0):
                    raise ValueError("gain-medium Domain contains cells outside the executable assembly")
                selectedValues.extend(values[indices].tolist())
            domains.append(component.domain)
            domainValues.append(np.asarray(selectedValues, dtype=np.float64))
        existing = getattr(self, "excitationState", None)
        if existing is None:
            self.excitationState = ExcitationState(tuple(domains), tuple(domainValues))
        else:
            if existing.domains != tuple(domains):
                raise RuntimeError("excitation domains cannot change after transport composition")
            object.__setattr__(existing, "values", tuple(domainValues))

    def _run_init_callbacks(self):
        if self._initialized:
            return
        self._initialized = True
        for callback, args, kwargs in self._init_callbacks:
            callback(self, *args, **kwargs)
