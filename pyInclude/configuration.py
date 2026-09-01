# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Schema-v3 construction of the composable physical object graph."""

from __future__ import annotations

from pathlib import Path

from hase_units import units
from material_library import Material
from .geometry import SurfaceOptics, VolumeTopology
from .laser import (
    PlanarPumpRelay,
    Pump,
    PumpAngularDistribution,
    PumpSpectrum,
    SuperGaussianPumpProfile,
    SurfacePumpInjector,
    UniformPumpProfile,
)
from .physical import Domain, GainMedium, OpticalComponent
from .simulation import PhiASE
from .timeIntegration import (
    ExplicitEuler,
    ExponentialEuler,
    FrozenPhiAseRungeKutta4,
    FrozenSourcesRungeKutta4,
    Heun,
    ImplicitEuler,
    Midpoint,
    RungeKutta4,
)


_REGISTRY_KEYS = {
    Material: "materials",
    Domain: "domains",
    OpticalComponent: "optical_components",
    GainMedium: "gain_media",
}


def _mapping(value, path):
    if not isinstance(value, dict):
        raise ValueError(f"{path} must be a mapping")
    return dict(value)


def _rejectUnknown(mapping, allowed, path):
    unknown = sorted(set(mapping) - set(allowed))
    if unknown:
        raise ValueError(f"unsupported {path} options: {unknown}")


def _oneOf(mapping, names, path):
    selected = [name for name in names if name in mapping]
    if len(selected) != 1:
        raise ValueError(f"{path} requires exactly one of: {', '.join(names)}")
    return selected[0]


def _loadYaml(filename):
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("schema-v3 YAML requires PyYAML") from exc
    path = Path(filename).expanduser().resolve()
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"configuration {path} must contain a mapping")
    if data.get("schema_version") != 3:
        raise ValueError("configuration requires schema_version: 3")
    allowed = {
        "schema_version",
        "materials",
        "topologies",
        "domains",
        "optical_components",
        "gain_media",
        "simulation",
    }
    _rejectUnknown(data, allowed, "top-level")
    return path, data


class _YamlContext:
    def __init__(self, filename, **injected):
        self.path, self.data = _loadYaml(filename)
        self.root = self.path.parent
        self.cache = {
            "materials": dict(injected.get("materials", {})),
            "topologies": dict(injected.get("topologies", {})),
            "domains": dict(injected.get("domains", {})),
            "optical_components": dict(
                injected.get("opticalComponents", injected.get("optical_components", {}))
            ),
            "gain_media": dict(injected.get("gainMedia", injected.get("gain_media", {}))),
        }
        unknown = sorted(
            set(injected)
            - {"materials", "topologies", "domains", "opticalComponents", "optical_components", "gainMedia", "gain_media"}
        )
        if unknown:
            raise TypeError(f"unknown injected object registries: {unknown}")
        self._resolving = []

    def resolve(self, registry, name):
        name = str(name)
        if name in self.cache[registry]:
            return self.cache[registry][name]
        marker = (registry, name)
        if marker in self._resolving:
            chain = " -> ".join(f"{kind}.{item}" for kind, item in (*self._resolving, marker))
            raise ValueError(f"cyclic YAML reference: {chain}")
        section = _mapping(self.data.get(registry, {}), registry)
        if name not in section:
            raise ValueError(f"unknown {registry} reference: {name}")
        self._resolving.append(marker)
        try:
            builders = {
                "materials": self._material,
                "topologies": self._topology,
                "domains": self._domain,
                "optical_components": self._component,
                "gain_media": self._gainMedium,
            }
            result = builders[registry](name, section[name])
            self.cache[registry][name] = result
            return result
        finally:
            self._resolving.pop()

    def _path(self, value):
        path = Path(value).expanduser()
        return path if path.is_absolute() else self.root / path

    def _material(self, registryName, value):
        path = f"materials.{registryName}"
        spec = _mapping(value, path)
        allowed = {
            "from_hdf5",
            "temperature",
            "active_ion_density",
            "interpolation",
            "spectral_resolution",
            "name",
            "optical_axis",
            "refractive_index",
            "fluorescence_lifetime",
            "bulk_attenuation",
            "absorption_coefficient",
        }
        _rejectUnknown(spec, allowed, path)
        source = _mapping(spec.get("from_hdf5"), f"{path}.from_hdf5")
        _rejectUnknown(source, {"path", "key"}, f"{path}.from_hdf5")
        material = Material.fromHdf5(
            self._path(source["path"]),
            key=source.get("key"),
            temperature=(None if "temperature" not in spec else float(spec["temperature"]) * units.K),
            activeIonDensity=float(spec.get("active_ion_density", 0.0)) / units.m**3,
            spectralResolution=spec.get("spectral_resolution"),
            interpolation=spec.get("interpolation", "exact"),
            name=spec.get("name", registryName),
            opticalAxis=spec.get("optical_axis"),
        )
        if "refractive_index" in spec:
            material.refractiveIndex = float(spec["refractive_index"])
        if "fluorescence_lifetime" in spec:
            material.fluorescenceLifetime = float(spec["fluorescence_lifetime"]) * units.s
        attenuationNames = {
            name for name in ("bulk_attenuation", "absorption_coefficient") if name in spec
        }
        if len(attenuationNames) > 1:
            raise ValueError(
                f"{path} provides both bulk_attenuation and absorption_coefficient"
            )
        if attenuationNames:
            name = attenuationNames.pop()
            material.bulkAttenuation = float(spec[name]) / units.m
        return material.validate()

    def _topologyValue(self, value, path):
        if isinstance(value, str):
            return self.resolve("topologies", value)
        spec = _mapping(value, path)
        form = _oneOf(spec, ("from_file", "from_tetrahedra"), path)
        _rejectUnknown(spec, {form}, path)
        values = _mapping(spec[form], f"{path}.{form}")
        if form == "from_file":
            _rejectUnknown(values, {"path", "format", "boundary_default", "mesh_size"}, f"{path}.{form}")
            kwargs = {}
            if values.get("boundary_default") is not None:
                from .geometry.volume import BOUND_INTERNAL, BOUND_STOP

                kwargs["boundaryDefault"] = {
                    "stop": BOUND_STOP,
                    "internal": BOUND_INTERNAL,
                }.get(values["boundary_default"], values["boundary_default"])
            if values.get("mesh_size") is not None:
                kwargs["meshSize"] = values["mesh_size"]
            meshFormat = values.get("format")
            return VolumeTopology.fromFile(
                self._path(values["path"]),
                format=None if meshFormat == "auto" else meshFormat,
                **kwargs,
            )
        _rejectUnknown(
            values,
            {"points", "cell_point_indices", "cell_domains", "face_boundaries", "metadata"},
            f"{path}.{form}",
        )
        return VolumeTopology.fromTetrahedra(
            values["points"],
            values["cell_point_indices"],
            cellDomains=values.get("cell_domains"),
            faceBoundaries=values.get("face_boundaries"),
            metadata=values.get("metadata"),
        )

    def _topology(self, name, value):
        return self._topologyValue(value, f"topologies.{name}")

    def _domainOperand(self, value, path):
        if isinstance(value, str):
            return self.resolve("domains", value)
        return self._domain(path, value, anonymous=True)

    def _domain(self, name, value, anonymous=False):
        path = name if anonymous else f"domains.{name}"
        spec = _mapping(value, path)
        forms = (
            "from_gmsh",
            "where",
            "topology",
            "component",
            "exterior_cells",
            "exterior_tets",
            "union",
            "difference",
            "boundary",
        )
        form = _oneOf(spec, forms, path)
        _rejectUnknown(spec, {form}, path)
        data = spec[form]
        if form == "from_gmsh":
            data = _mapping(data, f"{path}.from_gmsh")
            _rejectUnknown(data, {"topology", "physical_group", "entity_kind"}, f"{path}.from_gmsh")
            return Domain.fromGmsh(
                self.resolve("topologies", data["topology"]),
                data["physical_group"],
                entityKind=data.get("entity_kind"),
            )
        if form == "where":
            data = _mapping(data, f"{path}.where")
            _rejectUnknown(data, {"topology", "selector", "entity_kind"}, f"{path}.where")
            return Domain.where(
                self.resolve("topologies", data["topology"]),
                data["selector"],
                entityKind=data.get("entity_kind", "surface"),
            )
        if form == "topology":
            if isinstance(data, str):
                topology, entityKind = self.resolve("topologies", data), "volume"
            else:
                data = _mapping(data, f"{path}.topology")
                _rejectUnknown(data, {"name", "entity_kind"}, f"{path}.topology")
                topology, entityKind = self.resolve("topologies", data["name"]), data.get("entity_kind", "volume")
            return Domain.fromTopology(topology, entityKind=entityKind)
        if form == "component":
            return self.resolve("optical_components", data).domain
        if form in {"exterior_cells", "exterior_tets"}:
            return self.resolve("optical_components", data).exteriorCells
        if form == "boundary":
            return self._domainOperand(data, f"{path}.boundary").boundary()
        if form == "union":
            if not isinstance(data, list) or not data:
                raise ValueError(f"{path}.union must be a non-empty sequence")
            return Domain(self._domainOperand(item, f"{path}.union") for item in data)
        if not isinstance(data, list) or len(data) != 2:
            raise ValueError(f"{path}.difference must contain exactly two operands")
        return self._domainOperand(data[0], f"{path}.difference") - self._domainOperand(
            data[1], f"{path}.difference"
        )

    def _component(self, registryName, value):
        path = f"optical_components.{registryName}"
        spec = _mapping(value, path)
        allowed = {"material", "domain", "name", "optical_role", "surface_optics"}
        _rejectUnknown(spec, allowed, path)
        if "domain" not in spec:
            raise ValueError(f"{path} requires domain")
        kwargs = {
            "material": self.resolve("materials", spec["material"]),
            "name": spec.get("name", registryName),
            "opticalRole": spec.get("optical_role"),
            "domain": self.resolve("domains", spec["domain"]),
        }
        component = OpticalComponent(**kwargs)
        optics = spec.get("surface_optics", [])
        if not isinstance(optics, list):
            raise ValueError(f"{path}.surface_optics must be a sequence")
        for index, assignment in enumerate(optics):
            assignment = _mapping(assignment, f"{path}.surface_optics[{index}]")
            _rejectUnknown(
                assignment,
                {
                    "domain",
                    "reflectivity",
                    "interior_refractive_index",
                    "exterior_refractive_index",
                },
                f"{path}.surface_optics[{index}]",
            )
            component.assignSurfaceOptics(
                self.resolve("domains", assignment["domain"]),
                SurfaceOptics(
                    reflectivity=float(assignment.get("reflectivity", 0.0)),
                    n_inside=float(
                        assignment.get("interior_refractive_index", component.material.refractiveIndex)
                    ),
                    n_outside=float(assignment.get("exterior_refractive_index", 1.0)),
                ),
            )
        return component

    def _gainMedium(self, registryName, value):
        path = f"gain_media.{registryName}"
        spec = _mapping(value, path)
        _rejectUnknown(spec, {"components", "name"}, path)
        components = spec.get("components")
        if not isinstance(components, list) or not components:
            raise ValueError(f"{path}.components must be a non-empty sequence")
        medium = GainMedium(
            [self.resolve("optical_components", name) for name in components],
            name=spec.get("name", registryName),
        )
        return medium


def _phiAse(spec):
    spec = _mapping(spec, "simulation.phi_ase")
    aliases = {
        "propagation_mode": "propagationMode",
        "min_rays": "minRays",
        "max_rays": "maxRays",
        "forward_ray_count": "forwardRayCount",
        "relative_standard_error_threshold": "relativeStandardErrorThreshold",
        "track_ray_visits": "trackRayVisits",
        "adaptive_steps": "adaptiveSteps",
        "use_reflections": "useReflections",
        "reflection_mode": "reflectionMode",
        "surface_reservoir_size": "surfaceReservoirSize",
        "srm_position_mode": "srmPositionMode",
        "reflection_max_iterations": "reflectionMaxIterations",
        "reflection_tolerance": "reflectionTolerance",
        "rng_seed": "rngSeed",
        "openpmd_backend": "openpmdBackend",
        "parallel_mode": "parallelMode",
        "num_devices": "numDevices",
        "n_per_node": "nPerNode",
        "min_sample_range": "minSampleRange",
        "max_sample_range": "maxSampleRange",
    }
    unchanged = {"repetitions", "monochromatic", "backend", "ase_steps"}
    _rejectUnknown(spec, set(aliases) | unchanged, "simulation.phi_ase")
    values = {aliases.get(name, name): value for name, value in spec.items()}
    return PhiASE(**values)


def _pumpSpectrum(spec):
    spec = _mapping(spec, "pump.spectrum")
    if "monochromatic" in spec:
        _rejectUnknown(spec, {"monochromatic"}, "pump.spectrum")
        return PumpSpectrum.monochromatic(spec["monochromatic"])
    _rejectUnknown(spec, {"wavelengths", "weights"}, "pump.spectrum")
    return PumpSpectrum(spec["wavelengths"], spec["weights"])


def _angularDistribution(spec):
    spec = _mapping(spec, "pump.angular_distribution")
    form = _oneOf(spec, ("collimated", "uniform_cone", "discrete"), "pump.angular_distribution")
    _rejectUnknown(spec, {form}, "pump.angular_distribution")
    values = _mapping(spec[form], f"pump.angular_distribution.{form}")
    if form == "collimated":
        _rejectUnknown(values, set(), "pump.angular_distribution.collimated")
        return PumpAngularDistribution.collimated()
    if form == "uniform_cone":
        _rejectUnknown(values, {"half_angle", "polar_samples", "azimuthal_samples"}, "pump.angular_distribution.uniform_cone")
        return PumpAngularDistribution.uniformCone(
            values["half_angle"],
            polarSamples=values.get("polar_samples", 8),
            azimuthalSamples=values.get("azimuthal_samples", 16),
        )
    _rejectUnknown(values, {"polar_angles", "azimuthal_angles", "weights"}, "pump.angular_distribution.discrete")
    return PumpAngularDistribution(**values)


def _profile(spec):
    spec = _mapping(spec, "pump.profile")
    kind = spec.pop("kind", None)
    if kind == "uniform":
        _rejectUnknown(spec, set(), "pump.profile")
        return UniformPumpProfile()
    if kind != "super_gaussian":
        raise ValueError("pump.profile.kind must be 'uniform' or 'super_gaussian'")
    _rejectUnknown(spec, {"radius_u", "radius_v", "exponent", "center", "axis_u", "axis_v"}, "pump.profile")
    return SuperGaussianPumpProfile(**spec)


def _pump(spec, context):
    spec = _mapping(spec, "simulation.pumps[]")
    allowed = {
        "name", "total_power", "ray_count", "pump_steps", "rng_seed", "spectrum",
        "angular_distribution", "profile", "injection", "relays",
    }
    _rejectUnknown(spec, allowed, "simulation.pumps[]")
    pump = Pump(
        name=spec.get("name"),
        total_power=spec["total_power"],
        ray_count=spec["ray_count"],
        pump_steps=spec.get("pump_steps"),
        rng_seed=spec.get("rng_seed", 5489),
        spectrum=_pumpSpectrum(spec["spectrum"]),
        angular_distribution=_angularDistribution(spec.get("angular_distribution", {"collimated": {}})),
        profile=_profile(spec.get("profile", {"kind": "uniform"})),
    )
    injection = _mapping(spec["injection"], "pump.injection")
    _rejectUnknown(injection, {"domain"}, "pump.injection")
    injector = SurfacePumpInjector(context.resolve("domains", injection["domain"]))
    relays = []
    relayKeys = {
        "exit_domain", "entry_domain", "flip_u", "flip_v", "rotation", "offset",
        "tilt", "magnification", "transmission",
    }
    for value in spec.get("relays", []):
        value = _mapping(value, "pump.relays[]")
        _rejectUnknown(value, relayKeys, "pump.relays[]")
        relays.append(
            PlanarPumpRelay(
                context.resolve("domains", value.pop("exit_domain")),
                context.resolve("domains", value.pop("entry_domain")),
                **value,
            )
        )
    return pump, injector, tuple(relays)


def _timeIntegrator(spec):
    spec = _mapping(spec, "simulation.time_integrator")
    method = spec.pop("method", None)
    classes = {
        "explicit_euler": ExplicitEuler,
        "heun": Heun,
        "midpoint": Midpoint,
        "runge_kutta4": RungeKutta4,
        "frozen_phi_ase_runge_kutta4": FrozenPhiAseRungeKutta4,
        "frozen_sources_runge_kutta4": FrozenSourcesRungeKutta4,
        "implicit_euler": ImplicitEuler,
        "exponential_euler": ExponentialEuler,
    }
    if method not in classes:
        raise ValueError("unsupported simulation.time_integrator.method")
    if method == "implicit_euler":
        _rejectUnknown(spec, {"iterations", "tolerance"}, "simulation.time_integrator")
        return ImplicitEuler(**spec)
    _rejectUnknown(spec, set(), "simulation.time_integrator")
    return classes[method]()


def objectFromYaml(cls, filename, name, **objects):
    """Resolve only one named major primitive from schema-v3 YAML."""
    try:
        registry = _REGISTRY_KEYS[cls]
    except KeyError as exc:
        raise TypeError(f"{cls.__name__}.fromYaml is not a schema-v3 major primitive") from exc
    return _YamlContext(filename, **objects).resolve(registry, name)


def simulationFromYaml(filename, *, simulationCls, **objects):
    """Construct the complete schema-v3 simulation graph without executing it."""
    context = _YamlContext(filename, **objects)
    spec = _mapping(context.data.get("simulation"), "simulation")
    allowed = {
        "optical_components",
        "gain_medium",
        "exterior_surface",
        "initial_excitation",
        "phi_ase",
        "pumps",
        "time_integrator",
        "time_step_size",
        "simulation_steps",
        "max_time",
        "pre_pump",
        "report_timings",
        "execution_mode",
        "output_steps",
        "output_fields",
        "control_fields",
    }
    _rejectUnknown(spec, allowed, "simulation")
    if "simulation_steps" in spec and "max_time" in spec:
        raise ValueError("simulation configures at most one of simulation_steps and max_time")
    componentNames = spec.get("optical_components")
    if not isinstance(componentNames, list) or not componentNames:
        raise ValueError("simulation.optical_components must be a non-empty sequence")
    components = [context.resolve("optical_components", name) for name in componentNames]
    medium = context.resolve("gain_media", spec["gain_medium"])
    excitationSpec = spec.get("initial_excitation", {"value": 0.0})
    excitationSpec = _mapping(excitationSpec, "simulation.initial_excitation")
    form = _oneOf(excitationSpec, ("value", "domains"), "simulation.initial_excitation")
    _rejectUnknown(excitationSpec, {form}, "simulation.initial_excitation")
    if form == "value":
        excitation = excitationSpec["value"]
    else:
        assignments = excitationSpec["domains"]
        if not isinstance(assignments, list) or not assignments:
            raise ValueError("simulation.initial_excitation.domains must be a non-empty sequence")
        excitation = {}
        for index, assignment in enumerate(assignments):
            assignment = _mapping(assignment, f"simulation.initial_excitation.domains[{index}]")
            _rejectUnknown(assignment, {"domain", "value"}, f"simulation.initial_excitation.domains[{index}]")
            excitation[context.resolve("domains", assignment["domain"])] = assignment["value"]
    simulation = simulationCls(
        opticalComponents=components,
        gainMedium=medium,
        exteriorSurface=(
            None
            if spec.get("exterior_surface") is None
            else context.resolve("domains", spec["exterior_surface"])
        ),
        initialExcitation=excitation,
        phiASE=_phiAse(spec["phi_ase"]),
        timeIntegrator=_timeIntegrator(spec["time_integrator"]),
        timeStepSize=spec["time_step_size"],
        simulationSteps=spec.get("simulation_steps"),
        maxTime=spec.get("max_time"),
        prePump=spec.get("pre_pump", False),
        reportTimings=spec.get("report_timings", False),
        executionMode=str(spec.get("execution_mode", "autonomous")).replace("_", "-"),
        outputSteps=spec.get("output_steps"),
        outputFields=spec.get("output_fields"),
        controlFields=spec.get("control_fields", ()),
    )
    for pumpSpec in spec.get("pumps", []):
        pump, injector, relays = _pump(pumpSpec, context)
        simulation.addPump(pump, injector, relays=relays)
    return simulation


__all__ = ["objectFromYaml", "simulationFromYaml"]
