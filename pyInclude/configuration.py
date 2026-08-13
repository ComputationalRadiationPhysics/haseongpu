# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Declarative schema-v2 construction of public HASEonGPU Python objects."""

from __future__ import annotations

from pathlib import Path

from .geometry import GainMedium, SurfaceOptics, VolumeTopology
from .geometry.volume import BOUND_INTERNAL, BOUND_STOP
from .laser import (
    CrossSectionData,
    PlanarPumpRelay,
    Pump,
    PumpAngularDistribution,
    PumpSpectrum,
    SuperGaussianPumpProfile,
    SurfacePumpInjector,
    UniformPumpProfile,
)
from .simulation import PhiASE
from .timeIntegration import (
    ExplicitEuler,
    ExponentialEuler,
    FrozenPhiAseRungeKutta4,
    Heun,
    ImplicitEuler,
    Midpoint,
    RungeKutta4,
)


_REMOVED_OPTICS_FIELDS = {
    "refractive_indices",
    "reflectivities",
    "surface_reflectivity",
    "surface_refractive_index_inside",
    "surface_refractive_index_outside",
}


def _mapping(value, path):
    if not isinstance(value, dict):
        raise ValueError(f"{path} must be a mapping")
    return dict(value)


def _reject_unknown(mapping, allowed, path):
    unknown = sorted(set(mapping) - set(allowed))
    if unknown:
        raise ValueError(f"unsupported {path} options: {unknown}")


def _one_of(mapping, names, path):
    selected = [name for name in names if name in mapping]
    if len(selected) != 1:
        raise ValueError(f"{path} requires exactly one of: {', '.join(names)}")
    return selected[0]


def _resolve_path(root, value):
    path = Path(value).expanduser()
    return path if path.is_absolute() else root / path


def _load_yaml(filename):
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("simulation YAML requires PyYAML") from exc
    path = Path(filename).expanduser().resolve()
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"simulation config '{path}' must contain a mapping")
    if data.get("schema_version") != 2:
        raise ValueError("simulation YAML requires schema_version: 2")
    _reject_unknown(data, {"schema_version", "cross_sections", "simulation"}, "top-level")
    return path, data


def _cross_section_data(spec, root, path):
    spec = _mapping(spec, path)
    form = _one_of(spec, ("from_directory", "monochromatic", "inline"), path)
    _reject_unknown(spec, {form}, path)
    values = spec[form]
    if form == "from_directory":
        if isinstance(values, (str, Path)):
            directory = values
            resolution = 1000
        else:
            values = _mapping(values, f"{path}.from_directory")
            _reject_unknown(values, {"path", "resolution"}, f"{path}.from_directory")
            directory = values["path"]
            resolution = values.get("resolution", 1000)
        return CrossSectionData.fromDirectory(_resolve_path(root, directory), resolution=resolution)
    values = _mapping(values, f"{path}.{form}")
    if form == "monochromatic":
        _reject_unknown(
            values,
            {"wavelength", "cross_section_absorption", "cross_section_emission"},
            f"{path}.monochromatic",
        )
        return CrossSectionData.monochromatic(
            wavelength=values["wavelength"],
            crossSectionAbsorption=values["cross_section_absorption"],
            crossSectionEmission=values["cross_section_emission"],
        )
    _reject_unknown(
        values,
        {
            "wavelengths_absorption",
            "cross_section_absorption",
            "wavelengths_emission",
            "cross_section_emission",
            "resolution",
        },
        f"{path}.inline",
    )
    return CrossSectionData(
        wavelengthsAbsorption=values["wavelengths_absorption"],
        crossSectionAbsorption=values["cross_section_absorption"],
        wavelengthsEmission=values["wavelengths_emission"],
        crossSectionEmission=values["cross_section_emission"],
        resolution=values.get("resolution", 1),
    )


def _load_cross_sections(spec, root):
    spec = _mapping(spec, "cross_sections")
    if not spec:
        raise ValueError("cross_sections must define at least one CrossSectionData object")
    return {
        str(name): _cross_section_data(value, root, f"cross_sections.{name}")
        for name, value in spec.items()
    }


def _topology(spec, root):
    spec = _mapping(spec, "simulation.gain_medium.topology")
    form = _one_of(spec, ("from_file", "from_tetrahedra"), "simulation.gain_medium.topology")
    _reject_unknown(spec, {form}, "simulation.gain_medium.topology")
    values = _mapping(spec[form], f"simulation.gain_medium.topology.{form}")
    if form == "from_tetrahedra":
        _reject_unknown(
            values,
            {"points", "cell_point_indices", "cell_domains", "face_boundaries", "metadata"},
            "simulation.gain_medium.topology.from_tetrahedra",
        )
        return VolumeTopology.fromTetrahedra(
            values["points"],
            values["cell_point_indices"],
            cellDomains=values.get("cell_domains"),
            faceBoundaries=values.get("face_boundaries"),
            metadata=values.get("metadata"),
        )

    _reject_unknown(
        values,
        {"path", "format", "boundary_default", "mesh_size"},
        "simulation.gain_medium.topology.from_file",
    )
    path = _resolve_path(root, values["path"])
    mesh_format = values.get("format")
    if mesh_format == "auto":
        mesh_format = None
    kwargs = {}
    if "boundary_default" in values:
        boundary = values["boundary_default"]
        boundary = {"stop": BOUND_STOP, "internal": BOUND_INTERNAL}.get(boundary, boundary)
        kwargs["boundaryDefault"] = boundary
    if "mesh_size" in values and values["mesh_size"] is not None:
        kwargs["meshSize"] = values["mesh_size"]
    return VolumeTopology.fromFile(path, format=mesh_format, **kwargs)


def _assignment(spec, *, surface):
    spec = _mapping(spec, "surface domain assignment" if surface else "cell domain assignment")
    aliases = {
        "cell_indices": "cellIndices",
        "face_indices": "faceIndices",
        "gmsh_name": "gmshName",
        "gmsh_tag": "gmshTag",
        "allow_internal": "allowInternal",
    }
    allowed = {"domain", "name", "where", "cell_indices", "gmsh_name", "gmsh_tag"}
    if surface:
        allowed = {"domain", "name", "where", "face_indices", "gmsh_name", "gmsh_tag", "allow_internal"}
    _reject_unknown(spec, allowed, "surface domain assignment" if surface else "cell domain assignment")
    return {aliases.get(name, name): value for name, value in spec.items()}


def _gain_medium(spec, root):
    spec = _mapping(spec, "simulation.gain_medium")
    _reject_unknown(
        spec,
        {"from_vtk", "topology", "cell_domains", "surface_domains", "properties", "surface_optics", "custom_fields"},
        "simulation.gain_medium",
    )
    if ("from_vtk" in spec) == ("topology" in spec):
        raise ValueError("simulation.gain_medium requires exactly one of from_vtk or topology")
    medium = (
        GainMedium.fromVtk(_resolve_path(root, spec["from_vtk"]))
        if "from_vtk" in spec
        else GainMedium(_topology(spec["topology"], root))
    )
    if "cell_domains" in spec:
        assignments = [_assignment(value, surface=False) for value in spec["cell_domains"]]
        medium.topology = medium.topology.withCellDomains(assignments)
    if "surface_domains" in spec:
        assignments = [_assignment(value, surface=True) for value in spec["surface_domains"]]
        medium.topology = medium.topology.withSurfaceDomains(assignments)

    properties = _mapping(spec.get("properties", {}), "simulation.gain_medium.properties")
    removed = sorted(set(properties) & _REMOVED_OPTICS_FIELDS)
    if removed:
        raise ValueError(
            "raw optics fields are not public schema-v2 properties; configure surface_optics instead: "
            + ", ".join(removed)
        )
    property_names = {
        "beta_volume": "betaVolume",
        "cladding_cell_types": "claddingCellTypes",
        "n_tot": "nTot",
        "fluorescence_lifetime": "crystalTFluo",
        "cladding_number": "claddingNumber",
        "cladding_absorption": "claddingAbsorption",
    }
    _reject_unknown(properties, property_names, "simulation.gain_medium.properties")
    medium.withPhysicalProperties(**{property_names[name]: value for name, value in properties.items()})

    optics = _mapping(spec.get("surface_optics", {}), "simulation.gain_medium.surface_optics")
    if optics:
        optics_by_domain = {}
        for domain, value in optics.items():
            value = _mapping(value, f"simulation.gain_medium.surface_optics.{domain}")
            _reject_unknown(
                value,
                {"reflectivity", "n_inside", "n_outside"},
                f"simulation.gain_medium.surface_optics.{domain}",
            )
            optics_by_domain[domain] = SurfaceOptics(**value)
        medium.with_surface_optics(
            optics_by_domain
        )
    custom_fields = spec.get("custom_fields", [])
    if not isinstance(custom_fields, list):
        raise ValueError("simulation.gain_medium.custom_fields must be a sequence")
    aliases = {
        "unit_si": "unitSI",
        "unit_dimension": "unitDimension",
        "backend_required": "backendRequired",
    }
    for index, field in enumerate(custom_fields):
        field = _mapping(field, f"simulation.gain_medium.custom_fields[{index}]")
        _reject_unknown(
            field,
            {"name", "entity", "values", "dtype", "unit", "unit_si", "unit_dimension", "dynamic", "backend_required"},
            f"simulation.gain_medium.custom_fields[{index}]",
        )
        name = field.pop("name")
        medium.defineField(name, **{aliases.get(key, key): value for key, value in field.items()})
    return medium


def _phi_ase(spec, cross_sections):
    spec = _mapping(spec, "simulation.phi_ase")
    reference = spec.pop("cross_sections", None)
    if reference is None:
        raise ValueError("simulation.phi_ase.cross_sections is required")
    try:
        spectra = cross_sections[str(reference)]
    except KeyError as exc:
        raise ValueError(f"unknown cross_sections reference: {reference}") from exc
    aliases = {
        "propagation_mode": "propagationMode",
        "min_rays": "minRays",
        "max_rays": "maxRays",
        "forward_ray_count": "forwardRayCount",
        "relative_standard_error_threshold": "relativeStandardErrorThreshold",
        "adaptive_steps": "adaptiveSteps",
        "use_reflections": "useReflections",
        "reflection_max_iterations": "reflectionMaxIterations",
        "reflection_tolerance": "reflectionTolerance",
        "surface_reservoir_size": "surfaceReservoirSize",
        "rng_seed": "rngSeed",
        "openpmd_backend": "openpmdBackend",
        "parallel_mode": "parallelMode",
        "num_devices": "numDevices",
        "n_per_node": "nPerNode",
        "min_sample_range": "minSampleRange",
        "max_sample_range": "maxSampleRange",
    }
    unchanged = {"repetitions", "monochromatic", "backend", "ase_steps"}
    _reject_unknown(spec, set(aliases) | unchanged, "simulation.phi_ase")
    values = {aliases.get(name, name): value for name, value in spec.items()}
    return PhiASE(crossSections=spectra, spectralProperties=spectra, **values), spectra


def _pump_spectrum(spec):
    spec = _mapping(spec, "pump.spectrum")
    if "monochromatic" in spec:
        _reject_unknown(spec, {"monochromatic"}, "pump.spectrum")
        return PumpSpectrum.monochromatic(spec["monochromatic"])
    _reject_unknown(spec, {"wavelengths", "weights"}, "pump.spectrum")
    return PumpSpectrum(spec["wavelengths"], spec["weights"])


def _angular_distribution(spec):
    spec = _mapping(spec, "pump.angular_distribution")
    form = _one_of(spec, ("collimated", "uniform_cone", "discrete"), "pump.angular_distribution")
    _reject_unknown(spec, {form}, "pump.angular_distribution")
    values = _mapping(spec[form], f"pump.angular_distribution.{form}")
    if form == "collimated":
        _reject_unknown(values, set(), "pump.angular_distribution.collimated")
        return PumpAngularDistribution.collimated()
    if form == "uniform_cone":
        _reject_unknown(values, {"half_angle", "polar_samples", "azimuthal_samples"}, "pump.angular_distribution.uniform_cone")
        return PumpAngularDistribution.uniform_cone(**values)
    _reject_unknown(values, {"polar_angles", "azimuthal_angles", "weights"}, "pump.angular_distribution.discrete")
    return PumpAngularDistribution(**values)


def _profile(spec):
    spec = _mapping(spec, "pump.profile")
    kind = spec.pop("kind", None)
    if kind == "uniform":
        _reject_unknown(spec, set(), "pump.profile")
        return UniformPumpProfile()
    if kind != "super_gaussian":
        raise ValueError("pump.profile.kind must be 'uniform' or 'super_gaussian'")
    _reject_unknown(spec, {"radius_u", "radius_v", "exponent", "center", "axis_u", "axis_v"}, "pump.profile")
    return SuperGaussianPumpProfile(**spec)


def _pump(spec, cross_sections):
    spec = _mapping(spec, "simulation.pumps[]")
    _reject_unknown(
        spec,
        {
            "name", "total_power", "cross_sections", "ray_count", "pump_steps", "rng_seed",
            "spectrum", "angular_distribution", "profile", "injection", "relays",
        },
        "simulation.pumps[]",
    )
    reference = str(spec["cross_sections"])
    if reference not in cross_sections:
        raise ValueError(f"unknown cross_sections reference: {reference}")
    pump = Pump(
        name=spec.get("name"),
        total_power=spec["total_power"],
        ray_count=spec["ray_count"],
        pump_steps=spec.get("pump_steps"),
        rng_seed=spec.get("rng_seed", 5489),
        cross_sections=cross_sections[reference],
        spectrum=_pump_spectrum(spec["spectrum"]),
        angular_distribution=_angular_distribution(
            spec.get("angular_distribution", {"collimated": {}})
        ),
        profile=_profile(spec.get("profile", {"kind": "uniform"})),
    )
    injection = _mapping(spec["injection"], "pump.injection")
    _reject_unknown(injection, {"surface_domains"}, "pump.injection")
    relays = []
    relay_keys = {
        "exit_domains",
        "entry_domains",
        "flip_u",
        "flip_v",
        "rotation",
        "offset",
        "tilt",
        "magnification",
        "transmission",
    }
    for relay in spec.get("relays", []):
        relay = _mapping(relay, "pump.relays[]")
        _reject_unknown(relay, relay_keys, "pump.relays[]")
        relays.append(PlanarPumpRelay(**relay))
    return pump, SurfacePumpInjector(injection["surface_domains"]), tuple(relays)


def _time_integrator(spec):
    spec = _mapping(spec, "simulation.time_integrator")
    method = spec.pop("method", None)
    classes = {
        "explicit_euler": ExplicitEuler,
        "heun": Heun,
        "midpoint": Midpoint,
        "runge_kutta4": RungeKutta4,
        "frozen_phi_ase_runge_kutta4": FrozenPhiAseRungeKutta4,
        "implicit_euler": ImplicitEuler,
        "exponential_euler": ExponentialEuler,
    }
    if method not in classes:
        raise ValueError("unsupported simulation.time_integrator.method")
    if method == "implicit_euler":
        _reject_unknown(spec, {"iterations", "tolerance"}, "simulation.time_integrator")
        return ImplicitEuler(**spec)
    _reject_unknown(spec, set(), "simulation.time_integrator")
    return classes[method]()


def simulation_from_yaml(filename, *, simulation_cls):
    """Construct one ``Simulation`` from schema-v2 YAML without executing it."""
    path, data = _load_yaml(filename)
    root = path.parent
    cross_sections = _load_cross_sections(data.get("cross_sections", {}), root)
    spec = _mapping(data.get("simulation"), "simulation")
    allowed = {
        "gain_medium",
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
    _reject_unknown(spec, allowed, "simulation")
    if "simulation_steps" in spec and "max_time" in spec:
        raise ValueError("simulation configures at most one of: simulation_steps, max_time")
    medium = _gain_medium(spec["gain_medium"], root)
    phi_ase, spectra = _phi_ase(spec["phi_ase"], cross_sections)
    mode = str(spec.get("execution_mode", "autonomous")).replace("_", "-")
    simulation = simulation_cls(
        gain_medium=medium,
        phi_ase=phi_ase,
        time_integrator=_time_integrator(spec["time_integrator"]),
        time_step_size=spec["time_step_size"],
        cross_sections=spectra,
        simulation_steps=spec.get("simulation_steps"),
        max_time=spec.get("max_time"),
        pre_pump=spec.get("pre_pump", False),
        report_timings=spec.get("report_timings", False),
        execution_mode=mode,
        output_steps=spec.get("output_steps"),
        output_fields=spec.get("output_fields"),
        control_fields=spec.get("control_fields", ()),
    )
    pumps = spec.get("pumps", [])
    if not pumps:
        raise ValueError("simulation.pumps must contain at least one Pump registration")
    for pump_spec in pumps:
        pump, injection, relays = _pump(pump_spec, cross_sections)
        simulation.add_pump(pump, injection_method=injection, relays=relays)
    return simulation
