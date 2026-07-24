from __future__ import annotations

import contextlib
import os
import queue
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ..geometry import OpenPmdComponentField, OpenPmdScalarField
from .._runtime import runtime_config, runtime_executable_candidates, runtime_root
from .backends import _clean_backend_names, _load_backend_names
from . import HASE_TRANSPORT_VERSION, FieldSpec, backendFlatArray, fieldSpec, flatEntityLabel, haseTransportAttributes, resultFieldSpecs, simulationAttributeSpecs, spectralContext, unitDimension
from ..structures import Result


CANONICAL_POINTS_SPEC = FieldSpec(
    "canonicalPoints",
    "points",
    ("coordinate", "mesh_point"),
    np.float64,
    lambda context: (3, context.numberOfMeshPoints),
    unit="m",
    unitDimension=unitDimension.canonicalPoints,
    backendRequired=False,
)
CANONICAL_CONNECTIVITY_SPEC = FieldSpec(
    "canonicalConnectivity",
    "cells_connectivity",
    ("cell", "local_vertex"),
    np.uint32,
    lambda context: (context.numberOfPrisms, 6),
    unitDimension=unitDimension.canonicalConnectivity,
    backendRequired=False,
)
CANONICAL_OFFSETS_SPEC = FieldSpec(
    "canonicalOffsets",
    "cells_offsets",
    ("cell_offset",),
    np.uint32,
    lambda context: (context.numberOfPrisms + 1,),
    unitDimension=unitDimension.canonicalOffsets,
    backendRequired=False,
)
CANONICAL_CELL_TYPES_SPEC = FieldSpec(
    "canonicalCellTypes",
    "cells_types",
    ("cell",),
    np.uint32,
    lambda context: (context.numberOfPrisms,),
    unitDimension=unitDimension.canonicalCellTypes,
    backendRequired=False,
)
DYNAMIC_FIELD_NAMES = {"betaVolume", "pointBeta"}


def _env_flag(name):
    value = os.environ.get(name)
    if value is None:
        return None
    return value.strip().lower() in {"1", "true", "on", "yes"}


def _runtime_config():
    return runtime_config()


def _forward_backend_logging_enabled():
    override = _env_flag("HASE_FORWARD_LOGGING")
    if override is not None:
        return override
    return bool(getattr(_runtime_config(), "HASE_FORWARD_LOGGING", False))


def _forward_backend_logging(stdout="", stderr=""):
    if not _forward_backend_logging_enabled():
        return
    if stdout:
        sys.stdout.write(stdout)
        sys.stdout.flush()
    if stderr:
        sys.stderr.write(stderr)
        sys.stderr.flush()


def _backend_failure_detail(stdout="", stderr=""):
    sections = []
    if stdout and stdout.strip():
        sections.append("calcPhiASE stdout:\n" + stdout.strip())
    if stderr and stderr.strip():
        sections.append("calcPhiASE stderr:\n" + stderr.strip())
    return ": " + "\n".join(sections) if sections else ""


@dataclass(frozen=True)
class _AttributeField:
    name: str
    value: object


@dataclass(frozen=True)
class _ScalarArrayField:
    spec: FieldSpec
    values: object
    context: object
    prefix: str = "core_"


@dataclass(frozen=True)
class _ComponentArrayField:
    recordName: str
    spec: FieldSpec
    components: dict[str, object]
    axisLabels: list[str]
    context: object
    prefix: str = "core_"


@dataclass(frozen=True)
class _BackendSpec:
    name: str
    suffix: str
    config: dict
    streaming: bool = False


ADIOS2_CONFIG = {"backend": "adios2"}
HDF5_CONFIG = {"backend": "hdf5"}
SST_CONFIG = {
    "backend": "adios2",
    "adios2": {
        "engine": {
            "type": "sst",
            "parameters": {
                "DataTransport": "WAN",
                "OpenTimeoutSecs": "600",
            }
        }
    },
}

_STREAMING_RESULT_EOF = object()
_STREAMING_RESULT_POLL_SECONDS = 0.1
_STREAMING_THREAD_JOIN_TIMEOUT_SECONDS = 10.0

OPENPMD_BACKENDS = {
    "adios": _BackendSpec("adios", ".bp", ADIOS2_CONFIG),
    "adios-sst": _BackendSpec("adios-sst", ".sst", SST_CONFIG, streaming=True),
    "hdf5": _BackendSpec("hdf5", ".h5", HDF5_CONFIG),
}
DEFAULT_OPENPMD_BACKEND = "auto"
OPENPMD_BACKEND_PRIORITY = ("adios-sst", "adios", "hdf5")
HASE_CONFIGURE_HINT = "Run `hase-configure` to generate a matching backend/openPMD setup."
_OPENPMD_BACKEND_PROBE_CACHE = {}


def _normalize_backend(backend=None):
    value = backend if backend is not None else DEFAULT_OPENPMD_BACKEND
    normalized = str(value).strip().lower()
    if normalized != "auto" and normalized not in OPENPMD_BACKENDS:
        allowed = ", ".join(("auto", *OPENPMD_BACKEND_PRIORITY))
        raise ValueError(f"unsupported openPMD backend '{value}'; expected one of: {allowed}. {HASE_CONFIGURE_HINT}")
    return normalized


def _backend_spec(backend=None):
    normalized = _normalize_backend(backend)
    if normalized == "auto":
        raise RuntimeError("openPMD backend 'auto' must be resolved against the installed runtime")
    return OPENPMD_BACKENDS[normalized]


def _probed_openpmd_backends(executable):
    executable = Path(executable).resolve()
    cache_key = str(executable)
    if cache_key in _OPENPMD_BACKEND_PROBE_CACHE:
        return _OPENPMD_BACKEND_PROBE_CACHE[cache_key]

    backends, probe_library = _load_backend_names(
        (
            executable.parent,
            executable.parent / "python" / "pyInclude" / "_runtime",
        )
    )
    backends = _clean_backend_names(backends)
    if not backends:
        raise RuntimeError(
            "openPMD backend-probe library did not report any supported backends. " + HASE_CONFIGURE_HINT
        )
    result = (backends, probe_library)
    _OPENPMD_BACKEND_PROBE_CACHE[cache_key] = result
    return result


def _ensure_compiled_backend_available(spec, executable):
    backends, probe_library = _probed_openpmd_backends(executable)
    if spec.name not in backends:
        available = ", ".join(backends)
        raise RuntimeError(
            f"openPMD backend '{spec.name}' is selected by PhiASE.openpmdBackend/YAML "
            f"'openpmd_backend', but the runtime openPMD backend-probe library reports "
            f"available backends: {available}. Probe library: {probe_library}. Change "
            f"openpmd_backend in the YAML or fix the runtime openPMD provider environment. {HASE_CONFIGURE_HINT}"
        )


def _python_supported_backend_names(io):
    variants = getattr(io, "variants", {})
    extensions = set(getattr(io, "file_extensions", []))
    supported = []
    for name in OPENPMD_BACKEND_PRIORITY:
        spec = OPENPMD_BACKENDS[name]
        provider_enabled = variants.get("hdf5", False) if name == "hdf5" else variants.get("adios2", False)
        if provider_enabled and spec.suffix.lstrip(".") in extensions:
            supported.append(name)
    return tuple(supported)


def _resolve_backend(backend, executable):
    requested = _normalize_backend(backend)
    executable = Path(executable)
    io = _io(executable)
    compiled, probe_library = _probed_openpmd_backends(executable)
    python_backends = _python_supported_backend_names(io)
    supported = tuple(name for name in OPENPMD_BACKEND_PRIORITY if name in compiled and name in python_backends)
    if requested == "auto":
        if supported:
            return OPENPMD_BACKENDS[supported[0]]
        raise RuntimeError(
            "No compatible openPMD runtime backend was found. "
            f"Compiled provider supports: {', '.join(compiled)}; Python openPMD-api supports: "
            f"{', '.join(python_backends)}; probe library: {probe_library}. {HASE_CONFIGURE_HINT}"
        )
    spec = OPENPMD_BACKENDS[requested]
    if requested not in compiled:
        _ensure_compiled_backend_available(spec, executable)
    if requested not in python_backends:
        raise RuntimeError(
            f"openPMD backend '{requested}' is selected by PhiASE.openpmdBackend/YAML "
            f"'openpmd_backend', but the active Python openPMD-api supports: "
            f"{', '.join(python_backends)}. {HASE_CONFIGURE_HINT}"
        )
    return spec


def _truthy(value):
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _streaming_thread_join_timeout(value=None):
    if value is None:
        value = os.environ.get(
            "HASE_OPENPMD_THREAD_JOIN_TIMEOUT",
            str(_STREAMING_THREAD_JOIN_TIMEOUT_SECONDS),
        )
    try:
        seconds = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("HASE_OPENPMD_THREAD_JOIN_TIMEOUT must be a positive number of seconds") from exc
    if seconds <= 0.0:
        raise ValueError("HASE_OPENPMD_THREAD_JOIN_TIMEOUT must be a positive number of seconds")
    return seconds


def _watchdog_interval(value=None):
    if value is None:
        value = os.environ.get("HASE_OPENPMD_WATCHDOG_INTERVAL", "30")
    text = str(value).strip().lower()
    if text in {"0", "none", "off", "false", "no"}:
        return None
    try:
        seconds = float(text)
    except ValueError as exc:
        raise ValueError("HASE_OPENPMD_WATCHDOG_INTERVAL must be a positive number of seconds, 0, or 'none'") from exc
    if seconds <= 0.0:
        raise ValueError("HASE_OPENPMD_WATCHDOG_INTERVAL must be a positive number of seconds, 0, or 'none'")
    return seconds


def _artifact_root():
    explicit = os.environ.get("HASE_OPENPMD_ARTIFACT_DIR")
    if explicit:
        return Path(explicit)
    if _truthy(os.environ.get("HASE_OPENPMD_KEEP_ARTIFACTS", "")):
        return Path.cwd() / "hase-openpmd-artifacts"
    return None


def _safe_artifact_name(value):
    allowed = []
    for char in str(value):
        allowed.append(char if char.isalnum() or char in {"-", "_", "."} else "-")
    return "".join(allowed).strip(".-_") or "transport"


def _artifact_run_id():
    explicit = os.environ.get("HASE_OPENPMD_ARTIFACT_RUN_ID")
    if explicit:
        return _safe_artifact_name(explicit)
    prefix = _safe_artifact_name(os.environ.get("HASE_OPENPMD_ARTIFACT_PREFIX", "transport"))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return f"{prefix}-{stamp}-{os.getpid()}"


def _write_openpmd_handle(handle_path: Path, series_path: Path):
    handle_path.write_text(series_path.name + "\n", encoding="utf-8")


def _write_artifact_manifest(path: Path, *, backend, input_path, output_path, input_handle, output_handle, status, return_code=None):
    lines = [
        f"backend={backend}",
        f"status={status}",
        f"input={input_path}",
        f"inputHandle={input_handle}",
        f"output={output_path}",
        f"outputHandle={output_handle}",
    ]
    if return_code is not None:
        lines.append(f"returnCode={return_code}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _fieldContext(gainMedium):
    topology = gainMedium.topology
    topology._require_levels()
    if topology.thickness is None:
        raise ValueError("topology thickness is required before running a simulation")
    return SimpleNamespace(
        numberOfPoints=topology.numberOfPoints,
        numberOfTriangles=topology.numberOfTriangles,
        numberOfLevels=int(topology.levels),
    )


def _validatePhiAseTransportOptions(phiAse):
    if bool(phiAse.writeVtk):
        raise ValueError("PhiASE.writeVtk is not supported by the openPMD transport")
    if getattr(phiAse, "devices", None):
        raise ValueError("PhiASE.devices is not supported by the openPMD transport")


def _attributeFields(phiAse, gainMedium, crossSections):
    values = _attributeValues(phiAse, gainMedium, crossSections)
    for spec in simulationAttributeSpecs:
        if spec.name not in values:
            if spec.name == "rngSeed":
                continue
            raise KeyError(spec.name)
        yield _AttributeField(spec.attribute, spec.cast(values[spec.name]))


def _attributeValues(phiAse, gainMedium, crossSections):
    _validatePhiAseTransportOptions(phiAse)
    context = _fieldContext(gainMedium)
    number_of_samples = context.numberOfPoints * context.numberOfLevels
    values = {}
    values.update(gainMedium.topology.openPmdAttributes(context))
    values.update(gainMedium.openPmdAttributes(context))
    values.update(crossSections.openPmdAttributes())
    values.update(phiAse.openPmdAttributes(numberOfSamples=number_of_samples))
    return values

def _arrayFields(gainMedium, crossSections, *, include_static=True):
    context = _fieldContext(gainMedium)
    for field in _fieldsFromDomain(gainMedium.openPmdFields(context)):
        if include_static or field.spec.name in DYNAMIC_FIELD_NAMES:
            yield field
    if include_static:
        yield from _fieldsFromDomain(crossSections.openPmdFields(spectralContext))


def _fieldsFromDomain(fields):
    for field in fields:
        if isinstance(field, OpenPmdComponentField):
            spec = fieldSpec(field.name)
            yield _ComponentArrayField(
                recordName=field.recordName or spec.recordName,
                spec=spec,
                components=field.components,
                axisLabels=field.axisLabels,
                context=field.context,
                prefix=field.prefix,
            )
            continue
        if isinstance(field, OpenPmdScalarField):
            yield _ScalarArrayField(
                field.spec if field.spec is not None else fieldSpec(field.name),
                field.values,
                field.context,
                prefix=field.prefix,
            )
            continue
        name, values, context = field
        yield _ScalarArrayField(fieldSpec(name), values, context)



def _unit_dimension(io, exponents):
    labels = (
        io.Unit_Dimension.L,
        io.Unit_Dimension.M,
        io.Unit_Dimension.T,
        io.Unit_Dimension.I,
        io.Unit_Dimension.theta,
        io.Unit_Dimension.N,
        io.Unit_Dimension.J,
    )
    return {label: float(exponent) for label, exponent in zip(labels, exponents) if exponent != 0.0}


def _dimensionless_dimension():
    return {}


def _series_config(path: Path, backend=None):
    if backend is not None:
        return _backend_spec(backend).config
    if path.suffix == ".sst":
        return SST_CONFIG
    if path.suffix == ".h5":
        return HDF5_CONFIG
    return {}


def _io(executable=None):
    _prefer_matching_openpmd_api(findCalcPhiAse() if executable is None else Path(executable))
    try:
        import openpmd_api as io
    except ImportError as exc:
        raise ImportError(
            "The openPMD transport requires an openpmd_api Python module matching "
            "the calcPhiASE/openPMD C++ stack. Install openpmd-api in this Python "
            "environment with the same backend/MPI options as the openPMD C++ "
            "package found by CMake, or build HASEonGPU with "
            "HASE_BUILD_OPENPMD_FROM_SOURCE=ON. " + HASE_CONFIGURE_HINT
        ) from exc
    return io


def _ensure_backend_available(backend, executable=None):
    executable = findCalcPhiAse() if executable is None else Path(executable)
    return _resolve_backend(backend, executable)


def _openpmd_python_package_parent(path):
    path = Path(path)
    return path.parent if path.name == "openpmd_api" else path


def _env_openpmd_python_paths():
    for name in ("HASE_OPENPMD_PYTHONPATH", "HASE_OPENPMD_PYTHON_PACKAGE_DIR"):
        value = os.environ.get(name)
        if not value:
            continue
        for entry in value.split(os.pathsep):
            if entry:
                yield _openpmd_python_package_parent(Path(entry))


def _configured_openpmd_python_paths():
    configured = getattr(_runtime_config(), "HASE_OPENPMD_PYTHON_PACKAGE_DIR", "")
    if configured:
        yield _openpmd_python_package_parent(Path(configured))


def _using_external_openpmd():
    return bool(getattr(_runtime_config(), "HASE_USE_SYSTEM_OPENPMD", False))


def _candidate_python_paths(executable: Path):
    yield from _env_openpmd_python_paths()

    build_dir = executable.parent
    if (build_dir / "openpmd_api").is_dir():
        yield build_dir
    if (build_dir.parent / "openpmd_api").is_dir():
        yield build_dir.parent
    if (build_dir / "site-packages" / "openpmd_api").is_dir():
        yield build_dir / "site-packages"

    yield from _configured_openpmd_python_paths()

    yield from build_dir.glob("_deps/openpmd-build/lib/python*/site-packages")
    for parent in build_dir.parents:
        yield from parent.glob("_deps/openpmd-build/lib/python*/site-packages")


def _unique_existing_directories(paths):
    seen = set()
    for path in paths:
        resolved = Path(path).resolve()
        if resolved in seen or not resolved.is_dir():
            continue
        seen.add(resolved)
        yield resolved


def _prefer_matching_openpmd_api(executable: Path):
    candidates = list(_unique_existing_directories(_candidate_python_paths(executable)))
    if not candidates:
        if _using_external_openpmd():
            return
        raise RuntimeError(
            "The openPMD transport requires an openpmd_api Python module "
            "compatible with the openPMD C++ provider used by calcPhiASE. "
            "Install/load a matching provider, set "
            "-DHASE_OPENPMD_PYTHON_PACKAGE_DIR=<site-packages directory>, or set "
            "HASE_OPENPMD_PYTHONPATH at runtime. " + HASE_CONFIGURE_HINT
        )
    if "openpmd_api" in sys.modules:
        active = Path(getattr(sys.modules["openpmd_api"], "__file__", "")).resolve()
        for candidate in candidates:
            try:
                active.relative_to(candidate.resolve())
                return
            except ValueError:
                pass
        raise RuntimeError(
            "The openPMD transport requires the Python writer and C++ reader to use the same "
            "openPMD-api build/provider. Restart Python with the CMake-selected "
            "openpmd_api package first on PYTHONPATH, e.g. "
            f"PYTHONPATH={candidates[0]}:$PYTHONPATH. {HASE_CONFIGURE_HINT}"
        )

    for candidate in candidates:
        sys.path.insert(0, str(candidate))
        return


def _access(name):
    io = _io()
    if hasattr(io, "Access_Type"):
        return getattr(io.Access_Type, name)
    return getattr(io.Access, name)


def _length_dimension():
    io = _io()
    return _unit_dimension(io, (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))


def _as_array(values, dtype, shape=None, order="C"):
    arr = np.asarray(values, dtype=dtype)
    if shape is not None:
        arr = arr.reshape(shape, order=order)
    return np.ascontiguousarray(arr)


def _reset_scalar_record(
    record,
    data,
    axis_labels,
    unit_dimension=None,
    unit_si=1.0,
    grid_unit_si=1.0,
    grid_spacing=None,
    grid_global_offset=None,
    geometry_parameters="topology=unstructured_triangular_prism",
):
    io = _io()
    record.set_attribute("geometry", "other")
    record.set_attribute("geometryParameters", geometry_parameters)
    record.set_attribute("dataOrder", "C")
    record.axis_labels = axis_labels
    # ADIOS2 SST may not preserve openPMD axisLabels on streamed mesh records.
    # Keep the canonical axisLabels property, plus a scalar fallback for readers.
    record.set_attribute("haseAxisLabelsString", ",".join(axis_labels))
    record.grid_spacing = [1.0] * data.ndim if grid_spacing is None else list(grid_spacing)
    record.grid_global_offset = [0.0] * data.ndim if grid_global_offset is None else list(grid_global_offset)
    record.grid_unit_SI = float(grid_unit_si)
    record.unit_dimension = _dimensionless_dimension() if unit_dimension is None else unit_dimension
    component = record[io.Mesh_Record_Component.SCALAR]
    component.unit_SI = float(unit_si)
    component.position = [0.0] * data.ndim
    component.reset_dataset(io.Dataset(data.dtype, data.shape))
    component.store_chunk(data)


def _record_metadata(record, spec: FieldSpec):
    record.set_attribute("haseTransportVersion", HASE_TRANSPORT_VERSION)
    record.set_attribute("haseEntity", spec.entity)
    record.set_attribute("haseAxes", list(spec.axes))
    # ADIOS2 SST has been observed to stream string-list attributes as an empty
    # scalar string. Keep haseAxes canonical, plus a scalar fallback.
    record.set_attribute("haseAxesString", ",".join(spec.axes))
    record.set_attribute("haseLayoutOrder", "backendFlat")
    record.set_attribute("haseStatic", not spec.dynamic)
    record.set_attribute("haseDynamic", spec.dynamic)
    record.set_attribute("haseBackendRequired", spec.backendRequired)
    record.set_attribute("haseUnit", spec.unit)
    record.set_attribute("haseUserDefined", spec.userDefined)
    if spec.userDefined:
        record.set_attribute("haseUserFieldName", spec.name)


def _resetFlatField(record, spec: FieldSpec, values, context):
    io = _io()
    data = np.ascontiguousarray(backendFlatArray(values, spec, context, layoutOrder="backendFlat"))
    _reset_scalar_record(
        record,
        data,
        [flatEntityLabel(spec)],
        _unit_dimension(io, spec.unitDimension),
        spec.unitSI,
    )
    _record_metadata(record, spec)
    record.set_attribute("hasePrimitiveShape", list(spec.expectedShape(context)))


def _resetComponent(record, component_name, data, axis_labels, unit_dimension, unit_si=1.0):
    io = _io()
    record.set_attribute("geometry", "other")
    record.set_attribute("geometryParameters", "topology=unstructured_triangular_prism")
    record.set_attribute("dataOrder", "C")
    record.axis_labels = axis_labels
    # ADIOS2 SST may not preserve openPMD axisLabels on streamed mesh records.
    # Keep the canonical axisLabels property, plus a scalar fallback for readers.
    record.set_attribute("haseAxisLabelsString", ",".join(axis_labels))
    record.grid_spacing = [1.0] * data.ndim
    record.grid_global_offset = [0.0] * data.ndim
    record.grid_unit_SI = 1.0
    record.unit_dimension = unit_dimension
    component = record[component_name]
    component.unit_SI = float(unit_si)
    component.position = [0.0] * data.ndim
    component.reset_dataset(io.Dataset(data.dtype, data.shape))
    component.store_chunk(data)


def _canonical_topology_context(topology):
    topology._require_levels()
    return SimpleNamespace(
        numberOfMeshPoints=topology.numberOfPoints * int(topology.levels),
        numberOfPrisms=topology.numberOfTriangles * (int(topology.levels) - 1),
    )


def _canonical_point_components(topology):
    points = np.asarray(topology.points, dtype=np.float64)
    z_values = topology.levelCoordinates()
    return {
        "x": np.tile(points[:, 0], z_values.size),
        "y": np.tile(points[:, 1], z_values.size),
        "z": np.repeat(z_values, topology.numberOfPoints),
    }


def _canonical_cell_connectivity(topology):
    triangles = np.asarray(topology.trianglePointIndices, dtype=np.uint32)
    rows = []
    for level in range(int(topology.levels) - 1):
        lower = level * topology.numberOfPoints
        upper = (level + 1) * topology.numberOfPoints
        for tri in triangles:
            ids = [int(vertex) for vertex in tri]
            rows.append([
                ids[0] + lower,
                ids[1] + lower,
                ids[2] + lower,
                ids[0] + upper,
                ids[1] + upper,
                ids[2] + upper,
            ])
    return np.asarray(rows, dtype=np.uint32).reshape(-1)


def _write_canonical_static_topology(iteration, gainMedium):
    topology = gainMedium.topology
    context = _canonical_topology_context(topology)
    record = iteration.meshes["core_" + CANONICAL_POINTS_SPEC.recordName]
    for component_name, values in _canonical_point_components(topology).items():
        _resetComponent(
            record,
            component_name,
            np.ascontiguousarray(values),
            ["mesh_point"],
            _unit_dimension(_io(), CANONICAL_POINTS_SPEC.unitDimension),
            CANONICAL_POINTS_SPEC.unitSI,
        )
    _record_metadata(record, CANONICAL_POINTS_SPEC)
    record.set_attribute("hasePrimitiveShape", list(CANONICAL_POINTS_SPEC.expectedShape(context)))

    connectivity = _canonical_cell_connectivity(topology)
    _resetFlatField(
        iteration.meshes["core_" + CANONICAL_CONNECTIVITY_SPEC.recordName],
        CANONICAL_CONNECTIVITY_SPEC,
        connectivity,
        context,
    )
    _resetFlatField(
        iteration.meshes["core_" + CANONICAL_OFFSETS_SPEC.recordName],
        CANONICAL_OFFSETS_SPEC,
        np.arange(context.numberOfPrisms + 1, dtype=np.uint32) * np.uint32(6),
        context,
    )
    _resetFlatField(
        iteration.meshes["core_" + CANONICAL_CELL_TYPES_SPEC.recordName],
        CANONICAL_CELL_TYPES_SPEC,
        np.full(context.numberOfPrisms, 13, dtype=np.uint32),
        context,
    )

def _loadScalar(series, iteration, name, dtype):
    io = _io()
    component = iteration.meshes[name][io.Mesh_Record_Component.SCALAR]
    chunk = component.load_chunk()
    series.flush()
    return np.array(chunk, dtype=dtype, copy=True).reshape(-1)


def _build_dir_for_executable(executable: Path):
    path = Path(executable).resolve()
    for parent in [path.parent, *path.parents]:
        if (parent / "CMakeCache.txt").is_file():
            return parent
        if parent == Path.cwd().resolve():
            break
    return None


def _target_uses_openpmd_main(build_dir):
    if build_dir is None:
        return True
    manifests = [build_dir / name for name in ("build.ninja", "Makefile", "compile_commands.json")]
    existing = [path for path in manifests if path.is_file()]
    if not existing:
        return True
    return any("src/openpmd_main.cpp" in path.read_text(encoding="utf-8", errors="ignore") for path in existing)


def _is_openpmd_calc_phi_ase(executable: Path):
    return executable.is_file() and _target_uses_openpmd_main(_build_dir_for_executable(executable))


def findCalcPhiAse():
    env = os.environ.get("HASE_CALCPHIASE")
    if env:
        path = Path(env)
        if _is_openpmd_calc_phi_ase(path):
            return path
        raise RuntimeError(f"HASE_CALCPHIASE does not point to an openPMD calcPhiASE binary: {path}")

    root = runtime_root()
    for candidate in runtime_executable_candidates(("calcPhiASE",)):
        if _is_openpmd_calc_phi_ase(candidate):
            return candidate

    raise FileNotFoundError(
        f"Could not find an openPMD calcPhiASE binary in the selected HASE runtime '{root}'. "
        "Build that runtime, set HASE_RUNTIME_DIR, or override the executable with "
        "HASE_CALCPHIASE. " + HASE_CONFIGURE_HINT
    )



def _open_input_series(path, *, backend=None):
    series = _io().Series(str(path), _access("create_linear"), _series_config(path, backend))
    series.set_software("HASEonGPU-openPMD-python-frontend")
    for name, value in haseTransportAttributes.items():
        series.set_attribute(name, value)
    series.set_attribute("haseTransportVersion", HASE_TRANSPORT_VERSION)
    return series


def _write_input_iteration(series, iteration_index, phiAse, gainMedium, crossSections, *, include_static=True, runControl=None):
    iteration = series.snapshots()[int(iteration_index)]
    iteration.time = 0.0
    iteration.dt = 1.0
    iteration.time_unit_SI = 1.0

    for field in _attributeFields(phiAse, gainMedium, crossSections):
        iteration.set_attribute(field.name, field.value)
    if runControl:
        for name, value in runControl.items():
            if value is not None:
                iteration.set_attribute(name, value)

    if include_static:
        iteration.set_attribute("haseStaticUpdate", True)
        _write_canonical_static_topology(iteration, gainMedium)
    else:
        iteration.set_attribute("haseStaticUpdate", False)

    for field in _arrayFields(gainMedium, crossSections, include_static=include_static):
        _writeArrayField(iteration, field)

    iteration.close()


class OpenPmdInputSeries:
    """Context manager for writing HASE input iterations to one openPMD series."""

    def __init__(self, path, *, backend=None):
        self.path = Path(path)
        self.backend = backend
        self._series = None
        self._next_iteration = 0

    def __enter__(self):
        if self.backend is not None:
            self.backend = _ensure_backend_available(self.backend).name
        self._series = _open_input_series(self.path, backend=self.backend)
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()
        return False

    def write(self, phiAse, gainMedium, crossSections, *, iteration_index=None, include_static=None, runControl=None):
        if self._series is None:
            raise RuntimeError("OpenPmdInputSeries must be used as a context manager before writing")
        index = self._next_iteration if iteration_index is None else int(iteration_index)
        write_static = (index == 0) if include_static is None else bool(include_static)
        _write_input_iteration(
            self._series,
            index,
            phiAse,
            gainMedium,
            crossSections,
            include_static=write_static,
            runControl=runControl,
        )
        self._series.flush()
        self._next_iteration = max(self._next_iteration, index + 1)
        return index

    def close(self):
        if self._series is not None:
            self._series.close()
            self._series = None




def _writeArrayField(iteration, field):
    if isinstance(field, _ComponentArrayField):
        record = iteration.meshes[field.prefix + field.recordName]
        for component_name, values in field.components.items():
            data = np.ascontiguousarray(values)
            _resetComponent(
                record,
                component_name,
                data,
                field.axisLabels,
                _unit_dimension(_io(), field.spec.unitDimension),
                field.spec.unitSI,
            )
        _record_metadata(record, field.spec)
        record.set_attribute("hasePrimitiveShape", list(field.spec.expectedShape(field.context)))
        return

    _resetFlatField(
        iteration.meshes[field.prefix + field.spec.recordName],
        field.spec,
        field.values,
        field.context,
    )


def _iteration_index(iteration, fallback=None):
    for name in ("iteration_index", "iterationIndex"):
        if hasattr(iteration, name):
            return int(getattr(iteration, name))
    return fallback


def _read_result_iteration(series, iteration, *, fallback_index=None) -> tuple[int | None, Result]:
    iteration_index = _iteration_index(iteration, fallback_index)
    prefix = "core_result_"
    values = {
        spec.name: _loadScalar(series, iteration, prefix + spec.recordName, spec.dtypeObject)
        for spec in resultFieldSpecs()
    }
    iteration.close()
    return iteration_index, Result(**values)


def read_result(path, *, expected_iteration_index=0) -> Result:
    path = Path(path)
    series = _io().Series(str(path), _access("read_linear"), _series_config(path))
    for fallback_index, iteration in enumerate(series.read_iterations()):
        iteration_index, result = _read_result_iteration(series, iteration, fallback_index=fallback_index)
        series.close()
        if iteration_index is not None and iteration_index != expected_iteration_index:
            raise RuntimeError(
                f"Expected result iteration {expected_iteration_index} in {path}, got {iteration_index}"
            )
        return result
    series.close()
    raise RuntimeError(f"No result iteration was available in {path}")


def _read_optional_scalar(series, iteration, name, dtype, default_size=None):
    if name not in iteration.meshes:
        if default_size is None:
            return None
        return np.zeros(default_size, dtype=dtype)
    return _loadScalar(series, iteration, name, dtype)


def _has_attribute(obj, name):
    if hasattr(obj, "contains_attribute"):
        return obj.contains_attribute(name)
    try:
        obj.get_attribute(name)
        return True
    except Exception:
        return False


def read_simulation_output(path):
    """Read C++ time-stepped simulation snapshots from an openPMD output series."""
    path = Path(path)
    series = _io().Series(str(path), _access("read_linear"), _series_config(path))
    states = []
    for fallback_index, iteration in enumerate(series.read_iterations()):
        iteration_index = _iteration_index(iteration, fallback_index)
        number_of_points = int(iteration.get_attribute("number_of_points"))
        number_of_levels = int(iteration.get_attribute("number_of_levels"))
        number_of_cells = int(iteration.get_attribute("number_of_cells"))
        point_count = number_of_points * number_of_levels
        prism_count = number_of_cells * (number_of_levels - 1)
        states.append(SimpleNamespace(
            iterationIndex=iteration_index,
            step=int(iteration.get_attribute("step_index")) if _has_attribute(iteration, "step_index") else iteration_index + 1,
            time=float(iteration.get_attribute("time")) if _has_attribute(iteration, "time") else float(iteration.time),
            betaCells=_loadScalar(series, iteration, "core_point_beta", np.float64).reshape((number_of_points, number_of_levels), order="F"),
            betaVolume=_loadScalar(series, iteration, "core_beta_volume", np.float64).reshape((number_of_cells, number_of_levels - 1), order="F"),
            phiAse=_loadScalar(series, iteration, "core_result_phi_ase", np.float32).reshape((number_of_points, number_of_levels), order="F"),
            dndtAse=_loadScalar(series, iteration, "core_result_dndt_ase", np.float64).reshape((number_of_points, number_of_levels), order="F"),
            dndtPump=_read_optional_scalar(
                series,
                iteration,
                "core_result_dndt_pump",
                np.float64,
                point_count,
            ).reshape((number_of_points, number_of_levels), order="F"),
            aseResult=Result(
                phiAse=_read_optional_scalar(series, iteration, "core_result_phi_ase", np.float32, point_count),
                mse=_read_optional_scalar(series, iteration, "core_result_mse", np.float64, point_count),
                totalRays=_read_optional_scalar(series, iteration, "core_result_total_rays", np.uint32, point_count),
                dndtAse=_read_optional_scalar(series, iteration, "core_result_dndt_ase", np.float64, point_count),
            ),
            staticUpdate=bool(iteration.get_attribute("haseStaticUpdate")) if _has_attribute(iteration, "haseStaticUpdate") else iteration_index == 0,
        ))
        iteration.close()
    series.close()
    if not states:
        raise RuntimeError(f"No simulation iterations were available in {path}")
    return states


class OpenPmdPhiAseSession:
    """Run PhiASE requests through openPMD and wait for matching result iterations."""

    def __init__(self, *, transport=None, watchdog_interval=None, command_prefix=None, workspace_dir=None):
        self.requested_backend = _normalize_backend(transport)
        self.spec = None
        self.watchdog_interval = _watchdog_interval(watchdog_interval)
        self.command_prefix = [] if command_prefix is None else list(command_prefix)
        self.workspace_dir = None if workspace_dir is None else Path(workspace_dir)
        self._workspace = None
        self._tmp_path = None
        self._manifest_path = None
        self._input_handle = None
        self._output_handle = None
        self._input_path = None
        self._output_path = None
        self._executable = None
        self._proc = None
        self._input_series = None
        self._result_queue = None
        self._send_queue = None
        self._sender_errors = None
        self._watchdog_events = None
        self._watchdog_stop = None
        self._reader = None
        self._sender = None
        self._watchdog = None
        self._pending_results = {}
        self._result_reader_done = False
        self._next_iteration = 0
        self._entered = False

    def __enter__(self):
        artifact_root = _artifact_root()
        if artifact_root is None and self.workspace_dir is not None:
            self.workspace_dir.mkdir(parents=True, exist_ok=True)
        self._workspace = (
            tempfile.TemporaryDirectory(prefix="hase-openpmd-", dir=self.workspace_dir)
            if artifact_root is None
            else contextlib.nullcontext(artifact_root)
        )
        tmp = self._workspace.__enter__()
        self._tmp_path = Path(tmp)
        self._tmp_path.mkdir(parents=True, exist_ok=True)
        self._executable = findCalcPhiAse()
        self.spec = _ensure_backend_available(self.requested_backend, self._executable)

        artifact_id = _artifact_run_id() if artifact_root is not None else None
        if self.spec.streaming:
            stem = f"{artifact_id}-" if artifact_id else ""
            self._input_path = self._tmp_path / f"{stem}input{self.spec.suffix}"
            self._output_path = self._tmp_path / f"{stem}output{self.spec.suffix}"
            self._manifest_path = None if artifact_root is None else self._tmp_path / f"{artifact_id}-manifest.txt"
            self._input_handle = None if artifact_root is None else self._tmp_path / f"{artifact_id}-input.pmd"
            self._output_handle = None if artifact_root is None else self._tmp_path / f"{artifact_id}-output.pmd"
            self._write_handles_and_manifest(status="created")
            self._start_streaming_backend()

        self._entered = True
        return self

    def _calc_phi_ase_command(self, input_path, output_path):
        return [
            *self.command_prefix,
            str(self._executable),
            f"--input-path={input_path}",
            f"--output-path={output_path}",
        ]

    def __exit__(self, exc_type, exc, traceback):
        close_error = None
        try:
            self.close()
        except BaseException as error:
            close_error = error
        finally:
            self._entered = False
            if self._workspace is not None:
                self._workspace.__exit__(exc_type, exc, traceback)
                self._workspace = None
        if exc_type is None and close_error is not None:
            raise close_error
        return False

    def run(self, phiAse, gainMedium, crossSections):
        if not self._entered:
            raise RuntimeError("OpenPmdPhiAseSession must be used as a context manager before running")
        iteration_index = self._next_iteration
        if self.spec.streaming:
            result = self._run_streaming_iteration(iteration_index, phiAse, gainMedium, crossSections)
        else:
            result = self._run_file_iteration(iteration_index, phiAse, gainMedium, crossSections)
        self._next_iteration += 1
        return result

    def close(self):
        if self.spec.streaming:
            self._close_streaming()
            return
        if self._proc is None:
            return

        return_code, stderr = self._finish_backend_process()
        self._write_handles_and_manifest(status="completed" if return_code == 0 else "failed", return_code=return_code)
        if return_code != 0:
            detail = _backend_failure_detail(stderr=stderr)
            raise RuntimeError(f"calcPhiASE failed with return code {return_code}{detail}")

    def _finish_backend_process(self):
        return_code = self._proc.wait()
        stderr = "" if self._proc.stderr is None else self._proc.stderr.read()
        self._proc = None
        return return_code, stderr

    def _join_streaming_thread(self, thread, description):
        if thread is None:
            return None
        timeout = _streaming_thread_join_timeout()
        thread.join(timeout=timeout)
        if thread.is_alive():
            return RuntimeError(f"openPMD {description} thread did not stop within {timeout:g} seconds")
        return None

    def _pop_sender_error(self):
        if self._sender_errors is None:
            return None
        try:
            return self._sender_errors.get_nowait()
        except queue.Empty:
            return None

    def _pop_watchdog_error(self):
        if self._watchdog_events is None:
            return None
        while True:
            try:
                ok, payload = self._watchdog_events.get_nowait()
            except queue.Empty:
                return None
            if not ok:
                return payload

    def _close_streaming(self):
        close_error = None
        if self._send_queue is not None:
            self._send_queue.put(None)

        sender_error = self._join_streaming_thread(self._sender, "input sender")
        self._send_queue = None
        if sender_error is not None:
            close_error = sender_error
            if self._proc is not None and self._proc.poll() is None:
                self._proc.kill()
        elif self._sender is not None:
            self._sender = None

        return_code = None
        stderr = ""
        if self._proc is not None:
            return_code, stderr = self._finish_backend_process()
            self._write_handles_and_manifest(
                status="completed" if return_code == 0 else "failed",
                return_code=return_code,
            )

        reader_error = self._join_streaming_thread(self._reader, "result receiver")
        if reader_error is not None and close_error is None:
            close_error = reader_error
        elif self._reader is not None:
            self._reader = None

        pending_sender_error = self._pop_sender_error()
        if pending_sender_error is not None and close_error is None:
            _, close_error = pending_sender_error

        if return_code not in (None, 0) and close_error is None:
            detail = _backend_failure_detail(stderr=stderr)
            close_error = RuntimeError(f"calcPhiASE failed with return code {return_code}{detail}")

        if self._watchdog_stop is not None:
            self._watchdog_stop.set()
        watchdog_error = self._join_streaming_thread(self._watchdog, "backend watchdog")
        if watchdog_error is not None and close_error is None:
            close_error = watchdog_error
        self._watchdog = None
        self._watchdog_stop = None

        if close_error is not None:
            raise close_error

    def _paths_for_file_iteration(self, iteration_index):
        artifact_root = _artifact_root()
        artifact_id = _artifact_run_id() if artifact_root is not None else None
        if artifact_id is None:
            if iteration_index == 0:
                return self._tmp_path / ("input" + self.spec.suffix), self._tmp_path / ("output" + self.spec.suffix), None
            return (
                self._tmp_path / f"input-{iteration_index}{self.spec.suffix}",
                self._tmp_path / f"output-{iteration_index}{self.spec.suffix}",
                None,
            )
        stem = f"{artifact_id}-{iteration_index}"
        return (
            self._tmp_path / f"{stem}-input{self.spec.suffix}",
            self._tmp_path / f"{stem}-output{self.spec.suffix}",
            self._tmp_path / f"{stem}-manifest.txt",
        )

    def _run_file_iteration(self, iteration_index, phiAse, gainMedium, crossSections):
        input_path, output_path, manifest_path = self._paths_for_file_iteration(iteration_index)
        input_handle = None
        output_handle = None
        if manifest_path is not None:
            input_handle = manifest_path.with_name(manifest_path.stem + "-input.pmd")
            output_handle = manifest_path.with_name(manifest_path.stem + "-output.pmd")
            _write_openpmd_handle(input_handle, input_path)
            _write_openpmd_handle(output_handle, output_path)
            _write_artifact_manifest(
                manifest_path,
                backend=self.spec.name,
                input_path=input_path,
                output_path=output_path,
                input_handle=input_handle,
                output_handle=output_handle,
                status="created",
            )

        with OpenPmdInputSeries(input_path, backend=self.spec.name) as writer:
            writer.write(phiAse, gainMedium, crossSections, iteration_index=iteration_index, include_static=True)
        completed = subprocess.run(
            self._calc_phi_ase_command(input_path, output_path),
            check=False,
            text=True,
            capture_output=True,
        )
        _forward_backend_logging(stdout=completed.stdout, stderr=completed.stderr)
        if manifest_path is not None:
            _write_artifact_manifest(
                manifest_path,
                backend=self.spec.name,
                input_path=input_path,
                output_path=output_path,
                input_handle=input_handle,
                output_handle=output_handle,
                status="completed" if completed.returncode == 0 else "failed",
                return_code=completed.returncode,
            )
        if completed.returncode != 0:
            detail = _backend_failure_detail(stdout=completed.stdout, stderr=completed.stderr)
            raise RuntimeError(f"calcPhiASE failed with return code {completed.returncode}{detail}")
        return read_result(output_path, expected_iteration_index=iteration_index)

    def _write_handles_and_manifest(self, *, status, return_code=None):
        if self._manifest_path is None:
            return
        _write_openpmd_handle(self._input_handle, self._input_path)
        _write_openpmd_handle(self._output_handle, self._output_path)
        _write_artifact_manifest(
            self._manifest_path,
            backend=self.spec.name,
            input_path=self._input_path,
            output_path=self._output_path,
            input_handle=self._input_handle,
            output_handle=self._output_handle,
            status=status,
            return_code=return_code,
        )

    def _start_streaming_backend(self):
        self._result_reader_done = False
        self._result_queue = queue.Queue()
        self._send_queue = queue.Queue()
        self._sender_errors = queue.Queue()
        self._watchdog_events = queue.Queue()
        self._watchdog_stop = threading.Event()
        self._proc = subprocess.Popen(
            self._calc_phi_ase_command(self._input_path, self._output_path),
            stderr=subprocess.PIPE,
            text=True,
        )
        self._start_result_reader()
        self._start_input_sender()
        self._start_watchdog()

    def _start_input_sender(self):
        if self._sender is not None:
            return
        self._sender = threading.Thread(
            target=self._send_streaming_inputs,
            name="HASE openPMD input sender",
            daemon=True,
        )
        self._sender.start()

    def _start_result_reader(self):
        if self._reader is not None:
            return
        self._reader = threading.Thread(
            target=self._read_streaming_results,
            name="HASE openPMD result receiver",
            daemon=True,
        )
        self._reader.start()

    def _start_watchdog(self):
        if self._watchdog is not None or self.watchdog_interval is None:
            return
        self._watchdog = threading.Thread(
            target=self._watch_streaming_backend,
            name="HASE openPMD backend watchdog",
            daemon=True,
        )
        self._watchdog.start()

    def _watch_streaming_backend(self):
        try:
            while self._watchdog_stop is not None and not self._watchdog_stop.wait(self.watchdog_interval):
                proc = self._proc
                if proc is None:
                    return
                return_code = proc.poll()
                if return_code is not None:
                    self._watchdog_events.put((False, RuntimeError(f"calcPhiASE exited with return code {return_code}")))
                    return
                os.kill(proc.pid, 0)
                self._watchdog_events.put((True, None))
        except BaseException as exc:
            if self._watchdog_events is not None:
                self._watchdog_events.put((False, exc))

    def _queue_streaming_result(self, item):
        if self._result_queue is not None:
            self._result_queue.put(item)

    def _read_streaming_results(self):
        try:
            series = _io().Series(str(self._output_path), _access("read_linear"), _series_config(self._output_path))
            try:
                for fallback_index, iteration in enumerate(series.read_iterations()):
                    result = _read_result_iteration(series, iteration, fallback_index=fallback_index)
                    self._queue_streaming_result((True, result))
            finally:
                series.close()
            self._queue_streaming_result((True, _STREAMING_RESULT_EOF))
        except BaseException as exc:
            self._queue_streaming_result((False, exc))

    def _send_streaming_inputs(self):
        series = None
        try:
            series = _open_input_series(self._input_path, backend=self.spec.name)
            self._input_series = series
            while True:
                request = self._send_queue.get()
                if request is None:
                    return
                iteration_index, phiAse, gainMedium, crossSections = request
                _write_input_iteration(
                    series,
                    iteration_index,
                    phiAse,
                    gainMedium,
                    crossSections,
                    include_static=(iteration_index == 0),
                )
                series.flush()
        except BaseException as exc:
            self._sender_errors.put((None, exc))
        finally:
            if series is not None:
                try:
                    series.close()
                except BaseException as exc:
                    self._sender_errors.put((None, exc))
            self._input_series = None

    def _run_streaming_iteration(self, iteration_index, phiAse, gainMedium, crossSections):
        if self._send_queue is None:
            raise RuntimeError("openPMD input sender thread is not running")
        self._send_queue.put((iteration_index, phiAse, gainMedium, crossSections))
        return self._wait_for_result(iteration_index)

    def _raise_if_streaming_finished_without_result(self, expected_iteration_index):
        pending = ""
        if self._pending_results:
            pending = f"; buffered iterations: {sorted(self._pending_results)}"
        if self._result_reader_done:
            raise RuntimeError(
                f"Expected result iteration {expected_iteration_index} was not received "
                f"before the result stream ended{pending}"
            )
        if self._reader is not None and not self._reader.is_alive():
            raise RuntimeError(
                f"Expected result iteration {expected_iteration_index} was not received "
                f"before the result receiver thread stopped{pending}"
            )
        if self._proc is not None and self._proc.poll() == 0:
            raise RuntimeError(
                f"calcPhiASE completed before result iteration {expected_iteration_index} was received{pending}"
            )

    def _wait_for_result(self, expected_iteration_index):
        if expected_iteration_index in self._pending_results:
            return self._pending_results.pop(expected_iteration_index)

        while True:
            sender_error = self._pop_sender_error()
            if sender_error is not None:
                sender_iteration, exc = sender_error
                suffix = "" if sender_iteration is None else f" for iteration {sender_iteration}"
                raise RuntimeError(f"openPMD input sender failed{suffix}") from exc

            watchdog_error = self._pop_watchdog_error()
            if watchdog_error is not None:
                stderr = ""
                if self._proc is not None and self._proc.poll() is not None and self._proc.stderr is not None:
                    stderr = self._proc.stderr.read()
                detail = _backend_failure_detail(stderr=stderr)
                raise RuntimeError(f"openPMD backend watchdog failed{detail}") from watchdog_error

            if self._proc is not None and self._proc.poll() not in (None, 0):
                stderr = "" if self._proc.stderr is None else self._proc.stderr.read()
                detail = _backend_failure_detail(stderr=stderr)
                raise RuntimeError(f"calcPhiASE failed with return code {self._proc.returncode}{detail}")

            try:
                ok, payload = self._result_queue.get(timeout=_STREAMING_RESULT_POLL_SECONDS)
            except queue.Empty:
                self._raise_if_streaming_finished_without_result(expected_iteration_index)
                continue
            if not ok:
                raise payload
            if payload is _STREAMING_RESULT_EOF:
                self._result_reader_done = True
                self._raise_if_streaming_finished_without_result(expected_iteration_index)
                continue
            iteration_index, result = payload
            if iteration_index is None or iteration_index == expected_iteration_index:
                return result
            self._pending_results[iteration_index] = result


def _runOpenPmdAndExecuteHaseBinary(
    phiAse,
    gainMedium,
    crossSections,
    *,
    transport=None,
    command_prefix=None,
    workspace_dir=None,
    watchdog_interval=None,
    openpmdSession=None,
):
    if openpmdSession is not None:
        return openpmdSession.run(phiAse, gainMedium, crossSections)

    kwargs = {"transport": transport}
    if command_prefix is not None:
        kwargs["command_prefix"] = command_prefix
    if workspace_dir is not None:
        kwargs["workspace_dir"] = workspace_dir
    if watchdog_interval is not None:
        kwargs["watchdog_interval"] = watchdog_interval
    with OpenPmdPhiAseSession(**kwargs) as session:
        return session.run(phiAse, gainMedium, crossSections)

def runPhiASE(
    phiAse,
    gainMedium,
    crossSections,
    *,
    transport=None,
    command_prefix=None,
    workspace_dir=None,
    watchdog_interval=None,
    openpmdSession=None,
):
    return _runOpenPmdAndExecuteHaseBinary(
        phiAse,
        gainMedium,
        crossSections,
        transport=transport,
        command_prefix=command_prefix,
        workspace_dir=workspace_dir,
        watchdog_interval=watchdog_interval,
        openpmdSession=openpmdSession,
    )


_TIME_INTEGRATORS = {
    "explicit-euler",
    "heun",
    "midpoint",
    "runge-kutta-4",
    "frozen-phi-ase-runge-kutta-4",
    "implicit-euler",
    "exponential-euler",
}


def _time_integrator_name(solver):
    if isinstance(solver, str):
        name = solver
    else:
        name = getattr(solver, "name", None)
    if name not in _TIME_INTEGRATORS:
        raise ValueError(
            "compiled Simulation supports time integrators: "
            + ", ".join(sorted(_TIME_INTEGRATORS))
        )
    return name


def _simulation_run_control(simulation, *, steps, pumpSteps):
    pump = simulation.pump
    if pump.solver is not None:
        raise ValueError(
            "compiled Simulation does not support custom Python pump solvers; "
            "configure the built-in one-dimensional-z-traversal pump instead"
        )
    pump_dict = pump.toDict(timeFrame=simulation.timeStep)
    if pumpSteps is None:
        pumpSteps = pump.getProperty("pumpSteps")
    pump_steps_value = (2**32 - 1) if pumpSteps is None else int(pumpSteps)
    if pump_steps_value < 0:
        raise ValueError("pumpSteps must be non-negative")
    solver = simulation.timeIntegrationSolver
    control = {
        "time_step": float(simulation.timeStep),
        "number_of_steps": int(steps),
        "enable_ase": bool(getattr(simulation, "enableASE", True)),
        "pump_steps": pump_steps_value,
        "time_integrator": _time_integrator_name(solver),
        "pump_routine": "one-dimensional-z-traversal",
        "pump_intensity": float(pump.intensity),
        "pump_wavelength": float(pump.wavelength),
        "pump_radius_x": float(pump.radiusX),
        "pump_radius_y": float(pump.radiusY),
        "pump_exponent": float(pump.exponent),
        "pump_duration": float(pump.activeDuration(simulation.timeStep)),
        "pump_substeps": int(pump.pumpSubsteps),
        "pump_sigma_absorption": float(np.asarray(pump_dict["s_abs"], dtype=np.float64).reshape(-1)[0]),
        "pump_sigma_emission": float(np.asarray(pump_dict["s_ems"], dtype=np.float64).reshape(-1)[0]),
        "pump_back_reflection": bool(pump.backReflection),
        "pump_reflectivity": float(pump.reflectivity),
        "pump_extraction": bool(pump.extraction),
        "pump_temporary_fluorescence": 0.0 if pump.temporaryFluorescence is None else float(pump.temporaryFluorescence),
    }
    if hasattr(solver, "iterations"):
        control["implicit_iterations"] = int(solver.iterations)
    if hasattr(solver, "tolerance"):
        control["implicit_tolerance"] = float(solver.tolerance)
    return control


def _write_simulation_input(input_path, spec, simulation, run_control):
    with OpenPmdInputSeries(input_path, backend=spec.name) as writer:
        writer.write(
            simulation.phiASE,
            simulation.gainMedium,
            simulation.crossSections,
            iteration_index=0,
            include_static=True,
            runControl=run_control,
        )


def _run_streaming_simulation(command, input_path, output_path, spec, simulation, run_control):
    """Run compiled simulation while a Python receiver thread drains snapshots."""
    result_queue = queue.Queue(maxsize=1)

    def read_output():
        try:
            result_queue.put((True, read_simulation_output(output_path)))
        except BaseException as exc:
            result_queue.put((False, exc))

    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    reader = threading.Thread(
        target=read_output,
        name="HASE compiled simulation snapshot receiver",
        daemon=True,
    )
    reader.start()

    input_error = None
    try:
        _write_simulation_input(input_path, spec, simulation, run_control)
    except BaseException as exc:
        input_error = exc
        if proc.poll() is None:
            proc.kill()

    stdout, stderr = proc.communicate()
    _forward_backend_logging(stdout=stdout, stderr=stderr)

    timeout = _streaming_thread_join_timeout()
    reader.join(timeout=timeout)
    if reader.is_alive():
        raise RuntimeError(f"openPMD simulation snapshot receiver thread did not stop within {timeout:g} seconds")

    if input_error is not None:
        raise RuntimeError("openPMD simulation input writer failed") from input_error
    if proc.returncode != 0:
        detail = _backend_failure_detail(stdout=stdout, stderr=stderr)
        raise RuntimeError(f"calcPhiASE failed with return code {proc.returncode}{detail}")

    try:
        ok, payload = result_queue.get_nowait()
    except queue.Empty as exc:
        raise RuntimeError("openPMD simulation snapshot receiver stopped without returning snapshots") from exc
    if not ok:
        raise payload
    return payload


def runSimulation(simulation, *, steps, pumpSteps=None, transport=None, command_prefix=None, workspace_dir=None):
    """Run the full time-stepped Simulation in the compiled C++/alpaka backend."""
    steps = int(steps)
    if steps <= 0:
        raise ValueError("steps must be positive")
    spec = _backend_spec(transport)
    _ensure_backend_available(spec.name)
    executable = findCalcPhiAse()
    run_control = _simulation_run_control(simulation, steps=steps, pumpSteps=pumpSteps)

    artifact_root = _artifact_root()
    workspace_context = (
        tempfile.TemporaryDirectory(prefix="hase-openpmd-simulation-", dir=workspace_dir)
        if artifact_root is None
        else contextlib.nullcontext(artifact_root)
    )
    with workspace_context as tmp:
        tmp_path = Path(tmp)
        tmp_path.mkdir(parents=True, exist_ok=True)
        stem = _artifact_run_id() if artifact_root is not None else "simulation"
        input_path = tmp_path / f"{stem}-input{spec.suffix}"
        output_path = tmp_path / f"{stem}-output{spec.suffix}"

        command = [
            *(command_prefix or []),
            str(executable),
            f"--input-path={input_path}",
            f"--output-path={output_path}",
            "--run-simulation",
        ]

        if spec.streaming:
            return _run_streaming_simulation(command, input_path, output_path, spec, simulation, run_control)

        _write_simulation_input(input_path, spec, simulation, run_control)
        completed = subprocess.run(command, check=False, text=True, capture_output=True)
        _forward_backend_logging(stdout=completed.stdout, stderr=completed.stderr)
        if completed.returncode != 0:
            detail = _backend_failure_detail(stdout=completed.stdout, stderr=completed.stderr)
            raise RuntimeError(f"calcPhiASE failed with return code {completed.returncode}{detail}")
        return read_simulation_output(output_path)


def openStream(*, transport=None, command_prefix=None, workspace_dir=None, watchdog_interval=None):
    session = OpenPmdPhiAseSession(
        transport=transport,
        command_prefix=command_prefix,
        workspace_dir=workspace_dir,
        watchdog_interval=watchdog_interval,
    )
    return session.__enter__()


def closeStream(openpmdSession):
    if openpmdSession is None:
        return None
    return openpmdSession.__exit__(None, None, None)
