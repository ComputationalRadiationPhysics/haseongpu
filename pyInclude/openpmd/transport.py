from __future__ import annotations

import contextlib
import os
import queue
import re
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

from hase_transport import TransportComposer
from .._progress import _ProgressBar
from .._runtime import runtime_config, runtime_executable_candidates, runtime_root
from .backends import _clean_backend_names, _load_backend_names
from .graph import TRANSPORT_VERSION, attributeName, recordName, writeGraph
from ..structures import Result


_STEP_COMPLETE_RE = re.compile(r"\[HASE_STEP_COMPLETE\]\s+step=(\d+)/(\d+)(?:\s|$)")


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
    return True


def _frontend_progress_enabled():
    if not _forward_backend_logging_enabled():
        return False
    return not bool(getattr(_runtime_config(), "HASE_DEBUG_LOGGING", True))


def _backend_log_path():
    value = os.environ.get("HASE_BACKEND_LOG_FILE", "").strip()
    return None if not value else Path(value).expanduser().resolve()


def _backend_failure_detail(log_path=None):
    return "" if log_path is None else f": backend stdout/stderr log: {log_path}"


class _BackendProcess:
    """Run calcPhiASE while continuously draining its diagnostic streams."""

    def __init__(self, command, *, progress=False):
        self.log_path = _backend_log_path()
        self._log = None
        self._log_lock = threading.Lock()
        self._console_lock = threading.Lock()
        self._stream_errors = queue.Queue()
        self._threads = []
        self._finished = False
        self.backend_error = None
        if self.log_path is not None:
            self._log = self.log_path.open("a", encoding="utf-8")
        try:
            self._proc = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
            )
        except BaseException:
            if self._log is not None:
                self._log.close()
            raise
        self._progress = _ProgressBar(sys.stdout) if progress else None
        stdout = getattr(self._proc, "stdout", None)
        stderr = getattr(self._proc, "stderr", None)
        if stdout is not None:
            self._start_reader(stdout, sys.stdout, "stdout")
        if stderr is not None:
            self._start_reader(stderr, sys.stderr, "stderr")

    @property
    def returncode(self):
        return self._proc.returncode

    @property
    def pid(self):
        return self._proc.pid

    def poll(self):
        return self._proc.poll()

    def kill(self):
        return self._proc.kill()

    def _start_reader(self, source, destination, stream_name):
        thread = threading.Thread(
            target=self._drain_stream,
            args=(source, destination, stream_name),
            name=f"HASE calcPhiASE {stream_name} reader",
            daemon=True,
        )
        self._threads.append(thread)
        thread.start()

    def _drain_stream(self, source, destination, stream_name):
        forward = _forward_backend_logging_enabled()
        try:
            for line in source:
                if stream_name == "stderr" and line.startswith("calcPhiASE failed: "):
                    self.backend_error = line.removeprefix("calcPhiASE failed: ").strip()
                if forward:
                    try:
                        self._forward_line(line, destination, stream_name)
                    except BaseException as exc:
                        self._stream_errors.put((stream_name, exc))
                        forward = False
                try:
                    with self._log_lock:
                        if self._log is not None:
                            self._log.write(f"[{stream_name}] {line}")
                            self._log.flush()
                except BaseException as exc:
                    self._stream_errors.put(("debug log", exc))
                    with self._log_lock:
                        if self._log is not None:
                            self._log.close()
                            self._log = None
        except BaseException as exc:
            self._stream_errors.put((stream_name, exc))
        finally:
            source.close()

    def _forward_line(self, line, destination, stream_name):
        with self._console_lock:
            match = (
                _STEP_COMPLETE_RE.search(line)
                if self._progress is not None and stream_name == "stdout"
                else None
            )
            if match is not None:
                self._progress.update(int(match.group(1)), int(match.group(2)))
                return
            if self._progress is not None:
                self._progress.clear()
            destination.write(line)
            destination.flush()
            if self._progress is not None:
                self._progress.redraw()

    def wait(self):
        if hasattr(self._proc, "wait"):
            return_code = self._proc.wait()
        else:
            stdout, stderr = self._proc.communicate()
            self._forward_completed_text(stdout, sys.stdout, "stdout")
            self._forward_completed_text(stderr, sys.stderr, "stderr")
            return_code = self._proc.returncode
        if not self._finished:
            timeout = _streaming_thread_join_timeout()
            for thread in self._threads:
                thread.join(timeout=timeout)
                if thread.is_alive():
                    raise RuntimeError(f"{thread.name} did not stop within {timeout:g} seconds")
            if self._progress is not None:
                with self._console_lock:
                    self._progress.finish()
            if self._log is not None:
                with self._log_lock:
                    self._log.close()
                    self._log = None
            self._finished = True
            try:
                stream_name, error = self._stream_errors.get_nowait()
            except queue.Empty:
                pass
            else:
                raise RuntimeError(f"failed to forward calcPhiASE {stream_name}") from error
        return return_code

    def _forward_completed_text(self, value, destination, stream_name):
        if not value:
            return
        if stream_name == "stderr":
            for line in value.splitlines():
                if line.startswith("calcPhiASE failed: "):
                    self.backend_error = line.removeprefix("calcPhiASE failed: ").strip()
        if _forward_backend_logging_enabled():
            for line in value.splitlines(keepends=True):
                self._forward_line(line, destination, stream_name)
        with self._log_lock:
            if self._log is not None:
                for line in value.splitlines(keepends=True):
                    self._log.write(f"[{stream_name}] {line}")
                self._log.flush()


def _run_backend_process(command, *, progress=False):
    process = _BackendProcess(command, progress=progress)
    return SimpleNamespace(returncode=process.wait(), log_path=process.log_path, backend_error=process.backend_error)


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
                # Preserve every submitted iteration when a reader lags. SST
                # defines zero as an unlimited writer-side step queue.
                "QueueLimit": "0",
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
OPENPMD_BACKEND_PRIORITY = ("adios", "adios-sst", "hdf5")
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
    series.set_attribute("haseTransportVersion", TRANSPORT_VERSION)
    return series


def _write_input_iteration(
    series,
    iteration_index,
    root,
    *,
    dynamic_only=False,
    graph=None,
):
    refresh = getattr(root, "_refreshExcitationState", None)
    if refresh is not None:
        refresh()
    iteration = series.snapshots()[int(iteration_index)]
    iteration.time = float(getattr(root, "currentTime", 0.0))
    iteration.dt = float(getattr(root, "timeStep", 1.0))
    iteration.time_unit_SI = 1.0
    writeGraph(
        iteration,
        TransportComposer().compose(root) if graph is None else graph,
        _io(),
        dynamicOnly=dynamic_only,
    )
    iteration.close()


class OpenPmdInputSeries:
    """Context manager for writing HASE input iterations to one openPMD series."""

    def __init__(self, path, *, backend=None):
        self.path = Path(path)
        self.backend = backend
        self._series = None
        self._next_iteration = 0
        self._graph = None
        self._root = None

    def __enter__(self):
        if self.backend is not None:
            self.backend = _ensure_backend_available(self.backend).name
        self._series = _open_input_series(self.path, backend=self.backend)
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()
        return False

    def write(
        self,
        root,
        *,
        iteration_index=None,
        dynamic_only=False,
    ):
        if self._series is None:
            raise RuntimeError("OpenPmdInputSeries must be used as a context manager before writing")
        if self._graph is None:
            self._graph = TransportComposer().compose(root)
            self._root = root
        elif root is not self._root:
            raise ValueError("one OpenPmdInputSeries cannot change its transport root")
        index = self._next_iteration if iteration_index is None else int(iteration_index)
        _write_input_iteration(
            self._series,
            index,
            root,
            dynamic_only=dynamic_only,
            graph=self._graph,
        )
        self._next_iteration = max(self._next_iteration, index + 1)
        return index

    def close(self):
        if self._series is not None:
            self._series.close()
            self._series = None
            self._graph = None
            self._root = None




def _iteration_index(iteration, fallback=None):
    for name in ("iteration_index", "iterationIndex"):
        if hasattr(iteration, name):
            return int(getattr(iteration, name))
    return fallback


def _read_result_iteration(series, iteration, *, fallback_index=None) -> tuple[int | None, Result]:
    iteration_index = _iteration_index(iteration, fallback_index)
    values = {
        "phiAse": _loadScalar(series, iteration, recordName("phiAseResult/phiAse"), np.float32),
        "standardError": _loadScalar(series, iteration, recordName("phiAseResult/standardError"), np.float64),
        "relativeStandardError": _loadScalar(series, iteration, recordName("phiAseResult/relativeStandardError"), np.float64),
        "totalRays": _loadScalar(series, iteration, recordName("phiAseResult/totalRays"), np.uint32),
        "dndtAse": _loadScalar(series, iteration, recordName("phiAseResult/dndtAse"), np.float64),
    }
    values.update(_result_status_values(iteration, "phiAseResult"))
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


def _read_optional_scalar(series, iteration, name, dtype):
    if name not in iteration.meshes:
        return None
    return _loadScalar(series, iteration, name, dtype)


def _reshape_optional(value, shape):
    return None if value is None else value.reshape(shape, order="F")


def _has_attribute(obj, name):
    if hasattr(obj, "contains_attribute"):
        return obj.contains_attribute(name)
    try:
        obj.get_attribute(name)
        return True
    except Exception:
        return False


def _result_status_values(iteration, root="phiAseResult"):
    fields = {
        "boundaryStatus": ("disabled", str),
        "boundaryPasses": (0, int),
        "boundaryRemainingFraction": (0.0, float),
        "boundaryMaxPasses": (0, int),
        "boundaryDivergenceStreak": (0, int),
        "boundaryGamma": (0.0, float),
        "boundaryGammaStandardError": (0.0, float),
        "boundaryTailFactor": (0.0, float),
        "boundaryTailClosure": (0.0, float),
    }
    result = {}
    for field_name, (default, cast) in fields.items():
        name = attributeName(f"{root}/{field_name}")
        result[field_name] = cast(iteration.get_attribute(name)) if _has_attribute(iteration, name) else default
    return result


class _NoSimulationIterationsError(RuntimeError):
    pass


def read_simulation_output(path, *, on_state=None):
    """Read C++ time-stepped simulation snapshots from an openPMD output series."""
    path = Path(path)
    series = _io().Series(str(path), _access("read_linear"), _series_config(path))
    states = []
    for fallback_index, iteration in enumerate(series.read_iterations()):
        iteration_index = _iteration_index(iteration, fallback_index)
        beta_volume = _read_optional_scalar(series, iteration, recordName("simulationSnapshot/betaVolume"), np.float64)
        phi_ase = _read_optional_scalar(series, iteration, recordName("simulationSnapshot/phiAse"), np.float32)
        dndt_ase = _read_optional_scalar(series, iteration, recordName("simulationSnapshot/dndtAse"), np.float64)
        standard_error = _read_optional_scalar(series, iteration, recordName("simulationSnapshot/standardError"), np.float64)
        relative_standard_error = _read_optional_scalar(
            series, iteration, recordName("simulationSnapshot/relativeStandardError"), np.float64
        )
        total_rays = _read_optional_scalar(series, iteration, recordName("simulationSnapshot/totalRays"), np.uint32)
        present = next(
            (value for value in (beta_volume, phi_ase, dndt_ase, standard_error, relative_standard_error, total_rays) if value is not None),
            np.empty(0),
        )
        cell_shape = (present.size,)
        ase_result = None
        if any(value is not None for value in (phi_ase, standard_error, relative_standard_error, total_rays, dndt_ase)):
            ase_result = Result(
                phiAse=[] if phi_ase is None else phi_ase,
                standardError=[] if standard_error is None else standard_error,
                relativeStandardError=[] if relative_standard_error is None else relative_standard_error,
                totalRays=[] if total_rays is None else total_rays,
                dndtAse=[] if dndt_ase is None else dndt_ase,
                **_result_status_values(iteration, "simulationSnapshot"),
            )
        step_name = attributeName("simulationSnapshot/step")
        time_name = attributeName("simulationSnapshot/time")
        state = SimpleNamespace(
            iterationIndex=iteration_index,
            step=int(iteration.get_attribute(step_name)) if _has_attribute(iteration, step_name) else iteration_index + 1,
            time=float(iteration.get_attribute(time_name)) if _has_attribute(iteration, time_name) else float(iteration.time),
            betaVolume=_reshape_optional(beta_volume, cell_shape),
            phiAse=_reshape_optional(phi_ase, cell_shape),
            standardError=_reshape_optional(standard_error, cell_shape),
            relativeStandardError=_reshape_optional(relative_standard_error, cell_shape),
            totalRays=_reshape_optional(total_rays, cell_shape),
            dndtAse=_reshape_optional(dndt_ase, cell_shape),
            dndtPump=_reshape_optional(
                _read_optional_scalar(series, iteration, recordName("simulationSnapshot/dndtPump"), np.float64),
                cell_shape,
            ),
            aseResult=ase_result,
            staticUpdate=True,
        )
        states.append(state)
        if on_state is not None:
            on_state(state)
        iteration.close()
    series.close()
    if not states:
        raise _NoSimulationIterationsError(f"No simulation iterations were available in {path}")
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

    def run(self, root):
        if not self._entered:
            raise RuntimeError("OpenPmdPhiAseSession must be used as a context manager before running")
        iteration_index = self._next_iteration
        if self.spec.streaming:
            result = self._run_streaming_iteration(iteration_index, root)
        else:
            result = self._run_file_iteration(iteration_index, root)
        self._next_iteration += 1
        return result

    def close(self):
        if self.spec.streaming:
            self._close_streaming()
            return
        if self._proc is None:
            return

        return_code, log_path = self._finish_backend_process()
        self._write_handles_and_manifest(status="completed" if return_code == 0 else "failed", return_code=return_code)
        if return_code != 0:
            detail = _backend_failure_detail(log_path)
            raise RuntimeError(f"calcPhiASE failed with return code {return_code}{detail}")

    def _finish_backend_process(self):
        return_code = self._proc.wait()
        log_path = self._proc.log_path
        self._proc = None
        return return_code, log_path

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
        log_path = None
        if self._proc is not None:
            return_code, log_path = self._finish_backend_process()
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
            detail = _backend_failure_detail(log_path)
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

    def _run_file_iteration(self, iteration_index, root):
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
            writer.write(root, iteration_index=iteration_index)
        completed = _run_backend_process(self._calc_phi_ase_command(input_path, output_path))
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
            detail = _backend_failure_detail(completed.log_path)
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
        self._proc = _BackendProcess(self._calc_phi_ase_command(self._input_path, self._output_path))
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
                iteration_index, root = request
                _write_input_iteration(series, iteration_index, root)
        except BaseException as exc:
            self._sender_errors.put((None, exc))
        finally:
            if series is not None:
                try:
                    series.close()
                except BaseException as exc:
                    self._sender_errors.put((None, exc))
            self._input_series = None

    def _run_streaming_iteration(self, iteration_index, root):
        if self._send_queue is None:
            raise RuntimeError("openPMD input sender thread is not running")
        self._send_queue.put((iteration_index, root))
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
                log_path = None if self._proc is None else getattr(self._proc, "log_path", None)
                detail = _backend_failure_detail(log_path)
                raise RuntimeError(f"openPMD backend watchdog failed{detail}") from watchdog_error

            if self._proc is not None and self._proc.poll() not in (None, 0):
                self._proc.wait()
                detail = _backend_failure_detail(self._proc.log_path)
                if self._proc.backend_error:
                    detail = f": {self._proc.backend_error}{detail}"
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
    root,
    *,
    transport=None,
    command_prefix=None,
    workspace_dir=None,
    watchdog_interval=None,
    openpmdSession=None,
):
    if openpmdSession is not None:
        return openpmdSession.run(root)

    kwargs = {"transport": transport}
    if command_prefix is not None:
        kwargs["command_prefix"] = command_prefix
    if workspace_dir is not None:
        kwargs["workspace_dir"] = workspace_dir
    if watchdog_interval is not None:
        kwargs["watchdog_interval"] = watchdog_interval
    with OpenPmdPhiAseSession(**kwargs) as session:
        return session.run(root)

def runPhiASE(
    root,
    *,
    transport=None,
    command_prefix=None,
    workspace_dir=None,
    watchdog_interval=None,
    openpmdSession=None,
):
    return _runOpenPmdAndExecuteHaseBinary(
        root,
        transport=transport,
        command_prefix=command_prefix,
        workspace_dir=workspace_dir,
        watchdog_interval=watchdog_interval,
        openpmdSession=openpmdSession,
    )


def _write_simulation_input(input_path, spec, simulation, *, close_after=None):
    with OpenPmdInputSeries(input_path, backend=spec.name) as writer:
        writer.write(simulation, iteration_index=0)
        if close_after is not None:
            close_after.wait()


def _run_streaming_simulation(
    command,
    input_path,
    output_path,
    spec,
    simulation,
    steps,
    *,
    on_state=None,
    progress=False,
):
    """Run compiled simulation while Python threads exchange SST snapshots."""
    result_queue = queue.Queue(maxsize=1)
    input_queue = queue.Queue(maxsize=1)
    control_queue = queue.Queue(maxsize=1)
    control_ack_queue = queue.Queue(maxsize=1)
    backend_finished = threading.Event()
    synchronized_debug = getattr(simulation, "executionMode", "autonomous") == "synchronized-debug"

    def read_output():
        try:
            def receive_state(state):
                if on_state is not None:
                    on_state(state)
                if synchronized_debug and int(state.step) < int(steps):
                    completed_step = int(state.step)
                    control_queue.put(completed_step)
                    acknowledged_step = control_ack_queue.get()
                    if acknowledged_step != completed_step:
                        raise RuntimeError(
                            f"synchronized-debug control iteration {completed_step} was not published"
                        )

            while True:
                try:
                    if on_state is None and not synchronized_debug:
                        snapshots = read_simulation_output(output_path)
                    else:
                        snapshots = read_simulation_output(output_path, on_state=receive_state)
                    break
                except _NoSimulationIterationsError:
                    if proc.poll() is not None:
                        raise
                    backend_finished.wait(timeout=_STREAMING_RESULT_POLL_SECONDS)
            result_queue.put((True, snapshots))
        except BaseException as exc:
            result_queue.put((False, exc))
            if proc.poll() is None:
                proc.kill()
        finally:
            if synchronized_debug:
                try:
                    control_queue.put_nowait(None)
                except queue.Full:
                    pass

    def write_input():
        try:
            if not synchronized_debug:
                _write_simulation_input(
                    input_path,
                    spec,
                    simulation,
                    close_after=backend_finished,
                )
            else:
                with OpenPmdInputSeries(input_path, backend=spec.name) as writer:
                    writer.write(simulation, iteration_index=0)
                    for expected_step in range(1, int(steps)):
                        completed_step = control_queue.get()
                        if completed_step is None:
                            raise RuntimeError(
                                "synchronized-debug output ended before the next control exchange"
                            )
                        if completed_step != expected_step:
                            raise RuntimeError(
                                f"synchronized-debug expected completed step {expected_step}, "
                                f"received {completed_step}"
                            )
                        writer.write(
                            simulation,
                            iteration_index=completed_step,
                            dynamic_only=True,
                        )
                        control_ack_queue.put(completed_step)
                    backend_finished.wait()
            input_queue.put((True, None))
        except BaseException as exc:
            if synchronized_debug:
                try:
                    control_ack_queue.put_nowait(None)
                except queue.Full:
                    pass
            input_queue.put((False, exc))
            if proc.poll() is None:
                proc.kill()

    proc = _BackendProcess(command, progress=progress)
    reader = threading.Thread(
        target=read_output,
        name="HASE compiled simulation snapshot receiver",
        daemon=True,
    )
    reader.start()
    writer = threading.Thread(
        target=write_input,
        name="HASE compiled simulation input sender",
        daemon=True,
    )
    writer.start()

    try:
        proc.wait()
    finally:
        backend_finished.set()

    timeout = _streaming_thread_join_timeout()
    deadline = time.monotonic() + timeout
    writer.join(timeout=max(0.0, deadline - time.monotonic()))
    reader.join(timeout=max(0.0, deadline - time.monotonic()))

    input_result = None
    if not writer.is_alive():
        try:
            input_result = input_queue.get_nowait()
        except queue.Empty:
            if proc.returncode == 0:
                raise RuntimeError("openPMD simulation input sender stopped without reporting completion")
        if input_result is not None:
            input_ok, input_payload = input_result
            if not input_ok:
                raise RuntimeError("openPMD simulation input writer failed") from input_payload

    if proc.returncode != 0:
        detail = _backend_failure_detail(proc.log_path)
        if proc.backend_error:
            detail = f": {proc.backend_error}{detail}"
        raise RuntimeError(f"calcPhiASE failed with return code {proc.returncode}{detail}")
    if writer.is_alive():
        raise RuntimeError(f"openPMD simulation input sender thread did not stop within {timeout:g} seconds")
    if reader.is_alive():
        raise RuntimeError(f"openPMD simulation snapshot receiver thread did not stop within {timeout:g} seconds")

    try:
        ok, payload = result_queue.get_nowait()
    except queue.Empty as exc:
        raise RuntimeError("openPMD simulation snapshot receiver stopped without returning snapshots") from exc
    if not ok:
        raise payload
    return payload


def runSimulation(
    simulation,
    *,
    steps,
    transport=None,
    command_prefix=None,
    workspace_dir=None,
    on_state=None,
):
    """Run the full time-stepped Simulation in the compiled C++/alpaka backend."""
    steps = int(steps)
    if steps <= 0:
        raise ValueError("steps must be positive")
    executable = findCalcPhiAse()
    spec = _ensure_backend_available(transport, executable)
    progress = _frontend_progress_enabled()
    if simulation.executionMode == "synchronized-debug" and not spec.streaming:
        raise ValueError("synchronized-debug requires an openPMD streaming backend")

    artifact_root = _artifact_root()
    if artifact_root is None and workspace_dir is not None:
        Path(workspace_dir).mkdir(parents=True, exist_ok=True)
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
            "--cpp-control",
        ]

        if spec.streaming:
            return _run_streaming_simulation(
                command,
                input_path,
                output_path,
                spec,
                simulation,
                steps,
                on_state=on_state,
                progress=progress,
            )

        _write_simulation_input(input_path, spec, simulation)
        completed = _run_backend_process(command, progress=progress)
        if completed.returncode != 0:
            detail = _backend_failure_detail(completed.log_path)
            if completed.backend_error:
                detail = f": {completed.backend_error}{detail}"
            raise RuntimeError(f"calcPhiASE failed with return code {completed.returncode}{detail}")
        return read_simulation_output(output_path, on_state=on_state)


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
