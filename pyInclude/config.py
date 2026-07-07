# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Guided setup helper for HASEonGPU runtime/backend configuration."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from .openpmd import preflight


OPENPMD_BACKEND_PRIORITY = ("adios-sst", "adios", "hdf5")
PROVIDER_AUTO = "auto"
PROVIDER_BUNDLED = "bundled"
PROVIDER_SYSTEM = "system"
BUNDLED_ADIOS2_FETCH = "fetch"
BUNDLED_ADIOS2_OFF = "off"
BUNDLED_ADIOS2_SYSTEM = "system"
BUNDLED_HDF5_FETCH = "fetch"
BUNDLED_HDF5_OFF = "off"
BUNDLED_HDF5_SYSTEM = "system"
MPI_MODE_AUTO = "auto"
MPI_MODE_ON = "on"
MPI_MODE_OFF = "off"
OPENPMD_ADIOS2_AUTO = "auto"
OPENPMD_ADIOS2_YES = "yes"
OPENPMD_ADIOS2_NO = "no"
DEFAULT_PHI_ASE_CONFIG_PATH = Path("config") / "hase-phiase.yaml"
LAST_INSTALL_STATE_PATH = Path("build") / "hase-configure-last-install.json"
CCACHE_LAUNCHER_DEFINES = (
    "CMAKE_C_COMPILER_LAUNCHER",
    "CMAKE_CXX_COMPILER_LAUNCHER",
    "CMAKE_CUDA_COMPILER_LAUNCHER",
)
REINSTALL_CACHE_DEFINES = (
    "HASE_OPENPMD_PROVIDER",
    "HASE_OPENPMD_USE_ADIOS2",
    "HASE_OPENPMD_FETCH_ADIOS2",
    "HASE_OPENPMD_USE_HDF5",
    "HASE_OPENPMD_FETCH_HDF5",
    "HASE_OPENPMD_USE_SST",
    "HASE_OPENPMD_BUILD_PYTHON_BINDINGS",
    "HASE_USE_SYSTEM_ALPAKA",
    "HASE_CUDA_ARCHITECTURES",
    "DISABLE_MPI",
    "HASE_NATIVE_OPTIMIZATIONS",
    "CMAKE_PREFIX_PATH",
    "openPMD_DIR",
    "ADIOS2_DIR",
    "HDF5_DIR",
    *CCACHE_LAUNCHER_DEFINES,
)


@dataclass(frozen=True)
class WizardSelection:
    """Resolved first-run setup choices used to render commands and YAML."""

    provider: str
    openpmd_backend: str
    compute_backend: str
    supported_openpmd_backends: tuple[str, ...] = ()
    parallel_mode: str = "single"
    num_devices: int = 1
    n_per_node: int = 1
    cmake_prefix_path: str | None = None
    openpmd_dir: str | None = None
    bundled_adios2: str = BUNDLED_ADIOS2_FETCH
    bundled_hdf5: str = BUNDLED_HDF5_OFF
    openpmd_adios2: str = OPENPMD_ADIOS2_AUTO
    bundled_python_bindings: bool = True
    native_optimizations: bool = False
    mpi_mode: str = MPI_MODE_AUTO
    adios2_prefix: str | None = None
    adios2_dir: str | None = None
    hdf5_prefix: str | None = None
    hdf5_dir: str | None = None


def _csv(values):
    return ", ".join(values) if values else "<none>"


def _clean_supported_backends(backends):
    cleaned = []
    for backend in backends:
        name = str(backend).split(" ", 1)[0]
        if name in OPENPMD_BACKEND_PRIORITY and name not in cleaned:
            cleaned.append(name)
    return cleaned


def _combined_cmake_prefix_path(openpmd_prefix=None, adios2_prefix=None, hdf5_prefix=None):
    entries = [entry for entry in (hdf5_prefix, adios2_prefix, openpmd_prefix) if entry]
    return ";".join(entries) if entries else None


def supported_external_openpmd_backends(python_info, cmake_info):
    """Return runtime backends confirmed by both Python and CMake providers."""
    python_backends = set(preflight._python_supported_backends(python_info))
    cmake_backends = set(_clean_supported_backends(preflight._cmake_supported_backends(cmake_info)))
    return [backend for backend in OPENPMD_BACKEND_PRIORITY if backend in python_backends and backend in cmake_backends]


def preferred_openpmd_backend(supported, requested=None):
    """Choose a runtime openPMD backend from confirmed support."""
    supported = list(supported)
    if requested is not None:
        normalized = requested.strip().lower()
        if normalized not in OPENPMD_BACKEND_PRIORITY:
            raise ValueError(f"unsupported openPMD backend '{requested}'")
        if normalized not in supported:
            raise ValueError(
                f"openPMD backend '{requested}' was requested but supported backends are: {_csv(supported)}"
            )
        return normalized
    for backend in OPENPMD_BACKEND_PRIORITY:
        if backend in supported:
            return backend
    raise ValueError("no compatible openPMD runtime backend was confirmed")


def bundled_supported_openpmd_backends(adios2_mode, hdf5_mode=BUNDLED_HDF5_OFF):
    """Return runtime backends implied by the HASE-managed provider choices."""
    supported = []
    if adios2_mode != BUNDLED_ADIOS2_OFF:
        supported.extend(["adios-sst", "adios"])
    if hdf5_mode != BUNDLED_HDF5_OFF:
        supported.append("hdf5")
    return [backend for backend in OPENPMD_BACKEND_PRIORITY if backend in supported]


def available_alpaka_backends():
    """Return runtime Alpaka backends, or an explanatory error string."""
    try:
        from .alpakaUtils import AlpakaBackends

        return AlpakaBackends.all(), None
    except Exception as exc:  # pragma: no cover - exact failure is platform/build dependent
        return [], str(exc)


def preferred_compute_backend(backends, requested=None):
    """Choose a compute backend, preferring CPU host backends for first validation."""
    backends = list(backends)
    if requested:
        if backends and requested not in backends:
            raise ValueError(f"Alpaka backend '{requested}' was requested but available backends are: {_csv(backends)}")
        return requested
    for marker in ("Host_Cpu_CpuSerial", "CpuSerial"):
        for backend in backends:
            if marker in backend:
                return backend
    if backends:
        return backends[0]
    return "Host_Cpu_CpuSerial"


def _set_cmake_define(args, name, value):
    prefix = f"-D{name}="
    filtered = [arg for arg in args if not arg.startswith(prefix)]
    filtered.append(f"{prefix}{value}")
    return filtered


def _apply_ccache_args(args, use_ccache=False):
    if not use_ccache:
        return list(args)
    result = list(args)
    for define in CCACHE_LAUNCHER_DEFINES:
        result = _set_cmake_define(result, define, "ccache")
    return result


def cmake_args(selection: WizardSelection, *, use_ccache=False):
    """Return CMake arguments for the selected provider/install path."""
    args = []
    if selection.provider == PROVIDER_BUNDLED:
        args.append("-DHASE_OPENPMD_PROVIDER=bundled")
        if selection.bundled_python_bindings:
            args.append("-DHASE_OPENPMD_BUILD_PYTHON_BINDINGS=ON")
        if selection.bundled_adios2 == BUNDLED_ADIOS2_OFF:
            args.extend(["-DHASE_OPENPMD_USE_ADIOS2=OFF", "-DHASE_OPENPMD_USE_SST=OFF"])
        elif selection.bundled_adios2 == BUNDLED_ADIOS2_FETCH:
            args.extend(["-DHASE_OPENPMD_USE_ADIOS2=ON", "-DHASE_OPENPMD_FETCH_ADIOS2=ON"])
        elif selection.bundled_adios2 == BUNDLED_ADIOS2_SYSTEM:
            args.extend(["-DHASE_OPENPMD_USE_ADIOS2=ON", "-DHASE_OPENPMD_FETCH_ADIOS2=OFF"])
        if selection.bundled_hdf5 == BUNDLED_HDF5_OFF:
            args.append("-DHASE_OPENPMD_USE_HDF5=OFF")
        elif selection.bundled_hdf5 == BUNDLED_HDF5_FETCH:
            args.extend(["-DHASE_OPENPMD_USE_HDF5=ON", "-DHASE_OPENPMD_FETCH_HDF5=ON"])
        elif selection.bundled_hdf5 == BUNDLED_HDF5_SYSTEM:
            args.extend(["-DHASE_OPENPMD_USE_HDF5=ON", "-DHASE_OPENPMD_FETCH_HDF5=OFF"])
        prefix_path = _combined_cmake_prefix_path(
            adios2_prefix=selection.adios2_prefix,
            hdf5_prefix=selection.hdf5_prefix,
        )
        if prefix_path:
            args.append(f"-DCMAKE_PREFIX_PATH={prefix_path}")
        if selection.adios2_dir:
            args.append(f"-DADIOS2_DIR={selection.adios2_dir}")
        if selection.hdf5_dir:
            args.append(f"-DHDF5_DIR={selection.hdf5_dir}")
    elif selection.provider == PROVIDER_AUTO:
        args.append("-DHASE_OPENPMD_PROVIDER=auto")
    else:
        args.append("-DHASE_OPENPMD_PROVIDER=system")
        prefix_path = _combined_cmake_prefix_path(
            selection.cmake_prefix_path,
            selection.adios2_prefix,
            selection.hdf5_prefix,
        )
        if prefix_path:
            args.append(f"-DCMAKE_PREFIX_PATH={prefix_path}")
        if selection.openpmd_dir:
            args.append(f"-DopenPMD_DIR={selection.openpmd_dir}")
        if selection.adios2_dir:
            args.append(f"-DADIOS2_DIR={selection.adios2_dir}")
        if selection.hdf5_dir:
            args.append(f"-DHDF5_DIR={selection.hdf5_dir}")

    disable_mpi = {
        MPI_MODE_ON: "OFF",
        MPI_MODE_OFF: "ON",
        MPI_MODE_AUTO: "AUTO",
    }[selection.mpi_mode]
    args.append(f"-DDISABLE_MPI={disable_mpi}")
    args.append(f"-DHASE_NATIVE_OPTIMIZATIONS={'ON' if selection.native_optimizations else 'OFF'}")
    return _apply_ccache_args(args, use_ccache=use_ccache)


def _cmake_args_string(args):
    return " ".join(args)


def cmake_args_string(selection: WizardSelection, *, use_ccache=False):
    """Render CMake arguments for CMAKE_ARGS."""
    return _cmake_args_string(cmake_args(selection, use_ccache=use_ccache))


def _pip_install_args(*, break_system_packages=False):
    args = ["-m", "pip", "install", "-v"]
    if break_system_packages:
        args.append("--break-system-packages")
    args.append(".")
    return args


def _install_command_from_cmake_args(cmake_arg_list, *, break_system_packages=False):
    args = _cmake_args_string(cmake_arg_list)
    pip_args = " ".join(_pip_install_args(break_system_packages=break_system_packages))
    return (
        f'HASE_CONFIGURE_CMAKE_ARGS="{args}"\n'
        f'CMAKE_ARGS="$HASE_CONFIGURE_CMAKE_ARGS" python3 {pip_args}'
    )


def install_command(selection: WizardSelection, *, break_system_packages=False, use_ccache=False):
    """Render the recommended pip install command for the selected provider path."""
    return _install_command_from_cmake_args(
        cmake_args(selection, use_ccache=use_ccache),
        break_system_packages=break_system_packages,
    )


def _write_last_install_state(cmake_arg_list):
    LAST_INSTALL_STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    LAST_INSTALL_STATE_PATH.write_text(
        json.dumps({"version": 1, "cmake_args": list(cmake_arg_list)}, indent=2) + "\n",
        encoding="utf-8",
    )


def _read_last_install_state():
    if not LAST_INSTALL_STATE_PATH.is_file():
        return None
    state = json.loads(LAST_INSTALL_STATE_PATH.read_text(encoding="utf-8"))
    cmake_arg_list = state.get("cmake_args")
    if not isinstance(cmake_arg_list, list) or not all(isinstance(arg, str) for arg in cmake_arg_list):
        raise RuntimeError(f"invalid reinstall state in {LAST_INSTALL_STATE_PATH}")
    return cmake_arg_list


def _latest_cmake_cache():
    caches = [path for path in Path("build").glob("*/CMakeCache.txt") if path.is_file()]
    if not caches:
        return None
    return max(caches, key=lambda path: path.stat().st_mtime)


def _cmake_args_from_cache(cache_path):
    values = {}
    for line in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        key_type, value = line.split("=", 1)
        key = key_type.split(":", 1)[0]
        if key in REINSTALL_CACHE_DEFINES and value and not value.endswith("-NOTFOUND"):
            values[key] = value
    return [f"-D{name}={values[name]}" for name in REINSTALL_CACHE_DEFINES if name in values]


def previous_install_cmake_args():
    """Return CMake arguments from the last script install or CMake cache."""
    cmake_arg_list = _read_last_install_state()
    if cmake_arg_list is not None:
        return cmake_arg_list

    cache_path = _latest_cmake_cache()
    if cache_path is None:
        raise RuntimeError(
            "--reinstall requires a previous HASE CMake configure/install; "
            "no build/*/CMakeCache.txt or hase-configure install state was found"
        )
    cmake_arg_list = _cmake_args_from_cache(cache_path)
    if not cmake_arg_list:
        raise RuntimeError(f"could not recover HASE install settings from {cache_path}")
    return cmake_arg_list


def _check_ccache_available(cmake_arg_list):
    if any(arg in {f"-D{name}=ccache" for name in CCACHE_LAUNCHER_DEFINES} for arg in cmake_arg_list):
        if shutil.which("ccache") is None:
            raise RuntimeError("ccache compiler launcher was requested, but 'ccache' is not on PATH")


def _run_install_with_cmake_args(cmake_arg_list, *, break_system_packages=False, record_state=False):
    _check_ccache_available(cmake_arg_list)
    env = os.environ.copy()
    env["HASE_CONFIGURE_CMAKE_ARGS"] = _cmake_args_string(cmake_arg_list)
    env["CMAKE_ARGS"] = env["HASE_CONFIGURE_CMAKE_ARGS"]
    returncode = subprocess.run(
        [sys.executable, *_pip_install_args(break_system_packages=break_system_packages)],
        env=env,
        check=False,
    ).returncode
    if record_state and returncode == 0:
        _write_last_install_state(cmake_arg_list)
    return returncode


def run_install(selection: WizardSelection, *, break_system_packages=False, use_ccache=False):
    """Run pip install with the selected CMake arguments."""
    return _run_install_with_cmake_args(
        cmake_args(selection, use_ccache=use_ccache),
        break_system_packages=break_system_packages,
        record_state=True,
    )


def run_reinstall(*, break_system_packages=False, use_ccache=False):
    """Run pip install with settings from the previous HASE configure/install."""
    cmake_arg_list = _apply_ccache_args(previous_install_cmake_args(), use_ccache=use_ccache)
    print("Reinstall command:")
    print(_install_command_from_cmake_args(cmake_arg_list, break_system_packages=break_system_packages))
    return _run_install_with_cmake_args(
        cmake_arg_list,
        break_system_packages=break_system_packages,
        record_state=True,
    )


def yaml_config(selection: WizardSelection):
    """Render a PhiASE-compatible YAML run-control snippet."""
    lines = [
        "compute:",
        f"  backend: {selection.compute_backend}",
        f"  openpmd_backend: {selection.openpmd_backend}",
        f"  parallel_mode: {selection.parallel_mode}",
        f"  numDevices: {int(selection.num_devices)}",
        f"  n_per_node: {int(selection.n_per_node)}",
        "  write_vtk: false",
    ]
    return "\n".join(lines) + "\n"


def _format_alpaka_backends(alpaka_backends):
    if not alpaka_backends:
        return "\n\n    <none reported>"
    return "\n\n" + "\n".join(f"    - {backend}" for backend in alpaka_backends)


def alpaka_backend_guidance(selection: WizardSelection, *, alpaka_backends=(), alpaka_error=None):
    """Render backend-selection guidance grounded in HASE's CMake detection."""
    text = (
        "alpaka backend guidance:\n\n"
        "    backend is the alpaka compute backend, while openpmd_backend is the "
        "storage/streaming backend.\n\n"
        "    Currently available alpaka backends are the hardware configurations "
        "HASEonGPU can run on:"
        f"{_format_alpaka_backends(alpaka_backends)}"
    )
    backend = selection.compute_backend.lower()
    available_text = " ".join(alpaka_backends).lower()
    has_gpu_backend = "cuda" in available_text or "hip" in available_text or "gpu" in available_text

    if "cuda" in backend:
        text += (
            "\n\n    You selected a CUDA backend. If it does not work, make sure CMake/alpaka "
            "can find CUDA via CUDACXX or nvcc/CUDAToolkit_ROOT/CUDA_HOME/CUDA_PATH/CUDA_ROOT, "
            "and that an NVIDIA device is visible at runtime."
        )
    elif "hip" in backend:
        text += (
            "\n\n    You selected a HIP backend. If it does not work, make sure CMake/alpaka "
            "can find ROCm/HIP via HIPCXX, hipcc, or ROCM_PATH, and that an AMD device "
            "is visible at runtime."
        )
    elif not has_gpu_backend:
        text += (
            "\n\n    No GPU alpaka backend is listed right now; that usually means alpaka did "
            "not find a usable GPU toolchain/backend, or no matching device is available."
            "\n\n    To get NVIDIA/CUDA backends, provide CUDA through CUDACXX or "
            "nvcc/CUDAToolkit_ROOT/CUDA_HOME/CUDA_PATH/CUDA_ROOT and make an NVIDIA device "
            "visible.\n\n    To get AMD/HIP backends, provide ROCm/HIP through "
            "HIPCXX, hipcc, or ROCM_PATH and make an AMD device visible."
        )
    else:
        text += "\n\n    Use a CPU host backend for first validation, or select one of the listed GPU backends."

    if alpaka_error:
        text += " Could not query installed alpaka backends: " + alpaka_error
    return text


def guidance_items(selection: WizardSelection, yaml_path, *, alpaka_backends=(), alpaka_error=None):
    """Render concise final guidance for users."""
    if yaml_path == "<stdout>":
        yaml_text = (
            "The configuration file was printed to <stdout>. Write it to a file with "
            "--output and modify it afterwards."
        )
    else:
        yaml_text = (
            f"The configuration file is present under {yaml_path} and can be modified afterwards. "
            "Edit compute.backend to change the alpaka compute backend, "
            "compute.openpmd_backend to change the openPMD storage/streaming backend, "
            "and compute.parallel_mode / compute.n_per_node for MPI runs."
        )
    supported = _csv(list(selection.supported_openpmd_backends) or [selection.openpmd_backend])
    if selection.mpi_mode == MPI_MODE_ON:
        mpi_text = "MPI build support is enabled (-DDISABLE_MPI=OFF)."
    elif selection.mpi_mode == MPI_MODE_AUTO:
        mpi_text = "MPI build support is auto-detected (-DDISABLE_MPI=AUTO)."
    else:
        mpi_text = "MPI build support is disabled (-DDISABLE_MPI=ON)."
    provider_text = (
        "Bundled provider installs openPMD/ADIOS2/HDF5 into a build-local CMake prefix "
        "and reuses it on later installs when dependency settings are unchanged."
        if selection.provider == PROVIDER_BUNDLED
        else "System provider uses an existing openPMD-api C++ installation."
    )
    return [
        yaml_text,
        f"Supported openpmd_backends for this choice: {supported}.",
        mpi_text,
        provider_text,
        "Switch YAML parallel_mode to mpi and set n_per_node when running multi-rank jobs.",
        alpaka_backend_guidance(selection, alpaka_backends=alpaka_backends, alpaka_error=alpaka_error),
    ]


def provider_notes(selection: WizardSelection, *, alpaka_error=None):
    """Backward-compatible concise guidance helper."""
    return guidance_items(selection, "<configured YAML path>", alpaka_error=alpaka_error)


def _ask(prompt, default=None):
    suffix = f" [{default}]" if default not in (None, "") else ""
    value = input(f"{prompt}{suffix}: ").strip()
    return default if value == "" else value


def _ask_yes_no(prompt, default=False):
    default_text = "y" if default else "n"
    value = _ask(prompt + " (y/n)", default_text)
    return str(value).strip().lower() in {"y", "yes", "1", "true", "on"}


def _interactive_provider():
    print("Choose an openPMD provider setup:")
    print("  1) auto      use an existing provider if detected, otherwise build a bundled provider")
    print("  2) bundled   fetch/build/install openPMD dependencies with HASEonGPU")
    print("  3) system    use an existing openPMD-api C++ provider")
    choice = _ask("Provider choice", "1")
    return {"1": PROVIDER_AUTO, "2": PROVIDER_BUNDLED, "3": PROVIDER_SYSTEM}.get(str(choice).strip(), choice)


def _interactive_bundled_adios2():
    print("Choose ADIOS2 handling for the bundled openPMD provider:")
    print("  1) fetch   fetch/build pinned ADIOS2; enables adios-sst/adios")
    print("  2) off     do not use ADIOS2; hdf5-only")
    print("  3) system  use an existing ADIOS2 installation")
    choice = _ask("Bundled ADIOS2 choice", "1")
    return {
        "1": BUNDLED_ADIOS2_FETCH,
        "2": BUNDLED_ADIOS2_OFF,
        "3": BUNDLED_ADIOS2_SYSTEM,
    }.get(str(choice).strip(), choice)


def _interactive_bundled_hdf5():
    print("Choose HDF5 handling for the bundled openPMD provider:")
    print("  1) fetch   fetch/build pinned HDF5")
    print("  2) off     do not use HDF5; ADIOS2-only")
    print("  3) system  use an existing HDF5 installation")
    choice = _ask("Bundled HDF5 choice", "2")
    return {
        "1": BUNDLED_HDF5_FETCH,
        "2": BUNDLED_HDF5_OFF,
        "3": BUNDLED_HDF5_SYSTEM,
    }.get(str(choice).strip(), choice)


def _probe_external(args, *, backend=None):
    probe_args = SimpleNamespace(
        backend=backend or args.openpmd_backend or "adios-sst",
        cmake_prefix_path=_combined_cmake_prefix_path(args.cmake_prefix_path, args.adios2_prefix, args.hdf5_prefix),
        openpmd_dir=args.openpmd_dir,
        cmake=args.cmake,
        cmake_generator=args.cmake_generator,
    )
    errors: list[str] = []
    warnings: list[str] = []
    python_info = preflight._python_info(errors)
    cmake_info = preflight._cmake_probe(probe_args, errors)
    preflight._check_versions(python_info, cmake_info, warnings)
    return probe_args, python_info, cmake_info, errors, warnings


def _select_bundled_openpmd_backend(args, bundled_adios2, bundled_hdf5):
    supported = bundled_supported_openpmd_backends(bundled_adios2, bundled_hdf5)
    return preferred_openpmd_backend(supported, args.openpmd_backend)


def _external_openpmd_backend_result(args, openpmd_adios2):
    if openpmd_adios2 == OPENPMD_ADIOS2_NO:
        requested_backend = args.openpmd_backend or "hdf5"
    else:
        requested_backend = args.openpmd_backend
    probe_args, python_info, cmake_info, errors, warnings = _probe_external(args, backend=requested_backend)
    selected_backend = None
    supported = []
    if not errors:
        supported = supported_external_openpmd_backends(python_info, cmake_info)
        if openpmd_adios2 == OPENPMD_ADIOS2_NO:
            supported = [backend for backend in supported if backend == "hdf5"]
        selected_backend = preferred_openpmd_backend(supported, requested_backend)
        probe_args.backend = selected_backend
        preflight._check_python_backend(python_info, selected_backend, errors)
        preflight._check_cmake_backend(cmake_info, selected_backend, errors, warnings)
    return probe_args, selected_backend, tuple(supported), (python_info, cmake_info, warnings), errors, warnings


def _select_external_openpmd_backend(args, openpmd_adios2):
    probe_args, selected_backend, supported, probe_report, errors, warnings = _external_openpmd_backend_result(
        args,
        openpmd_adios2,
    )
    python_info, cmake_info, _probe_warnings = probe_report
    preflight.print_report(probe_args, python_info, cmake_info, errors, warnings)
    if errors:
        raise RuntimeError("openPMD provider preflight failed; see errors above")
    return selected_backend, supported, probe_report


def _try_external_openpmd_backend(args, openpmd_adios2):
    try:
        _probe_args, selected_backend, supported, probe_report, errors, _warnings = _external_openpmd_backend_result(
            args,
            openpmd_adios2,
        )
    except ValueError as exc:
        return None, f"system openPMD provider does not support the requested backend: {exc}"
    if errors:
        return None, "system openPMD provider probe failed: " + "; ".join(errors)
    return (selected_backend, supported, probe_report), None


def _build_selection(args):
    provider = args.provider
    if provider == PROVIDER_AUTO and args.interactive:
        provider = _interactive_provider()

    cmake_prefix_path = args.cmake_prefix_path
    openpmd_dir = args.openpmd_dir
    adios2_prefix = args.adios2_prefix
    adios2_dir = args.adios2_dir
    hdf5_prefix = args.hdf5_prefix
    hdf5_dir = args.hdf5_dir
    bundled_adios2 = args.bundled_adios2
    bundled_hdf5 = args.bundled_hdf5
    openpmd_adios2 = args.openpmd_adios2
    bundled_python_bindings = not args.no_bundled_python_bindings
    native_optimizations = args.native_optimizations == "on"
    mpi_mode = args.mpi
    probe_report = None

    if provider == PROVIDER_AUTO:
        probe_args = argparse.Namespace(**vars(args))
        probe_args.cmake_prefix_path = cmake_prefix_path
        probe_args.openpmd_dir = openpmd_dir
        probe_args.adios2_prefix = adios2_prefix
        probe_args.adios2_dir = adios2_dir
        probe_args.hdf5_prefix = hdf5_prefix
        probe_args.hdf5_dir = hdf5_dir
        external_result, fallback_reason = _try_external_openpmd_backend(probe_args, openpmd_adios2)
        if external_result is None:
            print(
                "hase-configure: auto provider falling back to bundled because "
                f"{fallback_reason}",
                file=sys.stderr,
            )
            provider = PROVIDER_BUNDLED
        else:
            provider = PROVIDER_SYSTEM
            selected_openpmd_backend, supported_openpmd_backends, probe_report = external_result

    if provider == PROVIDER_BUNDLED:
        if args.interactive:
            bundled_adios2 = _interactive_bundled_adios2()
            if bundled_adios2 == BUNDLED_ADIOS2_SYSTEM:
                adios2_prefix = _ask("ADIOS2 prefix for CMAKE_PREFIX_PATH (optional)", adios2_prefix or "") or None
                adios2_dir = _ask("ADIOS2_DIR CMake config directory (optional)", adios2_dir or "") or None
            bundled_hdf5 = _interactive_bundled_hdf5()
            if bundled_hdf5 == BUNDLED_HDF5_SYSTEM:
                hdf5_prefix = _ask("HDF5 prefix for CMAKE_PREFIX_PATH (optional)", hdf5_prefix or "") or None
                hdf5_dir = _ask("HDF5_DIR CMake config directory (optional)", hdf5_dir or "") or None
            bundled_python_bindings = _ask_yes_no(
                "Build matching bundled openPMD Python bindings",
                bundled_python_bindings,
            )
        supported_openpmd_backends = tuple(bundled_supported_openpmd_backends(bundled_adios2, bundled_hdf5))
        selected_openpmd_backend = _select_bundled_openpmd_backend(args, bundled_adios2, bundled_hdf5)
    elif provider == PROVIDER_SYSTEM and probe_report is None:
        if args.interactive:
            cmake_prefix_path = _ask("openPMD CMake prefix path", cmake_prefix_path or "") or None
            openpmd_dir = _ask("openPMD_DIR config directory (optional)", openpmd_dir or "") or None
            has_adios2 = _ask_yes_no("Was this openPMD provider built with ADIOS2 support", True)
            openpmd_adios2 = OPENPMD_ADIOS2_YES if has_adios2 else OPENPMD_ADIOS2_NO
            if has_adios2 and _ask_yes_no("Add an ADIOS2 prefix or ADIOS2_DIR for CMake discovery", False):
                adios2_prefix = _ask("ADIOS2 prefix for CMAKE_PREFIX_PATH (optional)", adios2_prefix or "") or None
                adios2_dir = _ask("ADIOS2_DIR CMake config directory (optional)", adios2_dir or "") or None
        probe_args = argparse.Namespace(**vars(args))
        probe_args.cmake_prefix_path = cmake_prefix_path
        probe_args.openpmd_dir = openpmd_dir
        probe_args.adios2_prefix = adios2_prefix
        probe_args.adios2_dir = adios2_dir
        probe_args.hdf5_prefix = hdf5_prefix
        probe_args.hdf5_dir = hdf5_dir
        selected_openpmd_backend, supported_openpmd_backends, probe_report = _select_external_openpmd_backend(
            probe_args,
            openpmd_adios2,
        )

    alpaka_backends, alpaka_error = available_alpaka_backends()
    if args.interactive and alpaka_backends:
        print("Currently available alpaka compute backends:")
        print()
        print(
            "These are the hardware configurations HASEonGPU can run on "
            "right now."
        )
        print()
        for index, backend in enumerate(alpaka_backends, start=1):
            print(f"  {index}) {backend}")
        print()
        print(
            "If no NVIDIA/CUDA or AMD/HIP backend is listed, this build probably did "
            "not find a usable GPU toolchain/backend or matching device.\n\n"
            "To make NVIDIA/CUDA backends available, provide CUDA via CUDACXX or "
            "nvcc/CUDAToolkit_ROOT/CUDA_HOME/CUDA_PATH/CUDA_ROOT.\n\n"
            "To make AMD/HIP backends available, provide ROCm/HIP via HIPCXX, hipcc, "
            "or ROCM_PATH."
        )
        backend_choice = _ask("Compute backend name or number", args.compute_backend or "1")
        if str(backend_choice).isdigit() and 1 <= int(backend_choice) <= len(alpaka_backends):
            compute_backend = alpaka_backends[int(backend_choice) - 1]
        else:
            compute_backend = str(backend_choice)
    else:
        compute_backend = preferred_compute_backend(alpaka_backends, args.compute_backend)

    parallel_mode = args.parallel_mode
    if args.interactive and _ask_yes_no(
        "Enable MPI support and use MPI for multi-rank runs now",
        parallel_mode == "mpi" or mpi_mode == MPI_MODE_ON,
    ):
        parallel_mode = "mpi"
        mpi_mode = MPI_MODE_ON
    if parallel_mode not in {"single", "mpi"}:
        raise ValueError("parallel mode must be 'single' or 'mpi'")
    if parallel_mode == "mpi" and mpi_mode == MPI_MODE_OFF:
        mpi_mode = MPI_MODE_ON

    if args.interactive:
        native_optimizations = _ask_yes_no(
            "Enable native CPU optimizations? This tunes for this machine and is not suitable for redistributable wheels",
            False,
        )

    selection = WizardSelection(
        provider=provider,
        openpmd_backend=selected_openpmd_backend,
        compute_backend=compute_backend,
        supported_openpmd_backends=tuple(supported_openpmd_backends),
        parallel_mode=parallel_mode,
        num_devices=int(args.num_devices),
        n_per_node=int(args.n_per_node),
        cmake_prefix_path=cmake_prefix_path,
        openpmd_dir=openpmd_dir,
        bundled_adios2=bundled_adios2,
        bundled_hdf5=bundled_hdf5,
        openpmd_adios2=openpmd_adios2,
        bundled_python_bindings=bundled_python_bindings,
        native_optimizations=native_optimizations,
        mpi_mode=mpi_mode,
        adios2_prefix=adios2_prefix,
        adios2_dir=adios2_dir,
        hdf5_prefix=hdf5_prefix,
        hdf5_dir=hdf5_dir,
    )
    return selection, alpaka_backends, alpaka_error, probe_report


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Guide HASEonGPU openPMD provider, MPI, and Alpaka backend setup."
    )
    parser.add_argument(
        "--provider",
        choices=(PROVIDER_AUTO, PROVIDER_BUNDLED, PROVIDER_SYSTEM),
        default="auto",
        help="Provider setup to guide. Defaults to auto-detecting a system provider before bundled.",
    )
    parser.add_argument(
        "--bundled-adios2",
        choices=(BUNDLED_ADIOS2_FETCH, BUNDLED_ADIOS2_OFF, BUNDLED_ADIOS2_SYSTEM),
        default=BUNDLED_ADIOS2_FETCH,
        help="ADIOS2 handling when --provider bundled is selected.",
    )
    parser.add_argument(
        "--bundled-hdf5",
        choices=(BUNDLED_HDF5_FETCH, BUNDLED_HDF5_OFF, BUNDLED_HDF5_SYSTEM),
        default=BUNDLED_HDF5_OFF,
        help="HDF5 handling when --provider bundled is selected. Defaults to off.",
    )
    parser.add_argument(
        "--openpmd-adios2",
        choices=(OPENPMD_ADIOS2_AUTO, OPENPMD_ADIOS2_YES, OPENPMD_ADIOS2_NO),
        default=OPENPMD_ADIOS2_AUTO,
        help="Whether a system openPMD provider was built with ADIOS2 support.",
    )
    parser.add_argument(
        "--no-bundled-python-bindings",
        action="store_true",
        help="Do not build matching openPMD Python bindings for bundled-provider installs.",
    )
    parser.add_argument("--openpmd-backend", choices=OPENPMD_BACKEND_PRIORITY, default=None)
    parser.add_argument("--compute-backend", default=None)
    parser.add_argument("--parallel-mode", choices=("single", "mpi"), default="single")
    parser.add_argument(
        "--mpi",
        choices=(MPI_MODE_AUTO, MPI_MODE_ON, MPI_MODE_OFF),
        default=MPI_MODE_AUTO,
        help="MPI build support: on -> -DDISABLE_MPI=OFF, off -> ON, auto -> AUTO. Defaults to auto.",
    )
    parser.add_argument(
        "--native-optimizations",
        choices=("on", "off"),
        default="off",
        help="Enable host-specific CPU tuning. Defaults to off for redistributable wheels.",
    )
    parser.add_argument("--num-devices", type=int, default=1)
    parser.add_argument("--n-per-node", type=int, default=1)
    parser.add_argument("--cmake-prefix-path", default=None)
    parser.add_argument("--openpmd-dir", default=None)
    parser.add_argument("--adios2-prefix", default=None)
    parser.add_argument("--adios2-dir", default=None)
    parser.add_argument("--hdf5-prefix", default=None)
    parser.add_argument("--hdf5-dir", default=None)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--cmake-generator", default=preflight._default_cmake_generator())
    parser.add_argument(
        "--output",
        default=str(DEFAULT_PHI_ASE_CONFIG_PATH),
        help="PhiASE YAML path to write, or '-' for stdout only. Defaults to config/hase-phiase.yaml.",
    )
    parser.add_argument("--yes", action="store_true", help="Use defaults and do not prompt.")
    parser.add_argument("--install", action="store_true", help="Run pip install after writing configuration.")
    parser.add_argument(
        "--autoinstall",
        action="store_true",
        help="Use all defaults non-interactively, write the YAML, and run pip install.",
    )
    parser.add_argument(
        "--reinstall",
        action="store_true",
        help="Run pip install with the previous HASE CMake settings. Requires a previous configure/install.",
    )
    parser.add_argument(
        "--use-ccache",
        action="store_true",
        help="Use ccache as C/C++/CUDA CMake compiler launcher for the install.",
    )
    parser.add_argument(
        "--break-system-packages",
        action="store_true",
        help=(
            "Pass pip --break-system-packages when installing. Use only when you "
            "intentionally install into an externally managed Python environment."
        ),
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    if args.autoinstall:
        args.yes = True
        args.install = True
    if args.reinstall:
        try:
            return run_reinstall(
                break_system_packages=args.break_system_packages,
                use_ccache=args.use_ccache,
            )
        except Exception as exc:
            print(f"hase-configure: error: {exc}", file=sys.stderr)
            return 1

    args.interactive = (not args.yes) and sys.stdin.isatty()
    try:
        selection, alpaka_backends, alpaka_error, _probe_report = _build_selection(args)
    except Exception as exc:
        print(f"hase-configure: error: {exc}", file=sys.stderr)
        return 1

    config_text = yaml_config(selection)
    print("\nRecommended install command:")
    print(install_command(selection, break_system_packages=args.break_system_packages, use_ccache=args.use_ccache))

    if alpaka_backends:
        print("\nCurrently available alpaka compute backends:")
        print()
        print("These are the hardware configurations HASEonGPU can run on right now.")
        print()
        for backend in alpaka_backends:
            print(f"  - {backend}")

    print("\nGenerated PhiASE YAML:")
    print(config_text, end="")

    if args.output != "-":
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(config_text, encoding="utf-8")
        print(f"Wrote {output_path}")

    yaml_path = args.output if args.output != "-" else "<stdout>"
    print("\nGuidance:")
    print()
    for note in guidance_items(selection, yaml_path, alpaka_backends=alpaka_backends, alpaka_error=alpaka_error):
        print(f"  - {note}\n")

    install_now = args.install
    if args.interactive and not install_now:
        install_now = _ask_yes_no("Install HASEonGPU now with these settings", True)
    if install_now:
        return run_install(
            selection,
            break_system_packages=args.break_system_packages,
            use_ccache=args.use_ccache,
        )

    print("\nInstall skipped. To install later, run:")
    print(install_command(selection, break_system_packages=args.break_system_packages, use_ccache=args.use_ccache))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
