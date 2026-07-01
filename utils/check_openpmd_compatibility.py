#!/usr/bin/env python3
"""Check that Python and CMake openPMD providers match HASEonGPU needs."""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path


BACKEND_EXTENSIONS = {
    "adios": "bp",
    "adios-sst": "sst",
    "hdf5": "h5",
}
HASE_SOURCE_BUILD_BACKENDS = ("adios-sst", "adios", "hdf5")


def _default_cmake_generator() -> str | None:
    if "CMAKE_GENERATOR" in os.environ:
        return None
    if shutil.which("ninja"):
        return "Ninja"
    return None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify the openpmd_api Python module and CMake openPMD::openPMD provider."
    )
    parser.add_argument(
        "--backend",
        choices=tuple(BACKEND_EXTENSIONS),
        default="adios-sst",
        help="Runtime openPMD backend to check.",
    )
    parser.add_argument(
        "--cmake-prefix-path",
        default=None,
        help="Prefix containing openPMDConfig.cmake, passed to CMAKE_PREFIX_PATH.",
    )
    parser.add_argument(
        "--openpmd-dir",
        default=None,
        help="Directory containing openPMDConfig.cmake, passed to openPMD_DIR.",
    )
    parser.add_argument(
        "--cmake",
        default="cmake",
        help="CMake executable.",
    )
    parser.add_argument(
        "--cmake-generator",
        default=_default_cmake_generator(),
        help=(
            "CMake generator for the provider probe. Defaults to Ninja when "
            "ninja is available and CMAKE_GENERATOR is unset."
        ),
    )
    return parser.parse_args()


def _version_tuple(value: str | None) -> tuple[int, ...] | None:
    if not value:
        return None
    match = re.search(r"(\d+(?:\.\d+)+)", value)
    if not match:
        return None
    return tuple(int(part) for part in match.group(1).split("."))


def _python_info(errors: list[str]) -> dict[str, object]:
    try:
        import openpmd_api as io
    except ImportError as exc:
        errors.append(
            "Python cannot import openpmd_api. Install the matching openpmd-api "
            "Python package in this environment."
        )
        return {"error": exc}

    return {
        "version": getattr(io, "__version__", None),
        "variants": dict(getattr(io, "variants", {})),
        "file_extensions": sorted(str(item) for item in getattr(io, "file_extensions", [])),
        "path": getattr(io, "__file__", None),
    }


def _check_python_backend(info: dict[str, object], backend: str, errors: list[str]) -> None:
    if "error" in info:
        return
    variants = info.get("variants", {})
    extensions = set(info.get("file_extensions", []))
    if backend == "hdf5":
        if not variants.get("hdf5", False):
            errors.append("Python openpmd_api does not report HDF5 support.")
    else:
        if not variants.get("adios2", False):
            errors.append("Python openpmd_api does not report ADIOS2 support.")

    extension = BACKEND_EXTENSIONS[backend]
    if extension not in extensions:
        errors.append(
            f"Python openpmd_api does not report file extension '{extension}'. "
            f"Reported extensions: {sorted(extensions)}"
        )


def _cmake_probe(args: argparse.Namespace, errors: list[str]) -> dict[str, str]:
    with tempfile.TemporaryDirectory(prefix="hase-openpmd-probe-") as tmp:
        source = Path(tmp) / "src"
        build = Path(tmp) / "build"
        source.mkdir()
        (source / "CMakeLists.txt").write_text(
            textwrap.dedent(
                """
                cmake_minimum_required(VERSION 3.24)
                project(HaseOpenPmdProbe LANGUAGES CXX)
                find_package(openPMD CONFIG REQUIRED)
                if(NOT TARGET openPMD::openPMD)
                    message(FATAL_ERROR "openPMD::openPMD target is missing")
                endif()
                foreach(name IN ITEMS
                    openPMD_VERSION
                    openPMD_HAVE_ADIOS2
                    openPMD_HAVE_HDF5
                    ADIOS2_HAVE_SST
                    ADIOS2_VERSION
                    HDF5_VERSION
                )
                    if(DEFINED ${name})
                        message(STATUS "HASE_PREFLIGHT_${name}=${${name}}")
                    else()
                        message(STATUS "HASE_PREFLIGHT_${name}=<undefined>")
                    endif()
                endforeach()
                """
            ).strip()
            + "\n",
            encoding="utf-8",
        )

        command = [args.cmake, "-S", str(source), "-B", str(build)]
        if args.cmake_generator:
            command.extend(["-G", args.cmake_generator])
        if args.cmake_prefix_path:
            command.append(f"-DCMAKE_PREFIX_PATH={args.cmake_prefix_path}")
        if args.openpmd_dir:
            command.append(f"-DopenPMD_DIR={args.openpmd_dir}")

        env = os.environ.copy()
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=env,
            check=False,
        )
        if result.returncode != 0:
            errors.append(
                "CMake cannot find a usable openPMD::openPMD provider. "
                "Pass --cmake-prefix-path or --openpmd-dir for the external C++ installation."
            )
            return {"log": result.stdout, "command": shlex.join(command)}

    values: dict[str, str] = {"log": result.stdout, "command": shlex.join(command)}
    for line in result.stdout.splitlines():
        marker = "HASE_PREFLIGHT_"
        if marker in line:
            name_value = line.split(marker, 1)[1]
            if "=" in name_value:
                name, value = name_value.split("=", 1)
                values[name.strip()] = value.strip()
    return values


def _cmake_bool(value: str | None) -> bool | None:
    if value is None or value == "<undefined>":
        return None
    return value.upper() in {"1", "ON", "TRUE", "YES"}


def _check_cmake_backend(info: dict[str, str], backend: str, errors: list[str], warnings: list[str]) -> None:
    if "command" not in info or "openPMD_VERSION" not in info:
        return
    if backend == "hdf5":
        hdf5 = _cmake_bool(info.get("openPMD_HAVE_HDF5"))
        if hdf5 is False:
            errors.append("CMake openPMD provider reports openPMD_HAVE_HDF5=OFF.")
        elif hdf5 is None:
            warnings.append("CMake openPMD provider did not expose openPMD_HAVE_HDF5.")
    else:
        adios2 = _cmake_bool(info.get("openPMD_HAVE_ADIOS2"))
        if adios2 is False:
            errors.append("CMake openPMD provider reports openPMD_HAVE_ADIOS2=OFF.")
        elif adios2 is None:
            warnings.append("CMake openPMD provider did not expose openPMD_HAVE_ADIOS2.")
        if backend == "adios-sst":
            sst = _cmake_bool(info.get("ADIOS2_HAVE_SST"))
            if sst is False:
                errors.append("CMake ADIOS2 provider reports ADIOS2_HAVE_SST=OFF.")
            elif sst is None:
                warnings.append(
                    "CMake probe could not confirm ADIOS2 SST support from CMake variables."
                )


def _check_versions(python_info: dict[str, object], cmake_info: dict[str, str], warnings: list[str]) -> None:
    python_version = _version_tuple(str(python_info.get("version") or ""))
    cmake_version = _version_tuple(cmake_info.get("openPMD_VERSION"))
    if python_version and cmake_version and python_version[:2] != cmake_version[:2]:
        warnings.append(
            "Python openpmd_api and CMake openPMD provider have different major/minor versions: "
            f"Python={python_info.get('version')}, CMake={cmake_info.get('openPMD_VERSION')}."
        )


def _python_supported_backends(info: dict[str, object]) -> list[str]:
    if "error" in info:
        return []
    variants = info.get("variants", {})
    extensions = set(info.get("file_extensions", []))
    backends: list[str] = []
    if variants.get("adios2", False) and "bp" in extensions:
        backends.append("adios")
    if variants.get("adios2", False) and "sst" in extensions:
        backends.append("adios-sst")
    if variants.get("hdf5", False) and "h5" in extensions:
        backends.append("hdf5")
    return backends


def _cmake_supported_backends(info: dict[str, str]) -> list[str]:
    if "openPMD_VERSION" not in info:
        return []
    backends: list[str] = []
    adios2 = _cmake_bool(info.get("openPMD_HAVE_ADIOS2"))
    hdf5 = _cmake_bool(info.get("openPMD_HAVE_HDF5"))
    sst = _cmake_bool(info.get("ADIOS2_HAVE_SST"))
    if adios2 is True:
        backends.append("adios")
        if sst is True:
            backends.append("adios-sst")
        elif sst is None:
            backends.append("adios-sst (SST unconfirmed)")
    if hdf5 is True:
        backends.append("hdf5")
    return backends


def _format_backend_list(backends: list[str]) -> str:
    return ", ".join(backends) if backends else "<none confirmed>"


def _print_backend_summary(
    python_info: dict[str, object],
    cmake_info: dict[str, str],
    selected_backend: str,
) -> None:
    source_build_backends = ", ".join(HASE_SOURCE_BUILD_BACKENDS[1:])
    print("Backend support summary:")
    print(
        "  HASE source-build provider: "
        f"{HASE_SOURCE_BUILD_BACKENDS[0]} (default), {source_build_backends}"
    )
    if "error" in python_info:
        print("  Python provider: unavailable")
    else:
        print(
            "  Python provider: "
            f"{_format_backend_list(_python_supported_backends(python_info))}"
        )
    print(
        "  CMake provider: "
        f"{_format_backend_list(_cmake_supported_backends(cmake_info))}"
    )
    print(f"  Selected backend: {selected_backend}")


def main() -> int:
    args = _parse_args()
    errors: list[str] = []
    warnings: list[str] = []

    python_info = _python_info(errors)
    _check_python_backend(python_info, args.backend, errors)
    cmake_info = _cmake_probe(args, errors)
    _check_cmake_backend(cmake_info, args.backend, errors, warnings)
    _check_versions(python_info, cmake_info, warnings)

    print("Python openPMD provider:")
    if "error" in python_info:
        print("  unavailable")
    else:
        print(f"  path: {python_info.get('path')}")
        print(f"  version: {python_info.get('version') or '<unknown>'}")
        print(f"  variants: {python_info.get('variants')}")
        print(f"  file extensions: {python_info.get('file_extensions')}")

    print("CMake openPMD provider:")
    print(f"  command: {cmake_info.get('command', '<not run>')}")
    if "openPMD_VERSION" in cmake_info:
        print(f"  version: {cmake_info.get('openPMD_VERSION')}")
        print(f"  ADIOS2: {cmake_info.get('openPMD_HAVE_ADIOS2', '<unknown>')}")
        print(f"  HDF5: {cmake_info.get('openPMD_HAVE_HDF5', '<unknown>')}")
        print(f"  SST: {cmake_info.get('ADIOS2_HAVE_SST', '<unknown>')}")

    _print_backend_summary(python_info, cmake_info, args.backend)

    for warning in warnings:
        print(f"warning: {warning}", file=sys.stderr)
    if errors:
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
        if "log" in cmake_info:
            print("CMake probe output:", file=sys.stderr)
            print(cmake_info["log"], file=sys.stderr)
        return 1

    print(f"openPMD compatibility check passed for backend '{args.backend}'.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
