import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


repoRoot = Path(__file__).resolve().parents[3]


def _cmake_generator_arguments():
    configured_generator = os.environ.get("CMAKE_GENERATOR")
    if configured_generator:
        return ["-G", configured_generator]
    if shutil.which("ninja"):
        return ["-G", "Ninja"]
    if shutil.which("make"):
        return ["-G", "Unix Makefiles"]
    pytest.skip("The CMake contract tests require Ninja or Make")


def _write_fake_system_dependencies(prefix):
    """Create lightweight system packages for a real top-level HASE configure."""
    alpaka_config = prefix / "lib" / "cmake" / "alpaka" / "alpakaConfig.cmake"
    alpaka_config.parent.mkdir(parents=True)
    alpaka_config.write_text(
        """
if(NOT TARGET alpaka::alpaka)
    add_library(alpaka::alpaka INTERFACE IMPORTED)
endif()
function(alpaka_finalize targetName)
endfunction()
""".lstrip(),
        encoding="utf-8",
    )

    adios2_config = prefix / "lib" / "cmake" / "ADIOS2" / "ADIOS2Config.cmake"
    adios2_config.parent.mkdir(parents=True)
    adios2_config.write_text(
        """
set(ADIOS2_FOUND TRUE)
set(ADIOS2_VERSION "2.12.1")
""".lstrip(),
        encoding="utf-8",
    )


def _write_fake_system_openpmd(prefix):
    openpmd_config = (
        prefix / "lib" / "cmake" / "openPMD" / "openPMDConfig.cmake"
    )
    openpmd_config.parent.mkdir(parents=True)
    openpmd_config.write_text(
        """
if(NOT TARGET openPMD::openPMD)
    add_library(openPMD::openPMD INTERFACE IMPORTED)
endif()
set(openPMD_FOUND TRUE)
set(openPMD_VERSION "0.17.1")
set(openPMD_HAVE_ADIOS2 TRUE)
set(openPMD_HAVE_HDF5 FALSE)
""".lstrip(),
        encoding="utf-8",
    )
    return openpmd_config.parent


def _create_fake_openpmd_repository(
    source_dir,
    *,
    record_adios2_dir=False,
    record_python_executable=False,
    record_toolchain_state=False,
    require_toolchain_marker=False,
):
    """Create a buildable local openPMD repository that requires system ADIOS2."""
    source_dir.mkdir()
    toolchain_checks = ""
    if require_toolchain_marker:
        toolchain_checks = """
if(NOT HASE_CONTRACT_TOOLCHAIN_MARKER)
    message(FATAL_ERROR "The HASE toolchain file did not reach openPMD")
endif()
if(NOT CMAKE_CXX_FLAGS MATCHES "HASE_CONTRACT_ABI_MARKER")
    message(FATAL_ERROR "The HASE ABI flags did not reach openPMD")
endif()
"""
    provider_recording = ""
    if record_adios2_dir:
        provider_recording += """
file(
    WRITE
    "${CMAKE_CURRENT_BINARY_DIR}/provider-adios2-dir.txt"
    "${ADIOS2_DIR}\\n"
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/provider-adios2-dir.txt"
    DESTINATION .
)
"""
    if record_python_executable:
        provider_recording += """
file(
    WRITE
    "${CMAKE_CURRENT_BINARY_DIR}/provider-python-executable.txt"
    "${Python3_EXECUTABLE}\\n"
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/provider-python-executable.txt"
    DESTINATION .
)
"""
    if record_toolchain_state:
        provider_recording += """
file(
    WRITE
    "${CMAKE_CURRENT_BINARY_DIR}/provider-toolchain-state.txt"
    "${HASE_CONTRACT_TOOLCHAIN_MARKER}\\n${CMAKE_CXX_FLAGS}\\n"
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/provider-toolchain-state.txt"
    DESTINATION .
)
"""
    (source_dir / "openPMDConfig.cmake.in").write_text(
        """
include("${CMAKE_CURRENT_LIST_DIR}/openPMDTargets.cmake")
set(openPMD_FOUND TRUE)
set(openPMD_VERSION "0.17.1")
set(openPMD_HAVE_ADIOS2 TRUE)
set(openPMD_HAVE_HDF5 FALSE)
""".lstrip(),
        encoding="utf-8",
    )
    (source_dir / "__init__.py").write_text("", encoding="utf-8")
    (source_dir / "CMakeLists.txt").write_text(
        (
            """
cmake_minimum_required(VERSION 3.24)
project(openPMD VERSION 0.17.1 LANGUAGES CXX)

include(GNUInstallDirs)
find_package(ADIOS2 CONFIG REQUIRED)
"""
            + toolchain_checks
            + """
add_library(openPMD INTERFACE)
add_library(openPMD::openPMD ALIAS openPMD)
"""
            + provider_recording
            + """
install(TARGETS openPMD EXPORT openPMDTargets)
install(
    EXPORT openPMDTargets
    FILE openPMDTargets.cmake
    NAMESPACE openPMD::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/openPMD"
)
configure_file(
    openPMDConfig.cmake.in
    openPMDConfig.cmake
    @ONLY
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/openPMDConfig.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/openPMD"
)
install(
    FILES "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
    DESTINATION "lib/python3/site-packages/openpmd_api"
)
"""
        ).lstrip(),
        encoding="utf-8",
    )
    subprocess.run(["git", "init", "-q"], cwd=source_dir, check=True)
    subprocess.run(
        [
            "git",
            "add",
            "CMakeLists.txt",
            "openPMDConfig.cmake.in",
            "__init__.py",
        ],
        cwd=source_dir,
        check=True,
    )
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=HASE contract fixture",
            "-c",
            "user.email=fixture@example.invalid",
            "commit",
            "-qm",
            "Add fake openPMD provider",
        ],
        cwd=source_dir,
        check=True,
    )
    subprocess.run(
        ["git", "tag", "0.17.1"],
        cwd=source_dir,
        check=True,
    )


def _git_rewrite_environment(openpmd_source):
    return {
        **os.environ,
        "GIT_ALLOW_PROTOCOL": "file",
        "GIT_CONFIG_COUNT": "1",
        "GIT_CONFIG_KEY_0": (
            f"url.file://{openpmd_source.as_posix()}.insteadOf"
        ),
        "GIT_CONFIG_VALUE_0": "https://github.com/openPMD/openPMD-api.git",
    }


def _cmake_cache_values(cache_path):
    values = {}
    for line in cache_path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        key_type, value = line.split("=", 1)
        values[key_type.split(":", 1)[0]] = value
    return values


def _write_python_wrapper(path):
    path.write_text(
        f"#!/bin/sh\nexec {shlex.quote(sys.executable)} \"$@\"\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def _write_empty_compile_launcher(path):
    """Write a fast compiler launcher for contract-only native artifacts."""
    path.write_text(
        """
import subprocess
import sys

compiler = sys.argv[1]
arguments = sys.argv[2:]
output = arguments[arguments.index("-o") + 1]
source = next(
    (
        argument
        for argument in arguments
        if argument.endswith((".c", ".cc", ".cpp", ".cxx"))
    ),
    "",
)
is_executable = "calcPhiASE.dir" in output
language = "c++" if source.endswith((".cc", ".cpp", ".cxx")) else "c"
source_text = "int main(int, char**) { return 0; }\\n" if is_executable else ""
command = [compiler, "-x", language, "-c", "-o", output, "-"]
if "-fPIC" in arguments:
    command.insert(1, "-fPIC")
subprocess.run(command, input=source_text, text=True, check=True)
""".lstrip(),
        encoding="utf-8",
    )


def test_documentedConfiguratorPreservesRuntimePathWithSpaces(tmp_path):
    runtime_dir = (
        tmp_path / "runtime with spaces $dollar;semicolon'quote\"double"
    )
    configured = subprocess.run(
        [
            sys.executable,
            str(repoRoot / "utils" / "configure_hase.py"),
            "--provider",
            "bundled",
            "--openpmd-backend",
            "adios",
            "--runtime-dir",
            str(runtime_dir),
            "--mpi",
            "off",
            "--yes",
            "--output",
            "-",
        ],
        cwd=tmp_path,
        env={
            **os.environ,
            "PYTHONPATH": str(repoRoot),
        },
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assignment = next(
        line
        for line in configured.stdout.splitlines()
        if line.startswith("HASE_CONFIGURE_CMAKE_ARGS=")
    )
    decoded_assignment = shlex.split(assignment)
    assert len(decoded_assignment) == 1, assignment
    variable, serialized_arguments = decoded_assignment[0].split("=", 1)
    assert variable == "HASE_CONFIGURE_CMAKE_ARGS"
    decoded_arguments = shlex.split(serialized_arguments)

    assert f"-DHASE_RUNTIME_DIR={runtime_dir}" in decoded_arguments


@pytest.mark.skipif(
    os.name == "nt",
    reason="The lightweight compiler-launcher fixture uses GCC-style arguments",
)
def test_thinFrontendPropagatesBuildSettingsToResidentRuntime(tmp_path):
    dependency_prefix = tmp_path / "system dependencies"
    runtime_dir = tmp_path / "resident runtime"
    frontend_dir = tmp_path / "thin frontend"
    compiler_launcher = tmp_path / "emptyCompileLauncher.py"
    _write_fake_system_dependencies(dependency_prefix)
    system_openpmd_dir = _write_fake_system_openpmd(dependency_prefix)
    _write_empty_compile_launcher(compiler_launcher)

    configured_arguments = [
        "-DHASE_TESTING=OFF",
        "-DDISABLE_MPI=ON",
        "-DHASE_SELECT_BACKEND_ALPAKA=ON",
        "-DHASE_USE_SYSTEM_ALPAKA=ON",
        (
            "-Dalpaka_DIR="
            f"{dependency_prefix / 'lib' / 'cmake' / 'alpaka'}"
        ),
        "-DHASE_OPENPMD_PROVIDER=system",
        f"-DopenPMD_DIR={system_openpmd_dir}",
        "-DHASE_BUILD_RELEASE=OFF",
        "-DCMAKE_BUILD_TYPE=Release",
        (
            "-DCMAKE_CXX_COMPILER_LAUNCHER="
            f"{sys.executable};{compiler_launcher}"
        ),
    ]
    first_class_build_setting = (
        "-DCMAKE_CXX_FLAGS=-DHASE_RUNTIME_ABI_MARKER"
    )
    configured = subprocess.run(
        [
            "cmake",
            *_cmake_generator_arguments(),
            "-S",
            str(repoRoot),
            "-B",
            str(frontend_dir),
            "-DHASE_BUILD_RUNTIME=OFF",
            f"-DHASE_RUNTIME_DIR={runtime_dir}",
            *configured_arguments,
            first_class_build_setting,
        ],
        env={
            **os.environ,
            # This is the public pip/scikit-build option channel. The ABI flag
            # deliberately remains a first-class CMake setting outside it.
            "CMAKE_ARGS": shlex.join(configured_arguments),
        },
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    assert configured.returncode == 0, configured.stdout
    runtime_cache = _cmake_cache_values(runtime_dir / "CMakeCache.txt")
    assert runtime_cache["CMAKE_CXX_FLAGS"] == "-DHASE_RUNTIME_ABI_MARKER"
    assert Path(runtime_cache["openPMD_DIR"]) == system_openpmd_dir
    assert (
        runtime_dir / "python" / "pyInclude" / "_runtime" / "calcPhiASE"
    ).is_file()


def test_documentedBundledProviderFindsSystemAdios2Prefix(tmp_path):
    dependency_prefix = tmp_path / "system dependencies"
    openpmd_source = tmp_path / "fake-openpmd"
    build_dir = tmp_path / "build"
    _write_fake_system_dependencies(dependency_prefix)
    _create_fake_openpmd_repository(openpmd_source)

    configured = subprocess.run(
        [
            "cmake",
            *_cmake_generator_arguments(),
            "-S",
            str(repoRoot),
            "-B",
            str(build_dir),
            "-DHASE_ENABLE_PYTHON=OFF",
            "-DHASE_TESTING=OFF",
            "-DDISABLE_MPI=ON",
            "-DHASE_SELECT_BACKEND_ALPAKA=ON",
            "-DHASE_USE_SYSTEM_ALPAKA=ON",
            "-DHASE_OPENPMD_PROVIDER=bundled",
            "-DHASE_OPENPMD_USE_ADIOS2=ON",
            "-DHASE_OPENPMD_FETCH_ADIOS2=OFF",
            "-DHASE_OPENPMD_USE_HDF5=OFF",
            "-DHASE_OPENPMD_USE_SST=OFF",
            f"-DCMAKE_PREFIX_PATH={dependency_prefix}",
        ],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    assert configured.returncode == 0, configured.stdout


def test_bundledProviderRebuildsForChangedToolchainAndAbi(tmp_path):
    dependency_prefix = tmp_path / "system dependencies"
    openpmd_source = tmp_path / "fake-openpmd"
    build_dir = tmp_path / "build"
    toolchain_file = tmp_path / "contract-toolchain.cmake"
    _write_fake_system_dependencies(dependency_prefix)
    _create_fake_openpmd_repository(
        openpmd_source,
        record_toolchain_state=True,
        require_toolchain_marker=True,
    )
    toolchain_file.write_text(
        'set(HASE_CONTRACT_TOOLCHAIN_MARKER "first toolchain")\n',
        encoding="utf-8",
    )

    common_arguments = [
        "cmake",
        *_cmake_generator_arguments(),
        "-S",
        str(repoRoot),
        "-B",
        str(build_dir),
        "-DHASE_ENABLE_PYTHON=OFF",
        "-DHASE_TESTING=OFF",
        "-DDISABLE_MPI=ON",
        "-DHASE_SELECT_BACKEND_ALPAKA=ON",
        "-DHASE_USE_SYSTEM_ALPAKA=ON",
        (
            "-Dalpaka_DIR="
            f"{dependency_prefix / 'lib' / 'cmake' / 'alpaka'}"
        ),
        "-DHASE_OPENPMD_PROVIDER=bundled",
        "-DHASE_OPENPMD_USE_ADIOS2=ON",
        "-DHASE_OPENPMD_FETCH_ADIOS2=OFF",
        "-DHASE_OPENPMD_USE_HDF5=OFF",
        "-DHASE_OPENPMD_USE_SST=OFF",
        (
            "-DADIOS2_DIR="
            f"{dependency_prefix / 'lib' / 'cmake' / 'ADIOS2'}"
        ),
        f"-DCMAKE_TOOLCHAIN_FILE={toolchain_file}",
    ]
    first_abi = "-DCMAKE_CXX_FLAGS=-DHASE_CONTRACT_ABI_MARKER_FIRST"
    configured = subprocess.run(
        [*common_arguments, first_abi],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    assert configured.returncode == 0, configured.stdout
    provider_record = (
        build_dir
        / "hase-openpmd-provider"
        / "install"
        / "provider-toolchain-state.txt"
    )
    assert provider_record.read_text(encoding="utf-8").splitlines() == [
        "first toolchain",
        "-DHASE_CONTRACT_ABI_MARKER_FIRST",
    ]

    toolchain_file.write_text(
        'set(HASE_CONTRACT_TOOLCHAIN_MARKER "second toolchain")\n',
        encoding="utf-8",
    )
    changed_toolchain = subprocess.run(
        [*common_arguments, first_abi],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert changed_toolchain.returncode == 0, changed_toolchain.stdout
    assert provider_record.read_text(encoding="utf-8").splitlines() == [
        "second toolchain",
        "-DHASE_CONTRACT_ABI_MARKER_FIRST",
    ]

    changed_abi = subprocess.run(
        [
            *common_arguments,
            "-DCMAKE_CXX_FLAGS=-DHASE_CONTRACT_ABI_MARKER_SECOND",
        ],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert changed_abi.returncode == 0, changed_abi.stdout
    assert provider_record.read_text(encoding="utf-8").splitlines() == [
        "second toolchain",
        "-DHASE_CONTRACT_ABI_MARKER_SECOND",
    ]


def test_documentedProviderSwitchReplacesManagedBundledPaths(tmp_path):
    dependency_prefix = tmp_path / "system dependencies"
    openpmd_source = tmp_path / "fake-openpmd"
    build_dir = tmp_path / "build"
    _write_fake_system_dependencies(dependency_prefix)
    system_openpmd_dir = _write_fake_system_openpmd(dependency_prefix)
    _create_fake_openpmd_repository(openpmd_source)

    bundled = subprocess.run(
        [
            "cmake",
            *_cmake_generator_arguments(),
            "-S",
            str(repoRoot),
            "-B",
            str(build_dir),
            "-DHASE_ENABLE_PYTHON=ON",
            "-DHASE_TESTING=OFF",
            "-DDISABLE_MPI=ON",
            "-DHASE_SELECT_BACKEND_ALPAKA=ON",
            "-DHASE_USE_SYSTEM_ALPAKA=ON",
            (
                "-Dalpaka_DIR="
                f"{dependency_prefix / 'lib' / 'cmake' / 'alpaka'}"
            ),
            "-DHASE_OPENPMD_PROVIDER=bundled",
            "-DHASE_OPENPMD_USE_ADIOS2=ON",
            "-DHASE_OPENPMD_FETCH_ADIOS2=OFF",
            "-DHASE_OPENPMD_USE_HDF5=OFF",
            "-DHASE_OPENPMD_USE_SST=OFF",
            (
                "-DADIOS2_DIR="
                f"{dependency_prefix / 'lib' / 'cmake' / 'ADIOS2'}"
            ),
        ],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert bundled.returncode == 0, bundled.stdout

    bundled_cache = _cmake_cache_values(build_dir / "CMakeCache.txt")
    bundled_prefix = build_dir / "hase-openpmd-provider" / "install"
    replacement_bundled_prefix = build_dir / "replacement-provider"
    assert Path(bundled_cache["openPMD_DIR"]).is_relative_to(bundled_prefix)
    assert Path(
        bundled_cache["HASE_OPENPMD_PYTHON_PACKAGE_DIR"]
    ).is_relative_to(bundled_prefix)

    system = subprocess.run(
        [
            "cmake",
            *_cmake_generator_arguments(),
            "-S",
            str(repoRoot),
            "-B",
            str(build_dir),
            "-DHASE_OPENPMD_PROVIDER=system",
            f"-DHASE_OPENPMD_BUNDLED_PREFIX={replacement_bundled_prefix}",
            f"-DCMAKE_PREFIX_PATH={dependency_prefix}",
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert system.returncode == 0, system.stdout

    system_cache = _cmake_cache_values(build_dir / "CMakeCache.txt")
    assert Path(system_cache["openPMD_DIR"]) == system_openpmd_dir
    assert system_cache.get("HASE_OPENPMD_PYTHON_PACKAGE_DIR", "") == ""
    assert system_cache["HASE_OPENPMD_BUILD_PYTHON_BINDINGS"] == "ON"

    bundled_again = subprocess.run(
        [
            "cmake",
            *_cmake_generator_arguments(),
            "-S",
            str(repoRoot),
            "-B",
            str(build_dir),
            "-DHASE_OPENPMD_PROVIDER=bundled",
        ],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert bundled_again.returncode == 0, bundled_again.stdout

    bundled_again_cache = _cmake_cache_values(build_dir / "CMakeCache.txt")
    assert Path(bundled_again_cache["openPMD_DIR"]).is_relative_to(
        replacement_bundled_prefix
    )
    assert Path(
        bundled_again_cache["HASE_OPENPMD_PYTHON_PACKAGE_DIR"]
    ).is_relative_to(replacement_bundled_prefix)


def test_bundledProviderRebuildsForChangedSystemDependency(tmp_path):
    first_prefix = tmp_path / "first system dependency"
    second_prefix = tmp_path / "second system dependency"
    openpmd_source = tmp_path / "fake-openpmd"
    build_dir = tmp_path / "build"
    _write_fake_system_dependencies(first_prefix)
    _write_fake_system_dependencies(second_prefix)
    _create_fake_openpmd_repository(
        openpmd_source,
        record_adios2_dir=True,
    )
    first_adios2_dir = first_prefix / "lib" / "cmake" / "ADIOS2"
    second_adios2_dir = second_prefix / "lib" / "cmake" / "ADIOS2"

    first = subprocess.run(
        [
            "cmake",
            *_cmake_generator_arguments(),
            "-S",
            str(repoRoot),
            "-B",
            str(build_dir),
            "-DHASE_ENABLE_PYTHON=OFF",
            "-DHASE_TESTING=OFF",
            "-DDISABLE_MPI=ON",
            "-DHASE_SELECT_BACKEND_ALPAKA=ON",
            "-DHASE_USE_SYSTEM_ALPAKA=ON",
            (
                "-Dalpaka_DIR="
                f"{first_prefix / 'lib' / 'cmake' / 'alpaka'}"
            ),
            "-DHASE_OPENPMD_PROVIDER=bundled",
            "-DHASE_OPENPMD_USE_ADIOS2=ON",
            "-DHASE_OPENPMD_FETCH_ADIOS2=OFF",
            "-DHASE_OPENPMD_USE_HDF5=OFF",
            "-DHASE_OPENPMD_USE_SST=OFF",
            f"-DADIOS2_DIR={first_adios2_dir}",
        ],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert first.returncode == 0, first.stdout

    provider_record = (
        build_dir
        / "hase-openpmd-provider"
        / "install"
        / "provider-adios2-dir.txt"
    )
    assert provider_record.read_text(
        encoding="utf-8"
    ).strip() == str(first_adios2_dir)

    second = subprocess.run(
        [
            "cmake",
            "-S",
            str(repoRoot),
            "-B",
            str(build_dir),
            f"-DADIOS2_DIR={second_adios2_dir}",
        ],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert second.returncode == 0, second.stdout
    assert provider_record.read_text(
        encoding="utf-8"
    ).strip() == str(second_adios2_dir)


@pytest.mark.skipif(
    os.name == "nt",
    reason="The lightweight alternate-interpreter fixture uses a POSIX wrapper",
)
def test_bundledProviderRebuildsForChangedPythonInterpreter(tmp_path):
    dependency_prefix = tmp_path / "system dependency"
    openpmd_source = tmp_path / "fake-openpmd"
    build_dir = tmp_path / "build"
    first_python = tmp_path / "python-first"
    second_python = tmp_path / "python-second"
    _write_fake_system_dependencies(dependency_prefix)
    _create_fake_openpmd_repository(
        openpmd_source,
        record_python_executable=True,
    )
    _write_python_wrapper(first_python)
    _write_python_wrapper(second_python)

    common_arguments = [
        "cmake",
        *_cmake_generator_arguments(),
        "-S",
        str(repoRoot),
        "-B",
        str(build_dir),
        "-DHASE_ENABLE_PYTHON=ON",
        "-DHASE_TESTING=OFF",
        "-DDISABLE_MPI=ON",
        "-DHASE_SELECT_BACKEND_ALPAKA=ON",
        "-DHASE_USE_SYSTEM_ALPAKA=ON",
        (
            "-Dalpaka_DIR="
            f"{dependency_prefix / 'lib' / 'cmake' / 'alpaka'}"
        ),
        "-DHASE_OPENPMD_PROVIDER=bundled",
        "-DHASE_OPENPMD_USE_ADIOS2=ON",
        "-DHASE_OPENPMD_FETCH_ADIOS2=OFF",
        "-DHASE_OPENPMD_USE_HDF5=OFF",
        "-DHASE_OPENPMD_USE_SST=OFF",
        (
            "-DADIOS2_DIR="
            f"{dependency_prefix / 'lib' / 'cmake' / 'ADIOS2'}"
        ),
    ]
    first = subprocess.run(
        [*common_arguments, f"-DPython3_EXECUTABLE={first_python}"],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert first.returncode == 0, first.stdout

    provider_record = (
        build_dir
        / "hase-openpmd-provider"
        / "install"
        / "provider-python-executable.txt"
    )
    assert provider_record.read_text(
        encoding="utf-8"
    ).strip() == str(first_python)

    second = subprocess.run(
        [*common_arguments, f"-DPython3_EXECUTABLE={second_python}"],
        env=_git_rewrite_environment(openpmd_source),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    assert second.returncode == 0, second.stdout
    assert provider_record.read_text(
        encoding="utf-8"
    ).strip() == str(second_python)


def test_thinFrontendBuildsSelectedMultiConfigRuntime(tmp_path):
    capabilities = json.loads(
        subprocess.run(
            ["cmake", "-E", "capabilities"],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        ).stdout
    )
    generators = {
        generator["name"] for generator in capabilities["generators"]
    }
    if "Ninja Multi-Config" not in generators:
        pytest.skip("CMake does not provide Ninja Multi-Config")

    runtime_source = tmp_path / "runtime-source"
    runtime_build = tmp_path / "runtime-build"
    frontend_build = tmp_path / "frontend-build"
    runtime_source.mkdir()
    (runtime_source / "createRuntime.cmake").write_text(
        """
if(CONFIG STREQUAL "Release")
    file(MAKE_DIRECTORY "${OUTPUT_DIR}")
    file(TOUCH "${OUTPUT_DIR}/calcPhiASE${EXECUTABLE_SUFFIX}")
endif()
file(WRITE "${BINARY_DIR}/built-config.txt" "${CONFIG}\\n")
""".lstrip(),
        encoding="utf-8",
    )
    (runtime_source / "CMakeLists.txt").write_text(
        """
cmake_minimum_required(VERSION 3.24)
project(HaseContractRuntime LANGUAGES C CXX)
add_custom_target(
    runtime ALL
    COMMAND
        "${CMAKE_COMMAND}"
        "-DCONFIG=$<CONFIG>"
        "-DOUTPUT_DIR=${CMAKE_BINARY_DIR}/python/pyInclude/_runtime"
        "-DBINARY_DIR=${CMAKE_BINARY_DIR}"
        "-DEXECUTABLE_SUFFIX=${CMAKE_EXECUTABLE_SUFFIX}"
        -P
        "${CMAKE_CURRENT_SOURCE_DIR}/createRuntime.cmake"
)
""".lstrip(),
        encoding="utf-8",
    )
    subprocess.run(
        [
            "cmake",
            "-S",
            str(runtime_source),
            "-B",
            str(runtime_build),
            "-G",
            "Ninja Multi-Config",
        ],
        check=True,
    )

    configured = subprocess.run(
        [
            "cmake",
            "-S",
            str(repoRoot),
            "-B",
            str(frontend_build),
            "-G",
            "Ninja",
            "-DHASE_BUILD_RUNTIME=OFF",
            f"-DHASE_RUNTIME_DIR={runtime_build}",
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    assert configured.returncode == 0, configured.stdout
    assert (runtime_build / "built-config.txt").read_text(
        encoding="utf-8"
    ).strip() == "Release"
