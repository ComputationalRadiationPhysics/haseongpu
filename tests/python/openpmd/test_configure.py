from pathlib import Path

import yaml

from pyInclude import PhiASE
from pyInclude import config as configure_module
from pyInclude.config import (
    BUNDLED_ADIOS2_FETCH,
    BUNDLED_ADIOS2_OFF,
    BUNDLED_ADIOS2_SYSTEM,
    BUNDLED_HDF5_FETCH,
    BUNDLED_HDF5_OFF,
    BUNDLED_HDF5_SYSTEM,
    MPI_MODE_OFF,
    MPI_MODE_ON,
    OPENPMD_ADIOS2_NO,
    PROVIDER_BUNDLED,
    PROVIDER_SYSTEM,
    WizardSelection,
    alpaka_backend_guidance,
    bundled_supported_openpmd_backends,
    install_command,
    preferred_compute_backend,
    preferred_openpmd_backend,
    provider_notes,
    supported_external_openpmd_backends,
    yaml_config,
)


def test_externalBackendsIntersectSupport():
    python_info = {
        "variants": {"adios2": True, "hdf5": True},
        "file_extensions": ["bp", "sst", "h5"],
    }
    cmake_info = {
        "openPMD_VERSION": "0.17.0",
        "openPMD_HAVE_ADIOS2": "TRUE",
        "openPMD_HAVE_HDF5": "FALSE",
        "ADIOS2_HAVE_SST": "TRUE",
    }

    assert supported_external_openpmd_backends(python_info, cmake_info) == ["adios-sst", "adios"]


def test_preferredBackendValidatesRequest():
    assert preferred_openpmd_backend(["adios", "hdf5"]) == "adios"
    assert preferred_openpmd_backend(["adios", "hdf5"], "hdf5") == "hdf5"


def test_bundledDefaultsAreAdios2Only():
    assert bundled_supported_openpmd_backends(BUNDLED_ADIOS2_FETCH) == [
        "adios-sst",
        "adios",
    ]
    assert bundled_supported_openpmd_backends(BUNDLED_ADIOS2_FETCH, BUNDLED_HDF5_FETCH) == [
        "adios-sst",
        "adios",
        "hdf5",
    ]
    assert bundled_supported_openpmd_backends(BUNDLED_ADIOS2_FETCH, BUNDLED_HDF5_SYSTEM) == [
        "adios-sst",
        "adios",
        "hdf5",
    ]
    assert bundled_supported_openpmd_backends(BUNDLED_ADIOS2_OFF, BUNDLED_HDF5_SYSTEM) == ["hdf5"]


def test_preferredComputeBackend():
    assert preferred_compute_backend(["Cuda_Gpu_CudaRt", "Host_Cpu_CpuSerial"]) == "Host_Cpu_CpuSerial"
    assert preferred_compute_backend(["Cuda_Gpu_CudaRt"]) == "Cuda_Gpu_CudaRt"


def test_bundledFetchInstallCommandMentionsSourceBuild():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        bundled_adios2=BUNDLED_ADIOS2_FETCH,
        bundled_hdf5=BUNDLED_HDF5_FETCH,
    )

    command = install_command(selection)

    assert "-DHASE_OPENPMD_PROVIDER=bundled" in command
    assert "-DHASE_OPENPMD_BUILD_PYTHON_BINDINGS=ON" in command
    assert "-DHASE_OPENPMD_USE_ADIOS2=ON" in command
    assert "-DHASE_OPENPMD_FETCH_ADIOS2=ON" in command
    assert "-DHASE_OPENPMD_USE_HDF5=ON" in command
    assert "-DHASE_OPENPMD_FETCH_HDF5=ON" in command


def test_bundledInstallSkipsPythonBindings():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        bundled_adios2=BUNDLED_ADIOS2_FETCH,
        bundled_python_bindings=False,
    )

    assert "-DHASE_OPENPMD_BUILD_PYTHON_BINDINGS=ON" not in install_command(selection)


def test_bundledHdf5UsesSystemHdf5():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="hdf5",
        compute_backend="Host_Cpu_CpuSerial",
        bundled_adios2=BUNDLED_ADIOS2_OFF,
        bundled_hdf5=BUNDLED_HDF5_SYSTEM,
    )

    command = install_command(selection)

    assert "-DHASE_OPENPMD_USE_ADIOS2=OFF" in command
    assert "-DHASE_OPENPMD_USE_SST=OFF" in command
    assert "-DHASE_OPENPMD_USE_HDF5=ON" in command
    assert "-DHASE_OPENPMD_FETCH_HDF5=OFF" in command


def test_bundledAdiosOnlyInstallCommandDisablesHdf5():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios",
        compute_backend="Host_Cpu_CpuSerial",
        bundled_hdf5=BUNDLED_HDF5_OFF,
    )

    command = install_command(selection)

    assert "-DHASE_OPENPMD_USE_HDF5=OFF" in command
    assert bundled_supported_openpmd_backends(BUNDLED_ADIOS2_FETCH, BUNDLED_HDF5_OFF) == ["adios-sst", "adios"]


def test_bundledSystemCarriesAdios2Paths():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        bundled_adios2=BUNDLED_ADIOS2_SYSTEM,
        adios2_prefix="/opt/adios2",
        adios2_dir="/opt/adios2/lib/cmake/adios2",
    )

    command = install_command(selection)

    assert "-DHASE_OPENPMD_FETCH_ADIOS2=OFF" in command
    assert "-DCMAKE_PREFIX_PATH=/opt/adios2" in command
    assert "-DADIOS2_DIR=/opt/adios2/lib/cmake/adios2" in command


def test_bundledSystemCarriesHdf5Paths():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="hdf5",
        compute_backend="Host_Cpu_CpuSerial",
        bundled_adios2=BUNDLED_ADIOS2_OFF,
        bundled_hdf5=BUNDLED_HDF5_SYSTEM,
        hdf5_prefix="/opt/hdf5",
        hdf5_dir="/opt/hdf5/lib/cmake/hdf5",
    )

    command = install_command(selection)

    assert "-DHASE_OPENPMD_FETCH_HDF5=OFF" in command
    assert "-DCMAKE_PREFIX_PATH=/opt/hdf5" in command
    assert "-DHDF5_DIR=/opt/hdf5/lib/cmake/hdf5" in command


def test_systemOpenPmdCarriesProviderPaths():
    selection = WizardSelection(
        provider=PROVIDER_SYSTEM,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        cmake_prefix_path="/opt/openpmd",
        adios2_prefix="/opt/adios2",
        adios2_dir="/opt/adios2/lib/cmake/adios2",
    )

    command = install_command(selection)

    assert "-DHASE_OPENPMD_PROVIDER=system" in command
    assert "-DCMAKE_PREFIX_PATH=/opt/adios2;/opt/openpmd" in command
    assert "-DADIOS2_DIR=/opt/adios2/lib/cmake/adios2" in command


def test_installCommandCanEnableMpiBuildSupport():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        mpi_mode=MPI_MODE_ON,
    )

    assert "-DDISABLE_MPI=OFF" in install_command(selection)


def test_providerNotesIncludeNewGuidanceShape():
    selection = WizardSelection(
        provider=PROVIDER_SYSTEM,
        openpmd_backend="hdf5",
        compute_backend="Host_Cpu_CpuSerial",
        openpmd_adios2=OPENPMD_ADIOS2_NO,
        supported_openpmd_backends=("hdf5",),
    )

    notes = "\n".join(provider_notes(selection))

    assert "configuration file is present under <configured YAML path>" in notes
    assert "Supported openpmd_backends for this choice: hdf5" in notes
    assert "alpaka backend guidance" in notes


def test_alpakaGuidanceExplainsMissingGpuBackends():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
    )

    guidance = alpaka_backend_guidance(selection, alpaka_backends=["Host_Cpu_CpuSerial"])

    assert "Currently available alpaka backends" in guidance
    assert "No GPU alpaka backend is listed" in guidance
    assert "To get NVIDIA/CUDA backends" in guidance
    assert "To get AMD/HIP backends" in guidance


def test_installCommandUsesEnvVarAndSafeDefaults():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
    )

    command = install_command(selection)

    assert 'HASE_CONFIGURE_CMAKE_ARGS="' in command
    assert 'CMAKE_ARGS="$HASE_CONFIGURE_CMAKE_ARGS" python3 -m pip install -v .' in command
    assert "-DDISABLE_MPI=AUTO" in command
    assert "-DHASE_NATIVE_OPTIMIZATIONS=OFF" in command


def test_installCommandCanPassBreakSystemPackages():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
    )

    command = install_command(selection, break_system_packages=True)

    assert "python3 -m pip install -v --break-system-packages ." in command


def test_installCommandCanDisableNativeOptimizations():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        native_optimizations=False,
    )

    assert "-DHASE_NATIVE_OPTIMIZATIONS=OFF" in install_command(selection)


def test_mainWritesDefaultYamlUnderConfig(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        configure_module,
        "available_alpaka_backends",
        lambda: (["Host_Cpu_CpuSerial"], None),
    )

    assert configure_module.main(["--yes", "--provider", PROVIDER_BUNDLED]) == 0

    output_path = tmp_path / configure_module.DEFAULT_PHI_ASE_CONFIG_PATH
    assert output_path.is_file()
    generated = yaml.safe_load(output_path.read_text(encoding="utf-8"))
    assert generated["compute"]["backend"] == "Host_Cpu_CpuSerial"

    output = capsys.readouterr().out
    assert "Wrote config/hase-phiase.yaml" in output
    assert "configuration file is present under config/hase-phiase.yaml" in output


def test_autoProviderFallsBackToBundled(monkeypatch, capsys):
    monkeypatch.setattr(
        configure_module,
        "available_alpaka_backends",
        lambda: (["Host_Cpu_CpuSerial"], None),
    )

    assert configure_module.main(["--yes", "--cmake", "/definitely/missing/cmake", "--output", "-"]) == 0

    captured = capsys.readouterr()
    assert "-DHASE_OPENPMD_PROVIDER=bundled" in captured.out
    assert "auto provider falling back to bundled" in captured.err
    assert "CMake executable '/definitely/missing/cmake' was not found" in captured.err


def test_bundledProviderTemplateSetsFetchedDependencyVersions():
    hdf5_template = Path("cmake/HaseHdf5Provider.cmake.in").read_text(encoding="utf-8")
    adios2_template = Path("cmake/HaseAdios2Provider.cmake.in").read_text(encoding="utf-8")
    openpmd_template = Path("cmake/HaseOpenPmdProvider.cmake.in").read_text(encoding="utf-8")

    assert 'string(REGEX REPLACE "^v" "" ADIOS2_VERSION "@HASE_ADIOS2_GIT_TAG@")' in adios2_template
    assert 'set(ADIOS2_VERSION' in adios2_template
    assert 'string(REGEX REPLACE "^hdf5_" "" HDF5_VERSION "@HASE_HDF5_GIT_TAG@")' in hdf5_template
    assert 'set(HDF5_VERSION' in hdf5_template
    assert "HDF5_PREFER_PARALLEL" in openpmd_template
    assert "HDF5_IS_PARALLEL" in openpmd_template


def test_cmakeMpiDisableUsesCurrentDisableMpiChoice():
    cmake_lists = Path("CMakeLists.txt").read_text(encoding="utf-8")
    provider_script = Path("cmake/findOpenPmd.cmake").read_text(encoding="utf-8")

    assert "set(HASE_MPI_ENABLED OFF)" in cmake_lists
    assert "elseif(DISABLE_MPI STREQUAL \"OFF\")\n    find_package(MPI COMPONENTS CXX REQUIRED)\n    set(HASE_MPI_ENABLED ON)" in cmake_lists
    assert "if(HASE_MPI_ENABLED)\n    target_compile_definitions(hase_core INTERFACE MPI_FOUND)" in cmake_lists
    assert "if(HASE_MPI_ENABLED)\n        set(HASE_OPENPMD_PROVIDER_MPI ON)" in provider_script
    assert "if(HASE_MPI_ENABLED)\n            set(HASE_PROVIDER_USE_MPI ON)" in provider_script


def test_bundledProviderStagesFetchedDependenciesBeforeOpenPmd():
    script = Path("cmake/findOpenPmd.cmake").read_text(encoding="utf-8")
    normalized_script = " ".join(script.split())

    hdf5_stage = "hase_openpmd_run_provider_stage(hdf5 HaseHdf5Provider.cmake.in)"
    adios2_stage = "hase_openpmd_run_provider_stage(adios2 HaseAdios2Provider.cmake.in)"
    openpmd_stage = "hase_openpmd_run_provider_stage(openpmd HaseOpenPmdProvider.cmake.in)"

    assert "is not supported for the installable bundled openPMD provider" not in script
    assert "HASE_OPENPMD_USE_HDF5\n    \"Enable HDF5 support in the HASE-managed openPMD provider\"\n    OFF" in script
    assert "stage_name STREQUAL \"openpmd\" AND HASE_OPENPMD_USE_HDF5 AND HASE_OPENPMD_FETCH_HDF5" in normalized_script
    assert "-DHDF5_DIR=${HASE_OPENPMD_BUNDLED_PREFIX}/lib/cmake/hdf5" in script
    assert script.index(hdf5_stage) < script.index(openpmd_stage)
    assert script.index(adios2_stage) < script.index(openpmd_stage)


def test_hdf5PromptRequiresHdf5(monkeypatch, capsys):
    monkeypatch.setattr("builtins.input", lambda prompt: "")

    choice = configure_module._interactive_bundled_hdf5(BUNDLED_ADIOS2_OFF)

    output = capsys.readouterr().out
    assert choice == BUNDLED_HDF5_FETCH
    assert "  2) off     do not use HDF5" not in output
    assert "  2) system  use an existing HDF5 installation" in output


def test_interactiveMpiNoDisablesMpiBuild(monkeypatch):
    args = configure_module._parse_args(["--provider", PROVIDER_BUNDLED])
    args.interactive = True
    answers = iter(["", "", "", "n", "n"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    monkeypatch.setattr(configure_module, "available_alpaka_backends", lambda: ([], None))

    selection, _alpaka_backends, _alpaka_error, _probe_report = configure_module._build_selection(args)

    assert selection.mpi_mode == MPI_MODE_OFF
    assert "-DDISABLE_MPI=ON" in install_command(selection)


def test_interactiveInstallPromptDefaultsToYes(tmp_path, monkeypatch):
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
    )

    class TtyStdin:
        @staticmethod
        def isatty():
            return True

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(configure_module.sys, "stdin", TtyStdin())
    monkeypatch.setattr(
        configure_module,
        "_build_selection",
        lambda args: (selection, [], None, None),
    )
    monkeypatch.setattr("builtins.input", lambda prompt: "")
    seen = {}

    def fake_run_install(selected, **kwargs):
        seen.update(kwargs)
        return 17

    monkeypatch.setattr(configure_module, "run_install", fake_run_install)

    assert configure_module.main(["--break-system-packages"]) == 17
    assert seen == {"break_system_packages": True, "use_ccache": False}


def test_installCommandCanEnableCcacheLaunchers():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
    )

    command = install_command(selection, use_ccache=True)

    assert "-DCMAKE_C_COMPILER_LAUNCHER=ccache" in command
    assert "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache" in command
    assert "-DCMAKE_CUDA_COMPILER_LAUNCHER=ccache" in command


def test_reinstallFailsWithoutPreviousBuild(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)

    assert configure_module.main(["--reinstall"]) == 1

    assert "--reinstall requires a previous HASE CMake configure/install" in capsys.readouterr().err


def test_reinstallUsesPreviousCmakeCache(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cache_dir = tmp_path / "build" / "cp-test"
    cache_dir.mkdir(parents=True)
    (cache_dir / "CMakeCache.txt").write_text(
        "HASE_OPENPMD_PROVIDER:STRING=bundled\n"
        "HASE_OPENPMD_USE_HDF5:BOOL=OFF\n"
        "DISABLE_MPI:STRING=AUTO\n"
        "HASE_NATIVE_OPTIMIZATIONS:BOOL=OFF\n",
        encoding="utf-8",
    )
    seen = {}

    def fake_run(cmake_arg_list, **kwargs):
        seen["cmake_arg_list"] = cmake_arg_list
        seen.update(kwargs)
        return 13

    monkeypatch.setattr(configure_module, "_run_install_with_cmake_args", fake_run)

    assert configure_module.main(["--reinstall"]) == 13
    assert seen["cmake_arg_list"] == [
        "-DHASE_OPENPMD_PROVIDER=bundled",
        "-DHASE_OPENPMD_USE_HDF5=OFF",
        "-DDISABLE_MPI=AUTO",
        "-DHASE_NATIVE_OPTIMIZATIONS=OFF",
    ]
    assert seen["record_state"] is True


def test_autoinstallUsesDefaultsAndRunsInstall(tmp_path, monkeypatch):
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        configure_module,
        "_build_selection",
        lambda args: (selection, [], None, None),
    )
    seen = {}

    def fake_run_install(selected, **kwargs):
        seen["selected"] = selected
        seen.update(kwargs)
        return 19

    monkeypatch.setattr(configure_module, "run_install", fake_run_install)

    assert configure_module.main(["--autoinstall"]) == 19
    assert seen == {"selected": selection, "break_system_packages": False, "use_ccache": False}


def test_generatedYamlLoadsAsPhiAseConfig():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios",
        compute_backend="Host_Cpu_CpuSerial",
        parallel_mode="mpi",
        num_devices=2,
        n_per_node=2,
    )

    config = yaml.safe_load(yaml_config(selection))
    phi_ase = PhiASE(config)

    assert phi_ase.backend == "Host_Cpu_CpuSerial"
    assert phi_ase.openpmdBackend == "adios"
    assert phi_ase.parallelMode == "mpi"
    assert phi_ase.numDevices == 2
    assert phi_ase.nPerNode == 2
