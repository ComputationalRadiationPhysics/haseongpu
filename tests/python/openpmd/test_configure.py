import yaml

from pyInclude import PhiASE
from pyInclude import config as configure_module
from pyInclude.config import (
    BUNDLED_ADIOS2_FETCH,
    BUNDLED_ADIOS2_OFF,
    BUNDLED_ADIOS2_SYSTEM,
    OPENPMD_ADIOS2_NO,
    PROVIDER_BUNDLED,
    PROVIDER_EXTERNAL,
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


def test_supported_external_openpmd_backends_intersects_python_and_cmake_support():
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


def test_preferred_openpmd_backend_uses_priority_and_validates_request():
    assert preferred_openpmd_backend(["adios", "hdf5"]) == "adios"
    assert preferred_openpmd_backend(["adios", "hdf5"], "hdf5") == "hdf5"


def test_bundled_adios2_off_is_hdf5_only():
    assert bundled_supported_openpmd_backends(BUNDLED_ADIOS2_OFF) == ["hdf5"]
    assert bundled_supported_openpmd_backends(BUNDLED_ADIOS2_FETCH) == ["adios-sst", "adios", "hdf5"]


def test_preferred_compute_backend_uses_cpu_backend_for_first_validation():
    assert preferred_compute_backend(["Cuda_Gpu_CudaRt", "Host_Cpu_CpuSerial"]) == "Host_Cpu_CpuSerial"
    assert preferred_compute_backend(["Cuda_Gpu_CudaRt"]) == "Cuda_Gpu_CudaRt"


def test_bundled_fetch_install_command_mentions_source_build_and_long_compile_note():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        bundled_adios2=BUNDLED_ADIOS2_FETCH,
    )

    command = install_command(selection)

    assert "-DHASE_BUILD_OPENPMD_FROM_SOURCE=ON" in command
    assert "-DHASE_OPENPMD_BUILD_PYTHON_BINDINGS=ON" in command


def test_bundled_install_command_can_skip_openpmd_python_bindings():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        bundled_adios2=BUNDLED_ADIOS2_FETCH,
        bundled_python_bindings=False,
    )

    assert "-DHASE_OPENPMD_BUILD_PYTHON_BINDINGS=ON" not in install_command(selection)


def test_bundled_hdf5_only_install_command_disables_adios2():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="hdf5",
        compute_backend="Host_Cpu_CpuSerial",
        bundled_adios2=BUNDLED_ADIOS2_OFF,
    )

    command = install_command(selection)

    assert "-DHASE_OPENPMD_USE_ADIOS2=OFF" in command
    assert "-DHASE_OPENPMD_USE_HDF5=ON" in command


def test_bundled_system_adios2_install_command_carries_adios2_paths():
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


def test_external_openpmd_with_adios2_install_command_carries_provider_paths():
    selection = WizardSelection(
        provider=PROVIDER_EXTERNAL,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        cmake_prefix_path="/opt/openpmd",
        adios2_prefix="/opt/adios2",
        adios2_dir="/opt/adios2/lib/cmake/adios2",
    )

    command = install_command(selection)

    assert "-DCMAKE_PREFIX_PATH=/opt/adios2;/opt/openpmd" in command
    assert "-DADIOS2_DIR=/opt/adios2/lib/cmake/adios2" in command


def test_provider_notes_include_new_guidance_shape():
    selection = WizardSelection(
        provider=PROVIDER_EXTERNAL,
        openpmd_backend="hdf5",
        compute_backend="Host_Cpu_CpuSerial",
        openpmd_adios2=OPENPMD_ADIOS2_NO,
        supported_openpmd_backends=("hdf5",),
    )

    notes = "\n".join(provider_notes(selection))

    assert "configuration file is present under <configured YAML path>" in notes
    assert "Supported openpmd_backends for this choice: hdf5" in notes
    assert "alpaka backend guidance" in notes


def test_alpaka_guidance_explains_missing_gpu_backends():
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


def test_install_command_uses_env_var_and_enables_native_optimizations_by_default():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
    )

    command = install_command(selection)

    assert 'HASE_CONFIGURE_CMAKE_ARGS="' in command
    assert 'CMAKE_ARGS="$HASE_CONFIGURE_CMAKE_ARGS" python3 -m pip install .' in command
    assert "-DHASE_NATIVE_OPTIMIZATIONS=ON" in command


def test_install_command_can_pass_break_system_packages():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
    )

    command = install_command(selection, break_system_packages=True)

    assert "python3 -m pip install --break-system-packages ." in command


def test_install_command_can_disable_native_optimizations():
    selection = WizardSelection(
        provider=PROVIDER_BUNDLED,
        openpmd_backend="adios-sst",
        compute_backend="Host_Cpu_CpuSerial",
        native_optimizations=False,
    )

    assert "-DHASE_NATIVE_OPTIMIZATIONS=OFF" in install_command(selection)


def test_main_writes_default_yaml_under_config(tmp_path, monkeypatch, capsys):
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


def test_interactive_install_prompt_defaults_to_yes(tmp_path, monkeypatch):
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
    assert seen == {"break_system_packages": True}


def test_generated_yaml_loads_as_phi_ase_config():
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
