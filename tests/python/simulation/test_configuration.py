# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path

import pytest
import yaml

from HASEonGPU import ImplicitEuler, Simulation, SuperGaussianPumpProfile


def _config():
    return {
        "schema_version": 2,
        "cross_sections": {
            "ase": {
                "inline": {
                    "wavelengths_absorption": [900e-9, 910e-9],
                    "cross_section_absorption": [1.0e-21, 1.1e-21],
                    "wavelengths_emission": [1020e-9, 1030e-9],
                    "cross_section_emission": [2.0e-20, 2.1e-20],
                    "resolution": 2,
                }
            },
            "pump": {
                "monochromatic": {
                    "wavelength": 940e-9,
                    "cross_section_absorption": 1.2e-21,
                    "cross_section_emission": 0.0,
                }
            },
        },
        "simulation": {
            "gain_medium": {
                "topology": {
                    "from_tetrahedra": {
                        "points": [
                            [0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0],
                        ],
                        "cell_point_indices": [[0, 1, 2, 3]],
                        "face_boundaries": [[1, 1, 1, 1]],
                    }
                },
                "properties": {
                    "beta_volume": [0.0],
                    "cladding_cell_types": [0],
                    "n_tot": 2.76e20,
                    "fluorescence_lifetime": 9.5e-4,
                    "cladding_number": 1,
                    "cladding_absorption": 0.0,
                },
                "surface_optics": {
                    1: {"reflectivity": 0.75, "n_inside": 1.8, "n_outside": 1.0}
                },
                "custom_fields": [
                    {
                        "name": "temperature",
                        "entity": "cell",
                        "values": [300.0],
                        "dtype": "float64",
                        "unit": "K",
                    }
                ],
            },
            "phi_ase": {
                "cross_sections": "ase",
                "propagation_mode": "forward",
                "min_rays": 10,
                "max_rays": 100,
                "forward_ray_count": 12,
                "relative_standard_error_threshold": 0.05,
                "repetitions": 2,
                "adaptive_steps": 3,
                "use_reflections": True,
                "reflection_max_iterations": 7,
                "reflection_tolerance": 1.0e-5,
                "surface_reservoir_size": 9,
                "monochromatic": False,
                "backend": "Host_Cpu_CpuSerial",
                "openpmd_backend": "auto",
                "parallel_mode": "single",
                "num_devices": 1,
                "n_per_node": 1,
                "min_sample_range": 0,
                "max_sample_range": 0,
                "rng_seed": 123,
                "ase_steps": 3,
            },
            "pumps": [
                {
                    "name": "main",
                    "total_power": 10.0,
                    "ray_count": 1234,
                    "pump_steps": 2,
                    "rng_seed": 99,
                    "cross_sections": "pump",
                    "spectrum": {"monochromatic": 940e-9},
                    "angular_distribution": {
                        "uniform_cone": {
                            "half_angle": 0.1,
                            "polar_samples": 2,
                            "azimuthal_samples": 3,
                        }
                    },
                    "profile": {
                        "kind": "super_gaussian",
                        "radius_u": 1.5,
                        "radius_v": 1.25,
                        "exponent": 40.0,
                        "center": [0.0, 0.0, 0.0],
                        "axis_u": [1.0, 0.0, 0.0],
                        "axis_v": [0.0, 1.0, 0.0],
                    },
                    "injection": {"surface_domains": [1]},
                    "relays": [
                        {
                            "exit_domains": [1],
                            "entry_domains": [1],
                            "flip_u": True,
                            "flip_v": False,
                            "rotation": 0.2,
                            "offset": [0.1, 0.2],
                            "tilt": [0.01, 0.02],
                            "magnification": 1.5,
                            "transmission": 0.8,
                        }
                    ],
                }
            ],
            "time_integrator": {
                "method": "implicit_euler",
                "iterations": 3,
                "tolerance": 1.0e-8,
            },
            "time_step_size": 2.0e-5,
            "simulation_steps": 3,
            "pre_pump": True,
            "report_timings": True,
            "execution_mode": "autonomous",
            "output_steps": [2, 3],
            "output_fields": ["beta_volume", "phi_ase", "dndt_pump"],
            "control_fields": [],
        },
    }


def _write(path: Path, config):
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


def testSimulationFromYamlBuildsPublicObjectGraph(tmp_path):
    simulation = Simulation.from_yaml(_write(tmp_path / "simulation.yaml", _config()))

    assert simulation.simulation_steps == 3
    assert isinstance(simulation.time_integrator, ImplicitEuler)
    assert simulation.time_integrator.iterations == 3
    assert simulation.phi_ase.useReflections is True
    assert simulation.phi_ase.forwardRayCount == 12
    assert simulation.phi_ase.ase_steps == 3
    assert simulation.pre_pump is True
    assert simulation.gain_medium.surface_optics[1].reflectivity == pytest.approx(0.75)
    assert simulation.gain_medium.getField("temperature").value[0] == pytest.approx(300.0)
    assert len(simulation.pumps) == 1
    assert simulation.pumps[0].name == "main"
    assert simulation.pumps[0].ray_count == 1234
    assert simulation.pumps[0].pump_steps == 2
    assert simulation.pumps[0].rng_seed == 99
    assert isinstance(simulation.pumps[0].profile, SuperGaussianPumpProfile)
    assert simulation.outputSteps == (2, 3)
    assert simulation.outputFields == ("beta_volume", "phi_ase", "dndt_pump")


@pytest.mark.parametrize(
    "name",
    [
        "refractive_indices",
        "reflectivities",
        "surface_reflectivity",
        "surface_refractive_index_inside",
        "surface_refractive_index_outside",
    ],
)
def testSchemaV2RejectsRawOpticsProperties(tmp_path, name):
    config = _config()
    config["simulation"]["gain_medium"]["properties"][name] = []

    with pytest.raises(ValueError, match="raw optics fields are not public"):
        Simulation.from_yaml(_write(tmp_path / "simulation.yaml", config))


def testSchemaV2RejectsRemovedPumpSolverSection(tmp_path):
    config = _config()
    config["simulation"]["pump_solver"] = {"ray_count": 12}

    with pytest.raises(ValueError, match="unsupported simulation options"):
        Simulation.from_yaml(_write(tmp_path / "simulation.yaml", config))


def testSchemaV2DerivesRunLimitWhenSimulationStepsIsOmitted(tmp_path):
    config = _config()
    config["simulation"].pop("simulation_steps")
    simulation = Simulation.from_yaml(_write(tmp_path / "simulation.yaml", config))
    assert simulation._derived_simulation_steps() == 3


def testSchemaV2RejectsTwoRunLimits(tmp_path):
    config = _config()
    config["simulation"]["max_time"] = 1.0

    with pytest.raises(ValueError, match="at most one of: simulation_steps, max_time"):
        Simulation.from_yaml(_write(tmp_path / "simulation.yaml", config))


def testSchemaV2RejectsUnknownSurfaceOpticsOptions(tmp_path):
    config = _config()
    config["simulation"]["gain_medium"]["surface_optics"][1]["refractive_indices"] = [1.0]

    with pytest.raises(ValueError, match="unsupported simulation.gain_medium.surface_optics.1 options"):
        Simulation.from_yaml(_write(tmp_path / "simulation.yaml", config))
