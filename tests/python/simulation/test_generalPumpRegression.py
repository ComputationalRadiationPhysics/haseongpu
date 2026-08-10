# Copyright 2026 Tim Hanel
# SPDX-License-Identifier: GPL-3.0-or-later

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[3] / "example"))

import numpy as np
import pytest

from HASEonGPU import (
    CrossSectionData,
    FrozenPhiAseRungeKutta4,
    PhiASE,
    MonteCarloPumpSolver,
    PlanarPumpRelay,
    Pump,
    PumpSpectrum,
    Simulation,
    SuperGaussianPumpProfile,
    SurfacePumpInjector,
    integrate_pump_profile,
)
import laserPumpCladding as example


def _normalized_wasserstein_distance(first_coordinates, first_weights, second_coordinates, second_weights):
    first_coordinates = np.asarray(first_coordinates, dtype=np.float64).reshape(-1)
    second_coordinates = np.asarray(second_coordinates, dtype=np.float64).reshape(-1)
    first_weights = np.asarray(first_weights, dtype=np.float64).reshape(-1)
    second_weights = np.asarray(second_weights, dtype=np.float64).reshape(-1)
    first_weights = first_weights / np.sum(first_weights)
    second_weights = second_weights / np.sum(second_weights)

    coordinates = np.unique(np.concatenate((first_coordinates, second_coordinates)))
    first_order = np.argsort(first_coordinates)
    second_order = np.argsort(second_coordinates)
    first_cdf = np.concatenate(([0.0], np.cumsum(first_weights[first_order])))
    second_cdf = np.concatenate(([0.0], np.cumsum(second_weights[second_order])))
    first_indices = np.searchsorted(first_coordinates[first_order], coordinates[:-1], side="right")
    second_indices = np.searchsorted(second_coordinates[second_order], coordinates[:-1], side="right")
    distance = np.sum(
        np.abs(first_cdf[first_indices] - second_cdf[second_indices]) * np.diff(coordinates)
    )
    coordinate_span = coordinates[-1] - coordinates[0]
    return float(distance / coordinate_span) if coordinate_span > 0.0 else 0.0


def _deposition_diagnostics(topology, states, reference, legacy_lumped_volume):
    point_coordinates = np.asarray(topology.points, dtype=np.float64)
    cell_coordinates = np.asarray(topology.cellCenters, dtype=np.float64)
    cell_volumes = np.asarray(topology.cellVolumes, dtype=np.float64)
    point_volumes = np.asarray(legacy_lumped_volume, dtype=np.float64).reshape(-1, order="F")
    point_radius = np.linalg.norm(point_coordinates[:, :2], axis=1)
    cell_radius = np.linalg.norm(cell_coordinates[:, :2], axis=1)
    diagnostics = []
    for step, (state, legacy_rate) in enumerate(zip(states, reference["dndtPump"], strict=True), start=1):
        cell_measure = np.asarray(state.dndt_pump, dtype=np.float64).reshape(-1) * cell_volumes
        point_measure = np.asarray(legacy_rate, dtype=np.float64).reshape(-1, order="F") * point_volumes
        cell_total = np.sum(cell_measure)
        point_total = np.sum(point_measure)
        diagnostics.append(
            {
                "step": step,
                "total_relative_error": float(abs(cell_total - point_total) / point_total),
                "axial_profile_distance": _normalized_wasserstein_distance(
                    cell_coordinates[:, 2], cell_measure, point_coordinates[:, 2], point_measure
                ),
                "radial_profile_distance": _normalized_wasserstein_distance(
                    cell_radius, cell_measure, point_radius, point_measure
                ),
                "cell_mean_z": float(np.sum(cell_coordinates[:, 2] * cell_measure) / cell_total),
                "legacy_mean_z": float(np.sum(point_coordinates[:, 2] * point_measure) / point_total),
                "cell_mean_radius": float(np.sum(cell_radius * cell_measure) / cell_total),
                "legacy_mean_radius": float(np.sum(point_radius * point_measure) / point_total),
            }
        )
    return diagnostics


@pytest.mark.integration
def test_general_pump_reproduces_legacy_crystal_inversion(openPmdFileBackend, alpakaRuntimeBackend):
    reference = np.load(
        Path(__file__).parents[2] / "data" / "pump" / "legacy_one_dimensional_reference.npz"
    )
    wavelength = 940e-9
    lambda_a, sigma_a, lambda_e, sigma_e = example._loadLaserPumpCladdingRawSpectra()
    pump_cross_sections = CrossSectionData.monochromatic(
        wavelength=wavelength,
        crossSectionAbsorption=np.interp(wavelength * 1e9, lambda_a, sigma_a),
        crossSectionEmission=np.interp(wavelength * 1e9, lambda_e, sigma_e),
    )
    medium = example.laserPumpCladdingMedium(cladAbsorption=5.5)
    profile = SuperGaussianPumpProfile(radius_u=1.5, radius_v=1.5, exponent=40)
    pump = Pump(
        total_power=16e3 * integrate_pump_profile(medium.topology, "ase_bottom", profile),
        spectrum=PumpSpectrum.monochromatic(wavelength),
        cross_sections=pump_cross_sections,
        profile=profile,
    )
    spectral = example.laserPumpCladdingSpectralProperties(191)
    phi_ase = PhiASE.fromYaml(
        example.defaultPhiAseConfigPath,
        spectralProperties=spectral,
        backend=alpakaRuntimeBackend,
        openpmdBackend=openPmdFileBackend,
    )
    simulation = Simulation(
        gain_medium=medium,
        phi_ase=phi_ase,
        time_integrator=FrozenPhiAseRungeKutta4(),
        time_step_size=2e-5,
        pump_solver=MonteCarloPumpSolver(ray_count=50_000, seed=5489, max_steps=3),
        cross_sections=spectral,
        enable_ase=False,
        pre_pump=True,
    ).add_pump(
        pump,
        injection_method=SurfacePumpInjector(surface_domains="ase_bottom"),
        relays=(PlanarPumpRelay.retroreflect("ase_top"),),
    )
    states = []
    simulation.on_step(states.append).step(3)

    beta_volume = np.stack([np.asarray(state.beta_volume) for state in states])
    relative_field_error = np.linalg.norm(beta_volume - reference["betaVolume"]) / np.linalg.norm(
        reference["betaVolume"]
    )

    cell_volumes = np.asarray(medium.topology.cellVolumes)
    legacy_lumped_volume = np.bincount(
        np.asarray(medium.topology.cellPointIndices).reshape(-1),
        weights=np.repeat(cell_volumes / 4.0, 4),
        minlength=medium.topology.numberOfPoints,
    ).reshape(reference["dndtPump"].shape[1:], order="F")
    new_total = np.asarray([np.sum(np.asarray(state.dndt_pump) * cell_volumes) for state in states])
    old_total = np.asarray(
        [np.sum(values * legacy_lumped_volume) for values in reference["dndtPump"]]
    )
    np.testing.assert_allclose(new_total, old_total, rtol=0.01, atol=1e-12)
    diagnostics = _deposition_diagnostics(medium.topology, states, reference, legacy_lumped_volume)
    assert relative_field_error < 0.05, f"deposition diagnostics: {diagnostics}"
