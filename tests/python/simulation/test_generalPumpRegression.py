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
from example import laserPumpCladding as example


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
    assert relative_field_error < 0.05

    cell_points = np.asarray(medium.topology.cellPointIndices).reshape(-1)
    lumped_volume = np.bincount(
        cell_points,
        weights=np.repeat(np.asarray(medium.topology.cellVolumes) / 4.0, 4),
        minlength=medium.topology.numberOfSamplePoints,
    ).reshape(np.asarray(states[0].dndtPump).shape, order="F")
    new_total = np.asarray([np.sum(np.asarray(state.dndt_pump) * lumped_volume) for state in states])
    old_total = np.asarray([np.sum(values * lumped_volume) for values in reference["dndtPump"]])
    np.testing.assert_allclose(new_total, old_total, rtol=0.01, atol=1e-12)
