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
    PlanarPumpRelay,
    PumpProperties,
    PumpSource,
    PumpSpectrum,
    Simulation,
    SuperGaussianPumpProfile,
    integratePumpProfile,
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
    profile = SuperGaussianPumpProfile(radiusU=1.5, radiusV=1.5, exponent=40)
    source = PumpSource(
        surfaceDomains=("ase_bottom",),
        totalPower=16e3 * integratePumpProfile(medium.topology, "ase_bottom", profile),
        spectrum=PumpSpectrum.monochromatic(wavelength),
        crossSections=pump_cross_sections,
        profile=profile,
        relays=(PlanarPumpRelay.retroreflect("ase_top"),),
    )
    spectral = example.laserPumpCladdingSpectralProperties(191)
    phi_ase = PhiASE.fromYaml(
        example.defaultPhiAseConfigPath,
        spectralProperties=spectral,
        backend=alpakaRuntimeBackend,
        openpmdBackend=openPmdFileBackend,
    )
    simulation = Simulation(
        gainMedium=medium,
        pump=PumpProperties((source,), rayCount=50_000, rngSeed=5489, pumpSteps=3),
        phiASE=phi_ase,
        timeIntegrationSolver=FrozenPhiAseRungeKutta4(),
        timeStep=2e-5,
        crossSections=spectral,
        enableASE=False,
        prePump=True,
    )
    states = []
    simulation.onStep(states.append).runSteps(3)

    beta_volume = np.stack([np.asarray(state.betaVolume) for state in states])
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
    new_total = np.asarray([np.sum(np.asarray(state.dndtPump) * lumped_volume) for state in states])
    old_total = np.asarray([np.sum(values * lumped_volume) for values in reference["dndtPump"]])
    np.testing.assert_allclose(new_total, old_total, rtol=0.01, atol=1e-12)
