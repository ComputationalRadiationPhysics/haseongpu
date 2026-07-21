# Copyright 2026 Tim Hanel
# SPDX-License-Identifier: GPL-3.0-or-later

import inspect

import numpy as np
import pytest

from HASEonGPU import (
    CrossSectionData,
    GaussianPump,
    MonteCarloPumpSolver,
    PlanarPumpRelay,
    Pump,
    PumpAngularDistribution,
    PumpSpectrum,
    Simulation,
    SuperGaussianPumpProfile,
    SurfacePumpInjector,
    UniformPumpProfile,
)


def monochromatic_cross_sections():
    return CrossSectionData.monochromatic(
        wavelength=940e-9,
        crossSectionAbsorption=1e-22,
        crossSectionEmission=2e-22,
    )


def test_public_pump_and_simulation_signatures_use_snake_case():
    public_classes = (
        GaussianPump,
        MonteCarloPumpSolver,
        PlanarPumpRelay,
        Pump,
        PumpAngularDistribution,
        PumpSpectrum,
        Simulation,
        SuperGaussianPumpProfile,
        SurfacePumpInjector,
        UniformPumpProfile,
    )
    for cls in public_classes:
        for name in inspect.signature(cls).parameters:
            assert not any(character.isupper() for character in name), (cls.__name__, name)


def test_gaussian_pump_keeps_physics_separate_from_injection_and_solver():
    pump = GaussianPump(
        total_power=12.5,
        spectrum=PumpSpectrum.monochromatic(940e-9),
        cross_sections=monochromatic_cross_sections(),
        waist=(1.5, 1.25),
        exponent=40,
        angular_distribution=PumpAngularDistribution.collimated(),
        name="lower_pump",
    )
    injector = SurfacePumpInjector(surface_domains=("lower",))
    solver = MonteCarloPumpSolver(ray_count=1234, seed=99, max_steps=4)

    assert pump.total_power == 12.5
    assert pump.profile.radius_u == 1.5
    assert pump.profile.radius_v == 1.25
    assert pump.profile.weight_at([[0.0, 0.0, 0.0]])[0] == pytest.approx(1.0)
    np.testing.assert_array_equal(pump.spectrum.weights, [1.0])
    assert injector.surface_domains == ("lower",)
    assert solver == MonteCarloPumpSolver(ray_count=1234, seed=99, max_steps=4)


def test_uniform_cone_uses_snake_case_sampling_controls():
    distribution = PumpAngularDistribution.uniform_cone(
        np.pi / 6.0,
        polar_samples=2,
        azimuthal_samples=3,
    )
    assert distribution.weights.size == 6
    assert distribution.weights.sum() == pytest.approx(1.0)
    assert np.all(distribution.polar_angles < np.pi / 6.0)


def test_simulation_step_signature_matches_picmi_default():
    assert inspect.signature(Simulation.step).parameters["nsteps"].default == 1
