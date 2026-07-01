# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Compiled time-integration descriptors for C++/Alpaka simulation runs."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class TimeIntegrationSolver:
    """Descriptor naming a compiled beta time integrator.

    Python no longer performs time-step algorithms. ``Simulation`` serializes
    this descriptor to openPMD run-control attributes and the C++/Alpaka backend
    executes the selected integrator.
    """

    name: str


class ExplicitEuler(TimeIntegrationSolver):
    """First-order explicit Euler compiled integrator descriptor."""

    def __init__(self):
        super().__init__("explicit-euler")


class Heun(TimeIntegrationSolver):
    """Second-order predictor-corrector compiled integrator descriptor."""

    def __init__(self):
        super().__init__("heun")


class Midpoint(TimeIntegrationSolver):
    """Second-order midpoint compiled integrator descriptor."""

    def __init__(self):
        super().__init__("midpoint")


class RungeKutta4(TimeIntegrationSolver):
    """Classical fourth-order Runge-Kutta compiled integrator descriptor."""

    def __init__(self):
        super().__init__("runge-kutta-4")


class ImplicitEuler(TimeIntegrationSolver):
    """Fixed-iteration implicit Euler compiled integrator descriptor."""

    def __init__(self, iterations=8, tolerance=1e-10):
        super().__init__("implicit-euler")
        self.iterations = int(iterations)
        self.tolerance = float(tolerance)


class ExponentialEuler(TimeIntegrationSolver):
    """Compiled exponential Euler descriptor with analytical fluorescence decay."""

    def __init__(self):
        super().__init__("exponential-euler")
