# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TimeDerivative:
    betaCells: np.ndarray
    dndtPump: np.ndarray
    dndtAse: np.ndarray
    derivative: np.ndarray
    tau: float
    phiAse: np.ndarray | None = None
    aseResult: object | None = None


@dataclass(frozen=True)
class TimeIntegrationResult:
    betaCells: np.ndarray
    evaluation: TimeDerivative


class TimeIntegrationSolver:
    name = "base"

    def step(self, rhs, betaCells, time, timeStep):
        raise NotImplementedError


class ExplicitEuler(TimeIntegrationSolver):
    name = "explicit-euler"

    def step(self, rhs, betaCells, time, timeStep):
        evaluation = rhs(betaCells, time)
        return TimeIntegrationResult(betaCells + timeStep * evaluation.derivative, evaluation)


class Heun(TimeIntegrationSolver):
    name = "heun"

    def step(self, rhs, betaCells, time, timeStep):
        first = rhs(betaCells, time)
        second = rhs(betaCells + timeStep * first.derivative, time + timeStep)
        updated = betaCells + 0.5 * timeStep * (first.derivative + second.derivative)
        return TimeIntegrationResult(updated, second)


class Midpoint(TimeIntegrationSolver):
    name = "midpoint"

    def step(self, rhs, betaCells, time, timeStep):
        first = rhs(betaCells, time)
        middle = rhs(betaCells + 0.5 * timeStep * first.derivative, time + 0.5 * timeStep)
        return TimeIntegrationResult(betaCells + timeStep * middle.derivative, middle)


class RungeKutta4(TimeIntegrationSolver):
    name = "runge-kutta-4"

    def step(self, rhs, betaCells, time, timeStep):
        k1 = rhs(betaCells, time)
        k2 = rhs(betaCells + 0.5 * timeStep * k1.derivative, time + 0.5 * timeStep)
        k3 = rhs(betaCells + 0.5 * timeStep * k2.derivative, time + 0.5 * timeStep)
        k4 = rhs(betaCells + timeStep * k3.derivative, time + timeStep)
        updated = betaCells + (timeStep / 6.0) * (
            k1.derivative + 2.0 * k2.derivative + 2.0 * k3.derivative + k4.derivative
        )
        return TimeIntegrationResult(updated, k4)


class ImplicitEuler(TimeIntegrationSolver):
    name = "implicit-euler"

    def __init__(self, iterations=8, tolerance=1e-10):
        self.iterations = int(iterations)
        self.tolerance = float(tolerance)

    def step(self, rhs, betaCells, time, timeStep):
        guess = np.asarray(betaCells, dtype=np.float64).copy()
        evaluation = rhs(guess, time + timeStep)
        for _ in range(max(self.iterations, 1)):
            previous = guess
            evaluation = rhs(guess, time + timeStep)
            guess = betaCells + timeStep * evaluation.derivative
            if np.linalg.norm(guess - previous, ord=np.inf) <= self.tolerance:
                break
        return TimeIntegrationResult(guess, evaluation)


class ExponentialEuler(TimeIntegrationSolver):
    name = "exponential-euler"

    def step(self, rhs, betaCells, time, timeStep):
        evaluation = rhs(betaCells, time)
        decay = np.exp(-timeStep / evaluation.tau)
        source = evaluation.dndtPump - evaluation.dndtAse
        updated = evaluation.tau * source * (1.0 - decay) + betaCells * decay
        return TimeIntegrationResult(updated, evaluation)
