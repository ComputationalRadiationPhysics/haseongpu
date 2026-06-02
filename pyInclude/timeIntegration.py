# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Time-integration solvers for the point-level beta population arrays."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TimeDerivative:
    """Derivative evaluation produced during one time-integration stage.

    ``Simulation`` builds this object from pump gain, ASE depletion, and
    fluorescence decay. Custom solvers receive it from their ``rhs`` callback
    while trying candidate ``betaCells`` values.
    """

    betaCells: np.ndarray
    """Point-level excited-state fraction used for this RHS evaluation."""
    dndtPump: np.ndarray
    """Pump source contribution to ``d beta / dt`` in ``1 / s``."""
    dndtAse: np.ndarray
    """ASE depletion contribution to ``d beta / dt`` in ``1 / s``."""
    derivative: np.ndarray
    """Total ``d beta / dt`` advanced by the time solver."""
    tau: float
    """Fluorescence lifetime used for spontaneous decay, in seconds."""
    phiAse: np.ndarray | None = None
    """ASE flux :math:`\Phi_i` from the derivative evaluation, if available."""
    aseResult: object | None = None
    """Raw lower-level ASE result object, if available."""


@dataclass(frozen=True)
class TimeIntegrationResult:
    """Return object from ``TimeIntegrationSolver.step(...)``."""

    betaCells: np.ndarray
    """Updated point-level beta array proposed by the solver."""
    evaluation: TimeDerivative
    """Derivative evaluation to expose in the resulting ``TimeStepState``."""


class TimeIntegrationSolver:
    """Base protocol for custom beta time-integration solvers.

    Implement ``step(rhs, betaCells, time, timeStep)``. ``rhs`` is a callable
    that accepts ``(betaCells, time)`` and returns a ``TimeDerivative``. Higher
    order solvers may call ``rhs`` several times for one public simulation step.
    """

    name = "base"

    def step(self, rhs, betaCells, time, timeStep):
        """Return a ``TimeIntegrationResult`` for one physical time step.

        ``time`` and ``timeStep`` are in seconds. ``betaCells`` has the same
        shape as ``GainMedium.get("betaCells").expectedShape``.
        """
        raise NotImplementedError


class ExplicitEuler(TimeIntegrationSolver):
    """First-order explicit Euler beta update using one RHS evaluation."""

    name = "explicit-euler"

    def step(self, rhs, betaCells, time, timeStep):
        evaluation = rhs(betaCells, time)
        return TimeIntegrationResult(betaCells + timeStep * evaluation.derivative, evaluation)


class Heun(TimeIntegrationSolver):
    """Second-order predictor-corrector beta update using two RHS evaluations."""

    name = "heun"

    def step(self, rhs, betaCells, time, timeStep):
        first = rhs(betaCells, time)
        second = rhs(betaCells + timeStep * first.derivative, time + timeStep)
        updated = betaCells + 0.5 * timeStep * (first.derivative + second.derivative)
        return TimeIntegrationResult(updated, second)


class Midpoint(TimeIntegrationSolver):
    """Second-order midpoint beta update using a half-step RHS evaluation."""

    name = "midpoint"

    def step(self, rhs, betaCells, time, timeStep):
        first = rhs(betaCells, time)
        middle = rhs(betaCells + 0.5 * timeStep * first.derivative, time + 0.5 * timeStep)
        return TimeIntegrationResult(betaCells + timeStep * middle.derivative, middle)


class RungeKutta4(TimeIntegrationSolver):
    """Classical fourth-order Runge-Kutta beta update using four RHS evaluations."""

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
    """Fixed-point implicit Euler beta update for stiff beta dynamics."""

    name = "implicit-euler"

    def __init__(self, iterations=8, tolerance=1e-10):
        """Set maximum fixed-point iterations and infinity-norm tolerance."""
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
    """Euler source update with analytical fluorescence decay over ``timeStep``."""

    name = "exponential-euler"

    def step(self, rhs, betaCells, time, timeStep):
        evaluation = rhs(betaCells, time)
        decay = np.exp(-timeStep / evaluation.tau)
        source = evaluation.dndtPump - evaluation.dndtAse
        updated = evaluation.tau * source * (1.0 - decay) + betaCells * decay
        return TimeIntegrationResult(updated, evaluation)
