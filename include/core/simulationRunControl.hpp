/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <limits>
#include <string>

namespace hase::core
{
    struct PumpRoutine
    {
        static inline std::string const ONE_DIMENSIONAL_Z_TRAVERSAL = "one-dimensional-z-traversal";
    };

    struct PumpParameters
    {
        std::string routine = PumpRoutine::ONE_DIMENSIONAL_Z_TRAVERSAL;
        double intensity = 0.0;
        double wavelength = 0.0;
        double radiusX = 0.0;
        double radiusY = 0.0;
        double exponent = 40.0;
        double duration = 0.0;
        double sigmaAbsorption = 0.0;
        double sigmaEmission = 0.0;
        double temporaryFluorescence = 0.0;
        unsigned substeps = 100u;
        bool backReflection = true;
        bool extraction = false;
        double reflectivity = 1.0;
    };

    struct TimeIntegrator
    {
        static inline std::string const EXPLICIT_EULER = "explicit-euler";
        static inline std::string const HEUN = "heun";
        static inline std::string const MIDPOINT = "midpoint";
        static inline std::string const RUNGE_KUTTA_4 = "runge-kutta-4";
        static inline std::string const IMPLICIT_EULER = "implicit-euler";
        static inline std::string const EXPONENTIAL_EULER = "exponential-euler";
    };

    struct TimeIntegrationParameters
    {
        std::string method = TimeIntegrator::EXPLICIT_EULER;
        unsigned implicitIterations = 8u;
        double implicitTolerance = 1.0e-10;
    };

    struct SimulationRunControl
    {
        double timeStep = 0.0;
        unsigned numberOfSteps = 0u;
        unsigned pumpSteps = std::numeric_limits<unsigned>::max();
        TimeIntegrationParameters timeIntegration;
        PumpParameters pump;
    };
} // namespace hase::core
