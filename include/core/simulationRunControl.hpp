/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace hase::core
{
    struct PumpProfileParameters
    {
        unsigned kind = 0u;
        double radiusU = 1.0;
        double radiusV = 1.0;
        double exponent = 2.0;
        double center[3] = {0.0, 0.0, 0.0};
        double axisU[3] = {1.0, 0.0, 0.0};
        double axisV[3] = {0.0, 1.0, 0.0};
    };

    struct PumpRelayParameters
    {
        std::vector<int> exitSurfaces;
        std::vector<int> entrySurfaces;
        bool flipU = false;
        bool flipV = false;
        double rotation = 0.0;
        double offset[2] = {0.0, 0.0};
        double tilt[2] = {0.0, 0.0};
        double magnification = 1.0;
        double transmission = 1.0;
    };

    struct PumpSourceParameters
    {
        std::vector<int> surfaces;
        double totalPower = 0.0;
        std::vector<double> wavelengths;
        std::vector<double> spectralWeights;
        std::vector<double> sigmaAbsorption;
        std::vector<double> sigmaEmission;
        std::vector<double> polarAngles;
        std::vector<double> azimuthalAngles;
        std::vector<double> angularWeights;
        PumpProfileParameters profile;
        std::vector<PumpRelayParameters> relays;
    };

    struct PumpParameters
    {
        unsigned schemaVersion = 1u;
        unsigned rayCount = 100000u;
        std::uint32_t rngSeed = 5489u;
        std::vector<PumpSourceParameters> sources;
    };

    struct TimeIntegrator
    {
        static inline std::string const EXPLICIT_EULER = "explicit-euler";
        static inline std::string const HEUN = "heun";
        static inline std::string const MIDPOINT = "midpoint";
        static inline std::string const RUNGE_KUTTA_4 = "runge-kutta-4";
        static inline std::string const FROZEN_PHI_ASE_RUNGE_KUTTA_4 = "frozen-phi-ase-runge-kutta-4";
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
        bool enableAse = true;
        bool prePump = false;
        unsigned pumpSteps = std::numeric_limits<unsigned>::max();
        TimeIntegrationParameters timeIntegration;
        PumpParameters pump;
    };
} // namespace hase::core
