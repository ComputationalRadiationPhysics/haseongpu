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

    struct SimulationExecutionMode
    {
        static inline std::string const AUTONOMOUS = "autonomous";
        static inline std::string const SYNCHRONIZED_DEBUG = "synchronized-debug";
    };

    struct SimulationOutputField
    {
        static inline std::string const BETA_VOLUME = "beta_volume";
        static inline std::string const PHI_ASE = "phi_ase";
        static inline std::string const STANDARD_ERROR = "standard_error";
        static inline std::string const RELATIVE_STANDARD_ERROR = "relative_standard_error";
        static inline std::string const TOTAL_RAYS = "total_rays";
        static inline std::string const DNDT_ASE = "dndt_ase";
        static inline std::string const DNDT_PUMP = "dndt_pump";

        [[nodiscard]] static std::vector<std::string> all()
        {
            return {BETA_VOLUME, PHI_ASE, STANDARD_ERROR, RELATIVE_STANDARD_ERROR, TOTAL_RAYS, DNDT_ASE, DNDT_PUMP};
        }
    };

    struct SimulationControlField
    {
        static inline std::string const BETA_VOLUME = "beta_volume";

        [[nodiscard]] static std::vector<std::string> all()
        {
            return {BETA_VOLUME};
        }
    };

    struct SimulationRunControl
    {
        double timeStep = 0.0;
        unsigned numberOfSteps = 0u;
        bool enableAse = true;
        bool prePump = false;
        unsigned pumpSteps = std::numeric_limits<unsigned>::max();
        std::string executionMode = SimulationExecutionMode::AUTONOMOUS;
        std::vector<unsigned> outputSteps;
        std::vector<std::string> outputFields = SimulationOutputField::all();
        std::vector<std::string> controlFields;
        TimeIntegrationParameters timeIntegration;
        PumpParameters pump;
    };
} // namespace hase::core
