/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/utils.hpp>
#include <concepts/concepts.hpp>
#include <core/mesh.hpp>
#include <core/simulationRunControl.hpp>

#include <cmath>

namespace hase::kernels
{
    struct OneDimensionalPumpKernelParameters
    {
        double intensity = 0.0;
        double wavelength = 0.0;
        double radiusX = 0.0;
        double radiusY = 0.0;
        double exponent = 40.0;
        double duration = 0.0;
        double sigmaAbsorption = 0.0;
        double sigmaEmission = 0.0;
        double reflectivity = 1.0;
        unsigned substeps = 100u;
        bool backReflection = true;
        bool extraction = false;
    };

    struct OneDimensionalZTraversalPump
    {
        OneDimensionalPumpKernelParameters pump;
        double nTot = 0.0;
        double crystalLength = 0.0;
        double speedOfLight = 299792458.0;
        double planckConstant = 6.62607015e-34;

        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            auto betaBefore,
            auto betaAfter,
            auto pumpForward,
            auto pumpBackward) const
        {
            for(auto [point] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfPoints}))
            {
                for(unsigned level = 0u; level < mesh.numberOfLevels; ++level)
                {
                    betaAfter[point + level * mesh.numberOfPoints] = betaBefore[point + level * mesh.numberOfPoints];
                }

                double const x = mesh.points[point];
                double const y = mesh.points[point + mesh.numberOfPoints];
                double const rx = alpaka::math::max(pump.radiusX, 1.0e-300);
                double const ry = alpaka::math::max(pump.radiusY, 1.0e-300);
                double const radius = alpaka::math::sqrt((x * x) / (ry * ry) + (y * y) / (rx * rx));
                double const inlet
                    = pump.extraction ? 0.0
                                      : pump.intensity * alpaka::math::exp(-alpaka::math::pow(radius, pump.exponent));
                double const crystalStep = crystalLength / static_cast<double>(mesh.numberOfLevels - 1u);
                double const invPhotonEnergy = pump.wavelength / (planckConstant * speedOfLight);

                auto const offset = point * mesh.numberOfLevels;
                pumpForward[offset] = inlet;
                for(unsigned level = 0u; level + 1u < mesh.numberOfLevels; ++level)
                {
                    double const betaA = betaBefore[point + level * mesh.numberOfPoints];
                    double const betaB = betaBefore[point + (level + 1u) * mesh.numberOfPoints];
                    double const betaAverage = 0.5 * (betaA + betaB);
                    double const exponent
                        = -(pump.sigmaAbsorption - betaAverage * (pump.sigmaAbsorption + pump.sigmaEmission)) * nTot
                          * crystalStep;
                    pumpForward[offset + level + 1u] = pumpForward[offset + level] * alpaka::math::exp(exponent);
                }

                if(pump.backReflection)
                {
                    pumpBackward[offset + mesh.numberOfLevels - 1u]
                        = pumpForward[offset + mesh.numberOfLevels - 1u] * pump.reflectivity;
                    for(unsigned reverse = mesh.numberOfLevels - 1u; reverse > 0u; --reverse)
                    {
                        unsigned const level = reverse - 1u;
                        double const betaA = betaBefore[point + level * mesh.numberOfPoints];
                        double const betaB = betaBefore[point + (level + 1u) * mesh.numberOfPoints];
                        double const betaAverage = 0.5 * (betaA + betaB);
                        double const exponent
                            = -(pump.sigmaAbsorption - betaAverage * (pump.sigmaAbsorption + pump.sigmaEmission))
                              * nTot * crystalStep;
                        pumpBackward[offset + level] = pumpBackward[offset + level + 1u] * alpaka::math::exp(exponent);
                    }
                }
                else
                {
                    for(unsigned level = 0u; level < mesh.numberOfLevels; ++level)
                    {
                        pumpBackward[offset + level] = 0.0;
                    }
                }

                for(unsigned level = 0u; level < mesh.numberOfLevels; ++level)
                {
                    double const localPump = pumpForward[offset + level] + pumpBackward[offset + level];
                    unsigned const sample = point + level * mesh.numberOfPoints;
                    double const beta = betaBefore[sample];
                    double const pumpRate = (pump.sigmaAbsorption - beta * (pump.sigmaAbsorption + pump.sigmaEmission))
                                            * localPump * invPhotonEnergy;
                    betaAfter[sample] = beta + pump.duration * pumpRate;
                }
            }
        }
    };

    void enqueueOneDimensionalPump(
        auto& devBundle,
        hase::concepts::Queue auto const& queue,
        auto const& mesh,
        hase::core::PumpParameters const& pump,
        auto& betaBefore,
        auto& betaAfter,
        auto& pumpForward,
        auto& pumpBackward)
    {
        auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfPoints});
        queue.enqueue(
            frameSpec,
            alpaka::KernelBundle{
                OneDimensionalZTraversalPump{
                    OneDimensionalPumpKernelParameters{
                        pump.intensity,
                        pump.wavelength,
                        pump.radiusX,
                        pump.radiusY,
                        pump.exponent,
                        pump.duration,
                        pump.sigmaAbsorption,
                        pump.sigmaEmission,
                        pump.reflectivity,
                        pump.substeps,
                        pump.backReflection,
                        pump.extraction},
                    static_cast<double>(mesh.nTot),
                    static_cast<double>(mesh.thickness) * static_cast<double>(mesh.numberOfLevels - 1u)},
                mesh,
                betaBefore,
                betaAfter,
                pumpForward,
                pumpBackward});
    }

} // namespace hase::kernels
