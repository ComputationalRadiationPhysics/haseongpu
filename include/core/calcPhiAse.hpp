/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once

#include <core/mesh.hpp>
#include <core/types.hpp>

/**
 * @brief Calculates Phi ASE. With minRaysPerSample < maxRaysPerSample
 *        adaptive sampling can be used to improve performance.
 *
 * @param minSample_i      Smallest Index of sample point to calculate.
 * @param maxSample_i      Biggest Index of sample point to calculate.
 * @param runtime          Reference to the needed runtime.
 *
 **/

#include <alpakaUtils/DevBundle.hpp>
#include <benchmark.hpp>
#include <core/logging.hpp>
#include <kernels/calcSampleGainSum.hpp>
#include <kernels/importanceSampling.hpp>
#include <kernels/mapRaysToPrisms.hpp>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <vector>
#define SEED 4321

namespace hase::core
{
    using alpakaUtils::Vec1D;

    struct Square
    {
        constexpr auto operator()(alpaka::concepts::Simd auto&& a) const
        {
            return a * a;
        }
    };

    template<std::floating_point T_Elem>
    T_Elem calcMSE(T_Elem const phiAse, T_Elem const phiAseSquare, unsigned const raysPerSample)
    {
        T_Elem a = phiAseSquare / raysPerSample;
        T_Elem b = (phiAse / raysPerSample) * (phiAse / raysPerSample);

        return alpaka::math::sqrt(alpaka::math::abs((a - b) / raysPerSample));
    }

    inline std::vector<int> generateRaysPerSampleExpList(int minRaysPerSample, int maxRaysPerSample, int steps)
    {
        std::vector<int> raysPerSample;

        if((minRaysPerSample == maxRaysPerSample) || steps < 2)
        {
            raysPerSample.push_back(minRaysPerSample);
            return raysPerSample;
        }

        for(int i = 0; i < steps; ++i)
        {
            int step_val = minRaysPerSample
                           * alpaka::math::pow((maxRaysPerSample / minRaysPerSample), (i / (float) (steps - 1)));
            raysPerSample.push_back(step_val);
        }

        return raysPerSample;
    }

    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec, typename... Args>
    constexpr float calcPhiAse(
        hase::alpakaUtils::DevBundle<T_Device, T_Exec>& devBundle,
        ExperimentParameters const& experiment,
        ComputeParameters const& compute,
        HostMesh const& hostMesh,
        DeviceMeshContainer<T_Device> const& meshContainer,
        Result& result,
        unsigned const minSampleIdx,
        unsigned const maxSampleIdx,
        float& runtime,
        Args&&... args)
    {
#ifdef HASE_ENABLE_BENCHMARK
        hase::benchmark::ScopedRunContext benchmarkContext{devBundle.device, devBundle.executor, compute, experiment};
        hase::benchmark::ScopedEvent benchmarkSimulation{
            "phiASE Simulation",
            "minSampleIdx=" + std::to_string(minSampleIdx) + ";maxSampleIdx=" + std::to_string(maxSampleIdx)};
#endif
        DeviceMeshView mesh = meshContainer.toView();
        time_t starttime = time(0);
        unsigned maxReflections = experiment.useReflections ? hostMesh.getMaxReflections() : 0;
        unsigned reflectionSlices = 1 + (2 * maxReflections);
        // Divide RaysPerSample range into steps
        std::vector<int> raysPerSampleList = generateRaysPerSampleExpList(
            experiment.minRaysPerSample,
            experiment.maxRaysPerSample,
            compute.adaptiveSteps);

        std::vector<int>::iterator raysPerSampleIter = raysPerSampleList.begin();

        // initialize our alpaka device queue
        auto queue = devBundle.device.makeQueue();
        unsigned threadLocalStridingRNG = 0;

        // stack allocated addresses for host-device result transfers
        double gainSumHost = 0.0;
        double gainSumSquareHost = 0.0;
        // view to the stack memory location
        auto hgainSumView = makeView(alpaka::api::host, &gainSumHost, alpaka::Vec{1u});
        auto hgainSumSquareView = makeView(alpaka::api::host, &gainSumSquareHost, alpaka::Vec{1u});
        //
        // Memory allocations
        auto dNumberOfReflectionSlices
            = alpaka::onHost::alloc<unsigned>(devBundle.device, experiment.maxRaysPerSample);

        auto dRaysPerPrism = alpaka::onHost::alloc<unsigned>(devBundle.device, mesh.numberOfPrisms * reflectionSlices);

        auto dPrefixSum = alpaka::onHost::alloc<unsigned>(devBundle.device, mesh.numberOfPrisms * reflectionSlices);

        auto dImportance = alpaka::onHost::alloc<double>(devBundle.device, mesh.numberOfPrisms * reflectionSlices);

        auto dPreImportance = alpaka::onHost::alloc<double>(devBundle.device, mesh.numberOfPrisms * reflectionSlices);

        auto dPreImpotanceReductionBuf = alpaka::onHost::alloc<double>(devBundle.device, 1);

        auto dGainOfRay = alpaka::onHost::alloc<double>(devBundle.device, experiment.maxRaysPerSample);

        auto dIndicesOfPrisms = alpaka::onHost::alloc<unsigned>(devBundle.device, experiment.maxRaysPerSample);

        // Memory allocations + copy
        auto dGainSum = hase::alpakaUtils::toDevice(queue, hgainSumView);
        auto dGainSumSquare = hase::alpakaUtils::toDevice(queue, hgainSumSquareView);

        auto dSigmaA = hase::alpakaUtils::toDevice(queue, experiment.sigmaA);
        auto dSigmaE = hase::alpakaUtils::toDevice(queue, experiment.sigmaE);

        // Memset + fill operations
        alpaka::onHost::memset(queue, dNumberOfReflectionSlices, 0);
        alpaka::onHost::memset(queue, dPrefixSum, 0);
        alpaka::onHost::memset(queue, dImportance, 0);
        alpaka::onHost::memset(queue, dPreImportance, 0);
        alpaka::onHost::memset(queue, dIndicesOfPrisms, 0);
        alpaka::onHost::fill(queue, dRaysPerPrism, 0u);
        // synchronize all previous operations
        alpaka::onHost::wait(queue);
        // Calculation for each sample point
        for(unsigned sampleIdx = minSampleIdx; sampleIdx <= maxSampleIdx; ++sampleIdx)
        {
            unsigned hRaysPerSampleDump = 0;
            raysPerSampleIter = raysPerSampleList.begin();
            bool mseTooHigh = true;
            alpaka::onHost::wait(queue);
            {
                BenchSync(queue, importanceSamplingPropagation)
                    // initial ray-trace to detect importance regions
                    hase::kernels::importanceSamplingPropagation(
                        devBundle,
                        sampleIdx,
                        reflectionSlices,
                        mesh,
                        experiment.maxSigmaA,
                        experiment.maxSigmaE,
                        dPreImportance);
                alpaka::onHost::reduce(
                    queue,
                    devBundle.executor,
                    0.0,
                    dPreImpotanceReductionBuf,
                    std::plus{},
                    dPreImportance);
            }
            alpaka::onHost::wait(queue);
            double hSumPhi = 0.0;

            alpaka::onHost::memcpy(
                queue,
                alpaka::makeView(alpaka::api::host, &hSumPhi, alpaka::Vec{1U}),
                dPreImpotanceReductionBuf);
            alpaka::onHost::wait(queue);
            while(mseTooHigh)
            {
                unsigned run = 0;

                while(run++ < compute.maxRepetitions && mseTooHigh)
                {
                    alpaka::onHost::wait(queue);
                    {
                        BenchSync(queue, importanceSamplingDistribution)
                            // map importance to raysPerPrism
                            hRaysPerSampleDump
                            = hase::kernels::importanceSamplingDistribution(
                                devBundle,
                                reflectionSlices,
                                mesh,
                                *raysPerSampleIter,
                                dPreImportance,
                                dImportance,
                                dRaysPerPrism,
                                hSumPhi,
                                threadLocalStridingRNG);
                    }
                    {
                        BenchSync(queue, mapRaysToPrisms)
                            // Prism scheduling for gpu threads
                            hase::kernels::mapRaysToPrisms(
                                devBundle,
                                dIndicesOfPrisms,
                                dNumberOfReflectionSlices,
                                dRaysPerPrism,
                                dPrefixSum,
                                reflectionSlices,
                                *raysPerSampleIter,
                                mesh.numberOfPrisms);
                    }

                    // Start Kernel
                    alpaka::onHost::fill(queue, dGainSum, double{0}, Vec1D{1});
                    alpaka::onHost::fill(queue, dGainSumSquare, double{0}, Vec1D{1});
                    alpaka::onHost::fill(queue, dGainOfRay, double{0}, alpaka::Vec{experiment.maxRaysPerSample});

                    auto frameSpec
                        = alpaka::onHost::getFrameSpec<unsigned>(devBundle.device, alpaka::Vec{*raysPerSampleIter});
                    auto const threadLocalStridingIndex = threadLocalStridingRNG;
                    if(experiment.useReflections)
                    {
                        BenchSync(queue, CalcSampleGainSumWithReflection)
                            // main ray propagation routine with reflection
                            queue.enqueue(
                                devBundle.executor,
                                frameSpec,
                                alpaka::KernelBundle{
                                    hase::kernels::CalcSampleGainSumWithReflection{},
                                    mesh,
                                    dIndicesOfPrisms,
                                    dNumberOfReflectionSlices,
                                    dImportance,
                                    hRaysPerSampleDump,
                                    dGainOfRay,
                                    sampleIdx,
                                    dSigmaA,
                                    dSigmaE,
                                    threadLocalStridingIndex});
                    }
                    else
                    {
                        BenchSync(queue, CalcSampleGainSum)
                            // main ray propagation routine
                            queue.enqueue(
                                devBundle.executor,
                                frameSpec,
                                alpaka::KernelBundle{
                                    hase::kernels::CalcSampleGainSum{},
                                    mesh,
                                    dIndicesOfPrisms,
                                    dImportance,
                                    hRaysPerSampleDump,
                                    dGainOfRay,
                                    sampleIdx,
                                    dSigmaA,
                                    dSigmaE,
                                    experiment.sigmaA.size(),
                                    threadLocalStridingIndex});
                    }
                    alpaka::onHost::reduce(queue, devBundle.executor, 0.0, dGainSum, std::plus{}, dGainOfRay);
                    alpaka::onHost::transformReduce(
                        queue,
                        devBundle.executor,
                        double{0},
                        dGainSumSquare,
                        std::plus{},
                        Square{},
                        dGainOfRay);
                    // add stride since rng seed has been used in CalcSampleGainSum
                    threadLocalStridingRNG
                        += (frameSpec.getNumFrames().product() * frameSpec.getFrameExtents().product());
                    alpaka::onHost::wait(queue);
                    alpaka::onHost::memcpy(queue, hgainSumView, dGainSum);
                    alpaka::onHost::memcpy(queue, hgainSumSquareView, dGainSumSquare);
                    alpaka::onHost::wait(queue);
                    float mseTmp = calcMSE(gainSumHost, gainSumSquareHost, hRaysPerSampleDump);

                    assert(!alpaka::math::isnan(hgainSumView[0]));
                    assert(!alpaka::math::isnan(hgainSumSquareView[0]));
                    assert(!alpaka::math::isnan(mseTmp));

                    if(result.mse.at(sampleIdx) > mseTmp)
                    {
                        result.mse.at(sampleIdx) = mseTmp;
                        result.phiAse.at(sampleIdx) = gainSumHost;
                        result.phiAse.at(sampleIdx) /= *raysPerSampleIter * 4.0f * M_PI;
                        result.totalRays.at(sampleIdx) = *raysPerSampleIter;
                    }
                    if(result.mse.at(sampleIdx) < experiment.mseThreshold)
                        mseTooHigh = false;
                }

                // Increase rays per sample or break, when mseThreshold was not met
                raysPerSampleIter++;
                if(raysPerSampleIter == raysPerSampleList.end())
                {
                    if(mseTooHigh)
                    {
                        dout(V_WARNING) << "For sample: " << sampleIdx
                                        << " the requested mse threshold: " << experiment.mseThreshold
                                        << " could not be reached given the maximum number of rays: "
                                        << experiment.maxRaysPerSample
                                        << " and the number of repetitions: " << compute.maxRepetitions;
                    }
                    break;
                }
            }
            if constexpr(sizeof...(Args) > 0)
            {
                auto& progressBar = std::get<0>(std::forward_as_tuple(args...));
                if(verbosity & V_PROGRESS)
                {
                    progressBar.printFancyProgressBar(mesh.numberOfSamples);
                }
            }
        }

        runtime = difftime(time(0), starttime);
        alpaka::onHost::wait(queue);
        return runtime;
    }

} // namespace hase::core
