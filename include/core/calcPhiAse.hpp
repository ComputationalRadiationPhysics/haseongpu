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
#include <core/cancellation.hpp>
#include <core/logging.hpp>
#include <kernels/calcSampleGainSum.hpp>
#include <kernels/importanceSampling.hpp>
#include <kernels/mapRaysToPrisms.hpp>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <limits>
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
    T_Elem calcMSE(T_Elem const gainSum, T_Elem const gainSumSquare, unsigned const raysPerSample)
    {
        if(raysPerSample < 2)
            return std::numeric_limits<T_Elem>::max();

        auto const n = static_cast<T_Elem>(raysPerSample);
        T_Elem const varianceOfMean = (gainSumSquare - gainSum * gainSum / n) / (n * (n - 1));

        return alpaka::math::sqrt(alpaka::math::max(T_Elem{0}, varianceOfMean)) / (T_Elem{4} * M_PI);
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
        unsigned threadLocalStridingRNG,
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

        // stack allocated addresses for host-device result transfers
        double gainSumHost = 0.0;
        double gainSumSquareHost = 0.0;
        unsigned droppedRaysHost = 0u;
        InfiniteRaySnapshot infiniteRaySnapshotHost{};
        // view to the stack memory location
        auto hgainSumView = makeView(alpaka::api::host, &gainSumHost, alpaka::Vec{1u});
        auto hgainSumSquareView = makeView(alpaka::api::host, &gainSumSquareHost, alpaka::Vec{1u});
        auto hDroppedRaysView = makeView(alpaka::api::host, &droppedRaysHost, alpaka::Vec{1u});
        auto hInfiniteRaySnapshotView = makeView(alpaka::api::host, &infiniteRaySnapshotHost, alpaka::Vec{1u});
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
        auto dDroppedRays = hase::alpakaUtils::toDevice(queue, hDroppedRaysView);
        auto dInfiniteRaySnapshots = hase::alpakaUtils::toDevice(queue, hInfiniteRaySnapshotView);

        auto dSigmaA = hase::alpakaUtils::toDevice(queue, experiment.sigmaA);
        auto dSigmaE = hase::alpakaUtils::toDevice(queue, experiment.sigmaE);

        // Memset + fill operations
        alpaka::onHost::memset(queue, dNumberOfReflectionSlices, 0);
        alpaka::onHost::memset(queue, dPrefixSum, 0);
        alpaka::onHost::memset(queue, dImportance, 0);
        alpaka::onHost::memset(queue, dPreImportance, 0);
        alpaka::onHost::memset(queue, dIndicesOfPrisms, 0);
        alpaka::onHost::memset(queue, dRaysPerPrism, 0);
        // synchronize all previous operations
        alpaka::onHost::wait(queue);
        // Calculation for each sample point
        for(unsigned sampleIdx = minSampleIdx; sampleIdx <= maxSampleIdx; ++sampleIdx)
        {
            // process potentially incoming Sigkill or Sigterm signals once per sample
            throwIfCancellationRequested();
            unsigned hRaysPerSampleDump = 0;
            raysPerSampleIter = raysPerSampleList.begin();
            bool mseTooHigh = true;
            unsigned observedDroppedRays = 0u;
            alpaka::onHost::wait(queue);
            alpaka::onHost::fill(queue, dDroppedRays, 0u, Vec1D{1});
            {
                BenchSync(queue, importanceSamplingPropagation);
                // initial ray-trace to detect importance regions
                kernels::importanceSamplingPropagation(
                    devBundle,
                    queue,
                    sampleIdx,
                    reflectionSlices,
                    mesh,
                    experiment.maxSigmaA,
                    experiment.maxSigmaE,
                    dPreImportance,
                    dDroppedRays,
                    dInfiniteRaySnapshots);
                alpaka::onHost::reduce(
                    queue,
                    devBundle.executor,
                    0.0,
                    dPreImpotanceReductionBuf,
                    std::plus{},
                    dPreImportance);
            }
            double hSumPhi = 0.0;

            alpaka::onHost::memcpy(
                queue,
                alpaka::makeView(alpaka::api::host, &hSumPhi, alpaka::Vec{1U}),
                dPreImpotanceReductionBuf);
            alpaka::onHost::memcpy(queue, hDroppedRaysView, dDroppedRays);
            alpaka::onHost::wait(queue);
            if(droppedRaysHost > 0u)
            {
                alpaka::onHost::memcpy(queue, hInfiniteRaySnapshotView, dInfiniteRaySnapshots);
                alpaka::onHost::wait(queue);
                auto const& snapshot = infiniteRaySnapshotHost;
                dout(V_WARNING) << "Non-finite importance ray: sample=" << sampleIdx << " prism=" << snapshot.prism
                                << " triangle=" << snapshot.triangle << " level=" << snapshot.level
                                << " gain=" << snapshot.gain << " start=(" << snapshot.start.x << ", "
                                << snapshot.start.y << ", " << snapshot.start.z << ") end=(" << snapshot.end.x << ", "
                                << snapshot.end.y << ", " << snapshot.end.z << ") direction=(" << snapshot.direction.x
                                << ", " << snapshot.direction.y << ", " << snapshot.direction.z
                                << ") accumulatedLength=" << snapshot.accumulatedLength
                                << " totalLength=" << snapshot.totalLength << std::endl;
            }
            while(mseTooHigh)
            {
                unsigned run = 0;

                while(run < compute.maxRepetitions && mseTooHigh)
                {
                    alpaka::onHost::wait(queue);
                    {
                        BenchSync(queue, importanceSamplingDistribution);
                        // map importance to raysPerPrism
                        hRaysPerSampleDump = kernels::importanceSamplingDistribution(
                            devBundle,
                            queue,
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
                        BenchSync(queue, mapRaysToPrisms);
                        // Prism scheduling for gpu threads
                        kernels::mapRaysToPrisms(
                            devBundle,
                            queue,
                            dIndicesOfPrisms,
                            dNumberOfReflectionSlices,
                            dRaysPerPrism,
                            dPrefixSum,
                            mesh.numberOfPrisms);
                    }

                    // Start Kernel
                    alpaka::onHost::fill(queue, dGainSum, double{0}, Vec1D{1});
                    alpaka::onHost::fill(queue, dGainSumSquare, double{0}, Vec1D{1});
                    alpaka::onHost::fill(queue, dDroppedRays, 0u, Vec1D{1});
                    alpaka::onHost::fill(queue, dGainOfRay, double{0}, alpaka::Vec{experiment.maxRaysPerSample});

                    auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                        devBundle.device,
                        devBundle.executor,
                        alpaka::Vec{static_cast<unsigned int>(*raysPerSampleIter)});
                    auto const threadLocalStridingIndex = threadLocalStridingRNG;
                    if(experiment.useReflections)
                    {
                        BenchSync(queue, CalcSampleGainSumWithReflection);
                        // main ray propagation routine with reflection
                        queue.enqueue(
                            frameSpec,
                            alpaka::KernelBundle{
                                hase::kernels::CalcSampleGainSumWithReflection{},
                                mesh,
                                dIndicesOfPrisms,
                                dNumberOfReflectionSlices,
                                dImportance,
                                hRaysPerSampleDump,
                                dGainOfRay,
                                dDroppedRays,
                                sampleIdx,
                                dSigmaA,
                                dSigmaE,
                                threadLocalStridingIndex});
                    }
                    else
                    {
                        BenchSync(queue, CalcSampleGainSum);
                        // main ray propagation routine
                        queue.enqueue(
                            frameSpec,
                            alpaka::KernelBundle{
                                hase::kernels::CalcSampleGainSum{},
                                mesh,
                                dIndicesOfPrisms,
                                dImportance,
                                hRaysPerSampleDump,
                                dGainOfRay,
                                dDroppedRays,
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
                    alpaka::onHost::memcpy(queue, hgainSumView, dGainSum);
                    alpaka::onHost::memcpy(queue, hgainSumSquareView, dGainSumSquare);
                    alpaka::onHost::memcpy(queue, hDroppedRaysView, dDroppedRays);
                    alpaka::onHost::wait(queue);
                    observedDroppedRays = std::max(observedDroppedRays, droppedRaysHost);
                    if(droppedRaysHost > 0u)
                    {
                        if(result.totalRays.at(sampleIdx) == 0u)
                        {
                            result.droppedRays.at(sampleIdx) = droppedRaysHost;
                        }
                        continue;
                    }
                    ++run;
                    double mseTmp = calcMSE(gainSumHost, gainSumSquareHost, hRaysPerSampleDump);

                    assert(!alpaka::math::isnan(hgainSumView[0]));
                    assert(!alpaka::math::isnan(hgainSumSquareView[0]));
                    assert(!alpaka::math::isnan(mseTmp));
                    auto updateSample = [&]()
                    {
                        result.mse.at(sampleIdx) = mseTmp;
                        result.phiAse.at(sampleIdx) = gainSumHost;
                        result.phiAse.at(sampleIdx) /= *raysPerSampleIter * 4.0f * M_PI;
                        result.totalRays.at(sampleIdx) = *raysPerSampleIter;
                        result.droppedRays.at(sampleIdx) = droppedRaysHost;
                    };

                    if(droppedRaysHost == 0u && result.mse.at(sampleIdx) > mseTmp)
                    {
                        updateSample();
                    }
                    else if(droppedRaysHost > 0u && result.totalRays.at(sampleIdx) == 0u)
                    {
                        result.droppedRays.at(sampleIdx) = droppedRaysHost;
                    }

                    if(result.mse.at(sampleIdx) < experiment.mseThreshold && result.droppedRays.at(sampleIdx) == 0u)
                        mseTooHigh = false;
                }

                // Increase rays per sample or break, when mseThreshold was not met
                raysPerSampleIter++;
                if(raysPerSampleIter == raysPerSampleList.end())
                {
                    if(mseTooHigh)
                    {
                        dout(V_WARNING)
                            << "For sample: " << sampleIdx
                            << " the requested mse threshold: " << experiment.mseThreshold
                            << " and zero dropped rays could not be reached given the maximum number of rays: "
                            << experiment.maxRaysPerSample
                            << " and the number of repetitions: " << compute.maxRepetitions
                            << ". Max dropped-ray count observed in an attempt: " << observedDroppedRays;
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
