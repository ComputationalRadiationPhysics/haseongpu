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
 * @licence GPLv3
 *
 */

#pragma once

#include <core/calcForwardPhiAse.hpp>
#include <core/forwardSrmDeviceState.hpp>
#include <core/logging.hpp>
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <random/random.hpp>
#include <utils/progressbar.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <memory>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>

/**
 * @brief Wrapper for forward PhiASE accumulation on pthread base.
 *        This function will spawn a thread for
 *        each function call and start calcForwardPhiAseRaw.
 *
 * @param dMesh            Explicit 3D cell mesh in device memory.
 * @param hMesh            Same as dMesh, but located in host memory.
 * @param sigmaA           Vector with Absorption values
 * @param sigmaE           Vector with Emission values
 * @param result           Reference to raw forward accumulators.
 * @param gpu_i            Number of device that should be used.
 * @param rayCount         Number of forward ray histories assigned to this thread.
 * @param runtime          Reference to the needed runtime.
 */
namespace hase::core
{

    static std::vector<std::thread> threadIds;
    static std::vector<std::exception_ptr> threadExceptions;
    static std::mutex threadExceptionsMutex;

    void joinAll();

    template<alpaka::onHost::concepts::Device T_Device>
    void calcForwardPhiAseThreaded(
        auto& devBundle,
        ExperimentParameters const& experiment,
        HostMesh const& hostMesh,
        DeviceMeshContainer<T_Device> const& mesh,
        ForwardPhiAseRawResult& result,
        unsigned const rayCount,
        unsigned rngSeed,
        float& runtime)
    {
        threadIds.emplace_back(
            std::thread(
                [&experiment, &hostMesh, &mesh, &result, &runtime, devBundle, rayCount, rngSeed]() mutable
                {
                    try
                    {
                        calcForwardPhiAseRaw(
                            devBundle,
                            experiment,
                            hostMesh,
                            mesh,
                            result,
                            runtime,
                            rayCount,
                            rngSeed);
                    }
                    catch(...)
                    {
                        std::lock_guard<std::mutex> lock(threadExceptionsMutex);
                        threadExceptions.emplace_back(std::current_exception());
                    }
                }));
    }

    template<typename T_Exec, alpaka::onHost::concepts::Device T_Device>
    ForwardPhiAseRawResult calcForwardPhiAseSrmOnDevices(
        T_Exec const& exec,
        ExperimentParameters const& experiment,
        HostMesh const& hostMesh,
        std::vector<DeviceMeshContainer<T_Device>> const& meshes,
        unsigned const firstDevice,
        unsigned const assignedDeviceCount,
        unsigned const rayCount,
        unsigned const baseSeed,
        unsigned const rank,
        std::vector<float>& runtimes)
    {
        if(experiment.surfaceReservoirSize == 0u)
        {
            throw std::runtime_error("Forward reflections require surfaceReservoirSize > 0.");
        }

        ForwardPhiAseRawResult combined = makeForwardRawResult(hostMesh.numberOfCells);
        SrmControls const controls = resolveSrmControls(experiment);
        bool const debugSrm = srmDebugLoggingEnabled();
        combined.srmMaxIterations = controls.maxIterations;
        combined.srmDivergenceStreak = controls.divergenceStreak;
        if(rayCount == 0u)
        {
            combined.srmStatus = SrmStatus::CONVERGED;
            return combined;
        }

        unsigned const activeDevices = std::min(assignedDeviceCount, rayCount);
        if(activeDevices == 0u)
        {
            return combined;
        }
        unsigned const raysPerDevice = rayCount / activeDevices;
        unsigned const remainder = rayCount % activeDevices;
        double const betaVolumeTotal = calcForwardBetaVolumeTotal(hostMesh);
        std::vector<std::unique_ptr<ForwardSrmDeviceState<T_Device, T_Exec>>> devices;
        devices.reserve(activeDevices);
        std::vector<unsigned> deviceIndices;
        deviceIndices.reserve(activeDevices);
        for(unsigned localDeviceIndex = 0u; localDeviceIndex < activeDevices; ++localDeviceIndex)
        {
            unsigned const deviceIndex = firstDevice + localDeviceIndex;
            unsigned const localRayCount
                = localDeviceIndex + 1u == activeDevices ? raysPerDevice + remainder : raysPerDevice;
            devices.emplace_back(
                std::make_unique<ForwardSrmDeviceState<T_Device, T_Exec>>(
                    meshes.at(deviceIndex).m_device,
                    exec,
                    meshes.at(deviceIndex),
                    experiment,
                    betaVolumeTotal,
                    localRayCount,
                    hase::random::seedForWorker(baseSeed, rank, deviceIndex)));
            deviceIndices.emplace_back(deviceIndex);
        }

        auto const started = std::chrono::steady_clock::now();
        auto const directStarted = started;
        if(debugSrm)
        {
            dout(V_INFO) << "SRM: direct pass: rays=" << rayCount << ", devices=" << activeDevices
                         << ", faces=" << hostMesh.numberOfCells * 4u << std::endl;
        }
        for(auto& device : devices)
        {
            device->traceDirect();
        }

        std::vector<double> sourceWeights(activeDevices, 0.0);
        for(unsigned deviceIndex = 0u; deviceIndex < activeDevices; ++deviceIndex)
        {
            sourceWeights[deviceIndex] = devices[deviceIndex]->prepareDirectSamplingCdf();
        }
        double const initialWeight = std::accumulate(sourceWeights.cbegin(), sourceWeights.cend(), 0.0);
        if(debugSrm)
        {
            double const directSeconds
                = std::chrono::duration<double>(std::chrono::steady_clock::now() - directStarted).count();
            dout(V_INFO) << "SRM: direct pass completed in " << directSeconds
                         << " s; reflected source weight=" << initialWeight << std::endl;
        }
        if(initialWeight == 0.0)
        {
            combined.srmStatus = SrmStatus::CONVERGED;
        }
        else
        {
            double previousWeight = initialWeight;
            unsigned growCount = 0u;
            combined.srmStatus = SrmStatus::MAX_ITERATIONS;
            combined.srmRemainingFraction = 1.0;
            for(unsigned pass = 1u; pass <= controls.maxIterations; ++pass)
            {
                auto const passStarted = std::chrono::steady_clock::now();
                if(debugSrm)
                {
                    dout(V_INFO) << "SRM: reflected pass " << pass << "/" << controls.maxIterations << ": launching "
                                 << rayCount << " rays from weight=" << previousWeight << " across " << activeDevices
                                 << " in-memory reservoir strata" << std::endl;
                }
                for(unsigned deviceIndex = 0u; deviceIndex < activeDevices; ++deviceIndex)
                {
                    double const sourceWeight
                        = sourceWeights[deviceIndex] / static_cast<double>(devices[deviceIndex]->rayCount());
                    devices[deviceIndex]->traceReflected(sourceWeight, pass);
                }

                for(unsigned deviceIndex = 0u; deviceIndex < activeDevices; ++deviceIndex)
                {
                    sourceWeights[deviceIndex] = devices[deviceIndex]->prepareReflectedSamplingCdf(pass);
                }
                double const currentWeight = std::accumulate(sourceWeights.cbegin(), sourceWeights.cend(), 0.0);
                combined.srmPasses = pass;
                combined.srmRemainingFraction = currentWeight / initialWeight;
                if(debugSrm)
                {
                    double const passSeconds
                        = std::chrono::duration<double>(std::chrono::steady_clock::now() - passStarted).count();
                    dout(V_INFO) << "SRM: reflected pass " << pass << " completed in " << passSeconds
                                 << " s; next weight=" << currentWeight << ", W/W0=" << combined.srmRemainingFraction
                                 << std::endl;
                }
                if(currentWeight > previousWeight)
                {
                    ++growCount;
                    if(growCount >= controls.divergenceStreak)
                    {
                        combined.srmStatus = SrmStatus::DIVERGED;
                        break;
                    }
                }
                else
                {
                    growCount = 0u;
                    if(std::abs(currentWeight - previousWeight) / std::max(currentWeight, 1.0e-30)
                       < experiment.reflectionTolerance)
                    {
                        combined.srmStatus = SrmStatus::STABLE;
                        break;
                    }
                }
                if(combined.srmRemainingFraction < experiment.reflectionTolerance)
                {
                    combined.srmStatus = SrmStatus::CONVERGED;
                    break;
                }
                previousWeight = currentWeight;
            }
        }

        for(auto const& device : devices)
        {
            ForwardPhiAseRawResult partial = makeForwardRawResult(hostMesh.numberOfCells);
            device->downloadAccumulation(partial);
            mergeForwardRawResult(combined, partial);
        }
        if(combined.rayCount != rayCount)
        {
            throw std::runtime_error("Forward ray partition accounting mismatch.");
        }
        float const runtime
            = static_cast<float>(std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count());
        for(unsigned const deviceIndex : deviceIndices)
        {
            runtimes.at(deviceIndex) = runtime;
        }
        return combined;
    }

    template<typename T_Exec, alpaka::onHost::concepts::Device T_Device>
    ForwardPhiAseRawResult calcForwardPhiAseOnDevices(
        T_Exec exec,
        ExperimentParameters const& experiment,
        HostMesh const& hostMesh,
        std::vector<DeviceMeshContainer<T_Device>> const& meshes,
        unsigned const firstDevice,
        unsigned const assignedDeviceCount,
        unsigned const rayCount,
        unsigned const baseSeed,
        unsigned const rank,
        std::vector<float>& runtimes)
    {
        unsigned const volumeCount = hostMesh.numberOfCells;
        // Each partial contributes the histories it actually launched; starting
        // this counter at the requested total would count every history twice.
        ForwardPhiAseRawResult combined = makeForwardRawResult(volumeCount);
        if(rayCount == 0u || assignedDeviceCount == 0u)
        {
            return combined;
        }

        if(experiment.useReflections)
        {
            return calcForwardPhiAseSrmOnDevices(
                exec,
                experiment,
                hostMesh,
                meshes,
                firstDevice,
                assignedDeviceCount,
                rayCount,
                baseSeed,
                rank,
                runtimes);
        }

        unsigned const activeDevices = std::min(assignedDeviceCount, rayCount);
        unsigned const raysPerDevice = rayCount / activeDevices;
        unsigned const remainder = rayCount % activeDevices;
        std::vector<ForwardPhiAseRawResult> partials(activeDevices, makeForwardRawResult(volumeCount));

        for(unsigned localDeviceIndex = 0u; localDeviceIndex < activeDevices; ++localDeviceIndex)
        {
            unsigned const deviceIndex = firstDevice + localDeviceIndex;
            unsigned const localRayCount
                = localDeviceIndex + 1u == activeDevices ? raysPerDevice + remainder : raysPerDevice;
            hase::alpakaUtils::DevBundle devBundle{meshes.at(deviceIndex).m_device, exec};
            unsigned const rngSeed = hase::random::seedForWorker(baseSeed, rank, deviceIndex);
            calcForwardPhiAseThreaded(
                devBundle,
                experiment,
                hostMesh,
                meshes.at(deviceIndex),
                partials.at(localDeviceIndex),
                localRayCount,
                rngSeed,
                runtimes.at(deviceIndex));
        }

        joinAll();
        for(auto const& partial : partials)
        {
            mergeForwardRawResult(combined, partial);
        }
        if(combined.rayCount != rayCount)
        {
            throw std::runtime_error("Forward ray partition accounting mismatch.");
        }
        return combined;
    }

    /**
     * @brief Wait for all threads to finish
     *
     */
    inline void joinAll()
    {
        for(auto& t : threadIds)
        {
            if(t.joinable())
            {
                t.join();
            }
        }
        threadIds.clear();

        if(!threadExceptions.empty())
        {
            auto exception = threadExceptions.front();
            threadExceptions.clear();
            std::rethrow_exception(exception);
        }
    }

} // namespace hase::core
