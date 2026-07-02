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
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <random/random.hpp>
#include <utils/progressbar.hpp>

#include <algorithm>
#include <exception>
#include <mutex>
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
        ForwardPhiAseRawResult combined = makeForwardRawResult(volumeCount);
        combined.rayCount = rayCount;
        if(rayCount == 0u || assignedDeviceCount == 0u)
        {
            return combined;
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
