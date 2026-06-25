/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpakaUtils/DevBundle.hpp>
#include <core/logging.hpp>
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <kernels/forward/accumulation.hpp>

#include <algorithm>
#include <ctime>
#include <stdexcept>
#include <vector>

namespace hase::core
{
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    float calcForwardPhiAse(
        hase::alpakaUtils::DevBundle<T_Device, T_Exec>& devBundle,
        ExperimentParameters const& experiment,
        HostMesh const& hostMesh,
        DeviceMeshContainer<T_Device> const& meshContainer,
        Result& result,
        float& runtime,
        unsigned threadLocalStridingRNG)
    {
        if(experiment.useReflections)
        {
            throw std::runtime_error("Forward volume propagation does not support useReflections=True.");
        }
        if(experiment.forwardRayLength <= 0.0)
        {
            throw std::runtime_error("Forward volume propagation requires forwardRayLength > 0.");
        }

        time_t starttime = time(0);
        auto queue = devBundle.device.makeQueue();
        DeviceMeshView mesh = meshContainer.toView();
        unsigned const volumeCount = mesh.numberOfCells;
        unsigned const rayCount = experiment.resolvedForwardRayCount();
        double const totalVolume = std::accumulate(hostMesh.cellVolumes.begin(), hostMesh.cellVolumes.end(), 0.0);

        auto dPhiAccumulator = alpaka::onHost::alloc<double>(devBundle.device, volumeCount);
        auto dVolumeRayVisits = alpaka::onHost::alloc<unsigned>(devBundle.device, volumeCount);
        auto dDroppedRays = alpaka::onHost::alloc<unsigned>(devBundle.device, volumeCount);
        auto dSigmaA = hase::alpakaUtils::toDevice(queue, experiment.sigmaA);
        auto dSigmaE = hase::alpakaUtils::toDevice(queue, experiment.sigmaE);

        alpaka::onHost::fill(queue, dPhiAccumulator, double{0}, alpaka::Vec{volumeCount});
        alpaka::onHost::fill(queue, dVolumeRayVisits, 0u, alpaka::Vec{volumeCount});
        alpaka::onHost::fill(queue, dDroppedRays, 0u, alpaka::Vec{volumeCount});
        alpaka::onHost::wait(queue);

        auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{static_cast<unsigned int>(rayCount)});
        queue.enqueue(
            frameSpec,
            alpaka::KernelBundle{
                hase::kernels::forward::AccumulateForwardPhiAse{},
                mesh,
                rayCount,
                experiment.forwardRayLength,
                totalVolume,
                dPhiAccumulator,
                dVolumeRayVisits,
                dDroppedRays,
                dSigmaA,
                dSigmaE,
                static_cast<unsigned>(experiment.sigmaA.size()),
                threadLocalStridingRNG});

        auto hPhiAccumulator = alpaka::onHost::allocHostLike(dPhiAccumulator);
        auto hVolumeRayVisits = alpaka::onHost::allocHostLike(dVolumeRayVisits);
        auto hDroppedRays = alpaka::onHost::allocHostLike(dDroppedRays);
        alpaka::onHost::memcpy(queue, hPhiAccumulator, dPhiAccumulator);
        alpaka::onHost::memcpy(queue, hVolumeRayVisits, dVolumeRayVisits);
        alpaka::onHost::memcpy(queue, hDroppedRays, dDroppedRays);
        alpaka::onHost::wait(queue);

        result = Result(
            std::vector<float>(volumeCount, 0.0f),
            std::vector<double>(volumeCount, 0.0),
            std::vector<unsigned>(volumeCount, 0u),
            std::vector<double>(volumeCount, 0.0),
            std::vector<unsigned>(volumeCount, 0u));
        auto* phiData = alpaka::onHost::data(hPhiAccumulator);
        auto* visitsData = alpaka::onHost::data(hVolumeRayVisits);
        auto* droppedData = alpaka::onHost::data(hDroppedRays);
        for(unsigned volume = 0u; volume < volumeCount; ++volume)
        {
            unsigned const visits = visitsData[volume];
            result.totalRays.at(volume) = visits;
            result.droppedRays.at(volume) = droppedData[volume];
            if(visits > 0u)
            {
                result.phiAse.at(volume) = static_cast<float>(phiData[volume] / static_cast<double>(visits));
            }
            else
            {
                result.mse.at(volume) = experiment.mseThreshold;
            }
        }

        runtime = difftime(time(0), starttime);
        return runtime;
    }
} // namespace hase::core
