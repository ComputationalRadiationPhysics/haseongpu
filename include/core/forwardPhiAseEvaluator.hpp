/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <alpakaUtils/memory.hpp>
#include <core/calcForwardPhiAse.hpp>
#include <core/calcPhiAseThreaded.hpp>
#include <core/forwardPhiAseUtilities.hpp>
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <random/random.hpp>

#if !defined(DISABLE_MPI) && defined(MPI_FOUND)
#    include <core/calcPhiAseMpi.hpp>
#endif

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hase::core
{
    namespace detail
    {
        template<typename T_Buffer>
        std::vector<typename T_Buffer::value_type> copyToVector(auto const& queue, T_Buffer const& buffer)
        {
            auto hostBuffer = alpaka::onHost::allocHostLike(buffer);
            alpaka::onHost::memcpy(queue, hostBuffer, buffer);
            alpaka::onHost::wait(queue);
            auto const* data = alpaka::onHost::data(hostBuffer);
            return {data, data + buffer.getExtents().product()};
        }

        template<typename T_Buffer, typename T>
        void copyVectorToBuffer(auto const& queue, std::vector<T> const& values, T_Buffer& buffer)
        {
            auto hostView = alpaka::makeView(alpaka::api::host, values.data(), alpaka::Vec{values.size()});
            alpaka::onHost::memcpy(queue, buffer, hostView);
        }
    } // namespace detail

    struct ForwardPhiAseEvaluation
    {
        bool deviceResidentPhi = false;
        float runtime = 0.0f;
        unsigned usedDevices = 0u;
        unsigned rayCount = 0u;
        unsigned adaptiveLaunches = 0u;
        RuntimeTopology topology;
        std::vector<unsigned> convergenceRayCounts;
    };

    template<typename T_Device, typename T_Executor>
    class ForwardPhiAseContext
    {
    public:
        ForwardPhiAseContext(
            std::vector<T_Device> devices,
            T_Executor executor,
            ExperimentParameters const& experiment,
            HostMesh& hostMesh)
            : m_executor(std::move(executor))
        {
            if(devices.empty())
                throw std::runtime_error("forward ASE context requires at least one device");
            m_meshes.reserve(devices.size());
            for(auto& device : devices)
                m_meshes.emplace_back(hostMesh.toDevice(device));
            m_deviceContexts.reserve(m_meshes.size());
            for(auto const& mesh : m_meshes)
            {
                m_deviceContexts.emplace_back(
                    std::make_unique<ForwardPhiAseDeviceContext<T_Device, T_Executor>>(
                        mesh.m_device,
                        m_executor,
                        experiment,
                        hostMesh.numberOfCells));
            }
        }

        [[nodiscard]] T_Device& primaryDevice()
        {
            return m_meshes.front().m_device;
        }

        [[nodiscard]] DeviceMeshContainer<T_Device>& primaryMesh()
        {
            return m_meshes.front();
        }

        [[nodiscard]] auto& primaryBetaVolume()
        {
            return m_meshes.front().betaVolume;
        }

        [[nodiscard]] bool requiresHostBetaVolume() const
        {
            return m_meshes.size() > 1u;
        }

        std::vector<double> downloadPrimaryVolumeDndtAse()
        {
            return m_deviceContexts.front()->downloadVolumeDndtAse();
        }

        Result downloadPrimaryResult(
            bool includePhiAse,
            bool includeStandardError,
            bool includeRelativeStandardError,
            bool includeTotalRays)
        {
            return m_deviceContexts.front()->downloadFinalizedResult(
                includePhiAse,
                includeStandardError,
                includeRelativeStandardError,
                includeTotalRays);
        }

        [[nodiscard]] auto& primaryVolumeDndtAse()
        {
            return m_deviceContexts.front()->volumeDndtAse();
        }

        [[nodiscard]] auto& primaryVolumePhiAse()
        {
            return m_deviceContexts.front()->volumePhiAse();
        }

        ForwardPhiAseEvaluation evaluate(
            ExperimentParameters& experiment,
            ComputeParameters& compute,
            HostMesh& hostMesh,
            auto const& betaVolume,
            Result& result,
            bool allowDeviceResident = true)
        {
            bool const mpiHostCombined = compute.parallelMode == ParallelMode::MPI;
            refreshDynamicMeshes(betaVolume, hostMesh, requiresHostBetaVolume() || mpiHostCombined, mpiHostCombined);
            if(!experiment.isForwardPropagation())
                throw std::runtime_error("Only forward volume propagation is supported by the openPMD backend.");

            float runtime = 0.0f;
            unsigned const maxDevices = static_cast<unsigned>(m_meshes.size());
            std::vector<float> runtimes(maxDevices, 0.0f);
            unsigned usedDevices = 0u;
            RuntimeTopology topology;
            ForwardPhiAseRawResult rawResult;
            unsigned adaptiveLaunches = 0u;
            std::vector<unsigned> convergenceRayCounts;
            bool hostResultAvailable = false;

            if(compute.parallelMode == ParallelMode::SINGLE)
            {
                unsigned const rngSeed = baseRngSeed(compute);
                for(unsigned completedIncreases = 0u;; ++completedIncreases)
                {
                    unsigned const targetRayCount = adaptiveRayTarget(experiment, compute, completedIncreases);
                    unsigned const batchRayCount = targetRayCount - rawResult.rayCount;
                    unsigned const activeDevices = std::min(maxDevices, batchRayCount);
                    if(batchRayCount == 0u)
                    {
                        if(targetRayCount == experiment.maxRays || experiment.forwardRayCount != 0u)
                            break;
                        continue;
                    }
                    if(activeDevices == 0u)
                        break;

                    std::fill(runtimes.begin(), runtimes.end(), 0.0f);
                    unsigned const launchSeed = random::seedForAdaptiveLaunch(rngSeed, adaptiveLaunches);
                    bool const terminalLaunch
                        = experiment.forwardRayCount != 0u || targetRayCount == experiment.maxRays;
                    bool const downloadAccumulators = !allowDeviceResident || activeDevices > 1u || !terminalLaunch;
                    auto const batchResult = evaluatePersistentBatch(
                        experiment,
                        hostMesh,
                        betaVolume,
                        activeDevices,
                        batchRayCount,
                        targetRayCount,
                        launchSeed,
                        runtimes,
                        adaptiveLaunches == 0u,
                        downloadAccumulators);
                    rawResult = batchResult;
                    runtime += *std::ranges::max_element(runtimes);
                    usedDevices = std::max(usedDevices, activeDevices);
                    ++adaptiveLaunches;
                    if(downloadAccumulators)
                    {
                        finalizeForwardPhiAse(hostMesh, rawResult, m_betaVolumeTotal, result);
                        hostResultAvailable = true;
                        recordAdaptiveRayConvergence(
                            result,
                            targetRayCount,
                            experiment.relativeStandardErrorThreshold,
                            convergenceRayCounts);
                    }
                    if(terminalLaunch
                       || forwardResultMeetsRelativeStandardError(result, experiment.relativeStandardErrorThreshold))
                        break;
                }
                topology.activeNodes = 1u;
                topology.activeRanks = 1u;
                topology.avgActiveRanksPerNode = 1.0;
                topology.minActiveRanksPerNode = 1u;
                topology.maxActiveRanksPerNode = 1u;
                topology.activeGpus = usedDevices;
                topology.avgGpusPerRank = static_cast<double>(usedDevices);
                topology.avgGpusPerNode = static_cast<double>(usedDevices);
                topology.minGpusPerNode = usedDevices;
                topology.maxGpusPerNode = usedDevices;
            }
            else if(compute.parallelMode == ParallelMode::MPI)
            {
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
                usedDevices = hase::core::calcForwardPhiAseMPI(
                    m_executor,
                    experiment,
                    compute,
                    hostMesh,
                    m_meshes,
                    rawResult,
                    topology,
                    runtime,
                    adaptiveLaunches,
                    convergenceRayCounts);
                hostResultAvailable = true;
#else
                throw std::runtime_error("MPI parallel mode is unavailable in this build");
#endif
            }
            else
                throw std::runtime_error("unsupported forward ASE parallel mode '" + compute.parallelMode + "'");

            if(usedDevices == 0u)
            {
                result = Result{};
                return ForwardPhiAseEvaluation{};
            }
            if(hostResultAvailable)
                finalizeForwardPhiAse(hostMesh, rawResult, m_betaVolumeTotal, result);
            else
            {
                result = Result{};
                result.srmStatus = rawResult.srmStatus;
                result.srmPasses = rawResult.srmPasses;
                result.srmRemainingFraction = rawResult.srmRemainingFraction;
                result.srmMaxIterations = rawResult.srmMaxIterations;
                result.srmDivergenceStreak = rawResult.srmDivergenceStreak;
            }
            // The persistent single-device path finalizes directly in
            // m_deviceContexts.  MPI currently combines host accumulators, so
            // its result must be copied back even when only one device exists.
            bool const deviceResidentPhi
                = allowDeviceResident && compute.parallelMode == ParallelMode::SINGLE && usedDevices == 1u;
            double const fluorescenceRate = hostMesh.nTot / hostMesh.crystalTFluo;
            for(unsigned volume = 0u;
                hostResultAvailable && volume < result.phiAse.size() && volume < hostMesh.betaVolume.size();
                ++volume)
            {
                result.phiAse[volume] *= fluorescenceRate;
                result.standardError[volume] *= fluorescenceRate;
                if(!deviceResidentPhi)
                {
                    result.dndtAse[volume] = calcVolumeDndtAse(
                        hostMesh,
                        experiment.maxSigmaA,
                        experiment.maxSigmaE,
                        result.phiAse[volume],
                        volume);
                }
            }
            if(deviceResidentPhi)
            {
                m_deviceContexts.front()->finalizeCellPhiAse(
                    primaryMeshView(betaVolume),
                    rawResult.rayCount,
                    m_betaVolumeTotal,
                    fluorescenceRate,
                    experiment.maxSigmaA,
                    experiment.maxSigmaE);
            }
            return ForwardPhiAseEvaluation{
                deviceResidentPhi,
                runtime,
                usedDevices,
                rawResult.rayCount,
                adaptiveLaunches,
                topology,
                std::move(convergenceRayCounts)};
        }

    private:
        ForwardPhiAseRawResult evaluatePersistentBatch(
            ExperimentParameters const& experiment,
            HostMesh const& hostMesh,
            auto const& betaVolume,
            unsigned activeDevices,
            unsigned rayCount,
            unsigned accumulatedRayCount,
            unsigned baseSeed,
            std::vector<float>& runtimes,
            bool resetAccumulators,
            bool downloadAccumulators)
        {
            ForwardPhiAseRawResult combined = makeForwardRawResult(hostMesh.numberOfCells);
            if(rayCount == 0u || activeDevices == 0u)
                return combined;

            unsigned const raysPerDevice = rayCount / activeDevices;
            unsigned const remainder = rayCount % activeDevices;
            double const betaVolumeTotal = m_betaVolumeTotal;
            double const sourceOffset = random::stratifiedUnitOffset(baseSeed);
            unsigned const spectrumPhase
                = random::stratifiedSpectrumPhase(baseSeed, static_cast<unsigned>(experiment.sigmaA.size()));
            std::vector<ForwardPhiAseRawResult> partials(activeDevices);
            unsigned rayOffset = 0u;
            for(unsigned deviceIndex = 0u; deviceIndex < activeDevices; ++deviceIndex)
            {
                unsigned const localRayCount
                    = deviceIndex + 1u == activeDevices ? raysPerDevice + remainder : raysPerDevice;
                m_deviceContexts[deviceIndex]->begin(
                    deviceIndex == 0u ? primaryMeshView(betaVolume) : m_meshes[deviceIndex].toView(),
                    localRayCount,
                    random::seedForWorker(baseSeed, 0u, deviceIndex),
                    rayOffset,
                    rayCount,
                    sourceOffset,
                    spectrumPhase,
                    betaVolumeTotal,
                    experiment,
                    resetAccumulators);
                rayOffset += localRayCount;
            }
            for(unsigned deviceIndex = 0u; deviceIndex < activeDevices; ++deviceIndex)
            {
                m_deviceContexts[deviceIndex]->finish(
                    partials[deviceIndex],
                    runtimes[deviceIndex],
                    downloadAccumulators);
                mergeForwardRawResult(combined, partials[deviceIndex]);
            }
            if(combined.rayCount != accumulatedRayCount)
                throw std::runtime_error("Forward ray partition accounting mismatch.");
            return combined;
        }

        [[nodiscard]] DeviceMeshView primaryMeshView(auto const& betaVolume) const
        {
            auto mesh = m_meshes.front().toView();
            mesh.betaVolume = std::span<double const>(betaVolume.data(), betaVolume.getMdSpan().getExtents().x());
            return mesh;
        }

        void refreshDynamicMeshes(
            auto const& betaVolume,
            HostMesh& hostMesh,
            bool requireHostValues,
            bool synchronizePrimaryMesh)
        {
            m_betaVolumeTotal = m_deviceContexts.front()->rebuildBetaVolumePrefix(m_meshes.front(), betaVolume);
            if(m_meshes.size() == 1u && !requireHostValues)
            {
                return;
            }

            auto queue = m_meshes.front().m_device.makeQueue(alpaka::queueKind::blocking);
            hostMesh.setBetaVolume(detail::copyToVector(queue, betaVolume));
            if(synchronizePrimaryMesh)
                detail::copyVectorToBuffer(queue, hostMesh.betaVolume, m_meshes.front().betaVolume);
            for(std::size_t index = 1u; index < m_meshes.size(); ++index)
            {
                auto& mesh = m_meshes[index];
                auto secondaryQueue = mesh.m_device.makeQueue(alpaka::queueKind::blocking);
                detail::copyVectorToBuffer(secondaryQueue, hostMesh.betaVolume, mesh.betaVolume);
                m_deviceContexts[index]->rebuildBetaVolumePrefix(mesh, mesh.betaVolume);
            }
        }

        T_Executor m_executor;
        std::vector<DeviceMeshContainer<T_Device>> m_meshes;
        std::vector<std::unique_ptr<ForwardPhiAseDeviceContext<T_Device, T_Executor>>> m_deviceContexts;
        double m_betaVolumeTotal = 0.0;
    };
} // namespace hase::core
