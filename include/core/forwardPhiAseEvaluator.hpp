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
#include <exception>
#include <memory>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <thread>
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

    /** @brief Metadata returned by one complete forward ASE evaluation. */
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

    /** @brief Policy-independent inputs for the adaptive forward simulation loop. */
    struct ForwardSimulationContext
    {
        ExperimentParameters const& experiment;
        ComputeParameters const& compute;
        HostMesh const& hostMesh;
        unsigned baseSeed;
        double betaVolumeTotal;
        unsigned batchCount;
    };

    /** @brief Raw and convergence results produced by the worker-group simulation loop. */
    struct ForwardSimulationResult
    {
        ForwardPhiAseRawResult raw;
        Result convergence;
        float runtime = 0.0f;
        unsigned adaptiveLaunches = 0u;
        std::vector<unsigned> convergenceRayCounts;
    };

    /** @brief Execute all complete batches mapped to one worker for one adaptive launch. */
    template<typename T_Worker>
    [[nodiscard]] ForwardRayBatchResults executeMappedBatches(
        T_Worker& worker,
        unsigned const launchRayCount,
        unsigned const launchSeed,
        unsigned const batchCount)
    {
        ForwardRayBatchResults localResults;
        for(auto [batch] : hase::mapIdx(worker, alpaka::IdxRange{batchCount}))
        {
            localResults.emplace_back(worker(
                ForwardRayBatch{
                    batch,
                    kernels::forward::rseBatchRayCount(0u, launchRayCount, batch, batchCount),
                    launchSeed}));
        }
        return localResults;
    }

    /** @brief Flatten gathered worker containers while retaining batch indices. */
    [[nodiscard]] inline ForwardRayBatchResults flattenBatchResults(
        std::vector<ForwardRayBatchResults> const& workerResults,
        unsigned const batchCount)
    {
        ForwardRayBatchResults batches;
        for(auto const& results : workerResults)
            batches.insert(batches.end(), results.begin(), results.end());
        std::ranges::sort(batches, {}, &ForwardRayBatchResult::index);
        if(batches.size() != batchCount)
            throw std::runtime_error("forward worker gather did not return every statistical batch");
        for(unsigned batch = 0u; batch < batches.size(); ++batch)
            if(batches[batch].index != batch)
                throw std::runtime_error("forward worker gather returned a missing or duplicate statistical batch");
        return batches;
    }

    /**
     * @brief Run the adaptive ASE simulation using a policy-selected worker group.
     *
     * Integration-stage orchestration calls this function once for the current
     * beta state. Complete ray batches are mapped to workers, traced on their
     * owned devices, gathered with their batch identities intact, and only then
     * combined for normalization and adaptive RSE evaluation.
     */
    template<typename T_WorkerPolicy>
    [[nodiscard]] ForwardSimulationResult runForwardSimulation(
        HaseWorker<T_WorkerPolicy>& worker,
        ForwardSimulationContext const& context)
    {
        ForwardSimulationResult simulation;
        simulation.raw = makeForwardRawResult(
            context.hostMesh.numberOfCells,
            context.hostMesh.numberOfMeshPoints,
            context.batchCount);
        simulation.convergenceRayCounts.assign(context.hostMesh.numberOfCells, 0u);
        unsigned const baseSeed = worker.scatter(context.baseSeed);
        // adaptive sampling loop
        for(unsigned completedIncreases = 0u;; ++completedIncreases)
        {
            unsigned const targetRayCount = adaptiveRayTarget(context.experiment, context.compute, completedIncreases);
            unsigned const launchRayCount = targetRayCount - simulation.raw.rayCount;
            unsigned const launchSeed = random::seedForAdaptiveLaunch(baseSeed, simulation.adaptiveLaunches);

            ForwardRayBatchResults localResults;
            // batch loop
            for(auto [batch] : hase::mapIdx(worker, alpaka::IdxRange{context.batchCount}))
            {
                localResults.emplace_back(worker(
                    ForwardRayBatch{
                        batch,
                        kernels::forward::rseBatchRayCount(0u, launchRayCount, batch, context.batchCount),
                        launchSeed}));
            }
            float const localRuntime = std::accumulate(
                localResults.cbegin(),
                localResults.cend(),
                0.0f,
                [](float const sum, ForwardRayBatchResult const& batch) { return sum + batch.runtime; });
            simulation.runtime
                += worker.reduce(localRuntime, [](float const lhs, float const rhs) { return std::max(lhs, rhs); });

            auto const gathered = worker.gather(std::move(localResults));
            auto batches = flattenBatchResults(*gathered, context.batchCount);
            for(auto const& batch : batches)
                mergeForwardRawResult(simulation.raw, batch.raw);
            if(simulation.raw.rayCount != targetRayCount)
                throw std::runtime_error("forward statistical batch accounting mismatch");

            ++simulation.adaptiveLaunches;
            simulation.convergence = worker(
                FinalizeForwardAse{
                    simulation.raw,
                    context.hostMesh.nTot / context.hostMesh.crystalTFluo,
                    context.experiment.maxSigmaA,
                    context.experiment.maxSigmaE});
            recordAdaptiveRayConvergence(
                simulation.convergence,
                targetRayCount,
                context.experiment.relativeStandardErrorThreshold,
                simulation.convergenceRayCounts);
            bool const stop = context.experiment.forwardRayCount != 0u || targetRayCount == context.experiment.maxRays
                              || forwardResultMeetsRelativeStandardError(
                                  simulation.convergence,
                                  context.experiment.relativeStandardErrorThreshold);
            if(stop)
                break;
        }
        return simulation;
    }

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
                m_deviceContexts.emplace_back(
                    std::make_unique<ForwardPhiAseDeviceContext<T_Device, T_Executor>>(
                        mesh.m_device,
                        m_executor,
                        experiment,
                        hostMesh));
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
            bool const includePhiAse,
            bool const includeStandardError,
            bool const includeRelativeStandardError,
            bool const includeTotalRays)
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
            bool const allowDeviceResident = true)
        {
            bool const mpiMode = compute.parallelMode == ParallelMode::MPI;
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
            if(mpiMode)
                detail::ensureMpiInitialized();
#endif
            unsigned workerCount = static_cast<unsigned>(m_meshes.size());
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
            if(mpiMode)
            {
                int mpiWorkerCount = 1;
                MPI_Comm_size(MPI_COMM_WORLD, &mpiWorkerCount);
                workerCount = static_cast<unsigned>(mpiWorkerCount);
            }
#endif
            unsigned const batchCount = hase::kernels::forward::forwardRseBatchCount(workerCount);
            for(auto& deviceContext : m_deviceContexts)
                deviceContext->configureBatchCount(batchCount);
            refreshDynamicMeshes(betaVolume, hostMesh, requiresHostBetaVolume() || mpiMode, mpiMode);
            if(!experiment.isForwardPropagation())
                throw std::runtime_error("Only forward volume propagation is supported by the openPMD backend.");

            unsigned seed = compute.rngSeed;
            if(seed == ComputeParameters::unspecifiedRngSeed)
            {
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
                int rank = 0;
                if(mpiMode)
                    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                seed = !mpiMode || rank == 0 ? random::SeedGenerator::get().getSeed() : 0u;
#else
                seed = random::SeedGenerator::get().getSeed();
#endif
            }
            ForwardSimulationContext
                simulationContext{experiment, compute, hostMesh, seed, m_betaVolumeTotal, batchCount};
            ForwardSimulationResult simulation;
            RuntimeTopology topology;
            unsigned usedDevices = 0u;
            unsigned residentDeviceIndex = 0u;

            if(compute.parallelMode == ParallelMode::SINGLE)
            {
                unsigned const threadWorkerCount = static_cast<unsigned>(m_meshes.size());
                detail::ThreadWorkerGroup group(threadWorkerCount);
                std::vector<ForwardSimulationResult> workerResults(threadWorkerCount);
                std::vector<std::exception_ptr> exceptions(threadWorkerCount);
                std::vector<std::thread> workers;
                workers.reserve(threadWorkerCount);
                for(unsigned workerIndex = 0u; workerIndex < threadWorkerCount; ++workerIndex)
                {
                    workers.emplace_back(
                        [&, workerIndex]
                        {
                            try
                            {
                                auto mesh
                                    = workerIndex == 0u ? primaryMeshView(betaVolume) : m_meshes[workerIndex].toView();
                                HaseWorker worker{ThreadOwnedDevices{
                                    workerIndex,
                                    threadWorkerCount,
                                    group,
                                    mesh,
                                    *m_deviceContexts[workerIndex],
                                    experiment,
                                    m_betaVolumeTotal}};
                                workerResults[workerIndex] = runForwardSimulation(worker, simulationContext);
                            }
                            catch(...)
                            {
                                exceptions[workerIndex] = std::current_exception();
                            }
                        });
                }
                for(auto& worker : workers)
                    worker.join();
                for(auto const& exception : exceptions)
                    if(exception)
                        std::rethrow_exception(exception);
                simulation = std::move(workerResults.front());
                usedDevices = threadWorkerCount;
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
                unsigned const deviceIndex = mpiRankDeviceIndex(static_cast<unsigned>(m_meshes.size()));
                residentDeviceIndex = deviceIndex;
                HaseWorker worker{MPIRank{
                    MPI_COMM_WORLD,
                    deviceIndex == 0u ? primaryMeshView(betaVolume) : m_meshes[deviceIndex].toView(),
                    *m_deviceContexts[deviceIndex],
                    experiment,
                    m_betaVolumeTotal,
                    hostMesh.numberOfCells,
                    hostMesh.numberOfMeshPoints,
                    batchCount}};
                simulation = runForwardSimulation(worker, simulationContext);
                topology = mpiWorkerTopology();
                usedDevices = topology.activeGpus;
#else
                throw std::runtime_error("MPI parallel mode is unavailable in this build");
#endif
            }
            else
                throw std::runtime_error("unsupported forward ASE parallel mode '" + compute.parallelMode + "'");

            result = std::move(simulation.convergence);
            if(allowDeviceResident && residentDeviceIndex != 0u)
            {
                m_deviceContexts.front()->uploadAndFinalize(
                    primaryMeshView(betaVolume),
                    simulation.raw,
                    m_betaVolumeTotal,
                    hostMesh.nTot / hostMesh.crystalTFluo,
                    experiment.maxSigmaA,
                    experiment.maxSigmaE);
            }

            return ForwardPhiAseEvaluation{
                allowDeviceResident,
                simulation.runtime,
                usedDevices,
                simulation.raw.rayCount,
                simulation.adaptiveLaunches,
                topology,
                std::move(simulation.convergenceRayCounts)};
        }

    private:
        [[nodiscard]] DeviceMeshView primaryMeshView(auto const& betaVolume) const
        {
            auto mesh = m_meshes.front().toView();
            mesh.betaVolume = std::span<double const>(betaVolume.data(), betaVolume.getMdSpan().getExtents().x());
            return mesh;
        }

        void refreshDynamicMeshes(
            auto const& betaVolume,
            HostMesh& hostMesh,
            bool const requireHostValues,
            bool const synchronizePrimaryMesh)
        {
            m_betaVolumeTotal = m_deviceContexts.front()->rebuildBetaVolumePrefix(m_meshes.front(), betaVolume);
            if(m_meshes.size() == 1u && !requireHostValues)
                return;

            auto queue = m_meshes.front().m_device.makeQueue(alpaka::queueKind::nonBlocking);
            hostMesh.setBetaVolume(detail::copyToVector(queue, betaVolume));
            if(synchronizePrimaryMesh)
            {
                detail::copyVectorToBuffer(queue, hostMesh.betaVolume, m_meshes.front().betaVolume);
                alpaka::onHost::wait(queue);
                m_deviceContexts.front()->rebuildBetaVolumePrefix(m_meshes.front(), m_meshes.front().betaVolume);
            }
            for(std::size_t index = 1u; index < m_meshes.size(); ++index)
            {
                auto& mesh = m_meshes[index];
                auto secondaryQueue = mesh.m_device.makeQueue(alpaka::queueKind::nonBlocking);
                detail::copyVectorToBuffer(secondaryQueue, hostMesh.betaVolume, mesh.betaVolume);
                alpaka::onHost::wait(secondaryQueue);
                m_deviceContexts[index]->rebuildBetaVolumePrefix(mesh, mesh.betaVolume);
            }
        }

        T_Executor m_executor;
        std::vector<DeviceMeshContainer<T_Device>> m_meshes;
        std::vector<std::unique_ptr<ForwardPhiAseDeviceContext<T_Device, T_Executor>>> m_deviceContexts;
        double m_betaVolumeTotal = 0.0;
    };
} // namespace hase::core
