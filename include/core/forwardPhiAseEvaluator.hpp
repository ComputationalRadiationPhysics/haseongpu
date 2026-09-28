/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <alpakaUtils/HybridBuffer.hpp>
#include <alpakaUtils/memory.hpp>
#include <benchmark.hpp>
#include <concepts/concepts.hpp>
#include <core/Runtime.hpp>
#include <core/calcForwardPhiAse.hpp>
#include <core/calcPhiAseThreaded.hpp>
#include <core/forwardPhiAseUtilities.hpp>
#include <core/forwardSamplingPlan.hpp>
#include <data/TraceData.hpp>
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
        /**
         * @brief Copy a device buffer into a newly allocated host vector.
         * @param queue Queue on the buffer's device; the copy waits for completion.
         * @param buffer Contiguous device buffer to download.
         * @return Host vector containing every buffer element.
         */
        auto copyToVector(hase::concepts::Queue auto const& queue, alpaka::concepts::IBuffer auto const& buffer)
        {
            using T_Value = alpaka::GetValueType_t<ALPAKA_TYPEOF(buffer)>;
            std::vector<T_Value> result(buffer.getExtents().product());
            auto hybridBuffer = hase::alpakaUtils::getHybridBuffer(result, buffer);
            hybridBuffer.toHost(queue);
            return result;
        }

        /**
         * @brief Enqueue a complete host-vector upload into an existing buffer.
         * @tparam T Element type shared by the vector and buffer.
         * @param queue Queue on the destination buffer's device.
         * @param values Host values to upload.
         * @param buffer Destination device buffer with a compatible extent.
         */
        template<typename T>
        void copyVectorToBuffer(
            hase::concepts::Queue auto const& queue,
            std::vector<T> const& values,
            alpaka::concepts::IBuffer<T> auto& buffer)
        {
            auto hybridBuffer = hase::alpakaUtils::getHybridBuffer(values, buffer);
            hybridBuffer.toDevice(queue);
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
    struct ForwardRunInputs
    {
        AseTraceControls const& experiment;
        ExecutionPolicy const& compute;
        hase::data::TraceData const& hostMesh;
        unsigned baseSeed;
        double betaVolumeTotal;
        unsigned numIndependentRayPopulations;
        std::vector<DomainCost> const& domainCosts;
        std::vector<DomainQuota> const& finalDomainQuotas;
        std::vector<hase::data::AseDomainInterface> const& interfaces;
    };

    /** Execute complete statistical populations through smaller, worker-owned transport batches. */
    template<typename T_WorkerPolicy>
    ForwardWorkerResult runForwardPopulationLaunch(
        HaseWorker<T_WorkerPolicy>& worker,
        ForwardRunInputs const& context,
        std::span<std::uint32_t const> const domainCounts,
        std::uint32_t const seed)
    {
        auto const started = std::chrono::steady_clock::now();
        auto const chunkSize = forwardExecutionChunkSize(
            std::accumulate(domainCounts.begin(), domainCounts.end(), 0u),
            worker.workerCount());
        auto const plan = makeForwardPopulationBatches(
            domainCounts,
            context.finalDomainQuotas,
            context.betaVolumeTotal,
            context.numIndependentRayPopulations,
            chunkSize);
        auto const descriptors = worker.gather(worker.descriptor());
        auto const schedule
            = makeForwardPopulationSchedule(*descriptors, context.domainCosts, plan, context.interfaces);
        std::vector<ForwardPopulationBatch> localWork;
        for(auto const& batch : plan)
            if(schedule.owner({batch.domainId, batch.rayPopulationId, batch.batchId}) == worker.workerIndex())
                localWork.push_back(batch);
        (void) worker(PrepareRayPopulationWork{{}, localWork, seed});
        // Every device has finished preparing its source rays before ANY worker traces.
        if(worker.workerCount() > 1u)
            (void) worker.gather(std::uint32_t{1u});
        ForwardPhiAseRawResult boundary;
        for(std::uint32_t population = 0u; population < context.numIndependentRayPopulations; ++population)
        {
            std::vector<ForwardPopulationRay> incoming;
            std::vector<double> fractions;
            double initialWeight = 0.0;
            double previousWeight = 0.0;
            std::uint32_t grows = 0u;
            auto status = !context.experiment.useReflections && context.experiment.domainCount <= 1u
                              ? data::BoundaryStatus::disabled
                              : data::BoundaryStatus::converged;
            std::uint32_t passes = 0u;
            double remaining = 0.0;
            bool const resident = worker.workerCount() == 1u;
            for(std::uint32_t pass = 0u;; ++pass)
            {
                double weight = 0.0;
                std::vector<ForwardPopulationRay> candidates;
                if(resident)
                    weight = worker(TraceResidentRayPopulationWork{{}, population, pass == 0u});
                else
                {
                    auto local = worker(TraceRayPopulationWork{{}, population, incoming, pass == 0u});
                    if(!context.experiment.useReflections && context.experiment.domainCount <= 1u)
                        break;
                    auto gathered = worker.gather(std::move(local));
                    if(worker.isRoot())
                    {
                        for(auto const& part : *gathered)
                            candidates.insert(candidates.end(), part.begin(), part.end());
                        std::ranges::sort(candidates, {}, &ForwardPopulationRay::ordinal);
                        for(auto const& ray : candidates)
                        {
                            if(!std::isfinite(ray.weight) || ray.weight < 0.0)
                                throw std::runtime_error("invalid boundary ray population weight");
                            weight += ray.weight;
                        }
                    }
                    weight = worker.scatter(weight);
                }
                if(!context.experiment.useReflections && context.experiment.domainCount <= 1u)
                    break;
                if(!std::isfinite(weight) || weight < 0.0)
                    throw std::runtime_error("non-finite boundary ray population total");
                if(pass == 0u)
                {
                    initialWeight = weight;
                    previousWeight = weight;
                    if(weight == 0.0)
                        break;
                    fractions.push_back(1.0);
                    remaining = 1.0;
                }
                else
                {
                    passes = pass;
                    remaining = weight / initialWeight;
                    fractions.push_back(remaining);
                    if(weight == 0.0 || remaining < context.experiment.reflectionTolerance)
                        break;
                    if(weight > previousWeight)
                    {
                        ++grows;
                        if(grows >= 3u && estimateBoundaryTail(fractions).divergent)
                        {
                            status = data::BoundaryStatus::diverged;
                            break;
                        }
                    }
                    else
                    {
                        grows = 0u;
                        if(std::abs(weight - previousWeight) / std::max(weight, 1.0e-30)
                           < context.experiment.reflectionTolerance)
                        {
                            status = data::BoundaryStatus::stable;
                            break;
                        }
                    }
                    previousWeight = weight;
                }
                if(pass >= context.experiment.resolvedBoundaryMaxPasses())
                {
                    status = data::BoundaryStatus::maxPasses;
                    break;
                }
                if(resident)
                    worker(
                        SelectResidentRayPopulationWork{
                            {},
                            kernels::forward::rayPopulationSeed(seed, population),
                            pass});
                else
                {
                    std::vector<std::vector<ForwardPopulationRay>> outgoing;
                    if(worker.isRoot())
                    {
                        incoming = worker(
                            SelectRayPopulationWork{
                                {},
                                std::move(candidates),
                                kernels::forward::rayPopulationSeed(seed, population),
                                pass});
                        std::vector<std::vector<std::size_t>> byDomain(context.domainCosts.size());
                        for(std::size_t index = 0u; index < incoming.size(); ++index)
                        {
                            auto const& ray = incoming[index];
                            if(ray.domainId >= byDomain.size())
                                throw std::runtime_error("boundary ray targets an unknown domain");
                            byDomain[ray.domainId].push_back(index);
                        }
                        std::vector<DomainWorkItem> routeWork;
                        for(std::uint32_t domain = 0u; domain < byDomain.size(); ++domain)
                            for(std::uint32_t begin = 0u, batch = 0u; begin < byDomain[domain].size(); ++batch)
                            {
                                auto const size = static_cast<std::uint32_t>(
                                    std::min<std::size_t>(chunkSize, byDomain[domain].size() - begin));
                                routeWork.push_back({{domain, population, batch}, size, 0u});
                                begin += size;
                            }
                        auto const routeSchedule
                            = makeDomainSchedule(*descriptors, context.domainCosts, routeWork, context.interfaces);
                        outgoing.resize(worker.workerCount());
                        for(auto const& item : routeWork)
                        {
                            auto const& indices = byDomain[item.id.domain];
                            auto const begin = static_cast<std::size_t>(item.id.batchId) * chunkSize;
                            auto& target = outgoing.at(routeSchedule.owner(item.id));
                            for(std::size_t offset = 0u; offset < item.rayCount; ++offset)
                                target.push_back(incoming[indices[begin + offset]]);
                        }
                    }
                    incoming = worker.distribute(std::move(outgoing));
                }
            }
            auto const tail = estimateBoundaryTail(fractions);
            if(tail.divergent)
                status = data::BoundaryStatus::diverged;
            if(boundaryStatusPriority(status) > boundaryStatusPriority(boundary.boundaryStatus))
                boundary.boundaryStatus = status;
            boundary.boundaryPasses = std::max(boundary.boundaryPasses, passes);
            boundary.boundaryRemainingFraction = std::max(boundary.boundaryRemainingFraction, remaining);
            if(status != data::BoundaryStatus::disabled)
                boundary.boundaryMaxPasses = context.experiment.resolvedBoundaryMaxPasses();
            boundary.boundaryDivergenceStreak = 3u;
            boundary.boundaryGamma = std::max(boundary.boundaryGamma, tail.gamma);
            boundary.boundaryGammaStandardError
                = std::max(boundary.boundaryGammaStandardError, tail.gammaStandardError);
            boundary.boundaryTailFactor = std::max(boundary.boundaryTailFactor, tail.tailFactor);
            boundary.boundaryTailClosure = std::max(boundary.boundaryTailClosure, tail.tailClosure);
        }
        auto raw = worker(CollectRayPopulationWork{});
        mergeForwardBoundaryResult(raw, boundary);
        auto const elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        return {std::move(raw), static_cast<float>(elapsed)};
    }

    /** Fixed logical SRM batches retain their own reservoirs, independent of device ownership. */
    template<typename T_WorkerPolicy>
    ForwardWorkerResult runForwardLogicalSrmLaunch(
        HaseWorker<T_WorkerPolicy>& worker,
        ForwardRunInputs const& context,
        std::span<std::uint32_t const> const domainCounts,
        std::uint32_t const seed)
    {
        auto const started = std::chrono::steady_clock::now();
        auto const plan = makeForwardPopulationBatches(
            domainCounts,
            context.finalDomainQuotas,
            context.betaVolumeTotal,
            context.numIndependentRayPopulations,
            maxLogicalSrmBatchRays);
        auto const descriptors = worker.gather(worker.descriptor());
        auto const schedule
            = makeForwardPopulationSchedule(*descriptors, context.domainCosts, plan, context.interfaces);
        std::vector<ForwardPopulationBatch> localWork;
        for(auto const& batch : plan)
            if(schedule.owner({batch.domainId, batch.rayPopulationId, batch.batchId}) == worker.workerIndex())
                localWork.push_back(batch);
        (void) worker(PrepareRayPopulationWork{{}, localWork, seed});
        // All populations' source/wavelength samples exist before transport starts.
        (void) worker.gather(std::uint32_t{1u});
        (void) worker(TraceLogicalSrmBatches{{}, seed});
        auto raw = worker(CollectRayPopulationWork{});
        auto const elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        return {std::move(raw), static_cast<float>(elapsed)};
    }

    /** @brief Raw and convergence results produced by the worker-group simulation loop. */
    struct ForwardSimulationResult
    {
        ForwardPhiAseRawResult raw;
        data::PhiAseResult convergence;
        float runtime = 0.0f;
        unsigned adaptiveLaunches = 0u;
        std::vector<unsigned> convergenceRayCounts;
    };

    /**
     * @brief Run the adaptive ASE simulation using a policy-selected worker group.
     *
     * Integration-stage orchestration calls this function once for the current
     * beta state. Domain/batch items are mapped to workers, grouped into
     * queue-efficient launches on their owned devices, gathered with statistical
     * batch identity intact, and only then combined for normalization and RSE.
     * @tparam T_WorkerPolicy Worker policy supplying mapping, collectives, and work dispatch.
     * @param worker Participating worker with a stable group identity.
     * @param context Trace controls and immutable launch inputs shared by the group.
     * @return Gathered raw accumulators, convergence result, runtime, and launch metadata.
     */
    template<typename T_WorkerPolicy>
    [[nodiscard]] ForwardSimulationResult runForwardSimulation(
        HaseWorker<T_WorkerPolicy>& worker,
        ForwardRunInputs const& context)
    {
        ForwardSimulationResult simulation;
        simulation.raw = makeForwardRawResult(
            context.hostMesh.numberOfCells,
            context.hostMesh.numberOfMaterials * context.hostMesh.numberOfMeshPoints,
            context.numIndependentRayPopulations);
        simulation.convergenceRayCounts.assign(context.hostMesh.numberOfCells, 0u);
        std::vector<std::uint64_t> previousDomainCounts(context.finalDomainQuotas.size(), 0u);
        unsigned const baseSeed = worker.scatter(context.baseSeed);
        // adaptive sampling loop
        for(unsigned completedIncreases = 0u;; ++completedIncreases)
        {
            auto const launch = planForwardLaunch(
                context.experiment,
                context.compute,
                context.finalDomainQuotas,
                previousDomainCounts,
                context.numIndependentRayPopulations,
                simulation.raw.rayCount,
                completedIncreases);
            unsigned const targetRayCount = launch.target;
            unsigned const launchSeed = random::seedForAdaptiveLaunch(baseSeed, simulation.adaptiveLaunches);
            auto const previousBatchRayCounts = simulation.raw.rayPopulationRayCounts;
            auto const& launchDomainCounts = launch.domainCounts;
            for(std::size_t domain = 0u; domain < launchDomainCounts.size(); ++domain)
                previousDomainCounts[domain] += launchDomainCounts[domain];

            auto localResult = context.experiment.reflectionMode == "direct"
                                   ? runForwardPopulationLaunch(worker, context, launchDomainCounts, launchSeed)
                                   : runForwardLogicalSrmLaunch(worker, context, launchDomainCounts, launchSeed);
            simulation.runtime += worker.reduce(
                localResult.runtime,
                [](float const lhs, float const rhs) { return std::max(lhs, rhs); });

            auto const gathered = worker.gather(std::move(localResult));
            for(auto const& workerResult : *gathered)
                mergeForwardRawResult(simulation.raw, workerResult.raw);
            if(simulation.raw.rayCount != targetRayCount)
                throw std::runtime_error("forward statistical batch accounting mismatch");
            for(unsigned batch = 0u; batch < context.numIndependentRayPopulations; ++batch)
            {
                unsigned expectedBatchRays = 0u;
                for(auto const domainCount : launchDomainCounts)
                    expectedBatchRays += domainCount / context.numIndependentRayPopulations
                                         + (batch < domainCount % context.numIndependentRayPopulations ? 1u : 0u);
                if(simulation.raw.rayPopulationRayCounts[batch] != previousBatchRayCounts[batch] + expectedBatchRays)
                    throw std::runtime_error("forward statistical batch accounting mismatch");
            }

            ++simulation.adaptiveLaunches;
            bool stop = false;
            if(worker.isRoot())
            {
                simulation.convergence = worker(FinalizeForwardAse{simulation.raw, true, false});
                recordAdaptiveRayConvergence(
                    simulation.convergence,
                    targetRayCount,
                    context.experiment.relativeStandardErrorThreshold,
                    simulation.convergenceRayCounts);
                stop = context.experiment.forwardRayCount != 0u || context.compute.adaptiveSteps == 0u
                       || targetRayCount == context.experiment.maxRays
                       || forwardResultMeetsRelativeStandardError(
                           simulation.convergence,
                           context.experiment.relativeStandardErrorThreshold);
            }
            stop = worker.scatter(stop);
            if(stop)
                break;
        }
        if(worker.isRoot())
            simulation.convergence = worker(FinalizeForwardAse{simulation.raw, false, true});
        else if(worker.requiresFinalizedDeviceState())
            simulation.convergence = worker(FinalizeForwardAse{simulation.raw, true, true});
        return simulation;
    }

    /**
     * @brief Owns persistent device state for one prepared multi-domain trace.
     *
     * The context keeps geometry, materials, routing tables, source CDFs, and
     * boundary workspaces resident across evaluations. The scheduler already
     * exposes each worker's required domain set; physical mesh allocation is
     * still replicated until domain-shard-only residency is introduced.
     */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Executor>
    class ForwardPhiAseContext
    {
    public:
        /**
         * @brief Allocate persistent trace and evaluator state on a non-empty device set.
         * @param devices Devices owned by this context, one per local worker.
         * @param executor Executor copied into each device context.
         * @param experiment Controls used to size optional reflection scratch storage.
         * @param hostMesh Host arrays that back and initialize the resident traces.
         */
        ForwardPhiAseContext(
            std::vector<T_Device> devices,
            T_Executor executor,
            AseTraceControls const& experiment,
            hase::data::TraceData& hostMesh,
            hase::data::AseDomainGraph& domains)
            : m_executor(std::move(executor))
            , m_domains(domains)
            , m_domainCosts(makeDomainCosts(domains))
        {
            if(devices.empty())
                throw std::runtime_error("forward ASE context requires at least one device");
            m_meshes.reserve(devices.size());
            m_interfaceMaps.reserve(devices.size());
            m_domainSources.reserve(devices.size());
            for(auto& device : devices)
            {
                m_meshes.emplace_back(hostMesh.makeResident(device));
                m_interfaceMaps.emplace_back(device, domains);
                m_domainSources.emplace_back(device, domains);
                auto queue = device.makeQueue(alpaka::queueKind::nonBlocking);
                m_meshes.back().toDevice(queue);
                m_interfaceMaps.back().toDevice(queue);
                m_domainSources.back().toDevice(queue);
                alpaka::onHost::wait(queue);
            }
            m_deviceContexts.reserve(m_meshes.size());
            for(auto const& mesh : m_meshes)
                m_deviceContexts.emplace_back(
                    std::make_unique<ForwardPhiAseDeviceContext<T_Device, T_Executor>>(
                        mesh.m_device,
                        m_executor,
                        experiment,
                        hostMesh));
        }

        /** @return First device, used as the time integrator's primary resident device. */
        [[nodiscard]] T_Device& primaryDevice()
        {
            return m_meshes.front().m_device;
        }

        /** @return Resident trace owned by the primary device. */
        [[nodiscard]] hase::data::ResidentTrace<T_Device>& primaryMesh()
        {
            return m_meshes.front();
        }

        /** @return Non-owning device view of the primary trace's excitation array. */
        [[nodiscard]] auto primaryBetaVolume()
        {
            return m_meshes.front().betaVolume.toDeviceView();
        }

        /** @return Whether excitation must be downloaded to distribute it to secondary devices. */
        [[nodiscard]] bool requiresHostBetaVolume() const
        {
            return m_meshes.size() > 1u;
        }

        /** @return Cell-ordered ASE population derivative downloaded from the primary device. */
        std::vector<double> downloadPrimaryVolumeDndtAse()
        {
            return m_deviceContexts.front()->downloadVolumeDndtAse();
        }

        /**
         * @param includePhiAse Whether to download cell ASE flux.
         * @param includeStandardError Whether to download absolute standard error.
         * @param includeRelativeStandardError Whether to download relative standard error.
         * @param includeTotalRays Whether to download per-cell ray visits.
         * @return Finalized result containing only the requested large arrays.
         */
        data::PhiAseResult downloadPrimaryResult(
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

        /** @return Owning device buffer containing the primary ASE population derivative. */
        [[nodiscard]] auto& primaryVolumeDndtAse()
        {
            return m_deviceContexts.front()->volumeDndtAse();
        }

        /** @return Owning device buffer containing the primary finalized ASE flux. */
        [[nodiscard]] auto& primaryVolumePhiAse()
        {
            return m_deviceContexts.front()->volumePhiAse();
        }

        /**
         * @brief Refresh only material and spectral buffers on every owned device.
         *
         * Geometry allocations are deliberately retained. The next evaluate
         * call rebuilds source-strength prefixes from its current beta buffer.
         * @param hostTrace Prepared host trace supplying replacement material arrays.
         */
        void refreshMaterials(hase::data::TraceData& hostTrace)
        {
            for(auto& resident : m_meshes)
            {
                auto queue = resident.m_device.makeQueue(alpaka::queueKind::nonBlocking);
                resident.refreshMaterials(hostTrace, queue);
                alpaka::onHost::wait(queue);
            }
        }

        /**
         * @brief Execute adaptive forward tracing for the supplied excitation state.
         * @param experiment Physical and statistical controls; retained for worker dispatch.
         * @param compute Backend, device, seed, and adaptive scheduling controls.
         * @param hostMesh Host trace used for multi-device or MPI synchronization.
         * @param betaVolume Cell excitation view resident on the primary device.
         * @param result Host convergence result replaced by the final adaptive result.
         * @param allowDeviceResident Whether finalized primary buffers may remain device-only.
         * @return Runtime, device-topology, ray-count, and convergence metadata.
         */
        ForwardPhiAseEvaluation evaluate(
            AseTraceControls& experiment,
            ExecutionPolicy& compute,
            hase::data::TraceData& hostMesh,
            alpaka::concepts::IView<double> auto const& betaVolume,
            data::PhiAseResult& result,
            bool const allowDeviceResident = true)
        {
#ifdef HASE_ENABLE_BENCHMARK
            hase::benchmark::ScopedRunContext benchmarkContext{primaryDevice(), m_executor, compute, experiment};
#endif
            BENCH(AseEvaluation);
            bool const mpiMode = compute.parallelMode == ParallelMode::MPI;
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
            if(mpiMode)
                detail::ensureMpiInitialized();
#endif
            refreshDynamicMeshes(betaVolume, hostMesh, requiresHostBetaVolume() || mpiMode, mpiMode);
            if(!experiment.isForwardPropagation())
                throw std::runtime_error("Only forward volume propagation is supported by the openPMD backend.");

            unsigned seed = compute.rngSeed;
            if(seed == ExecutionPolicy::unspecifiedRngSeed)
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
            if(experiment.forwardRayCount == 0u
               && (experiment.minRays == 0u || experiment.maxRays < experiment.minRays))
                throw std::invalid_argument("adaptive ASE requires 0 < minRays <= maxRays");
            auto const finalRayCount = experiment.forwardRayCount != 0u
                                           ? experiment.forwardRayCount
                                           : (compute.adaptiveSteps == 0u ? experiment.minRays : experiment.maxRays);
            auto const finalDomainQuotas = allocateDomainRays(m_domainCosts, finalRayCount);
            unsigned const numIndependentRayPopulations = experiment.numIndependentRayPopulations;
            if(numIndependentRayPopulations == 0u
               || domainRayPopulationCount(finalDomainQuotas, numIndependentRayPopulations)
                      != numIndependentRayPopulations)
                throw std::invalid_argument(
                    "ASE budget must represent every emitting domain in every independent ray population");
            for(auto& deviceContext : m_deviceContexts)
                deviceContext->configureRayPopulationCount(numIndependentRayPopulations);
            ForwardRunInputs simulationContext{
                experiment,
                compute,
                hostMesh,
                seed,
                m_betaVolumeTotal,
                numIndependentRayPopulations,
                m_domainCosts,
                finalDomainQuotas,
                m_domains.interfaces};
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
                std::vector<std::jthread> workers;
                workers.reserve(threadWorkerCount);
                try
                {
                    for(unsigned workerIndex = 0u; workerIndex < threadWorkerCount; ++workerIndex)
                    {
                        workers.emplace_back(
                            [&, workerIndex]
                            {
                                try
                                {
                                    auto mesh = workerIndex == 0u ? primaryMeshView(betaVolume)
                                                                  : m_meshes[workerIndex].view();
                                    HaseWorker worker{ThreadOwnedDevices{
                                        workerIndex,
                                        threadWorkerCount,
                                        group,
                                        mesh,
                                        *m_deviceContexts[workerIndex],
                                        experiment,
                                        m_betaVolumeTotal,
                                        m_interfaceMaps[workerIndex].view(),
                                        m_domainSources[workerIndex].view()}};
                                    workerResults[workerIndex] = runForwardSimulation(worker, simulationContext);
                                }
                                catch(...)
                                {
                                    exceptions[workerIndex] = std::current_exception();
                                    group.cancel(exceptions[workerIndex]);
                                }
                            });
                    }
                }
                catch(...)
                {
                    group.cancel(std::current_exception());
                    throw;
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
                    deviceIndex == 0u ? primaryMeshView(betaVolume) : m_meshes[deviceIndex].view(),
                    *m_deviceContexts[deviceIndex],
                    experiment,
                    m_betaVolumeTotal,
                    hostMesh.numberOfCells,
                    hostMesh.numberOfMaterials * hostMesh.numberOfMeshPoints,
                    numIndependentRayPopulations,
                    m_interfaceMaps[deviceIndex].view(),
                    m_domainSources[deviceIndex].view()}};
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
                    m_betaVolumeTotal);
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
        [[nodiscard]] hase::data::TraceView primaryMeshView(
            alpaka::concepts::IView<double> auto const& betaVolume) const
        {
            auto mesh = m_meshes.front().view();
            mesh.betaVolume = std::span<double const>(betaVolume.data(), betaVolume.getExtents().x());
            return mesh;
        }

        void refreshDynamicMeshes(
            alpaka::concepts::IView<double> auto const& betaVolume,
            hase::data::TraceData& hostMesh,
            bool const requireHostValues,
            bool const synchronizePrimaryMesh)
        {
            m_betaVolumeTotal = m_deviceContexts.front()->rebuildSourceStrengthPrefix(m_meshes.front(), betaVolume);
            auto primaryDevBundle = hase::alpakaUtils::DevBundle{m_meshes.front().m_device, m_executor};
            auto primaryQueue = m_meshes.front().m_device.makeQueue(alpaka::queueKind::nonBlocking);
            m_domainSources.front().rebuild(primaryDevBundle, primaryQueue, primaryMeshView(betaVolume));
            auto const domainSourceStrengths = m_domainSources.front().downloadSourceStrengthTotals(primaryQueue);
            if(domainSourceStrengths.size() != m_domainCosts.size())
                throw std::runtime_error("domain source totals do not match scheduling statistics");
            for(std::size_t domain = 0u; domain < m_domainCosts.size(); ++domain)
                m_domainCosts[domain].sourceStrength = domainSourceStrengths[domain];
            if(m_meshes.size() == 1u && !requireHostValues)
                return;

            auto queue = m_meshes.front().m_device.makeQueue(alpaka::queueKind::nonBlocking);
            auto synchronizedBetaVolume = hase::alpakaUtils::getHybridBuffer(hostMesh.betaVolume, betaVolume);
            synchronizedBetaVolume.toHost(queue);
            hostMesh.rebuildSourceStrengthPrefix();
            if(synchronizePrimaryMesh)
            {
                m_meshes.front().betaVolume.toDevice(queue);
                alpaka::onHost::wait(queue);
                m_deviceContexts.front()->rebuildSourceStrengthPrefix(
                    m_meshes.front(),
                    m_meshes.front().betaVolume.toDeviceView());
            }
            for(std::size_t index = 1u; index < m_meshes.size(); ++index)
            {
                auto& mesh = m_meshes[index];
                auto secondaryQueue = mesh.m_device.makeQueue(alpaka::queueKind::nonBlocking);
                mesh.betaVolume.toDevice(secondaryQueue);
                alpaka::onHost::wait(secondaryQueue);
                m_deviceContexts[index]->rebuildSourceStrengthPrefix(mesh, mesh.betaVolume.toDeviceView());
                auto secondaryDevBundle = hase::alpakaUtils::DevBundle{mesh.m_device, m_executor};
                m_domainSources[index].rebuild(secondaryDevBundle, secondaryQueue, mesh.view());
                alpaka::onHost::wait(secondaryQueue);
            }
        }

        T_Executor m_executor;
        hase::data::AseDomainGraph& m_domains;
        std::vector<DomainCost> m_domainCosts;
        std::vector<hase::data::ResidentTrace<T_Device>> m_meshes;
        std::vector<ResidentAseDomainInterfaces<T_Device>> m_interfaceMaps;
        std::vector<ResidentAseDomainSources<T_Device>> m_domainSources;
        std::vector<std::unique_ptr<ForwardPhiAseDeviceContext<T_Device, T_Executor>>> m_deviceContexts;
        double m_betaVolumeTotal = 0.0;
    };
} // namespace hase::core
