/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <core/types.hpp>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <sstream>
#include <string>
#include <utility>

namespace hase::benchmark
{
    struct RunContext
    {
        std::string api;
        std::string deviceKind;
        std::string executor;
        std::string deviceName;
        std::string additionalEnvironment;
        std::filesystem::path outputPath;
        unsigned gpuIndex = 0u;
        unsigned numDevices = 0u;
        unsigned maxRepetitions = 0u;
        unsigned adaptiveSteps = 0u;
        unsigned minSampleRange = 0u;
        unsigned maxSampleRange = 0u;
        unsigned minRaysPerSample = 0u;
        unsigned maxRaysPerSample = 0u;
        bool useReflections = false;
    };

    void setRunContext(RunContext context);
    RunContext currentRunContext();
    void flush();
    std::uint64_t nowUnixNs();

    class ScopedEvent
    {
    public:
        ScopedEvent(std::string kernel, std::string additionalArgs = {});
        ScopedEvent(ScopedEvent const&) = delete;
        ScopedEvent& operator=(ScopedEvent const&) = delete;
        ScopedEvent(ScopedEvent&&) = delete;
        ScopedEvent& operator=(ScopedEvent&&) = delete;
        ~ScopedEvent();

    private:
        std::string m_kernel;
        std::string m_additionalArgs;
        RunContext m_context;
        std::chrono::steady_clock::time_point m_start;
        std::uint64_t m_timestamp;
    };

    template<typename T_Queue>
    class ScopedSyncEvent
    {
    public:
        ScopedSyncEvent(T_Queue& queue, std::string kernel, std::string additionalArgs = {})
            : m_queue(queue)
            , m_kernel(std::move(kernel))
            , m_additionalArgs(std::move(additionalArgs))
        {
            alpaka::onHost::wait(m_queue);
            m_event.emplace(std::move(m_kernel), std::move(m_additionalArgs));
        }

        ScopedSyncEvent(ScopedSyncEvent const&) = delete;
        ScopedSyncEvent& operator=(ScopedSyncEvent const&) = delete;
        ScopedSyncEvent(ScopedSyncEvent&&) = delete;
        ScopedSyncEvent& operator=(ScopedSyncEvent&&) = delete;

        ~ScopedSyncEvent()
        {
            alpaka::onHost::wait(m_queue);
        }

    private:
        T_Queue& m_queue;
        std::string m_kernel;
        std::string m_additionalArgs;
        std::optional<ScopedEvent> m_event;
    };

    template<typename T_Device, typename T_Exec>
    RunContext makeRunContext(
        T_Device const& device,
        T_Exec const& executor,
        hase::core::ComputeParameters const& compute,
        hase::core::ExperimentParameters const& experiment)
    {
        auto const properties = device.getDeviceProperties();

        std::ostringstream environment;
        environment << "globalMemCapacityBytes=" << properties.globalMemCapacityBytes
                    << ";sharedMemPerBlockBytes=" << properties.sharedMemPerBlockBytes
                    << ";multiProcessorCount=" << properties.multiProcessorCount << ";warpSize=" << properties.warpSize
                    << ";maxThreadsPerBlock=" << properties.maxThreadsPerBlock;

        RunContext context;
        context.api = alpaka::onHost::getName(alpaka::getApi(device));
        context.deviceKind = alpaka::onHost::getName(alpaka::getDeviceKind(device));
        context.executor = alpaka::onHost::getName(executor);
        context.deviceName = properties.getName();
        context.additionalEnvironment = environment.str();
        context.outputPath = compute.outputPath;
        context.gpuIndex = compute.gpu_i;
        context.numDevices = compute.numDevices;
        context.maxRepetitions = compute.maxRepetitions;
        context.adaptiveSteps = compute.adaptiveSteps;
        context.minSampleRange = compute.minSampleRange;
        context.maxSampleRange = compute.maxSampleRange;
        context.minRaysPerSample = experiment.minRaysPerSample;
        context.maxRaysPerSample = experiment.maxRaysPerSample;
        context.useReflections = experiment.useReflections;
        return context;
    }

    class ScopedRunContext
    {
    public:
        template<typename T_Device, typename T_Exec>
        ScopedRunContext(
            T_Device const& device,
            T_Exec const& executor,
            hase::core::ComputeParameters const& compute,
            hase::core::ExperimentParameters const& experiment)
            : m_previous(currentRunContext())
        {
            setRunContext(makeRunContext(device, executor, compute, experiment));
        }

        ScopedRunContext(ScopedRunContext const&) = delete;
        ScopedRunContext& operator=(ScopedRunContext const&) = delete;

        ~ScopedRunContext()
        {
            flush();
            setRunContext(m_previous);
        }

    private:
        RunContext m_previous;
    };
} // namespace hase::benchmark

#ifdef HASE_ENABLE_BENCHMARK
#    define HASE_BENCHMARK_STRINGIFY_IMPL(...) #__VA_ARGS__
#    define HASE_BENCHMARK_STRINGIFY(...) HASE_BENCHMARK_STRINGIFY_IMPL(__VA_ARGS__)
#    define HASE_BENCHMARK_CONCAT_IMPL(lhs, rhs) lhs##rhs
#    define HASE_BENCHMARK_CONCAT(lhs, rhs) HASE_BENCHMARK_CONCAT_IMPL(lhs, rhs)
#    define BENCH(kernel, ...)                                                                                        \
        ::hase::benchmark::ScopedEvent HASE_BENCHMARK_CONCAT(haseBenchmarkScope_, __LINE__)(                          \
            #kernel,                                                                                                  \
            HASE_BENCHMARK_STRINGIFY(__VA_ARGS__));
#    define Bench(...) BENCH(__VA_ARGS__)
#    define BENCH_SYNC(queue, kernel, ...)                                                                            \
        ::hase::benchmark::ScopedSyncEvent HASE_BENCHMARK_CONCAT(                                                     \
            haseBenchmarkSyncScope_,                                                                                  \
            __LINE__)(queue, #kernel, HASE_BENCHMARK_STRINGIFY(__VA_ARGS__));
#    define BenchSync(...) BENCH_SYNC(__VA_ARGS__)
#else
#    define BENCH(...) static_cast<void>(0);
#    define Bench(...) static_cast<void>(0);
#    define BENCH_SYNC(...) static_cast<void>(0);
#    define BenchSync(...) static_cast<void>(0);
#endif
