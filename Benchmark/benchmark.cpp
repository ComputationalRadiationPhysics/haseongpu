/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "benchmark.hpp"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace hase::benchmark
{
    namespace
    {
        struct Row
        {
            RunContext context;
            std::string kernel;
            std::string additionalArgs;
            std::uint64_t timestamp;
            std::uint64_t endUnixNs;
            std::uint64_t durationNs;
        };

        std::mutex& benchmarkMutex()
        {
            static std::mutex mutex;
            return mutex;
        }

        std::vector<Row>& rows()
        {
            static std::vector<Row> storage;
            return storage;
        }

        RunContext& globalContext()
        {
            static thread_local RunContext context;
            return context;
        }

        std::string csvEscape(std::string value)
        {
            bool needsQuotes = false;
            std::string escaped;
            escaped.reserve(value.size());
            for(char const c : value)
            {
                if(c == '"')
                {
                    escaped += "\"\"";
                    needsQuotes = true;
                }
                else
                {
                    escaped += c;
                    needsQuotes = needsQuotes || c == ',' || c == '\n' || c == '\r';
                }
            }

            if(needsQuotes)
            {
                return "\"" + escaped + "\"";
            }
            return escaped;
        }

        std::filesystem::path csvPath(RunContext const& context)
        {
            if(context.outputPath.empty())
            {
                return "hase_benchmark.csv";
            }
            return context.outputPath / "hase_benchmark.csv";
        }

        void flushAtExit()
        {
            flush();
        }

        void registerFlushAtExit()
        {
            static bool registered = []()
            {
                std::atexit(flushAtExit);
                return true;
            }();
            static_cast<void>(registered);
        }
    } // namespace

    std::uint64_t nowUnixNs()
    {
        auto const now = std::chrono::system_clock::now().time_since_epoch();
        return static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(now).count());
    }

    void setRunContext(RunContext context)
    {
        std::lock_guard lock(benchmarkMutex());
        globalContext() = std::move(context);
        registerFlushAtExit();
    }

    RunContext currentRunContext()
    {
        std::lock_guard lock(benchmarkMutex());
        return globalContext();
    }

    ScopedEvent::ScopedEvent(std::string kernel, std::string additionalArgs)
        : m_kernel(std::move(kernel))
        , m_additionalArgs(std::move(additionalArgs))
        , m_context(currentRunContext())
        , m_start(std::chrono::steady_clock::now())
        , m_timestamp(nowUnixNs())
    {
        registerFlushAtExit();
    }

    ScopedEvent::~ScopedEvent()
    {
        auto const end = std::chrono::steady_clock::now();
        auto const endUnixNs = nowUnixNs();
        auto const durationNs
            = static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - m_start).count());

        std::lock_guard lock(benchmarkMutex());
        rows().push_back(
            Row{m_context, std::move(m_kernel), std::move(m_additionalArgs), m_timestamp, endUnixNs, durationNs});
    }

    void flush()
    {
        std::vector<Row> snapshot;
        std::filesystem::path path;
        {
            std::lock_guard lock(benchmarkMutex());
            if(rows().empty())
            {
                return;
            }
            snapshot = rows();
            rows().clear();
            path = csvPath(snapshot.front().context);
        }

        if(path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }

        bool const writeHeader = !std::filesystem::exists(path) || std::filesystem::is_empty(path);
        std::ofstream out(path, std::ios::app);
        if(!out)
        {
            return;
        }

        if(writeHeader)
        {
            out << "API,DeviceKind,Executor,Timestamp,DeviceName,DeviceIndex,Kernel,AdditionalArgs,EndUnixNs,"
                   "DurationNs,numDevices,MaxRepetitions,AdaptiveSteps,MinSampleRange,MaxSampleRange,MinRaysPerSample,"
                   "MaxRaysPerSample,UseReflections,Environment\n";
        }

        for(auto const& row : snapshot)
        {
            out << csvEscape(row.context.api) << ',' << csvEscape(row.context.deviceKind) << ','
                << csvEscape(row.context.executor) << ',' << row.timestamp << ',' << csvEscape(row.context.deviceName)
                << ',' << row.context.gpuIndex << ',' << csvEscape(row.kernel) << ',' << csvEscape(row.additionalArgs)
                << ',' << row.endUnixNs << ',' << row.durationNs << ',' << row.context.numDevices << ','
                << row.context.maxRepetitions << ',' << row.context.adaptiveSteps << ',' << row.context.minSampleRange
                << ',' << row.context.maxSampleRange << ',' << row.context.minRaysPerSample << ','
                << row.context.maxRaysPerSample << ',' << row.context.useReflections << ','
                << csvEscape(row.context.additionalEnvironment) << '\n';
        }
    }
} // namespace hase::benchmark
