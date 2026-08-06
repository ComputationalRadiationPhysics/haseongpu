/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#ifdef HASE_ENABLE_HOST_ROUTINE_TIMING

#    include <chrono>
#    include <cstdint>
#    include <cstdlib>
#    include <fstream>
#    include <map>
#    include <mutex>
#    include <string>
#    include <string_view>

namespace hase::core::timing
{
    struct HostRoutineAggregate
    {
        std::uint64_t calls = 0u;
        double totalSeconds = 0.0;
    };

    inline std::map<std::string, HostRoutineAggregate>& hostRoutineAggregates()
    {
        static std::map<std::string, HostRoutineAggregate> aggregates;
        return aggregates;
    }

    inline std::mutex& hostRoutineMutex()
    {
        static std::mutex mutex;
        return mutex;
    }

    class ScopedHostRoutine
    {
    public:
        explicit ScopedHostRoutine(std::string_view name) : m_name{name}, m_started{std::chrono::steady_clock::now()}
        {
        }

        ~ScopedHostRoutine()
        {
            try
            {
                std::chrono::duration<double> const elapsed = std::chrono::steady_clock::now() - m_started;
                std::scoped_lock lock{hostRoutineMutex()};
                auto& aggregate = hostRoutineAggregates()[m_name];
                ++aggregate.calls;
                aggregate.totalSeconds += elapsed.count();
            }
            catch(...)
            {
            }
        }

        ScopedHostRoutine(ScopedHostRoutine const&) = delete;
        ScopedHostRoutine& operator=(ScopedHostRoutine const&) = delete;

    private:
        std::string m_name;
        std::chrono::steady_clock::time_point m_started;
    };

    inline void writeHostRoutineTimingCsv()
    {
        auto const* path = std::getenv("HASE_HOST_ROUTINE_TIMING_CSV");
        if(path == nullptr || *path == '\0')
            return;
        auto const* revision = std::getenv("HASE_BENCHMARK_REVISION");
        auto const* backend = std::getenv("HASE_BENCHMARK_BACKEND");
        std::ofstream output{path};
        output << "revision,backend,routine,calls,total_seconds,mean_seconds\n";
        std::scoped_lock lock{hostRoutineMutex()};
        for(auto const& [name, aggregate] : hostRoutineAggregates())
        {
            double const mean = aggregate.calls == 0u ? 0.0 : aggregate.totalSeconds / aggregate.calls;
            output << (revision ? revision : "") << ',' << (backend ? backend : "") << ',' << name << ','
                   << aggregate.calls << ',' << aggregate.totalSeconds << ',' << mean << '\n';
        }
    }
} // namespace hase::core::timing

#    define HASE_DETAIL_JOIN_IMPL(lhs, rhs) lhs##rhs
#    define HASE_DETAIL_JOIN(lhs, rhs) HASE_DETAIL_JOIN_IMPL(lhs, rhs)
#    define HASE_HOST_ROUTINE_SCOPE(name)                                                                             \
        ::hase::core::timing::ScopedHostRoutine HASE_DETAIL_JOIN(haseHostRoutineScope_, __LINE__)                     \
        {                                                                                                             \
            name                                                                                                      \
        }

#else

namespace hase::core::timing
{
    inline void writeHostRoutineTimingCsv()
    {
    }
} // namespace hase::core::timing

#    define HASE_HOST_ROUTINE_SCOPE(name) static_cast<void>(0)

#endif
