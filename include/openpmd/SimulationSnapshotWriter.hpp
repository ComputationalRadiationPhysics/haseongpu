/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/simulationSnapshot.hpp>

#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <optional>
#include <queue>
#include <thread>
#include <utility>

namespace hase::openpmd
{
    class AsyncSimulationSnapshotWriter
    {
    public:
        using WriteSnapshot = std::function<void(core::SimulationSnapshot const&)>;

        AsyncSimulationSnapshotWriter(bool enabled, WriteSnapshot writeSnapshot)
            : m_enabled(enabled)
            , m_writeSnapshot(std::move(writeSnapshot))
        {
            if(m_enabled)
            {
                m_thread = std::thread([this] { drain(); });
            }
        }

        ~AsyncSimulationSnapshotWriter()
        {
            if(m_thread.joinable())
            {
                try
                {
                    finish();
                }
                catch(...)
                {
                }
            }
        }

        AsyncSimulationSnapshotWriter(AsyncSimulationSnapshotWriter const&) = delete;
        AsyncSimulationSnapshotWriter& operator=(AsyncSimulationSnapshotWriter const&) = delete;

        void enqueue(core::SimulationSnapshot const& snapshot)
        {
            if(!m_enabled)
            {
                return;
            }
            {
                std::scoped_lock lock{m_mutex};
                m_pending.push(snapshot);
            }
            m_ready.notify_one();
        }

        void finish()
        {
            if(!m_enabled || m_finished)
            {
                return;
            }
            {
                std::scoped_lock lock{m_mutex};
                m_pending.push(std::nullopt);
            }
            m_ready.notify_one();
            if(m_thread.joinable())
            {
                m_thread.join();
            }
            m_finished = true;
            if(m_error)
            {
                std::rethrow_exception(m_error);
            }
        }

    private:
        void drain()
        {
            try
            {
                while(true)
                {
                    std::optional<core::SimulationSnapshot> item;
                    {
                        std::unique_lock lock{m_mutex};
                        m_ready.wait(lock, [&] { return !m_pending.empty(); });
                        item = std::move(m_pending.front());
                        m_pending.pop();
                    }
                    if(!item)
                    {
                        break;
                    }
                    m_writeSnapshot(*item);
                }
            }
            catch(...)
            {
                m_error = std::current_exception();
            }
        }

        bool m_enabled = false;
        bool m_finished = false;
        WriteSnapshot m_writeSnapshot;
        std::mutex m_mutex;
        std::condition_variable m_ready;
        std::queue<std::optional<core::SimulationSnapshot>> m_pending;
        std::thread m_thread;
        std::exception_ptr m_error;
    };
} // namespace hase::openpmd
