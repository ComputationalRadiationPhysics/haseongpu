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
#include <queue>
#include <stdexcept>
#include <thread>
#include <utility>

namespace hase::openpmd
{
    class AsyncSimulationSnapshotWriter
    {
    public:
        using WriteSnapshot = std::function<void(core::SimulationSnapshot const&)>;

        AsyncSimulationSnapshotWriter(bool enabled, WriteSnapshot writeSnapshot, bool asynchronous = true)
            : m_enabled(enabled)
            , m_asynchronous(asynchronous)
            , m_writeSnapshot(std::move(writeSnapshot))
        {
            if(m_enabled && m_asynchronous)
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
            if(!m_asynchronous)
            {
                m_writeSnapshot(snapshot);
                return;
            }
            std::unique_lock lock{m_mutex};
            m_spaceAvailable.wait(
                lock,
                [&] { return m_pending.size() < maxPendingSnapshots || m_error || m_finishRequested; });
            if(m_error)
            {
                std::rethrow_exception(m_error);
            }
            if(m_finishRequested)
            {
                throw std::logic_error("cannot enqueue a simulation snapshot after finish was requested");
            }
            m_pending.push(snapshot);
            lock.unlock();
            m_ready.notify_one();
        }

        void finish()
        {
            if(!m_enabled || m_finished)
            {
                return;
            }
            if(!m_asynchronous)
            {
                m_finished = true;
                return;
            }
            {
                std::scoped_lock lock{m_mutex};
                m_finishRequested = true;
            }
            m_ready.notify_one();
            if(m_thread.joinable())
            {
                m_thread.join();
            }
            m_finished = true;
            std::exception_ptr error;
            {
                std::scoped_lock lock{m_mutex};
                error = m_error;
            }
            if(error)
            {
                std::rethrow_exception(error);
            }
        }

    private:
        void drain()
        {
            try
            {
                while(true)
                {
                    core::SimulationSnapshot item;
                    {
                        std::unique_lock lock{m_mutex};
                        m_ready.wait(lock, [&] { return !m_pending.empty() || m_finishRequested; });
                        if(m_pending.empty())
                        {
                            break;
                        }
                        item = std::move(m_pending.front());
                        m_pending.pop();
                    }
                    m_spaceAvailable.notify_one();
                    m_writeSnapshot(item);
                }
            }
            catch(...)
            {
                {
                    std::scoped_lock lock{m_mutex};
                    m_error = std::current_exception();
                    m_finishRequested = true;
                    std::queue<core::SimulationSnapshot> empty;
                    m_pending.swap(empty);
                }
                m_spaceAvailable.notify_all();
            }
        }

        static constexpr std::size_t maxPendingSnapshots = 1u;
        bool m_enabled = false;
        bool m_asynchronous = true;
        bool m_finished = false;
        bool m_finishRequested = false;
        WriteSnapshot m_writeSnapshot;
        std::mutex m_mutex;
        std::condition_variable m_ready;
        std::condition_variable m_spaceAvailable;
        std::queue<core::SimulationSnapshot> m_pending;
        std::thread m_thread;
        std::exception_ptr m_error;
    };
} // namespace hase::openpmd
