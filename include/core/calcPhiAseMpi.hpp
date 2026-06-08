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

#ifndef DISABLE_MPI
#    include <core/calcPhiAse.hpp>
#    include <core/calcPhiAseThreaded.hpp>
#    include <core/logging.hpp>
#    include <core/mesh.hpp>
#    include <core/types.hpp>
#    include <mpi.h>
#    include <random/random.hpp>
#    include <utils/progressbar.hpp>

#    include <algorithm>
#    include <array>
#    include <iostream>
#    include <limits>
#    include <ranges>
#    include <vector>

// Nodes
#    define HEAD_NODE 0

namespace hase::core
{

    /***
     * Performs ASE computation in MPI mode by statically partitioning the global
     * sample range across all MPI ranks. Each rank (including the HEAD-Node) computes its assigned subset
     * locally, and the partial results are gathered on the head node after all
     * ranks have finished.
     */
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    float calcPhiAseMPI(
        T_Exec exec,
        ExperimentParameters const& experiment,
        ComputeParameters const& compute,
        HostMesh const& hostMesh,
        std::vector<DeviceMeshContainer<T_Device>> const& meshes,
        Result& result,
        RuntimeTopology& topology,
        float& maxRankRuntime)
    {
        int mpiError = MPI_Init(nullptr, nullptr);
        if(mpiError != MPI_SUCCESS)
        {
            dout(V_ERROR) << "Error starting MPI program." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, mpiError);
            return 1.0f;
        }

        int rank = 0;
        int size = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        unsigned baseSeed = rank == HEAD_NODE ? (compute.rngSeed == ComputeParameters::unspecifiedRngSeed
                                                     ? hase::random::SeedGenerator::get().getSeed()
                                                     : compute.rngSeed)
                                              : 0u;
        MPI_Bcast(&baseSeed, 1, MPI_UNSIGNED, HEAD_NODE, MPI_COMM_WORLD);
        std::array<char, MPI_MAX_PROCESSOR_NAME> processorName{};
        int processorNameLength = 0;
        MPI_Get_processor_name(processorName.data(), &processorNameLength);
        processorName.at(std::min(processorNameLength, static_cast<int>(processorName.size()) - 1)) = '\0';

        MPI_Comm nodeComm;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &nodeComm);

        int localRank = 0;
        int localSize = 1;
        MPI_Comm_rank(nodeComm, &localRank);
        MPI_Comm_size(nodeComm, &localSize);

        if(rank != HEAD_NODE)
        {
            verbosity &= ~V_PROGRESS;
            verbosity &= ~V_STAT;
            verbosity &= ~V_INFO;
        }

        int const localDeviceCount = static_cast<int>(meshes.size());
        int const devicesPerRank = localDeviceCount / localSize;
        int const deviceRemainder = localDeviceCount % localSize;
        int const firstDevice = localRank * devicesPerRank + std::min(localRank, deviceRemainder);
        int const assignedDeviceCount = devicesPerRank + (localRank < deviceRemainder ? 1 : 0);

        MPI_Comm activeComm = MPI_COMM_NULL;
        int const activeColor = assignedDeviceCount > 0 ? 0 : MPI_UNDEFINED;
        MPI_Comm_split(MPI_COMM_WORLD, activeColor, rank, &activeComm);

        int activeRank = 0;
        int activeSize = 0;
        if(assignedDeviceCount > 0)
        {
            MPI_Comm_rank(activeComm, &activeRank);
            MPI_Comm_size(activeComm, &activeSize);
        }

        int const firstSample = static_cast<int>(compute.minSampleRange);
        int const sampleCount = static_cast<int>(compute.maxSampleRange - compute.minSampleRange + 1);

        int const base = activeSize > 0 ? sampleCount / activeSize : 0;
        int const rem = activeSize > 0 ? sampleCount % activeSize : 0;

        int const localBegin = firstSample + activeRank * base + std::min(activeRank, rem);
        int const localCount = assignedDeviceCount > 0 ? base + (activeRank < rem ? 1 : 0) : 0;
        int const actuallyUsedDevices = localCount > 0 ? std::min(assignedDeviceCount, localCount) : 0;
        int const activeRankContribution = actuallyUsedDevices > 0 ? 1 : 0;

        int activeRanksOnNode = 0;
        int activeDevicesOnNode = 0;
        MPI_Allreduce(&activeRankContribution, &activeRanksOnNode, 1, MPI_INT, MPI_SUM, nodeComm);
        MPI_Allreduce(&actuallyUsedDevices, &activeDevicesOnNode, 1, MPI_INT, MPI_SUM, nodeComm);

        int activeRanksTotal = 0;
        int activeDevicesTotal = 0;
        int activeNodesTotal = 0;
        int const activeNodeContribution = localRank == 0 && activeRanksOnNode > 0 ? 1 : 0;
        MPI_Allreduce(&activeRankContribution, &activeRanksTotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&actuallyUsedDevices, &activeDevicesTotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&activeNodeContribution, &activeNodesTotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        int const rankNodeMinContribution
            = activeNodeContribution ? activeRanksOnNode : std::numeric_limits<int>::max();
        int const rankNodeMaxContribution = activeNodeContribution ? activeRanksOnNode : 0;
        int const gpuNodeMinContribution
            = activeNodeContribution ? activeDevicesOnNode : std::numeric_limits<int>::max();
        int const gpuNodeMaxContribution = activeNodeContribution ? activeDevicesOnNode : 0;
        int minActiveRanksPerNode = 0;
        int maxActiveRanksPerNode = 0;
        int minGpusPerNode = 0;
        int maxGpusPerNode = 0;
        MPI_Allreduce(&rankNodeMinContribution, &minActiveRanksPerNode, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&rankNodeMaxContribution, &maxActiveRanksPerNode, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&gpuNodeMinContribution, &minGpusPerNode, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&gpuNodeMaxContribution, &maxGpusPerNode, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        if(rank == HEAD_NODE)
        {
            topology.activeNodes = static_cast<unsigned>(activeNodesTotal);
            topology.activeRanks = static_cast<unsigned>(activeRanksTotal);
            topology.activeGpus = static_cast<unsigned>(activeDevicesTotal);
            topology.avgActiveRanksPerNode
                = activeNodesTotal > 0 ? static_cast<double>(activeRanksTotal) / static_cast<double>(activeNodesTotal)
                                       : 0.0;
            topology.minActiveRanksPerNode = minActiveRanksPerNode == std::numeric_limits<int>::max()
                                                 ? 0u
                                                 : static_cast<unsigned>(minActiveRanksPerNode);
            topology.maxActiveRanksPerNode = static_cast<unsigned>(maxActiveRanksPerNode);
            topology.avgGpusPerRank = activeRanksTotal > 0 ? static_cast<double>(activeDevicesTotal)
                                                                 / static_cast<double>(activeRanksTotal)
                                                           : 0.0;
            topology.avgGpusPerNode = activeNodesTotal > 0 ? static_cast<double>(activeDevicesTotal)
                                                                 / static_cast<double>(activeNodesTotal)
                                                           : 0.0;
            topology.minGpusPerNode
                = minGpusPerNode == std::numeric_limits<int>::max() ? 0u : static_cast<unsigned>(minGpusPerNode);
            topology.maxGpusPerNode = static_cast<unsigned>(maxGpusPerNode);
        }

        std::array<int, 8> const rankWorkInfo{
            rank,
            localRank,
            localSize,
            firstDevice,
            assignedDeviceCount,
            localBegin,
            localCount,
            actuallyUsedDevices};
        std::vector<int> allRankWorkInfo;
        if(rank == HEAD_NODE)
        {
            allRankWorkInfo.resize(static_cast<std::size_t>(size) * rankWorkInfo.size());
        }
        MPI_Gather(
            rankWorkInfo.data(),
            static_cast<int>(rankWorkInfo.size()),
            MPI_INT,
            rank == HEAD_NODE ? allRankWorkInfo.data() : nullptr,
            static_cast<int>(rankWorkInfo.size()),
            MPI_INT,
            HEAD_NODE,
            MPI_COMM_WORLD);

        std::vector<char> allProcessorNames;
        if(rank == HEAD_NODE)
        {
            allProcessorNames.resize(static_cast<std::size_t>(size) * processorName.size());
        }
        MPI_Gather(
            processorName.data(),
            static_cast<int>(processorName.size()),
            MPI_CHAR,
            rank == HEAD_NODE ? allProcessorNames.data() : nullptr,
            static_cast<int>(processorName.size()),
            MPI_CHAR,
            HEAD_NODE,
            MPI_COMM_WORLD);

        if(rank == HEAD_NODE)
        {
            dout(V_INFO) << "MPI work distribution:" << std::endl;
            for(int r = 0; r < size; ++r)
            {
                auto const offset = static_cast<std::size_t>(r) * rankWorkInfo.size();
                int const worldRank = allRankWorkInfo.at(offset);
                int const nodeRank = allRankWorkInfo.at(offset + 1u);
                int const ranksOnNode = allRankWorkInfo.at(offset + 2u);
                int const firstAssignedDevice = allRankWorkInfo.at(offset + 3u);
                int const assignedDevices = allRankWorkInfo.at(offset + 4u);
                int const beginSample = allRankWorkInfo.at(offset + 5u);
                int const samples = allRankWorkInfo.at(offset + 6u);
                int const usedDevices = allRankWorkInfo.at(offset + 7u);
                char const* nodeName = allProcessorNames.data() + static_cast<std::size_t>(r) * processorName.size();

                dout(V_INFO) << "  rank " << worldRank << " node " << nodeName << " local " << nodeRank << "/"
                             << ranksOnNode << ": samples ";
                if(samples > 0)
                {
                    dout(V_INFO | V_NOLABEL)
                        << beginSample << "-" << (beginSample + samples - 1) << " (" << samples << ")";
                }
                else
                {
                    dout(V_INFO | V_NOLABEL) << "none";
                }

                dout(V_INFO | V_NOLABEL) << ", assigned GPUs ";
                if(assignedDevices > 0)
                {
                    dout(V_INFO | V_NOLABEL)
                        << firstAssignedDevice << "-" << (firstAssignedDevice + assignedDevices - 1);
                }
                else
                {
                    dout(V_INFO | V_NOLABEL) << "none";
                }

                dout(V_INFO | V_NOLABEL) << ", used GPUs " << usedDevices << std::endl;
            }
        }

        if(assignedDeviceCount == 0)
        {
            MPI_Comm_free(&nodeComm);
            MPI_Finalize();
            return 0.0f;
        }

        auto const& mesh = meshes.at(static_cast<std::size_t>(firstDevice));
        int const totalSamples = static_cast<int>(mesh.numberOfSamples);

        float runtime = 0.0f;

        if(localCount > 0)
        {
            unsigned const activeDevices
                = std::min(static_cast<unsigned>(assignedDeviceCount), static_cast<unsigned>(localCount));
            unsigned const samplesPerDevice = static_cast<unsigned>(localCount) / activeDevices;
            unsigned const sampleRemainder = static_cast<unsigned>(localCount) % activeDevices;

            std::vector<float> deviceRuntimes(meshes.size(), 0.0f);
            std::vector<ComputeParameters> rankComputes(meshes.size(), compute);
            hase::utils::ProgressBar bar;

            for(unsigned localDeviceIndex = 0; localDeviceIndex < activeDevices; ++localDeviceIndex)
            {
                unsigned const deviceIndex = static_cast<unsigned>(firstDevice) + localDeviceIndex;
                unsigned const workingSampleCount
                    = localDeviceIndex + 1u == activeDevices ? samplesPerDevice + sampleRemainder : samplesPerDevice;
                unsigned const minSampleIdx = static_cast<unsigned>(localBegin) + localDeviceIndex * samplesPerDevice;
                unsigned const maxSampleIdx = minSampleIdx + workingSampleCount - 1u;

                rankComputes.at(deviceIndex).gpu_i = compute.devices.at(deviceIndex);
                hase::alpakaUtils::DevBundle devBundle{meshes.at(deviceIndex).m_device, exec};
                unsigned const rngSeed
                    = hase::random::seedForWorker(baseSeed, static_cast<unsigned>(rank), deviceIndex);
                calcPhiAseThreaded(
                    devBundle,
                    experiment,
                    rankComputes.at(deviceIndex),
                    hostMesh,
                    meshes.at(deviceIndex),
                    result,
                    minSampleIdx,
                    maxSampleIdx,
                    rngSeed,
                    deviceRuntimes.at(deviceIndex),
                    bar);
            }

            joinAll();
            runtime = *std::ranges::max_element(deviceRuntimes);
        }

        std::vector<int> recvCounts;
        std::vector<int> displs;

        if(activeRank == HEAD_NODE)
        {
            recvCounts.resize(activeSize);
            displs.resize(activeSize);

            for(int r = 0; r < activeSize; ++r)
            {
                int rBegin = firstSample + r * base + std::min(r, rem);
                int rCount = base + (r < rem ? 1 : 0);
                recvCounts[r] = rCount;
                displs[r] = rBegin;
            }
        }

        MPI_Gatherv(
            localCount > 0 ? result.phiAse.data() + localBegin : nullptr,
            localCount,
            MPI_FLOAT,
            activeRank == HEAD_NODE ? result.phiAse.data() : nullptr,
            activeRank == HEAD_NODE ? recvCounts.data() : nullptr,
            activeRank == HEAD_NODE ? displs.data() : nullptr,
            MPI_FLOAT,
            HEAD_NODE,
            activeComm);

        MPI_Gatherv(
            localCount > 0 ? result.mse.data() + localBegin : nullptr,
            localCount,
            MPI_DOUBLE,
            activeRank == HEAD_NODE ? result.mse.data() : nullptr,
            activeRank == HEAD_NODE ? recvCounts.data() : nullptr,
            activeRank == HEAD_NODE ? displs.data() : nullptr,
            MPI_DOUBLE,
            HEAD_NODE,
            activeComm);

        MPI_Gatherv(
            localCount > 0 ? result.totalRays.data() + localBegin : nullptr,
            localCount,
            MPI_UNSIGNED,
            activeRank == HEAD_NODE ? result.totalRays.data() : nullptr,
            activeRank == HEAD_NODE ? recvCounts.data() : nullptr,
            activeRank == HEAD_NODE ? displs.data() : nullptr,
            MPI_UNSIGNED,
            HEAD_NODE,
            activeComm);

        MPI_Gatherv(
            localCount > 0 ? result.droppedRays.data() + localBegin : nullptr,
            localCount,
            MPI_UNSIGNED,
            activeRank == HEAD_NODE ? result.droppedRays.data() : nullptr,
            activeRank == HEAD_NODE ? recvCounts.data() : nullptr,
            activeRank == HEAD_NODE ? displs.data() : nullptr,
            MPI_UNSIGNED,
            HEAD_NODE,
            activeComm);

        MPI_Reduce(&runtime, &maxRankRuntime, 1, MPI_FLOAT, MPI_MAX, HEAD_NODE, activeComm);

        int totalUsedDevices = 0;
        MPI_Reduce(&actuallyUsedDevices, &totalUsedDevices, 1, MPI_INT, MPI_SUM, HEAD_NODE, activeComm);

        if(activeRank == HEAD_NODE)
        {
            hase::utils::ProgressBar bar;
            for(int i = 0; i < totalSamples; ++i)
            {
                bar.printFancyProgressBar(mesh.numberOfSamples);
            }
        }
        MPI_Comm_free(&activeComm);
        MPI_Comm_free(&nodeComm);
        MPI_Finalize();


        return activeRank == HEAD_NODE ? static_cast<float>(totalUsedDevices) : 0.0f;
    }
#endif

} // namespace hase::core
