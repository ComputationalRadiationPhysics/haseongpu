/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <algorithm> /* std::max */
#include <chrono> /* std::chrono::system_clock */
#include <cstdlib> /* getenv, strtoul */
#include <ctime> /* time */
#include <locale> /* std::locale */
#include <numeric> /* accumulate*/
#include <stdexcept>
#include <string> /* string */
#include <vector> /* vector */


// User header files
#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <core/SerialVersion.hpp>
#include <core/calcForwardPhiAse.hpp>
#include <core/calcPhiAseThreaded.hpp>
#include <core/logging.hpp>
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <random/random.hpp>
#include <utils/ray_histogram.hpp>
#include <utils/writeToVtk.hpp>
#if !defined(DISABLE_MPI) && defined(MPI_FOUND)
#    include <core/calcPhiAseMpi.hpp>
#endif

namespace hase::core
{

    /**
     * @brief Calculates dndt ASE from phi ASE values
     *
     * @param mesh needed for some constants
     * @param sigmaA absorption
     * @param sigmaE emission
     * @param phiAse results from calcPhiAse
     * @param sample_i index of sample point
     * @return dndtAse
     *
     */
    inline double calcDndtAse(
        HostMesh const& mesh,
        double const sigmaA,
        double const sigmaE,
        float const phiAse,
        unsigned const sample_i)
    {
        double const gainPerDensity = mesh.betaCells[sample_i] * (sigmaE + sigmaA) - sigmaA;
        return gainPerDensity * phiAse;
    }

    inline double calcVolumeDndtAse(
        HostMesh const& mesh,
        double const sigmaA,
        double const sigmaE,
        float const phiAse,
        unsigned const volume)
    {
        double const gainPerDensity = mesh.betaVolume[volume] * (sigmaE + sigmaA) - sigmaA;
        return gainPerDensity * phiAse;
    }

    inline unsigned baseRngSeed(ComputeParameters const& compute)
    {
        if(compute.rngSeed == ComputeParameters::unspecifiedRngSeed)
        {
            return random::SeedGenerator::get().getSeed();
        }
        return compute.rngSeed;
    }

    std::string getNameForBackend(auto const& backend, auto const& device)
    {
        std::string backendName;
        backendName += alpaka::onHost::getName(alpaka::getApi(device)) + "_";
        backendName += alpaka::onHost::getName(alpaka::getDeviceKind(device)) + "_";
        backendName += alpaka::onHost::getName(backend[alpaka::object::exec]);
        return backendName;
    }

    static std::vector<std::string> backendList()
    {
        auto backends = alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors);
        std::vector<std::string> list;
        alpaka::onHost::executeForEachIfHasDevice(
            [&](auto const& backend) -> int
            {
                auto devSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
                auto sampleDevice = devSelector.makeDevice(0);
                list.emplace_back(getNameForBackend(backend, sampleDevice));
                return 0;
            },
            backends);
        return list;
    }

    inline unsigned envUnsigned(char const* name, unsigned fallback = 0)
    {
        char const* value = std::getenv(name);
        if(value == nullptr || *value == '\0')
        {
            return fallback;
        }

        char* end = nullptr;
        unsigned long parsed = std::strtoul(value, &end, 10);
        return end == value ? fallback : static_cast<unsigned>(parsed);
    }

    inline RuntimeTopology detectRuntimeTopology()
    {
        RuntimeTopology topology;

        unsigned const worldSize = envUnsigned(
            "OMPI_COMM_WORLD_SIZE",
            envUnsigned("PMI_SIZE", envUnsigned("PMIX_SIZE", envUnsigned("SLURM_NTASKS", 1))));
        unsigned const localSize = std::max(
            1u,
            envUnsigned(
                "OMPI_COMM_WORLD_LOCAL_SIZE",
                envUnsigned("MPI_LOCALNRANKS", envUnsigned("MV2_COMM_WORLD_LOCAL_SIZE", 1))));
        unsigned const slurmNodes = envUnsigned("SLURM_JOB_NUM_NODES", 0);
        unsigned const activeNodes
            = slurmNodes > 0 ? slurmNodes : std::max(1u, (worldSize + localSize - 1u) / localSize);

        topology.activeNodes = activeNodes;
        topology.activeRanks = worldSize;
        topology.avgActiveRanksPerNode = static_cast<double>(worldSize) / static_cast<double>(activeNodes);
        topology.minActiveRanksPerNode = localSize;
        topology.maxActiveRanksPerNode = localSize;
        return topology;
    }

    bool isSelected(auto const& backend, auto const& device, ComputeParameters& compute)
    {
        if(getNameForBackend(backend, device) == compute.backend)
        {
            return true;
        }
        return false;
    }

    template<bool MATLAB>
    int startSimulation(
        ExperimentParameters& experiment,
        ComputeParameters& compute,
        Result& result,
        HostMesh& hostMesh)
    {
        auto backends = alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors);
        bool oneDidRun = false;
        auto i = alpaka::onHost::executeForEachIfHasDevice(
            [&](auto const& backend) -> int
            {
                auto deviceSpec = backend[alpaka::object::deviceSpec];
                auto exec = backend[alpaka::object::exec];

                auto devSelector = alpaka::onHost::makeDeviceSelector(deviceSpec);

                std::size_t deviceCount = devSelector.getDeviceCount();
                if(deviceCount == 0u)
                {
                    return 0;
                }
                compute.devices = std::vector<unsigned>(deviceCount);
                std::iota(compute.devices.begin(), compute.devices.end(), 0u);
                compute.gpu_i = compute.devices.front();

                if(compute.numDevices == 0)
                {
                    compute.numDevices = deviceCount;
                }
                using T_Device = ALPAKA_TYPEOF(devSelector.makeDevice(0));
                T_Device sampleDevice = devSelector.makeDevice(0);
                if(!isSelected(backend, sampleDevice, compute))
                {
                    return 0;
                }
                if(compute.numDevices > deviceCount)
                {
                    dout(V_WARNING) << "Requested number of GPUs via --numDevices (" << compute.numDevices
                                    << ") exceeds the number of available devices (" << deviceCount
                                    << ") on the current backend/node. "
                                       "HASEonGPU will use all available GPUs instead."
                                    << std::endl;
                    compute.numDevices = deviceCount;
                }
                compute.devices.resize(compute.numDevices);

                std::vector<DeviceMeshContainer<T_Device>> meshes;
                for(auto const& gpu_i : compute.devices)
                {
                    // use the first device
                    alpaka::onHost::Device device = devSelector.makeDevice(gpu_i);
                    meshes.emplace_back(hostMesh.toDevice(device));
                }

                oneDidRun = true;
                // Statistics data
                float runtime = 0.0;
                double maxRelativeStandardError = 0;
                double avgRelativeStandardError = 0;
                unsigned highRelativeStandardError = 0;
                unsigned definedRelativeStandardErrors = 0;
                time_t starttime = time(0);
                unsigned maxDevices = compute.devices.size();
                std::vector runtimes(maxDevices, 0.f);
                unsigned usedGPUs = 0;
                RuntimeTopology topology;
                if(!experiment.isForwardPropagation())
                {
                    throw std::runtime_error("Only forward volume propagation is supported by the openPMD backend.");
                }

                ForwardPhiAseRawResult rawResult;
                unsigned adaptiveLaunches = 0u;
                std::vector<unsigned> convergenceRayCounts;
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
                            {
                                break;
                            }
                            continue;
                        }
                        if(activeDevices == 0u)
                        {
                            break;
                        }

                        std::fill(runtimes.begin(), runtimes.end(), 0.0f);
                        ForwardPhiAseRawResult const batchResult = calcForwardPhiAseOnDevices(
                            exec,
                            experiment,
                            hostMesh,
                            meshes,
                            0u,
                            activeDevices,
                            batchRayCount,
                            0u,
                            batchRayCount,
                            random::seedForAdaptiveLaunch(rngSeed, adaptiveLaunches),
                            0u,
                            runtimes);
                        mergeForwardRawResult(rawResult, batchResult);
                        runtime += *std::ranges::max_element(runtimes);
                        usedGPUs = std::max(usedGPUs, activeDevices);
                        ++adaptiveLaunches;
                        finalizeForwardPhiAse(hostMesh, rawResult, result);
                        recordAdaptiveRayConvergence(
                            result,
                            targetRayCount,
                            experiment.relativeStandardErrorThreshold,
                            convergenceRayCounts);

                        if(experiment.forwardRayCount != 0u || targetRayCount == experiment.maxRays
                           || forwardResultMeetsRelativeStandardError(
                               result,
                               experiment.relativeStandardErrorThreshold))
                        {
                            break;
                        }
                    }
                    topology = RuntimeTopology{};
                    topology.activeNodes = 1u;
                    topology.activeRanks = 1u;
                    topology.avgActiveRanksPerNode = 1.0;
                    topology.minActiveRanksPerNode = 1u;
                    topology.maxActiveRanksPerNode = 1u;
                    topology.activeGpus = usedGPUs;
                    topology.avgGpusPerRank = static_cast<double>(usedGPUs);
                    topology.avgGpusPerNode = static_cast<double>(usedGPUs);
                    topology.minGpusPerNode = usedGPUs;
                    topology.maxGpusPerNode = usedGPUs;
                }
                else if(compute.parallelMode == ParallelMode::MPI)
                {
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
                    usedGPUs = hase::core::calcForwardPhiAseMPI(
                        exec,
                        experiment,
                        compute,
                        hostMesh,
                        meshes,
                        rawResult,
                        topology,
                        runtime,
                        adaptiveLaunches,
                        convergenceRayCounts);
#else
#    if !defined(MPI_FOUND)
                    dout(V_ERROR) << "Did not find MPI on your system!";
                    exit(1);
#    else
                    dout(V_ERROR) << "TURN 'DISABLE_MPI' to 'OFF' in order to run PhiASE on multiple nodes!";
                    exit(1);
#    endif
#endif
                }

                else
                {
                    dout(V_ERROR) << "No valid parallel-mode for GPU!" << std::endl;
                    exit(1);
                }

                if(usedGPUs == 0)
                {
                    return 0;
                }
                finalizeForwardPhiAse(hostMesh, rawResult, result);

                dout(V_INFO) << "Active nodes             : " << topology.activeNodes << std::endl;
                dout(V_INFO) << "Active ranks             : " << topology.activeRanks << std::endl;
                dout(V_INFO) << "Active ranks per node    : " << topology.avgActiveRanksPerNode
                             << " avg (min=" << topology.minActiveRanksPerNode
                             << ", max=" << topology.maxActiveRanksPerNode << ")" << std::endl;
                dout(V_INFO) << "Active GPUs              : " << topology.activeGpus << std::endl;
                dout(V_INFO) << "GPUs per active rank     : " << topology.avgGpusPerRank << " avg" << std::endl;
                dout(V_INFO) << "GPUs per active node     : " << topology.avgGpusPerNode
                             << " avg (min=" << topology.minGpusPerNode << ", max=" << topology.maxGpusPerNode << ")"
                             << std::endl;

                for(unsigned volume = 0u; volume < result.phiAse.size() && volume < hostMesh.betaVolume.size();
                    ++volume)
                {
                    double const fluorescenceRate = hostMesh.nTot / hostMesh.crystalTFluo;
                    result.phiAse.at(volume) *= fluorescenceRate;
                    result.standardError.at(volume) *= fluorescenceRate;
                    result.dndtAse.at(volume) = calcVolumeDndtAse(
                        hostMesh,
                        experiment.maxSigmaA,
                        experiment.maxSigmaE,
                        result.phiAse.at(volume),
                        volume);
                }
                /***************************************************************************
                 * PRINT SOLUTION
                 **************************************************************************/
                if(verbosity & V_DEBUG)
                {
                    unsigned const debugResultSize = static_cast<unsigned>(result.phiAse.size());
                    for(unsigned sample_i = 0; sample_i < debugResultSize; ++sample_i)
                    {
                        dout(V_DEBUG) << "PHI ASE[" << sample_i << "]: " << result.phiAse.at(sample_i) << " "
                                      << result.standardError.at(sample_i)
                                      << " (RSE=" << result.relativeStandardError.at(sample_i) << ")" << std::endl;
                        if(sample_i >= 10)
                            break;
                    }
                }

                /***************************************************************************
                 * WRITE VTK FILES
                 **************************************************************************/
                if(compute.writeVtk)
                {
                    std::vector<double> tmpPhiAse(result.phiAse.begin(), result.phiAse.end());
                    std::vector<double> tmpTotalRays(result.totalRays.begin(), result.totalRays.end());

                    fs::path vtkPath = compute.outputPath / fs::path("vtk");
                    if(fs::exists(compute.outputPath))
                    {
                        fs::create_directory(vtkPath);
                    }

                    std::string currentTime
                        = std::to_string(std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        result.dndtAse,
                        vtkPath / fs::path("dndt_" + currentTime + ".vtk"),
                        rawResult.rayCount,
                        experiment.maxRays,
                        experiment.relativeStandardErrorThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        tmpPhiAse,
                        vtkPath / fs::path("phiase_" + currentTime + ".vtk"),
                        rawResult.rayCount,
                        experiment.maxRays,
                        experiment.relativeStandardErrorThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        result.standardError,
                        vtkPath / fs::path("standard_error_" + currentTime + ".vtk"),
                        rawResult.rayCount,
                        experiment.maxRays,
                        experiment.relativeStandardErrorThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        tmpTotalRays,
                        vtkPath / fs::path("total_rays_" + currentTime + ".vtk"),
                        rawResult.rayCount,
                        experiment.maxRays,
                        experiment.relativeStandardErrorThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        result.relativeStandardError,
                        vtkPath / fs::path("relative_standard_error_" + currentTime + ".vtk"),
                        rawResult.rayCount,
                        experiment.maxRays,
                        experiment.relativeStandardErrorThreshold,
                        experiment.useReflections,
                        runtime);
                }

                /***************************************************************************
                 * PRINT STATISTICS
                 **************************************************************************/
                if(verbosity & V_STAT)
                {
                    unsigned numSamples = compute.maxSampleRange - compute.minSampleRange + 1;
                    for(unsigned sample_i = compute.minSampleRange;
                        sample_i <= compute.maxSampleRange && sample_i < result.relativeStandardError.size();
                        ++sample_i)
                    {
                        double const element = result.relativeStandardError[sample_i];
                        if(!std::isfinite(element))
                        {
                            continue;
                        }
                        maxRelativeStandardError = std::max(maxRelativeStandardError, element);
                        avgRelativeStandardError += element;
                        ++definedRelativeStandardErrors;
                        if(element >= experiment.relativeStandardErrorThreshold)
                        {
                            std::cout << " too high relative standard error at element: " << sample_i
                                      << " rse: " << element << std::endl;
                            highRelativeStandardError++;
                        }
                    }
                    if(definedRelativeStandardErrors > 0u)
                    {
                        avgRelativeStandardError /= static_cast<double>(definedRelativeStandardErrors);
                    }

                    try
                    {
                        std::cout.imbue(std::locale(""));
                    }
                    catch(std::runtime_error const&)
                    {
                    }

                    dout(V_STAT | V_NOLABEL) << std::endl;
                    dout(V_STAT) << "=== Statistics ===" << std::endl;
                    dout(V_STAT) << "Backend       : " << compute.backend << std::endl;
                    dout(V_STAT) << "RNG Seed      : " << baseRngSeed(compute) << std::endl;
                    dout(V_STAT) << "ParallelMode      : " << compute.parallelMode << std::endl;
                    dout(V_STAT) << "Prisms            : " << meshes[0].numberOfPrisms << std::endl;
                    dout(V_STAT) << "Samples           : "
                                 << std::min(static_cast<unsigned>(result.dndtAse.size()), numSamples) << std::endl;
                    if(experiment.forwardRayCount != 0u)
                    {
                        dout(V_STAT) << "Forward rays      : " << rawResult.rayCount << " (explicit)" << std::endl;
                    }
                    else if(experiment.maxRays > experiment.minRays)
                    {
                        dout(V_STAT) << "Forward rays      : " << rawResult.rayCount << " of " << experiment.minRays
                                     << " - " << experiment.maxRays << " (" << adaptiveLaunches << " launches)"
                                     << std::endl;
                    }
                    else
                    {
                        dout(V_STAT) << "Forward rays      : " << rawResult.rayCount << std::endl;
                    }
                    dout(V_STAT) << "sum(totalRays)    : "
                                 << std::accumulate(result.totalRays.begin(), result.totalRays.end(), 0.) << std::endl;
                    dout(V_STAT) << "RSE threshold     : " << experiment.relativeStandardErrorThreshold << std::endl;
                    dout(V_STAT) << "int. Wavelength   : " << experiment.sigmaA.size() << std::endl;
                    dout(V_STAT) << "max. RSE          : " << maxRelativeStandardError << std::endl;
                    dout(V_STAT) << "avg. RSE          : " << avgRelativeStandardError << std::endl;
                    dout(V_STAT) << "too high RSE      : " << highRelativeStandardError << " of "
                                 << definedRelativeStandardErrors << " defined" << std::endl;

                    if constexpr(alpaka::thisApi() == alpaka::api::cuda || alpaka::thisApi() == alpaka::api::hip)
                    {
                        dout(V_STAT) << "Nr of GPU's        : " << usedGPUs << std::endl;
                    }
                    else
                    {
                        dout(V_STAT) << "Nr of Device's   : " << usedGPUs << std::endl;
                    }
                    dout(V_STAT) << "Simulation runtime: " << runtime << "s" << std::endl;
                    dout(V_STAT) << "Total runtime     : " << difftime(time(0), starttime) << "s" << std::endl;
                    dout(V_STAT) << std::endl;
                    dout(V_STAT) << "=== Adaptive forward-ray convergence by cell (green: RSE target reached; red: "
                                    "budget exhausted) ==="
                                 << std::endl;
                    hase::utils::ray_histogram(convergenceRayCounts, rawResult.rayCount);
                    dout(V_STAT) << std::endl;
                }
                // Cleanup device memory
                // TODO: replace by smart pointer for device memory
                std::cout << " " << "\n";
                return 0;
            },
            backends);
        if(!oneDidRun)
        {
            std::cout << "\n------------------------------ HASEONGPU ERROR ------------------------------\n"
                      << std::endl;
            std::cout << " Backend did not match any available backend with available device! \n Available backend "
                         "specifications are: "
                      << std::endl;
            for(auto const& element : backendList())
            {
                std::cout << element << "\n";
            }
            std::cout << "Run hase-configure to generate a matching backend/openPMD setup." << std::endl;
            std::cout << "\n------------------------------------------------------------------------------\n"
                      << std::endl;
        }
        return i || !oneDidRun;
    }


} // namespace hase::core
