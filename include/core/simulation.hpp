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
#include <ctime> /* time */
#ifdef HASE_ENABLE_BACKEND_TIMING
#    include <cstdlib>
#    include <fstream>
#    include <iomanip>
#endif
#include <locale> /* std::locale */
#include <numeric> /* accumulate*/
#include <stdexcept>
#include <string> /* string */
#include <vector> /* vector */


// User header files
#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/backendNames.hpp>
#include <core/SerialVersion.hpp>
#include <core/forwardPhiAseEvaluator.hpp>
#include <core/forwardPhiAseUtilities.hpp>
#include <core/logging.hpp>
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <utils/ray_histogram.hpp>
#include <utils/writeToVtk.hpp>

namespace hase::core
{

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
                if(hase::alpakaUtils::getNameForBackend(backend, sampleDevice) != compute.backend)
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

                std::vector<T_Device> devices;
                devices.reserve(compute.devices.size());
                for(auto const& gpu_i : compute.devices)
                    devices.emplace_back(devSelector.makeDevice(gpu_i));

                oneDidRun = true;
                // Statistics data
                double maxRelativeStandardError = 0;
                double avgRelativeStandardError = 0;
                unsigned highRelativeStandardError = 0;
                unsigned definedRelativeStandardErrors = 0;
                time_t starttime = time(0);
                ForwardPhiAseContext context{std::move(devices), exec, experiment, hostMesh};
#ifdef HASE_ENABLE_BACKEND_TIMING
                auto const preciseBackendStarted = std::chrono::steady_clock::now();
#endif
                auto evaluation
                    = context.evaluate(experiment, compute, hostMesh, context.primaryBetaVolume(), result, false);
                float const runtime = evaluation.runtime;
                unsigned const usedGPUs = evaluation.usedDevices;
                RuntimeTopology const& topology = evaluation.topology;
                unsigned const rayCount = evaluation.rayCount;
                unsigned const adaptiveLaunches = evaluation.adaptiveLaunches;
                auto const& convergenceRayCounts = evaluation.convergenceRayCounts;

#ifdef HASE_ENABLE_BACKEND_TIMING
                if(auto const* timingPath = std::getenv("HASE_BACKEND_TIMING_CSV"))
                {
                    std::chrono::duration<double> const preciseBackendElapsed
                        = std::chrono::steady_clock::now() - preciseBackendStarted;
                    std::ofstream timing(timingPath, std::ios::app);
                    if(timing.tellp() == 0)
                    {
                        timing << "elapsed_seconds\n";
                    }
                    timing << std::setprecision(17) << preciseBackendElapsed.count() << '\n';
                }
#endif

                if(usedGPUs == 0)
                    throw std::runtime_error("forward ASE evaluator used no devices");

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
                        rayCount,
                        experiment.maxRays,
                        experiment.relativeStandardErrorThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        tmpPhiAse,
                        vtkPath / fs::path("phiase_" + currentTime + ".vtk"),
                        rayCount,
                        experiment.maxRays,
                        experiment.relativeStandardErrorThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        result.standardError,
                        vtkPath / fs::path("standard_error_" + currentTime + ".vtk"),
                        rayCount,
                        experiment.maxRays,
                        experiment.relativeStandardErrorThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        tmpTotalRays,
                        vtkPath / fs::path("total_rays_" + currentTime + ".vtk"),
                        rayCount,
                        experiment.maxRays,
                        experiment.relativeStandardErrorThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        result.relativeStandardError,
                        vtkPath / fs::path("relative_standard_error_" + currentTime + ".vtk"),
                        rayCount,
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
                    dout(V_STAT) << "Cells             : " << context.primaryMesh().numberOfCells << std::endl;
                    dout(V_STAT) << "Samples           : "
                                 << std::min(static_cast<unsigned>(result.dndtAse.size()), numSamples) << std::endl;
                    if(experiment.forwardRayCount != 0u)
                    {
                        dout(V_STAT) << "Forward rays      : " << rayCount << " (explicit)" << std::endl;
                    }
                    else if(experiment.maxRays > experiment.minRays)
                    {
                        dout(V_STAT) << "Forward rays      : " << rayCount << " of " << experiment.minRays << " - "
                                     << experiment.maxRays << " (" << adaptiveLaunches << " launches)" << std::endl;
                    }
                    else
                    {
                        dout(V_STAT) << "Forward rays      : " << rayCount << std::endl;
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
                    hase::utils::ray_histogram(convergenceRayCounts, rayCount);
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
            for(auto const& element : hase::alpakaUtils::availableBackendNames())
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
