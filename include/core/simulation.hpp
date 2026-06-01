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
#include <core/calcPhiAseThreaded.hpp>
#include <core/logging.hpp>
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <parse/parser.hpp> /* Backend, ParallelMode */
#include <utils/ray_histogram.hpp>
#include <utils/writeMatlabOutput.hpp>
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
            [&](auto const& backend)
            {
                auto devSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
                auto sampleDevice = devSelector.makeDevice(0);
                list.emplace_back(getNameForBackend(backend, sampleDevice));
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
            [&](auto const& backend)
            {
                auto deviceSpec = backend[alpaka::object::deviceSpec];
                auto exec = backend[alpaka::object::exec];

                auto devSelector = alpaka::onHost::makeDeviceSelector(deviceSpec);

                // this code block is not executed if deviceCount is zero
                std::size_t deviceCount = devSelector.getDeviceCount();
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
                double maxMSE = 0;
                float avgMSE = 0;
                unsigned highMSE = 0;
                time_t starttime = time(0);
                unsigned maxDevices = compute.devices.size();
                std::vector<float> runtimes(maxDevices, 0);
                unsigned usedGPUs = 0;
                RuntimeTopology topology;
                std::vector<ComputeParameters> computes(maxDevices, compute);
                hase::utils::ProgressBar bar;
                // Single Threaded @TODO use for Validation


                if(compute.parallelMode == ParallelMode::SINGLE)
                {
                    unsigned const samplesPerNode = compute.maxSampleRange - compute.minSampleRange + 1u;

                    unsigned const activeDevices = std::min(maxDevices, samplesPerNode);
                    unsigned const samplesPerGpu = samplesPerNode / activeDevices;
                    unsigned const remainder = samplesPerNode % activeDevices;

                    for(unsigned deviceIndex = 0u; deviceIndex < activeDevices; ++deviceIndex)
                    {
                        auto& device_i = compute.devices.at(deviceIndex);

                        // last device does the remaining work
                        unsigned const workingSampleCount
                            = deviceIndex + 1u == activeDevices ? samplesPerGpu + remainder : samplesPerGpu;

                        unsigned const beginOffset = deviceIndex * samplesPerGpu;

                        unsigned const minSampleIdx = compute.minSampleRange + beginOffset;
                        unsigned const maxSampleIdx = minSampleIdx + workingSampleCount - 1u;
                        hase::alpakaUtils::DevBundle devBundle{meshes[device_i].m_device, exec};
                        computes[device_i].gpu_i = compute.devices.at(device_i);
                        calcPhiAseThreaded(
                            devBundle,
                            experiment,
                            computes[device_i],
                            hostMesh,
                            meshes[device_i],
                            result,
                            minSampleIdx,
                            maxSampleIdx,
                            runtimes.at(device_i),
                            bar);
                    }

                    joinAll();
                    usedGPUs = activeDevices;
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
                    runtime = *std::ranges::max_element(runtimes);
                }

                else if(compute.parallelMode == ParallelMode::MPI)
                {
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
                    usedGPUs = hase::core::calcPhiAseMPI(
                        exec,
                        experiment,
                        compute,
                        hostMesh,
                        meshes,
                        result,
                        topology,
                        runtime);
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

                for(unsigned sample_i = compute.minSampleRange;
                    sample_i < meshes[0].numberOfSamples && sample_i <= compute.maxSampleRange;
                    ++sample_i)
                {
                    result.phiAse.at(sample_i) *= (hostMesh.nTot / hostMesh.crystalTFluo);
                    result.dndtAse.at(sample_i) = calcDndtAse(
                        hostMesh,
                        experiment.maxSigmaA,
                        experiment.maxSigmaE,
                        result.phiAse.at(sample_i),
                        sample_i);
                    if(sample_i <= 10)
                        dout(V_DEBUG) << " Dndt ASE[ " << sample_i << " ]: " << result.dndtAse.at(sample_i) << " "
                                      << result.mse.at(sample_i) << std::endl;
                }
                /***************************************************************************
                 * PRINT SOLUTION
                 **************************************************************************/
                if(verbosity & V_DEBUG)
                {
                    for(unsigned sample_i = 0; sample_i < meshes[0].numberOfSamples; ++sample_i)
                    {
                        dout(V_DEBUG) << "PHI ASE[" << sample_i << "]: " << result.phiAse.at(sample_i) << " "
                                      << result.mse.at(sample_i) << std::endl;
                        if(sample_i >= 10)
                            break;
                    }
                }

                if constexpr(MATLAB)
                {
                    /***************************************************************************
                     * WRITE MATLAB OUTPUT
                     **************************************************************************/
                    // output folder has to be the same as TMP_FOLDER in the calling MatLab script
                    hase::utils::writeMatlabOutput(
                        compute.outputPath,
                        result.phiAse,
                        result.totalRays,
                        result.mse,
                        meshes[0].numberOfSamples,
                        meshes[0].numberOfLevels);
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
                        experiment.minRaysPerSample,
                        experiment.maxRaysPerSample,
                        experiment.mseThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        tmpPhiAse,
                        vtkPath / fs::path("phiase_" + currentTime + ".vtk"),
                        experiment.minRaysPerSample,
                        experiment.maxRaysPerSample,
                        experiment.mseThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        result.mse,
                        vtkPath / fs::path("mse_" + currentTime + ".vtk"),
                        experiment.minRaysPerSample,
                        experiment.maxRaysPerSample,
                        experiment.mseThreshold,
                        experiment.useReflections,
                        runtime);

                    hase::utils::writePointsToVtk(
                        hostMesh,
                        tmpTotalRays,
                        vtkPath / fs::path("total_rays_" + currentTime + ".vtk"),
                        experiment.minRaysPerSample,
                        experiment.maxRaysPerSample,
                        experiment.mseThreshold,
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
                        sample_i <= compute.maxSampleRange && sample_i < result.mse.size();
                        ++sample_i)
                    {
                        auto& element = result.mse[sample_i];
                        maxMSE = std::max(maxMSE, element);
                        avgMSE += element;
                        if(element >= experiment.mseThreshold)
                        {
                            std::cout << " to high mse at for element: " << sample_i << " mse: " << element
                                      << std::endl;
                            highMSE++;
                        }
                    }
                    avgMSE /= numSamples;

                    try
                    {
                        std::cout.imbue(std::locale(""));
                    }
                    catch(std::runtime_error e)
                    {
                    }

                    dout(V_STAT | V_NOLABEL) << std::endl;
                    dout(V_STAT) << "=== Statistics ===" << std::endl;
                    dout(V_STAT) << "Backend       : " << compute.backend << std::endl;
                    dout(V_STAT) << "ParallelMode      : " << compute.parallelMode << std::endl;
                    dout(V_STAT) << "Prisms            : " << meshes[0].numberOfPrisms << std::endl;
                    dout(V_STAT) << "Samples           : "
                                 << std::min(static_cast<unsigned>(result.dndtAse.size()), numSamples) << std::endl;
                    dout(V_STAT) << "RaysPerSample     : " << experiment.minRaysPerSample;
                    if(experiment.maxRaysPerSample > experiment.minRaysPerSample)
                    {
                        dout(V_STAT | V_NOLABEL) << " - " << experiment.maxRaysPerSample << " (adaptive)";
                    }
                    dout(V_STAT | V_NOLABEL) << std::endl;
                    dout(V_STAT) << "sum(totalRays)    : "
                                 << std::accumulate(result.totalRays.begin(), result.totalRays.end(), 0.) << std::endl;
                    dout(V_STAT) << "MSE threshold     : " << experiment.mseThreshold << std::endl;
                    dout(V_STAT) << "int. Wavelength   : " << experiment.sigmaA.size() << std::endl;
                    dout(V_STAT) << "max. MSE          : " << maxMSE << std::endl;
                    dout(V_STAT) << "avg. MSE          : " << avgMSE << std::endl;
                    dout(V_STAT) << "too high MSE      : " << highMSE << " of " << numSamples << std::endl;

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
                    if(experiment.maxRaysPerSample > experiment.minRaysPerSample)
                    {
                        dout(V_STAT) << "=== Sampling resolution as Histogram ===" << std::endl;
                        hase::utils::ray_histogram(
                            result.totalRays,
                            experiment.maxRaysPerSample,
                            experiment.mseThreshold,
                            result.mse);
                    }
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
            std::cout << "\n------------------------------------------------------------------------------\n"
                      << std::endl;
        }
        return i || !oneDidRun;
    }

    int pythonEntry(ExperimentParameters experiment, ComputeParameters compute, Result& result, HostMesh host_mesh)
    {
        hase::parse::pythonParse(experiment, compute, host_mesh, result);
        return startSimulation<false>(experiment, compute, result, host_mesh);
    }

} // namespace hase::core
