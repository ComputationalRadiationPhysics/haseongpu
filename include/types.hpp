/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert
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

#pragma once

#include <boost/filesystem/path.hpp> /* fs::path */

#include <string>
#include <vector>

namespace fs = boost::filesystem;

struct DeviceMode
{
    static inline std::string const NONE = "no_device_mode";
    static inline std::string const GPU = "gpu";
    static inline std::string const CPU = "cpu";
};

struct ParallelMode
{
    static inline std::string const NONE = "no_parallel_mode";
    static inline std::string const THREADED = "threaded";
    static inline std::string const MPI = "mpi";
    static inline std::string const GRAYBAT = "graybat";
};

struct CompSwitch
{
    static inline std::string const parallel_mode = "parallel-mode";
    static inline std::string const device_mode = "device-mode";
    static inline std::string const ngpus = "ngpus";
    static inline std::string const repetitions = "repetitions";
    static inline std::string const adaptive_steps = "adaptive-steps";
    static inline std::string const min_sample_i = "min-sample-i";
    static inline std::string const max_sample_i = "max-sample-i";
    static inline std::string const write_vtk = "write-vtk";
};

struct ExpSwitch
{
    static inline std::string const input_path = "input-path";
    static inline std::string const output_path = "output-path";
    static inline std::string const min_rays = "min-rays";
    static inline std::string const max_rays = "max-rays";
    static inline std::string const mse = "mse-threshold";
    static inline std::string const reflection = "reflection";
    static inline std::string const spectral = "spectral-resolution";
};

struct ComputeParameters
{
    ComputeParameters()
    {
    }

    ComputeParameters(
        unsigned maxRepetitions,
        unsigned adaptiveSteps,
        unsigned gpu_i,
        std::string deviceMode,
        std::string parallelMode,
        bool writeVtk,
        fs::path inputPath,
        fs::path outputPath,
        std::vector<unsigned> devices,
        unsigned minSampleRange,
        unsigned maxSampleRange,
        unsigned maxGpus)
        : maxRepetitions(maxRepetitions)
        , adaptiveSteps(adaptiveSteps)
        , maxGpus(maxGpus)
        , gpu_i(gpu_i)
        , deviceMode(deviceMode)
        , parallelMode(parallelMode)
        , writeVtk(writeVtk)
        , inputPath(inputPath)
        , outputPath(outputPath)
        , devices(devices)
        , minSampleRange(minSampleRange)
        , maxSampleRange(maxSampleRange)
    {
    }

    unsigned maxRepetitions;
    unsigned adaptiveSteps;
    //user defined nr of gpus
    unsigned maxGpus;
    unsigned gpu_i;
    std::string deviceMode;
    std::string parallelMode;
    bool writeVtk;
    fs::path inputPath;
    fs::path outputPath;
    //gpu ids from cuda api
    std::vector<unsigned> devices;
    unsigned minSampleRange;
    unsigned maxSampleRange;
};

struct Result
{
    Result()
    {
    }

    Result(
        std::vector<float> phiAse,
        std::vector<double> mse,
        std::vector<unsigned> totalRays,
        std::vector<double> dndtAse)
        : phiAse(phiAse)
        , mse(mse)
        , totalRays(totalRays)
        , dndtAse(dndtAse)
    {
    }

    std::vector<float> phiAse;
    std::vector<double> mse;
    std::vector<unsigned> totalRays;
    std::vector<double> dndtAse;
};

struct ExperimentParameters
{
    ExperimentParameters()
    {
    }

    ExperimentParameters(
        unsigned minRaysPerSample,
        unsigned maxRaysPerSample,
        std::vector<double> sigmaA,
        std::vector<double> sigmaE,
        double maxSigmaA,
        double maxSigmaE,
        double mseThreshold,
        bool useReflections)
        : minRaysPerSample(minRaysPerSample)
        , maxRaysPerSample(maxRaysPerSample)
        , sigmaA(sigmaA)
        , sigmaE(sigmaE)
        , maxSigmaA(maxSigmaA)
        , maxSigmaE(maxSigmaE)
        , mseThreshold(mseThreshold)
        , useReflections(useReflections)
    {
    }

    unsigned minRaysPerSample;
    unsigned maxRaysPerSample;
    std::vector<double> lambdaA;
    std::vector<double> lambdaE;
    std::vector<double> sigmaA;
    std::vector<double> sigmaE;
    double maxSigmaA;
    double maxSigmaE;
    double mseThreshold;
    bool useReflections;
    unsigned spectral;
};
