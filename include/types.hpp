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


struct DeviceMode {
    static const inline std::string NONE  = "no_device_mode";
    static const inline std::string GPU   = "gpu";
    static const inline std::string CPU   = "cpu";
};


struct ParallelMode {
    static const inline std::string NONE        = "no_parallel_mode";
    static const inline std::string THREADED    = "threaded";
    static const inline std::string MPI         = "mpi";
    static const inline std::string GRAYBAT     = "graybat";

};

struct CompSwitch{
    static const inline std::string parallel_mode     = "parallel-mode";
    static const inline std::string device_mode       = "device-mode";
    static const inline std::string ngpus             = "ngpus";
    static const inline std::string repetitions       = "repetitions";
    static const inline std::string adaptive_steps    = "adaptive-steps";
    static const inline std::string min_sample_i      = "min-sample-i";
    static const inline std::string max_sample_i      = "max-sample-i";
    static const inline std::string write_vtk         = "write-vtk";
};

struct ExpSwitch{
    static const inline std::string input_path     = "input-path";
    static const inline std::string output_path    = "output-path";
    static const inline std::string min_rays       = "min-rays";
    static const inline std::string max_rays       = "max-rays";
    static const inline std::string mse            = "mse-threshold";
    static const inline std::string reflection     = "reflection";
    static const inline std::string spectral       = "spectral-resolution";
};

struct ComputeParameters {

    ComputeParameters() {}

    ComputeParameters(  unsigned maxRepetitions,
            unsigned adaptiveSteps,
            unsigned gpu_i,
            std::string deviceMode,
            std::string parallelMode,
            bool writeVtk,
            fs::path inputPath,
            fs::path outputPath,
            std::vector<unsigned> devices,
            unsigned minSampleRange,
            unsigned maxSampleRange) :
        maxRepetitions(maxRepetitions),
        adaptiveSteps(adaptiveSteps),
        gpu_i(gpu_i),
        deviceMode(deviceMode),
        parallelMode(parallelMode),
        writeVtk(writeVtk),
        inputPath(inputPath),
        outputPath(outputPath),
        devices(devices),
        minSampleRange(minSampleRange),
        maxSampleRange(maxSampleRange){ }

    unsigned maxRepetitions;
    unsigned adaptiveSteps;
    unsigned gpu_i;
    std::string deviceMode;
    std::string parallelMode;
    bool writeVtk;
    fs::path inputPath;
    fs::path outputPath;
    std::vector<unsigned> devices;
    unsigned minSampleRange;
    unsigned maxSampleRange;


};

struct Result {

    Result(){}

    Result( std::vector<float> phiAse,
            std::vector<double> mse,
            std::vector<unsigned> totalRays,
            std::vector<double>   dndtAse) :
        phiAse(phiAse),
        mse(mse),
        totalRays(totalRays),
        dndtAse(dndtAse){}

    std::vector<float> phiAse;
    std::vector<double> mse;
    std::vector<unsigned> totalRays;
    std::vector<double>   dndtAse;

};

struct ExperimentParameters {

    ExperimentParameters() {}

    ExperimentParameters(  unsigned minRaysPerSample,
            unsigned maxRaysPerSample,
            std::vector<double> sigmaA,
            std::vector<double> sigmaE,
            double maxSigmaA,
            double maxSigmaE,
            double mseThreshold,
            bool useReflections) :
        minRaysPerSample(minRaysPerSample),
        maxRaysPerSample(maxRaysPerSample),
        sigmaA(sigmaA),
        sigmaE(sigmaE),
        maxSigmaA(maxSigmaA),
        maxSigmaE(maxSigmaE),
        mseThreshold(mseThreshold),
        useReflections(useReflections) { }

    unsigned minRaysPerSample;
    unsigned maxRaysPerSample;
    std::vector<double> sigmaA;
    std::vector<double> sigmaE;
    double maxSigmaA;
    double maxSigmaE;
    double mseThreshold;
    bool useReflections;

};
