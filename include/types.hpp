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
    static const std::string NONE;
    static const std::string GPU;
    static const std::string CPU;
};


struct ParallelMode {
    static const std::string NONE;
    static const std::string THREADED;
#if defined(MPI_FOUND)
    static const std::string MPI;
#endif
#if defined(BOOST_MPI_FOUND) || defined(ZMQ_FOUND)
    static const std::string GRAYBAT;
#endif
};

struct CompSwitch{
    static const std::string parallel_mode;
    static const std::string device_mode;
    static const std::string ngpus;
    static const std::string repetitions;
    static const std::string adaptive_steps;
    static const std::string min_sample_i;
    static const std::string max_sample_i;
};

struct ExpSwitch{
    static const std::string input_path;
    static const std::string output_path;
    static const std::string min_rays;
    static const std::string max_rays;
    static const std::string mse;
    static const std::string reflection;
    static const std::string spectral;
};

struct ComputeParameters {

    ComputeParameters() {}

    ComputeParameters(  unsigned maxRepetitions,
			unsigned adaptiveSteps,
			std::string deviceMode,
			std::string parallelMode,
			bool writeVtk,
			fs::path inputPath,
			fs::path outputPath,
			unsigned minSampleRange,
			unsigned maxSampleRange) :
        maxRepetitions(maxRepetitions),
        adaptiveSteps(adaptiveSteps),
        deviceMode(deviceMode),
        parallelMode(parallelMode),
        writeVtk(writeVtk),
        inputPath(inputPath),
        outputPath(outputPath),
        minSampleRange(minSampleRange),
        maxSampleRange(maxSampleRange){ }

    unsigned maxRepetitions;
    unsigned adaptiveSteps;
    //unsigned gpu_i;
    std::string deviceMode;
    std::string parallelMode;
    bool writeVtk;
    fs::path inputPath;
    fs::path outputPath;
    //std::vector<unsigned> devices;
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
			   bool useReflections,
			   unsigned numberOfSamples,
			   unsigned numberOfLevels,
			   unsigned numberOfPrisms) :
        minRaysPerSample(minRaysPerSample),
        maxRaysPerSample(maxRaysPerSample),
        sigmaA(sigmaA),
        sigmaE(sigmaE),
        maxSigmaA(maxSigmaA),
        maxSigmaE(maxSigmaE),
        mseThreshold(mseThreshold),
        useReflections(useReflections),
	numberOfSamples(numberOfSamples),
	numberOfLevels(numberOfLevels),
	numberOfPrisms(numberOfPrisms){ }

    unsigned minRaysPerSample;
    unsigned maxRaysPerSample;
    std::vector<double> sigmaA;
    std::vector<double> sigmaE;
    double maxSigmaA;
    double maxSigmaE;
    double mseThreshold;
    bool useReflections;
    unsigned numberOfSamples;
    unsigned numberOfLevels;
    unsigned numberOfPrisms;

};
