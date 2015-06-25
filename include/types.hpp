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

#include <boost/filesystem.hpp> /* fs::path */
namespace fs = boost::filesystem;


enum DeviceMode { NO_DEVICE_MODE, GPU_DEVICE_MODE, CPU_DEVICE_MODE};
enum ParallelMode { NO_PARALLEL_MODE, THREADED_PARALLEL_MODE, MPI_PARALLEL_MODE, GRAYBAT_PARALLEL_MODE};

std::ostream& operator<<(std::ostream& out, const DeviceMode value);
std::ostream& operator<<(std::ostream& out, const ParallelMode value);

struct ComputeParameters {

    ComputeParameters() {}
    
    ComputeParameters(  unsigned maxRepetitions,
		        unsigned gpu_i,
			DeviceMode deviceMode,
			ParallelMode parallelMode,
			bool writeVtk,
			fs::path inputPath,
			fs::path outputPath,
			std::vector<unsigned> devices,
			int minSampleRange,
			int maxSampleRange) :
	maxRepetitions(maxRepetitions),
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
    unsigned gpu_i;
    DeviceMode deviceMode;
    ParallelMode parallelMode;
    bool writeVtk;
    fs::path inputPath;
    fs::path outputPath;
    std::vector<unsigned> devices;
    int minSampleRange;
    int maxSampleRange;

    
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

    
