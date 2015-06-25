#pragma once

#include <boost/filesystem.hpp> /* fs::path */
namespace fs = boost::filesystem;


enum DeviceMode { NO_DEVICE_MODE, GPU_DEVICE_MODE, CPU_DEVICE_MODE};
enum ParallelMode { NO_PARALLEL_MODE, THREADED_PARALLEL_MODE, MPI_PARALLEL_MODE, GRAYBAT_PARALLEL_MODE};

struct ComputeParameters {

    ComputeParameters() {}
    
    ComputeParameters(  unsigned maxRepetitions,
            unsigned adaptiveSteps,
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

    
