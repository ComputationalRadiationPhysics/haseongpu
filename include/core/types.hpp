/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert
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

#pragma once

#include <core/mesh.hpp>

#include <filesystem>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace hase::core
{
    namespace fs = std::filesystem;

    struct Backend
    {
        static inline std::string const NONE = "no_device_mode";
        static inline std::string const GPU = "gpu";
        static inline std::string const CPU = "cpu";
    };

    struct ParallelMode
    {
        static inline std::string const NONE = "no_parallel_mode";
        static inline std::string const SINGLE = "single";
        static inline std::string const MPI = "mpi";
    };

    struct CompSwitch
    {
        static inline std::string const parallel_mode = "parallel-mode";
        static inline std::string const backend = "backend";
        static inline std::string const numDevices = "numDevices";
        static inline std::string const repetitions = "repetitions";
        static inline std::string const adaptive_steps = "adaptive-steps";
        static inline std::string const min_sample_i = "min-sample-i";
        static inline std::string const max_sample_i = "max-sample-i";
        static inline std::string const write_vtk = "write-vtk";
        static inline std::string const rng_seed = "rng-seed";
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
        static inline std::string const monochromatic = "monochromatic";
    };

    struct ComputeParameters
    {
        static constexpr unsigned unspecifiedRngSeed = std::numeric_limits<unsigned>::max();

        ComputeParameters()
        {
        }

        ComputeParameters(
            unsigned maxRepetitions,
            unsigned adaptiveSteps,
            unsigned numDevices,
            unsigned gpu_i,
            std::string backend,
            std::string parallelMode,
            bool writeVtk,
            std::vector<unsigned> devices,
            unsigned minSampleRange,
            unsigned maxSampleRange,
            unsigned rngSeed = unspecifiedRngSeed)
            : maxRepetitions(maxRepetitions)
            , adaptiveSteps(adaptiveSteps)
            , numDevices(numDevices)
            , gpu_i(gpu_i)
            , backend(std::move(backend))
            , parallelMode(std::move(parallelMode))
            , writeVtk(writeVtk)
            , devices(std::move(devices))
            , minSampleRange(minSampleRange)
            , maxSampleRange(maxSampleRange)
            , rngSeed(rngSeed)
        {
        }

        ComputeParameters(
            unsigned maxRepetitions,
            unsigned adaptiveSteps,
            unsigned gpu_i,
            std::string backend,
            std::string parallelMode,
            bool writeVtk,
            fs::path inputPath,
            fs::path outputPath,
            std::vector<unsigned> devices,
            unsigned minSampleRange,
            unsigned maxSampleRange,
            unsigned numDevices,
            unsigned rngSeed = unspecifiedRngSeed)
            : maxRepetitions(maxRepetitions)
            , adaptiveSteps(adaptiveSteps)
            , numDevices(numDevices)
            , gpu_i(gpu_i)
            , backend(std::move(backend))
            , parallelMode(std::move(parallelMode))
            , writeVtk(writeVtk)
            , inputPath(std::move(inputPath))
            , outputPath(std::move(outputPath))
            , devices(std::move(devices))
            , minSampleRange(minSampleRange)
            , maxSampleRange(maxSampleRange)
            , rngSeed(rngSeed)
        {
        }

        unsigned maxRepetitions;
        unsigned adaptiveSteps;
        // user defined nr of gpus
        unsigned numDevices;
        unsigned gpu_i;
        std::string backend;
        std::string parallelMode;
        bool writeVtk;
        fs::path inputPath;
        fs::path outputPath;
        // gpu ids from cuda api
        std::vector<unsigned> devices;
        unsigned minSampleRange;
        unsigned maxSampleRange;
        unsigned rngSeed = unspecifiedRngSeed;
    };

    struct RuntimeTopology
    {
        unsigned activeNodes = 1;
        unsigned activeRanks = 1;
        unsigned activeGpus = 0;
        double avgActiveRanksPerNode = 1.0;
        unsigned minActiveRanksPerNode = 1;
        unsigned maxActiveRanksPerNode = 1;
        double avgGpusPerRank = 0.0;
        double avgGpusPerNode = 0.0;
        unsigned minGpusPerNode = 0;
        unsigned maxGpusPerNode = 0;
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
            std::vector<double> dndtAse,
            std::vector<unsigned> droppedRays = {})
            : phiAse(std::move(phiAse))
            , mse(std::move(mse))
            , totalRays(std::move(totalRays))
            , dndtAse(std::move(dndtAse))
            , droppedRays(std::move(droppedRays))
        {
            if(this->droppedRays.empty())
                this->droppedRays.assign(this->phiAse.size(), 0u);
        }

        std::vector<float> phiAse;
        std::vector<double> mse;
        std::vector<unsigned> totalRays;
        std::vector<double> dndtAse;
        std::vector<unsigned> droppedRays;
    };

    struct ExperimentParameters
    {
        ExperimentParameters()
        {
        }

        ExperimentParameters(
            unsigned minRaysPerSample,
            unsigned maxRaysPerSample,
            std::vector<double> lambdaA,
            std::vector<double> lambdaE,
            std::vector<double> sigmaA,
            std::vector<double> sigmaE,
            double maxSigmaA,
            double maxSigmaE,
            double mseThreshold,
            bool useReflections,
            unsigned spectral,
            bool monochromatic = false)
            : minRaysPerSample(minRaysPerSample)
            , maxRaysPerSample(maxRaysPerSample)
            , lambdaA(std::move(lambdaA))
            , lambdaE(std::move(lambdaE))
            , sigmaA(std::move(sigmaA))
            , sigmaE(std::move(sigmaE))
            , maxSigmaA(maxSigmaA)
            , maxSigmaE(maxSigmaE)
            , mseThreshold(mseThreshold)
            , useReflections(useReflections)
            , monochromatic(monochromatic)
            , spectral(spectral)
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
            bool useReflections,
            bool monochromatic = false)
            : minRaysPerSample(minRaysPerSample)
            , maxRaysPerSample(maxRaysPerSample)
            , sigmaA(std::move(sigmaA))
            , sigmaE(std::move(sigmaE))
            , maxSigmaA(maxSigmaA)
            , maxSigmaE(maxSigmaE)
            , mseThreshold(mseThreshold)
            , useReflections(useReflections)
            , monochromatic(monochromatic)
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
        bool monochromatic = false;
        unsigned spectral;
    };


} // namespace hase::core
