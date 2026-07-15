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

#include <alpaka/math.hpp>

#include <core/mesh.hpp>

#include <algorithm>
#include <cmath>
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
        static inline std::string const relativeStandardError = "relative-standard-error-threshold";
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

    enum class SrmStatus
    {
        DISABLED,
        CONVERGED,
        STABLE,
        DIVERGED,
        MAX_ITERATIONS
    };

    [[nodiscard]] inline char const* toString(SrmStatus const status)
    {
        switch(status)
        {
        case SrmStatus::DISABLED:
            return "disabled";
        case SrmStatus::CONVERGED:
            return "converged";
        case SrmStatus::STABLE:
            return "stable";
        case SrmStatus::DIVERGED:
            return "diverged";
        case SrmStatus::MAX_ITERATIONS:
            return "max_iterations";
        }
        return "max_iterations";
    }

    struct Result
    {
        Result()
        {
        }

        Result(
            std::vector<float> phiAse,
            std::vector<double> standardError,
            std::vector<double> relativeStandardError,
            std::vector<unsigned> totalRays,
            std::vector<double> dndtAse,
            std::vector<unsigned> droppedRays = {},
            SrmStatus srmStatus = SrmStatus::DISABLED,
            unsigned srmPasses = 0u,
            double srmRemainingFraction = 0.0,
            unsigned srmMaxIterations = 0u,
            unsigned srmDivergenceStreak = 0u)
            : phiAse(std::move(phiAse))
            , standardError(std::move(standardError))
            , relativeStandardError(std::move(relativeStandardError))
            , totalRays(std::move(totalRays))
            , dndtAse(std::move(dndtAse))
            , droppedRays(std::move(droppedRays))
            , srmStatus(srmStatus)
            , srmPasses(srmPasses)
            , srmRemainingFraction(srmRemainingFraction)
            , srmMaxIterations(srmMaxIterations)
            , srmDivergenceStreak(srmDivergenceStreak)
        {
            if(this->droppedRays.empty())
                this->droppedRays.assign(this->phiAse.size(), 0u);
        }

        std::vector<float> phiAse;
        std::vector<double> standardError;
        std::vector<double> relativeStandardError;
        std::vector<unsigned> totalRays;
        std::vector<double> dndtAse;
        std::vector<unsigned> droppedRays;
        SrmStatus srmStatus = SrmStatus::DISABLED;
        unsigned srmPasses = 0u;
        double srmRemainingFraction = 0.0;
        unsigned srmMaxIterations = 0u;
        unsigned srmDivergenceStreak = 0u;
    };

    struct ExperimentParameters
    {
        ExperimentParameters()
        {
        }

        ExperimentParameters(
            unsigned minRays,
            unsigned maxRays,
            std::vector<double> lambdaA,
            std::vector<double> lambdaE,
            std::vector<double> sigmaA,
            std::vector<double> sigmaE,
            double maxSigmaA,
            double maxSigmaE,
            double relativeStandardErrorThreshold,
            bool useReflections,
            unsigned spectral,
            bool monochromatic = false)
            : minRays(minRays)
            , maxRays(maxRays)
            , lambdaA(std::move(lambdaA))
            , lambdaE(std::move(lambdaE))
            , sigmaA(std::move(sigmaA))
            , sigmaE(std::move(sigmaE))
            , maxSigmaA(maxSigmaA)
            , maxSigmaE(maxSigmaE)
            , relativeStandardErrorThreshold(relativeStandardErrorThreshold)
            , useReflections(useReflections)
            , monochromatic(monochromatic)
            , spectral(spectral)
        {
        }

        ExperimentParameters(
            unsigned minRays,
            unsigned maxRays,
            std::vector<double> sigmaA,
            std::vector<double> sigmaE,
            double maxSigmaA,
            double maxSigmaE,
            double relativeStandardErrorThreshold,
            bool useReflections,
            bool monochromatic = false)
            : minRays(minRays)
            , maxRays(maxRays)
            , sigmaA(std::move(sigmaA))
            , sigmaE(std::move(sigmaE))
            , maxSigmaA(maxSigmaA)
            , maxSigmaE(maxSigmaE)
            , relativeStandardErrorThreshold(relativeStandardErrorThreshold)
            , useReflections(useReflections)
            , monochromatic(monochromatic)
        {
        }

        [[nodiscard]] bool isForwardPropagation() const
        {
            return propagationMode == "forward";
        }

        [[nodiscard]] unsigned resolvedForwardRayCount() const
        {
            return forwardRayCount == 0u ? minRays : forwardRayCount;
        }

        unsigned minRays;
        unsigned maxRays;
        unsigned forwardRayCount = 0u;
        std::string propagationMode = "forward";
        std::vector<double> lambdaA;
        std::vector<double> lambdaE;
        std::vector<double> sigmaA;
        std::vector<double> sigmaE;
        double maxSigmaA;
        double maxSigmaE;
        double relativeStandardErrorThreshold;
        bool useReflections;
        bool monochromatic = false;
        unsigned spectral;
        unsigned reflectionMaxIterations = 8u;
        double reflectionTolerance = 1.0e-4;
        unsigned surfaceReservoirSize = 32u;
    };

    [[nodiscard]] inline unsigned adaptiveRayTarget(
        ExperimentParameters const& experiment,
        ComputeParameters const& compute,
        unsigned const completedIncreases)
    {
        if(experiment.forwardRayCount != 0u || experiment.maxRays <= experiment.minRays || compute.adaptiveSteps == 0u)
        {
            return experiment.resolvedForwardRayCount();
        }
        if(completedIncreases >= compute.adaptiveSteps)
        {
            return experiment.maxRays;
        }

        double const growth = std::pow(
            static_cast<double>(experiment.maxRays) / static_cast<double>(experiment.minRays),
            1.0 / static_cast<double>(compute.adaptiveSteps));
        double const target
            = static_cast<double>(experiment.minRays) * std::pow(growth, static_cast<double>(completedIncreases));
        unsigned const rounded = static_cast<unsigned>(std::ceil(target));
        return std::clamp(rounded, experiment.minRays, experiment.maxRays);
    }

    [[nodiscard]] inline bool forwardResultMeetsRelativeStandardError(
        Result const& result,
        double const relativeStandardErrorThreshold)
    {
        return !result.relativeStandardError.empty()
               && std::all_of(
                   result.relativeStandardError.cbegin(),
                   result.relativeStandardError.cend(),
                   [relativeStandardErrorThreshold](double const relativeStandardError)
                   {
                       return alpaka::math::isfinite(relativeStandardError)
                              && relativeStandardError <= relativeStandardErrorThreshold;
                   });
    }

    inline void recordAdaptiveRayConvergence(
        Result const& result,
        unsigned const targetRayCount,
        double const relativeStandardErrorThreshold,
        std::vector<unsigned>& convergenceRayCounts)
    {
        if(convergenceRayCounts.empty())
        {
            convergenceRayCounts.assign(result.relativeStandardError.size(), 0u);
        }

        unsigned const volumeCount = std::min(
            static_cast<unsigned>(convergenceRayCounts.size()),
            static_cast<unsigned>(result.relativeStandardError.size()));
        for(unsigned volume = 0u; volume < volumeCount; ++volume)
        {
            bool const hasDroppedRays = !result.droppedRays.empty() && result.droppedRays.at(volume) != 0u;
            double const relativeStandardError = result.relativeStandardError.at(volume);
            if(convergenceRayCounts.at(volume) == 0u && !hasDroppedRays
               && alpaka::math::isfinite(relativeStandardError)
               && relativeStandardError <= relativeStandardErrorThreshold)
            {
                convergenceRayCounts.at(volume) = targetRayCount;
            }
        }
    }


} // namespace hase::core
