/**
 * Copyright 2013-2026 Erik Zenker, Carlchristian Eckert, Marius Melzer, Tim Hanel
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU : it contains device function definitions corresponding to importance sampling
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
#include <alpakaUtils/utils.hpp>
#include <core/geometry.hpp>
#include <core/mesh.hpp>
#include <kernels/propagateRay.hpp>

#include <cassert>
#include <cstdio>
#include <iostream>
#include <ostream>

namespace hase::kernels
{

    inline ALPAKA_FN_HOST_ACC void assertMeshPropagationInputs(
        core::DeviceMeshView const& mesh,
        unsigned cell,
        unsigned sample_i)
    {
        assert(cell < mesh.numberOfCells);
        if(sample_i >= mesh.numberOfSamples)
        {
            assert(false);
        }

        [[maybe_unused]] double const beta = mesh.getBetaVolume(cell);
        assert(alpaka::math::isfinite(beta));
    }

    struct ValidateMeshKernel
    {
        ALPAKA_FN_HOST_ACC void operator()(auto const& acc, core::DeviceMeshView const& mesh, unsigned sample_i) const
        {
            for(auto [prism] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfPrisms}))
            {
                assertMeshPropagationInputs(mesh, prism, sample_i);
            }
        }
    };

    /**
     * @brief calculates a first estimate on the importance of each prism, based on a single ray started in the center
     * of each prism
     *
     * @param *importance will contain the initial importance for each prism
     *
     * @param *sumPhi will contain the cumulative sum of the importance values
     *
     * For other parameters, see documentation of importanceSampling()
     *
     */
    struct PropagateFromCellCenter
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            core::DeviceMeshView mesh,
            unsigned const sample_i,
            double const sigmaA,
            double const sigmaE,
            alpaka::concepts::IMdSpan auto importance,
            alpaka::concepts::IMdSpan auto droppedRays,
            alpaka::concepts::IMdSpan auto infiniteRaySnapshots) const
        {
            core::Point samplePoint = mesh.getSamplePoint(sample_i);

            for(auto [prismIndex] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{importance.getExtents().product()}))
            {
                unsigned const reflection_i = prismIndex / mesh.numberOfPrisms;
                unsigned const startCell = prismIndex % mesh.numberOfPrisms;
                unsigned const reflections = (reflection_i + 1u) / 2u;
                core::Point startPoint = mesh.getCellCenterPoint(startCell);
                core::ReflectionPlane const reflectionPlane
                    = (reflection_i % 2u == 0u) ? core::BOTTOM_REFLECTION : core::TOP_REFLECTION;
                double gain = propagateRayWithReflection(
                    startPoint,
                    samplePoint,
                    reflections,
                    reflectionPlane,
                    startCell,
                    mesh,
                    sigmaA,
                    sigmaE);
                if(!alpaka::math::isfinite(gain))
                {
                    if(droppedRays[0] == 0u)
                    {
                        core::Ray ray = core::generateRay(startPoint, samplePoint);
                        infiniteRaySnapshots[0] = core::InfiniteRaySnapshot{
                            startPoint,
                            samplePoint,
                            ray.dir,
                            0.0,
                            ray.length,
                            startCell,
                            gain};
                    }
                    alpaka::onAcc::atomicAdd(acc, &droppedRays[0], 1u);
                    gain = 0.0;
                }
                importance[prismIndex] = mesh.getBetaVolume(startCell) * gain;
            }
        }
    };

    struct DistributeRaysByImportance
    {
        double const sumPhi;
        unsigned raysPerSample;

        ALPAKA_FN_HOST_ACC auto operator()(alpaka::concepts::Simd auto const importance) const
        {
            alpaka::concepts::Simd auto const raysForPrism = static_cast<double>(raysPerSample) * importance / sumPhi;
            return alpaka::pCast<unsigned>(raysForPrism);
        }
    };

    /**
     * @brief takes a number of rays and distributes them randomly over the available prisms
     *        Warning: Does not distribute to reflection slices !!!
     *
     * @param *raysPerPrism the number of rays for each prism (will be changed)
     * @param *raysDump the number of rays which were already distributed
     *
     * for other parameters, see documentation of importanceSampling()
     *
     */
    struct DistributeRemaingRaysRandomly
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            core::DeviceMeshView mesh,
            alpaka::concepts::IMdSpan auto raysPerPrism,
            unsigned raysLeft,
            alpaka::concepts::IMdSpan auto raysDump,
            auto const threadLocalStridingRNG) const
        {
            auto const tIdx = hase::alpakaUtils::getLinGlobalIdx(acc);
            alpaka::rand::engine::Philox4x32x10 engine(threadLocalStridingRNG + tIdx);
            for(auto id :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{raysLeft}))
            {
                unsigned randomPrism = alpaka::rand::distribution::UniformReal{
                    0.0f,
                    static_cast<float>(mesh.numberOfPrisms),
                    alpaka::rand::interval::co}(engine);
                alpaka::onAcc::atomicAdd(acc, &raysPerPrism[randomPrism], static_cast<unsigned>(1));
                alpaka::onAcc::atomicAdd(acc, &raysDump[0], static_cast<unsigned>(1));
            }
        }
    };

    struct DrawRaysByImportance
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            alpaka::concepts::IMdSpan auto cdfInclusive,
            alpaka::concepts::IMdSpan auto raysPerPrism,
            double sumPhi,
            unsigned numberOfWeightedPrisms,
            unsigned raysPerSample,
            auto const threadLocalStridingRNG) const
        {
            auto const tIdx = hase::alpakaUtils::getLinGlobalIdx(acc);
            alpaka::rand::engine::Philox4x32x10 engine(threadLocalStridingRNG + tIdx);

            for(auto ray :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{raysPerSample}))
            {
                if(ray >= raysPerSample)
                {
                    continue;
                }
                double const draw
                    = alpaka::rand::distribution::UniformReal<double>{0.0, sumPhi, alpaka::rand::interval::co}(engine);

                unsigned lower = 0u;
                unsigned upper = numberOfWeightedPrisms;
                while(lower < upper)
                {
                    unsigned const mid = lower + ((upper - lower) / 2u);
                    double const binEnd = cdfInclusive[mid];
                    if(draw < binEnd)
                    {
                        upper = mid;
                    }
                    else
                    {
                        lower = mid + 1u;
                    }
                }

                if(numberOfWeightedPrisms > 0u)
                {
                    if(lower >= numberOfWeightedPrisms)
                    {
                        lower = numberOfWeightedPrisms - 1u;
                    }
                    alpaka::onAcc::atomicAdd(acc, &raysPerPrism[lower], 1u);
                }
            }
        }
    };

    struct RecalculateUnbiasedImportance
    {
        mutable core::DeviceMeshView mesh;
        double sumPhi;

        ALPAKA_FN_HOST_ACC void operator()(
            alpaka::concepts::SimdPtr auto preImportance,
            alpaka::concepts::SimdPtr auto importance) const
        {
            alpaka::concepts::Vector auto packOffset = preImportance.getIdx();
            constexpr uint32_t width = ALPAKA_TYPEOF(preImportance)::width();
            alpaka::concepts::Simd auto weightedPrism = alpaka::Simd<unsigned, width>::fill(packOffset[0]);
            for(uint32_t laneIdx = 0; laneIdx < width; ++laneIdx)
            {
                weightedPrism[laneIdx] += laneIdx;
            }
            alpaka::concepts::Simd auto const prism = weightedPrism % mesh.numberOfPrisms;

            alpaka::concepts::Simd auto const preImportanceValues = preImportance.load();
            alpaka::concepts::Simd auto prismVolume = alpaka::Simd<double, width>::fill(0.0);
            // in order to avoid division by zero
            alpaka::concepts::Simd auto safePreImportanceValues = alpaka::Simd<double, width>::fill(1.0);

            for(uint32_t laneIdx = 0; laneIdx < width; ++laneIdx)
            {
                if(preImportanceValues[laneIdx] > 0.0)
                {
                    prismVolume[laneIdx] = mesh.getCellVolume(prism[laneIdx]);
                    safePreImportanceValues[laneIdx] = preImportanceValues[laneIdx];
                }
            }
            alpaka::concepts::Simd auto const samplingProbability = safePreImportanceValues / sumPhi;
            importance = prismVolume / samplingProbability;
        }
    };

} // namespace hase::kernels
