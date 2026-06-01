/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
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


/**
 * @brief Calculates the gain sum for the given
 *        sample point with or without reflections. This is done by a monte carlo
 *        simulation with randomly generated rays and
 *        wavelenghts.
 *
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 **/

#pragma once

#include <alpaka/alpaka.hpp>
#include <alpaka/core/common.hpp>

#include <core/geometry.hpp> /* generateRay */
#include <core/mesh.hpp>
#include <kernels/propagateRay.hpp> /* propagateRay */

#include <cassert> /* assert */
#include <iomanip>

namespace hase::kernels
{

    /**
     * @brief get a random number from [0..length)
     *
     * @param length the maximum number to return (exclusive)
     * @param globalState State for random number generation (mersenne twister).
     *                    The state need to be initialized before. See
     *                    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MTGP/
     *                    for more information.
     *
     * @return a random number
     *
     */
    struct GenRndSigmas
    {
        constexpr ALPAKA_FN_HOST_ACC unsigned operator()(
            unsigned const length,
            alpaka::rand::engine::Philox4x32x10& engine)
        {
            return alpaka::rand::distribution::UniformReal{
                0.0f,
                static_cast<float>(length),
                alpaka::rand::interval::oc}(engine);
        }
    };

    struct CalcSampleGainSumWithReflection
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            alpaka::concepts::IMdSpan auto const indicesOfPrisms,
            alpaka::concepts::IMdSpan auto const numberOfReflectionSlices,
            alpaka::concepts::IMdSpan auto const importance,
            unsigned const raysPerSample,
            alpaka::concepts::IMdSpan auto dGainOfRay,
            unsigned const sample_i,
            alpaka::concepts::IMdSpan auto const sigmaA,
            alpaka::concepts::IMdSpan auto const sigmaE,
            unsigned const threadLocalStridingRNG) const
        {
            using Vec1 = alpaka::Vec<uint32_t, 1>;
            hase::core::Point samplePoint = mesh.getSamplePoint(sample_i);
            auto& blockOffset
                = declareSharedVar<unsigned[4], alpaka::uniqueId()>(acc); // 4 in case of warp-based raynumber

            auto const tIdx = hase::alpakaUtils::getLinGlobalIdx(acc);
            auto rndEngine = alpaka::rand::engine::Philox4x32x10{threadLocalStridingRNG + tIdx};
            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{Vec1{raysPerSample}}))
            {
                // Get triangle/prism to start ray from
                unsigned startPrism = indicesOfPrisms[rayNumber];
                unsigned reflection_i = numberOfReflectionSlices[rayNumber];
                unsigned reflections = (reflection_i + 1) / 2;
                hase::core::ReflectionPlane reflectionPlane
                    = (reflection_i % 2 == 0) ? hase::core::BOTTOM_REFLECTION : hase::core::TOP_REFLECTION;
                unsigned startLevel = startPrism / mesh.numberOfTriangles;
                unsigned startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
                unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;
                hase::core::Point startPoint = mesh.genRndPoint(acc, startTriangle, startLevel, rndEngine);

                // get a random index in the wavelength array
                unsigned sigma_i = GenRndSigmas{}(sigmaA.getExtents().product(), rndEngine);

                // Calculate reflections as different ray propagations
                double gain = propagateRayWithReflection(
                    startPoint,
                    samplePoint,
                    reflections,
                    reflectionPlane,
                    startLevel,
                    startTriangle,
                    mesh,
                    sigmaA[sigma_i],
                    sigmaE[sigma_i]);

                // include the stimulus from the starting prism and the importance of that ray
                gain *= mesh.getBetaVolume(startPrism) * importance[startPrism + reflectionOffset];

                assert(!alpaka::math::isnan(mesh.getBetaVolume(startPrism)));
                assert(!alpaka::math::isnan(importance[startPrism + reflectionOffset]));
                assert(!alpaka::math::isnan(gain));

                dGainOfRay[rayNumber] = gain;
            }
        }
    };

    struct CalcSampleGainSum
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            alpaka::concepts::IMdSpan auto const indicesOfPrisms,

            alpaka::concepts::IMdSpan auto const importance,
            unsigned const raysPerSample,
            alpaka::concepts::IMdSpan auto dGainOfRay,
            unsigned const sample_i,
            alpaka::concepts::IMdSpan auto const sigmaA,
            alpaka::concepts::IMdSpan auto const sigmaE,
            unsigned const lambdaResolution,
            unsigned const threadLocalStridingRNG) const
        {
            using Vec1 = alpaka::Vec<uint32_t, 1>;
            hase::core::Point samplePoint = mesh.getSamplePoint(sample_i);

            auto const tIdx = hase::alpakaUtils::getLinGlobalIdx(acc);
            auto rndEngine = alpaka::rand::engine::Philox4x32x10{threadLocalStridingRNG + tIdx};

            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{Vec1{raysPerSample}}))
            {
                // Get triangle/prism to start ray from
                unsigned startPrism = indicesOfPrisms[rayNumber];

                unsigned startLevel = startPrism / mesh.numberOfTriangles;
                unsigned startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
                hase::core::Point startPoint = mesh.genRndPoint(acc, startTriangle, startLevel, rndEngine);
                hase::core::Ray ray = hase::core::generateRay(startPoint, samplePoint);

                // Get a random index in the wavelength array
                unsigned sigma_i = GenRndSigmas{}(lambdaResolution, rndEngine);
                assert(sigma_i < lambdaResolution);

                // Calculate the gain for the whole ray at once
                double gain = propagateRay(ray, &startLevel, &startTriangle, mesh, sigmaA[sigma_i], sigmaE[sigma_i]);
                gain /= ray.length * ray.length; // important, since usually done in the reflection device function

                gain *= mesh.getBetaVolume(startPrism) * importance[startPrism];
                dGainOfRay[rayNumber] = gain;
            }
        }
    };

} // namespace hase::kernels
