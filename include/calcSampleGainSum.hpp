/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
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
#include <geometry.hpp> /* generateRay */
#include <mesh.hpp>
#include <propagateRay.hpp> /* propagateRay */

#include <cassert> /* assert */
#include <iomanip>

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
    constexpr ALPAKA_FN_ACC unsigned operator()(unsigned const length, alpaka::rand::engine::Philox4x32x10& engine)
    {
        return alpaka::rand::distribution::UniformReal{0.0f, static_cast<float>(length), alpaka::rand::interval::oc}(
            engine);
    }
};

struct CalcSampleGainSumWithReflection
{
    ALPAKA_FN_ACC void operator()(
        auto const& acc,
        DeviceMeshView const mesh,
        unsigned const* indicesOfPrisms,
        unsigned const* numberOfReflectionSlices,
        double const* importance,
        unsigned const raysPerSample,
        float* gainSum,
        float* gainSumSquare,
        unsigned const sample_i,
        double const* sigmaA,
        double const* sigmaE,
        unsigned const maxInterpolation,
        unsigned const threadLocalStridingRNG) const
    {
        using Vec1 = alpaka::Vec<uint32_t, 1>;
        double gainSumTemp = 0.0;
        double gainSumSquareTemp = 0.0;
        Point samplePoint = mesh.getSamplePoint(sample_i);
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
            ReflectionPlane reflectionPlane = (reflection_i % 2 == 0) ? BOTTOM_REFLECTION : TOP_REFLECTION;
            unsigned startLevel = startPrism / mesh.numberOfTriangles;
            unsigned startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
            unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;
            Point startPoint = mesh.genRndPoint(acc, startTriangle, startLevel, rndEngine);

            // get a random index in the wavelength array
            unsigned sigma_i = GenRndSigmas{}(maxInterpolation, rndEngine);

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

            gainSumTemp += gain;
            gainSumSquareTemp += gain * gain;
        }

        alpaka::onAcc::atomicAdd(acc, &(gainSum[0]), float(gainSumTemp));
        alpaka::onAcc::atomicAdd(acc, &(gainSumSquare[0]), float(gainSumSquareTemp));
    }
};

struct CalcSampleGainSum
{
    ALPAKA_FN_ACC void operator()(
        auto const& acc,
        DeviceMeshView const mesh,
        unsigned const* indicesOfPrisms,
        double const* importance,
        unsigned const raysPerSample,
        float* gainSum,
        float* gainSumSquare,
        unsigned const sample_i,
        double const* sigmaA,
        double const* sigmaE,
        unsigned const lambdaResolution,
        unsigned const threadLocalStridingRNG) const
    {
        using Vec1 = alpaka::Vec<uint32_t, 1>;
        double gainSumTemp = 0.0;
        double gainSumSquareTemp = 0.0;
        Point samplePoint = mesh.getSamplePoint(sample_i);

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
            Point startPoint = mesh.genRndPoint(acc, startTriangle, startLevel, rndEngine);
            Ray ray = generateRay(startPoint, samplePoint);

            // Get a random index in the wavelength array
            unsigned sigma_i = GenRndSigmas{}(lambdaResolution, rndEngine);
            assert(sigma_i < lambdaResolution);

            // Calculate the gain for the whole ray at once
            double gain = propagateRay(ray, &startLevel, &startTriangle, mesh, sigmaA[sigma_i], sigmaE[sigma_i]);
            gain /= ray.length * ray.length; // important, since usually done in the reflection device function

            gain *= mesh.getBetaVolume(startPrism) * importance[startPrism];


            gainSumTemp += gain;
            gainSumSquareTemp += gain * gain;
        }
        alpaka::onAcc::atomicAdd(acc, &(gainSum[0]), float(gainSumTemp));
        alpaka::onAcc::atomicAdd(acc, &(gainSumSquare[0]), float(gainSumSquareTemp));
    }
};
