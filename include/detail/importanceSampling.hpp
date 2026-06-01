/**
 * Copyright 2013-2026 Erik Zenker, Carlchristian Eckert, Marius Melzer, Tim Hanel
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
#include <geometry.hpp>
#include <mesh.hpp>
#include <propagateRay.hpp>

#include <cassert>
#include <cstdio>
#include <iostream>
#include <ostream>

constexpr ALPAKA_FN_ACC bool deviceIsFinite(double x)
{
    return alpaka::math::isfinite(x);
}

constexpr ALPAKA_FN_ACC inline bool deviceIsFinite(float x)
{
    return alpaka::math::isfinite(x);
}

ALPAKA_FN_ACC void assertMeshPropagationInputs(DeviceMeshView const& mesh, unsigned prism, unsigned sample_i)
{
    assert(prism < mesh.numberOfPrisms);
    assert(sample_i < mesh.numberOfSamples);
    unsigned const level = prism / mesh.numberOfTriangles;
    [[maybe_unused]] unsigned const triangle = prism - mesh.numberOfTriangles * level;

    assert(level < mesh.numberOfLevels);
    assert(triangle < mesh.numberOfTriangles);

    [[maybe_unused]] double const beta = mesh.getBetaVolume(prism);
    assert(isfinite(beta));
}

struct ValidateMeshKernel
{
    ALPAKA_FN_ACC void operator()(auto const& acc, DeviceMeshView const& mesh, unsigned sample_i) const
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
 * @brief calculates a first estimate on the importance of each prism, based on a single ray started in the center of
 * each prism
 *
 * @param *importance will contain the initial importance for each prism
 *
 * @param *sumPhi will contain the cumulative sum of the importance values
 *
 * For other parameters, see documentation of importanceSampling()
 *
 */
struct PropagateFromTriangleCenter
{
    ALPAKA_FN_ACC void operator()(
        auto const& acc,
        DeviceMeshView const mesh,
        alpaka::concepts::IMdSpan auto importance,
        unsigned reflectionsSlices,
        unsigned raysPerSample,
        unsigned const sample_i,
        double const sigmaA,
        double const sigmaE) const
    {
        auto numFrames = acc[alpaka::frame::count];
        auto frameExtent = acc[alpaka::frame::extent];
        for(auto [reflection_i, _, startPrism] : alpaka::onAcc::makeIdxMap(
                acc,
                alpaka::onAcc::worker::threadsInGrid,
                alpaka::IdxRange{Vec3D{reflectionsSlices, 1, raysPerSample}}))
        {
            double gain = 0;
            unsigned reflections = (reflection_i + 1) / 2;
            ReflectionPlane reflectionPlane = (reflection_i % 2 == 0) ? BOTTOM_REFLECTION : TOP_REFLECTION;
            unsigned startLevel = startPrism / (mesh.numberOfTriangles);
            unsigned startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
            Point startPoint = mesh.getCenterPoint(startTriangle, startLevel);
            Point samplePoint = mesh.getSamplePoint(sample_i);
            unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;

            gain = propagateRayWithReflection(
                startPoint,
                samplePoint,
                reflections,
                reflectionPlane,
                startLevel,
                startTriangle,
                mesh,
                sigmaA,
                sigmaE);
            importance[startPrism + reflectionOffset] = mesh.getBetaVolume(startPrism) * gain;
            if(mesh.getBetaVolume(startPrism) < 0 || gain < 0 || importance[startPrism + reflectionOffset] < 0)
            {
                printf(
                    "beta: %f importance: %f gain: %f\n",
                    mesh.getBetaVolume(startPrism),
                    importance[startPrism + reflectionOffset],
                    gain);
            }
        }
    }
};

struct DistributeRaysByImportance
{
    ALPAKA_FN_ACC void operator()(
        auto const& acc,
        DeviceMeshView mesh,
        alpaka::concepts::IMdSpan auto raysPerPrism,
        alpaka::concepts::IMdSpan auto importance,

        double const sumPhi,
        unsigned reflectionsSlices,
        unsigned raysPerSample) const
    {
        for(auto [reflection_i, _, startPrism] : alpaka::onAcc::makeIdxMap(
                acc,
                alpaka::onAcc::worker::threadsInGrid,
                alpaka::IdxRange{Vec3D{reflectionsSlices, 1, mesh.numberOfPrisms}}))
        {
            unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;

            double const raysForPrism
                = static_cast<double>(raysPerSample) * importance[startPrism + reflectionOffset] / sumPhi;
            raysPerPrism[startPrism + reflectionOffset] = static_cast<unsigned>(alpaka::math::floor(raysForPrism));
            if(raysPerPrism[startPrism + reflectionOffset] > raysPerSample) // this condition is wrong
            {
                printf(
                    "importance: %f sumPhi: %f raysPerPrism[%d]: %d (max %d)\n",
                    importance[startPrism + reflectionOffset],
                    sumPhi,
                    startPrism + reflectionOffset,
                    raysPerPrism[startPrism + reflectionOffset],
                    raysPerSample);
            }
        }
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
    ALPAKA_FN_ACC void operator()(
        auto const& acc,
        DeviceMeshView mesh,
        alpaka::concepts::IMdSpan auto raysPerPrism,
        unsigned raysLeft,
        alpaka::concepts::IMdSpan auto raysDump,
        auto const threadLocalStridingRNG) const
    {
        auto const tIdx = hase::alpakaUtils::getLinGlobalIdx(acc);
        alpaka::rand::engine::Philox4x32x10 engine(threadLocalStridingRNG + tIdx);
        for(auto [id] : alpaka::onAcc::makeIdxMap(
                acc,
                alpaka::onAcc::worker::threadsInGrid,
                alpaka::IdxRange{alpaka::Vec{raysLeft}}))
        {
            // get random numbers floating for any vector size -> conversion float -> int does floor automatically
            unsigned rand_t = alpaka::rand::distribution::UniformReal{
                0.0f,
                static_cast<float>(mesh.numberOfTriangles),
                alpaka::rand::interval::co}(engine);
            unsigned rand_z = alpaka::rand::distribution::UniformReal{
                0.0f,
                static_cast<float>(mesh.numberOfLevels - 1),
                alpaka::rand::interval::co}(engine);
            unsigned randomPrism = rand_t + rand_z * mesh.numberOfTriangles;
            alpaka::onAcc::atomicAdd(acc, &raysPerPrism[randomPrism], static_cast<unsigned>(1));
            alpaka::onAcc::atomicAdd(acc, &raysDump[0], static_cast<unsigned>(1));
        }
    }
};

struct RecalculateImportance
{
    /**
     * @brief corrects the importance to match with the randomly distributed rays
     *
     * @param *raysPerPrism the number of rays to be launced for each prism
     *
     * @param *importance the importance for each prism (will be changed)
     *
     * for other parameters, see documentation of importanceSampling()
     */
    ALPAKA_FN_ACC void operator()(
        auto const& acc,
        DeviceMeshView mesh,
        alpaka::concepts::IMdSpan auto raysPerPrism,
        alpaka::concepts::IMdSpan auto dRaysDump,
        alpaka::concepts::IMdSpan auto importance,
        unsigned reflectionSlices) const
    {
        unsigned numberOfPrisms = mesh.numberOfPrisms;
        for(auto [reflection_i, y, startPrism] : alpaka::onAcc::makeIdxMap(
                acc,
                alpaka::onAcc::worker::threadsInGrid,
                alpaka::IdxRange{Vec3D{reflectionSlices, 1, numberOfPrisms}}))
        {
            unsigned reflectionOffset = reflection_i * numberOfPrisms;

            if(startPrism >= numberOfPrisms)
            {
                return;
            }
            int startLevel = startPrism / (mesh.numberOfTriangles);
            int startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);

            if(raysPerPrism[startPrism + reflectionOffset] > 0)
            {
                auto val = dRaysDump[0] * mesh.triangleSurfaces[startTriangle]
                           / (mesh.surfaceTotal * raysPerPrism[startPrism + reflectionOffset]);

                importance[startPrism + reflectionOffset] = val;
            }
            else
            {
                importance[startPrism + reflectionOffset] = 0;
            }
        }
    }
};
