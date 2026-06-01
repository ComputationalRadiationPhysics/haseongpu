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


#include <mapRaysToPrisms.hpp>

/**
 * @brief takes a prefix sum (obtained by an exclusive scan over raysPerPrism)
 *        and writes it into a unary representation. The value from the
 *        prefixSum at index i describes the offset where to start writing,
 *        whereas i itself is the new value to be stored in the output array.
 *
 *        example:
 *        raysPerPrism [3,0,2,1]
 *
 *        -> exclusive prefixSum [0,3,3,5]
 *
 *        beginning from place 0 in the output should be 0 (length 3 according to raysPerPrism[0])
 *        beginning from place 3 in the output should be 1 (EMPTY range at raysPerPrism[1])
 *        beginning from place 3 in the output should be 2 (length 2 according to raysPerPrism[2])
 *        beginning from place 5 in the output should be 3 (length 1 according to raysPerPrism[3])
 *
 *        resulting output arrays:
 *        [0 0 0 2 2 3] (indicesOfPrisms)
 *
 *        output numberOfReflections is handled in a similar way
 *
 *
 * @param numberOfPrisms the number of prisms. numberOfPrisms * reflectionSlices
 *                       must be equal to the the length of the prefixSum.
 * @param raysPerSample the size of indicesOfPrisms/numberOfReflections. Actually
 *                      identical to the sum of all values in raysPerPrism
 * @param reflectionSlices the number of reflectionSlices. see numberOfPrisms
 * @param raysPerPrism the input array from which prefixSum was generated
 * @param prefixSum the prefixSum generated from raysPerPrism
 * @param indicesOfPrisms a pointer to the OUTPUT generated like described in the example
 * @param numberOfReflections a pointer to the OUTPUT similar to indicesOfPrisms
 */
struct MapPrefixSumToPrisms
{
    //@TODO [performance] use alpaka scan or thrust scan instead
    ALPAKA_FN_ACC void operator()(
        auto const& acc,
        unsigned const numberOfPrisms,
        unsigned const raysPerSample,
        unsigned const reflectionSlices,
        alpaka::concepts::IMdSpan auto const raysPerPrism,
        alpaka::concepts::IMdSpan auto const prefixSum,
        alpaka::concepts::IMdSpan auto indicesOfPrisms,
        alpaka::concepts::IMdSpan auto numberOfReflections) const
    {
        for(auto [id] : alpaka::onAcc::makeIdxMap(
                acc,
                alpaka::onAcc::worker::threadsInGrid,
                alpaka::IdxRange{numberOfPrisms * reflectionSlices}))
        {
            unsigned const count = raysPerPrism[id];
            unsigned const startingPosition = prefixSum[id];
            unsigned const reflection_i = id / numberOfPrisms;
            unsigned const prism_i = id % numberOfPrisms;

            for(unsigned i = 0; i < count; ++i)
            {
                indicesOfPrisms[startingPosition + i] = prism_i;
                numberOfReflections[startingPosition + i] = reflection_i;
            }
        }
    }
};

#include <cassert>
#include <vector>

void mapRaysToPrisms(
    auto& devBundle,
    alpaka::concepts::IMdSpan auto& indicesOfPrisms,
    alpaka::concepts::IMdSpan auto& numberOfReflections,
    alpaka::concepts::IMdSpan auto const& raysPerPrism,
    alpaka::concepts::IMdSpan auto& prefixSum,
    unsigned const reflectionSlices,
    unsigned const raysPerSample,
    unsigned const numberOfPrisms)
{
    auto queue = devBundle.device.makeQueue();
    auto frameSpec = alpaka::onHost::getFrameSpec<Vec1D::type>(
        devBundle.device,
        alpaka::Vec{raysPerPrism.getExtents().product()});

    alpaka::onHost::exclusiveScan(queue, devBundle.executor, prefixSum, raysPerPrism);
    std::vector<unsigned> hRays(raysPerPrism.getExtents().product());
    std::vector<unsigned> hPrefix(prefixSum.getExtents().product());
    alpaka::onHost::memcpy(queue, hRays, raysPerPrism);
    alpaka::onHost::memcpy(queue, hPrefix, prefixSum);
    alpaka::onHost::wait(queue);
    queue.enqueue(
        devBundle.executor,
        frameSpec,
        alpaka::KernelBundle{
            MapPrefixSumToPrisms{},
            numberOfPrisms,
            raysPerSample,
            reflectionSlices,
            raysPerPrism,
            prefixSum,
            indicesOfPrisms,
            numberOfReflections});
    alpaka::onHost::wait(queue);
}
