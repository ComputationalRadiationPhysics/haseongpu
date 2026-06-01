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


#include <kernels/mapRaysToPrisms.hpp>

namespace hase::kernels
{

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
    template<alpaka::concepts::IMdSpan T1, alpaka::concepts::IMdSpan T2>
    struct MapPrefixSumToPrisms
    {
        uint32_t numberOfPrisms;
        mutable T1 indicesOfPrisms;
        mutable T2 numberOfReflections;

        //@TODO [performance] use alpaka scan or thrust scan instead
        ALPAKA_FN_HOST_ACC void operator()(
            alpaka::concepts::SimdPtr auto const raysPerPrism,
            alpaka::concepts::SimdPtr auto const prefixSum) const
        {
            alpaka::concepts::Vector auto packOffset = raysPerPrism.getIdx();
            alpaka::concepts::Simd auto const count = raysPerPrism.load();
            alpaka::concepts::Simd auto const startingPosition = prefixSum.load();
            constexpr uint32_t width = ALPAKA_TYPEOF(raysPerPrism)::width();
            for(uint32_t laneIdx = 0; laneIdx < width; ++laneIdx)
            {
                auto [elemIdx] = (packOffset + laneIdx);
                uint32_t const prism_i = elemIdx % numberOfPrisms;
                uint32_t const reflection_i = elemIdx / numberOfPrisms;
                for(unsigned i = 0; i < count[laneIdx]; ++i)
                {
                    indicesOfPrisms[startingPosition[laneIdx] + i] = prism_i;
                    numberOfReflections[startingPosition[laneIdx] + i] = reflection_i;
                }
            }
        }
    };

} // namespace hase::kernels

#include <cassert>
#include <vector>

namespace hase::kernels
{
    using hase::alpakaUtils::Vec1D;

    void mapRaysToPrisms(
        auto& devBundle,
        alpaka::concepts::IBuffer auto& indicesOfPrisms,
        alpaka::concepts::IBuffer auto& numberOfReflections,
        alpaka::concepts::IBuffer auto& raysPerPrism,
        alpaka::concepts::IBuffer auto& prefixSum,
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

        alpaka::concepts::IMdSpan auto indicesSpan = indicesOfPrisms.getMdSpan();
        alpaka::concepts::IMdSpan auto reflectionSpan = numberOfReflections.getMdSpan();
        alpaka::onHost::concurrent<typename ALPAKA_TYPEOF(raysPerPrism)::value_type>(
            queue,
            devBundle.executor,
            prefixSum.getExtents(),
            hase::kernels::MapPrefixSumToPrisms<ALPAKA_TYPEOF(indicesSpan), ALPAKA_TYPEOF(reflectionSpan)>{
                numberOfPrisms,
                indicesSpan,
                reflectionSpan},
            raysPerPrism,
            prefixSum);
        alpaka::onHost::wait(queue);
    }

} // namespace hase::kernels
