#pragma once

#include <core/forwardSamplingPlan.hpp>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace hase::core
{
    // Statistical subdivision, not a worker-dependent launch-tuning parameter.
    inline constexpr std::uint32_t maxLogicalSrmBatchRays = 65536u;

    [[nodiscard]] inline std::uint32_t forwardPopulationDomainCount(std::size_t const count)
    {
        if(count > std::numeric_limits<std::uint32_t>::max())
            throw std::overflow_error("forward domain count exceeds the transport index range");
        return static_cast<std::uint32_t>(std::max(std::size_t{1u}, count));
    }

    /** Execution metadata; batchId and ownership never enter the sampling RNG. */
    struct ForwardPopulationBatch
    {
        std::uint32_t rayPopulationId{};
        data::DomainId domainId{};
        std::uint32_t batchId{};
        std::uint32_t rayOffset{}; //!< Index within this source domain's population.
        std::uint32_t rayCount{};
        std::uint32_t populationOffset{}; //!< First ray of this domain in the complete population.
        std::uint32_t domainRayCount{};
        double sourceWeight{};
    };

    /** Fixed populations are resolved before subdividing work for device owners. */
    [[nodiscard]] inline std::vector<ForwardPopulationBatch> makeForwardPopulationBatches(
        std::span<std::uint32_t const> const domainCounts,
        std::span<DomainQuota const> const quotas,
        double const totalSource,
        std::uint32_t const populationCount,
        std::uint32_t const chunkSize)
    {
        if(populationCount == 0u || chunkSize == 0u || domainCounts.size() != quotas.size())
            throw std::invalid_argument("invalid forward execution chunk layout");
        std::vector<ForwardPopulationBatch> result;
        for(std::uint32_t rayPopulationId = 0u; rayPopulationId < populationCount; ++rayPopulationId)
        {
            std::uint64_t total = 0u;
            for(auto const count : domainCounts)
                total += domainPopulationRayCount(count, rayPopulationId, populationCount);
            if(total > std::numeric_limits<std::uint32_t>::max())
                throw std::overflow_error("forward population exceeds the ray index range");
            std::uint32_t offset = 0u;
            for(std::uint32_t domain = 0u; domain < quotas.size(); ++domain)
            {
                auto const count = domainPopulationRayCount(domainCounts[domain], rayPopulationId, populationCount);
                if(quotas[domain].sourceStrength > 0.0 && count == 0u)
                    throw std::invalid_argument("every emitting domain requires rays in every independent population");
                auto const weight = domainPopulationSourceWeight(
                    quotas[domain].sourceStrength,
                    totalSource,
                    count,
                    static_cast<std::uint32_t>(total));
                std::uint32_t batchId = 0u;
                for(std::uint32_t begin = 0u; begin < count;)
                {
                    auto const size = std::min(chunkSize, count - begin);
                    result.push_back(
                        {.rayPopulationId = rayPopulationId,
                         .domainId = domain,
                         .batchId = batchId++,
                         .rayOffset = begin,
                         .rayCount = size,
                         .populationOffset = offset,
                         .domainRayCount = count,
                         .sourceWeight = weight});
                    begin += size;
                }
                offset += count;
            }
        }
        return result;
    }

    /** Enough independent work for scheduling, with a bounded device launch size. */
    [[nodiscard]] inline std::uint32_t forwardExecutionChunkSize(std::uint32_t const rays, std::uint32_t const workers)
    {
        if(workers == 0u)
            throw std::invalid_argument("forward execution requires a worker");
        auto const partitions = static_cast<std::uint64_t>(workers) * 4u;
        return static_cast<std::uint32_t>(
            std::clamp<std::uint64_t>((static_cast<std::uint64_t>(rays) + partitions - 1u) / partitions, 1u, 65536u));
    }

    /** Map the concrete execution chunks while retaining their independent sampling identities. */
    [[nodiscard]] inline DomainSchedule makeForwardPopulationSchedule(
        std::vector<WorkerDescriptor> const& workers,
        std::vector<DomainCost> const& domains,
        std::span<ForwardPopulationBatch const> const batches,
        std::span<data::AseDomainInterface const> const interfaces)
    {
        std::vector<DomainWorkItem> work;
        work.reserve(batches.size());
        for(auto const& batch : batches)
            work.push_back({{batch.domainId, batch.rayPopulationId, batch.batchId}, batch.rayCount, 0u});
        return makeDomainSchedule(workers, domains, work, interfaces);
    }
} // namespace hase::core
