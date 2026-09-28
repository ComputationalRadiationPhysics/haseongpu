#include <core/domainSchedule.hpp>

#include <ranges>
#include <unordered_map>

namespace hase::core
{
    namespace
    {
        struct PendingAssignment
        {
            DomainRayPopulationId id;
            std::uint64_t rayCount{};
            double work{};
        };

        [[nodiscard]] double estimatedWork(DomainCost const& domain, std::uint64_t const rayCount)
        {
            double const trace
                = static_cast<double>(rayCount) * std::max(1.0, std::cbrt(static_cast<double>(domain.cellCount)));
            double const boundary
                = static_cast<double>(std::min<std::uint64_t>(2u * rayCount, domain.boundaryFaceCount + rayCount));
            return trace + boundary;
        }
    } // namespace

    DomainSchedule makeDomainSchedule(
        std::vector<WorkerDescriptor> const& workers,
        std::vector<DomainCost> const& domains,
        std::vector<DomainQuota> const& quotas,
        std::uint32_t const numIndependentRayPopulations,
        std::span<data::AseDomainInterface const> const interfaces)
    {
        if(workers.empty() || domains.empty() || numIndependentRayPopulations == 0u)
            throw std::invalid_argument("domain scheduling requires workers, domains, and batches");
        std::unordered_map<data::DomainId, DomainQuota const*> quotaById;
        for(auto const& quota : quotas)
            if(!quotaById.emplace(quota.id, &quota).second)
                throw std::invalid_argument("domain quota ids must be unique");
        std::vector<DomainWorkItem> work;
        work.reserve(domains.size() * numIndependentRayPopulations);
        for(auto const& domain : domains)
        {
            auto const found = quotaById.find(domain.id);
            if(found == quotaById.end())
                throw std::invalid_argument("each scheduled domain requires a ray quota");
            auto const total = found->second->rayCount;
            for(std::uint32_t batch = 0u; batch < numIndependentRayPopulations; ++batch)
            {
                auto const rays
                    = total / numIndependentRayPopulations + (batch < total % numIndependentRayPopulations ? 1u : 0u);
                work.push_back({{domain.id, batch, 0u}, rays, 0u});
            }
        }
        return makeDomainSchedule(workers, domains, work, interfaces);
    }

    DomainSchedule makeDomainSchedule(
        std::vector<WorkerDescriptor> const& workers,
        std::vector<DomainCost> const& domains,
        std::span<DomainWorkItem const> const work,
        std::span<data::AseDomainInterface const> const interfaces)
    {
        if(workers.empty() || domains.empty())
            throw std::invalid_argument("domain scheduling requires workers and domains");
        for(auto const& worker : workers)
            if(worker.relativeCapacity <= 0.0 || !std::isfinite(worker.relativeCapacity))
                throw std::invalid_argument("worker relative capacity must be finite and positive");

        std::unordered_map<data::DomainId, DomainCost const*> domainById;
        for(auto const& domain : domains)
            if(!domainById.emplace(domain.id, &domain).second)
                throw std::invalid_argument("domain ids must be unique");
        std::vector<PendingAssignment> pending;
        pending.reserve(work.size());
        for(auto const& item : work)
        {
            auto const found = domainById.find(item.id.domain);
            if(found == domainById.end())
                throw std::invalid_argument("scheduled work references an unknown domain");
            pending.push_back({item.id, item.rayCount, estimatedWork(*found->second, item.rayCount)});
        }
        std::ranges::sort(pending, {}, &PendingAssignment::id);
        if(std::ranges::adjacent_find(pending, {}, &PendingAssignment::id) != pending.end())
            throw std::invalid_argument("domain execution batch ids must be unique");
        std::ranges::sort(
            pending,
            [](auto const& left, auto const& right)
            {
                if(left.work != right.work)
                    return left.work > right.work;
                return left.id < right.id;
            });

        std::vector<double> workerLoad(workers.size(), 0.0);
        std::vector<std::uint64_t> workerMemory(workers.size(), 0u);
        std::vector<std::unordered_set<data::DomainId>> residentDomains(workers.size());
        std::unordered_map<WorkerId, std::size_t> workerIndex;
        for(std::size_t index = 0u; index < workers.size(); ++index)
            if(!workerIndex.emplace(workers[index].id, index).second)
                throw std::invalid_argument("worker ids must be unique");

        std::vector<DomainAssignment> assignments;
        assignments.reserve(pending.size());
        for(auto const& item : pending)
        {
            auto const& domain = *domainById.at(item.id.domain);
            double minimumProjectedLoad = std::numeric_limits<double>::max();
            std::vector<std::size_t> feasible;
            for(std::size_t index = 0u; index < workers.size(); ++index)
            {
                bool const newReplica = !residentDomains[index].contains(item.id.domain);
                auto const projectedMemory
                    = workerMemory[index] + (newReplica ? domain.residentBytes : std::uint64_t{0u});
                if(projectedMemory > workers[index].deviceMemoryBytes)
                    continue;
                feasible.push_back(index);
                minimumProjectedLoad = std::min(
                    minimumProjectedLoad,
                    (workerLoad[index] + item.work) / workers[index].relativeCapacity);
            }
            if(feasible.empty())
                throw std::runtime_error("no worker has enough device memory for a domain batch");

            std::optional<std::size_t> selected;
            auto selectedScore = std::tuple{
                std::numeric_limits<std::uint64_t>::max(),
                std::numeric_limits<std::uint64_t>::max(),
                std::numeric_limits<WorkerId>::max()};
            double const loadLimit = minimumProjectedLoad == 0.0 ? 0.0 : minimumProjectedLoad * 1.05;
            for(auto const index : feasible)
            {
                double const projectedLoad = (workerLoad[index] + item.work) / workers[index].relativeCapacity;
                if(projectedLoad > loadLimit)
                    continue;
                std::uint64_t crossNodeEdges = 0u;
                for(auto const& interface : interfaces)
                {
                    data::DomainId other{};
                    if(interface.sourceDomain == item.id.domain)
                        other = interface.targetDomain;
                    else if(interface.targetDomain == item.id.domain)
                        other = interface.sourceDomain;
                    else
                        continue;
                    bool colocated = false;
                    bool placed = false;
                    for(std::size_t otherWorker = 0u; otherWorker < workers.size(); ++otherWorker)
                        if(residentDomains[otherWorker].contains(other))
                        {
                            placed = true;
                            colocated |= workers[otherWorker].node == workers[index].node;
                        }
                    if(placed && !colocated)
                        ++crossNodeEdges;
                }
                bool const newReplica = !residentDomains[index].contains(item.id.domain);
                auto const score = std::tuple{
                    crossNodeEdges,
                    newReplica ? domain.residentBytes : std::uint64_t{0u},
                    workers[index].id};
                if(!selected || score < selectedScore)
                {
                    selected = index;
                    selectedScore = score;
                }
            }
            if(!selected)
                throw std::runtime_error("domain scheduler failed to select a feasible worker");
            auto const index = *selected;
            if(residentDomains[index].insert(item.id.domain).second)
                workerMemory[index] += domain.residentBytes;
            workerLoad[index] += item.work;
            assignments.push_back({item.id, workers[index].id, item.rayCount, item.work});
        }

        std::ranges::sort(assignments, {}, &DomainAssignment::id);
        return DomainSchedule{std::move(assignments)};
    }

    std::vector<DomainCost> makeDomainCosts(data::AseDomainGraph const& graph)
    {
        auto bytes = []<typename T>(std::vector<T> const& values)
        { return static_cast<std::uint64_t>(values.size()) * sizeof(T); };
        std::vector<DomainCost> result;
        result.reserve(graph.domains.size());
        for(auto const& domain : graph.domains)
        {
            auto const& trace = domain.trace;
            std::uint64_t residentBytes = 0u;
            residentBytes += bytes(trace.points) + bytes(trace.betaVolume) + bytes(trace.cellMaterialIds);
            residentBytes += bytes(trace.materialActive) + bytes(trace.materialRefractiveIndices);
            residentBytes += bytes(trace.materialActiveIonDensities) + bytes(trace.materialFluorescenceLifetimes);
            residentBytes += bytes(trace.materialBulkAttenuations) + bytes(trace.materialPeakAbsorption);
            residentBytes += bytes(trace.materialPeakEmission) + bytes(trace.materialCrossSectionOffsets);
            residentBytes += bytes(trace.crossSectionWavelengths) + bytes(trace.crossSectionAbsorption);
            residentBytes += bytes(trace.crossSectionEmission) + bytes(trace.surfaceReflectivities);
            residentBytes += bytes(trace.surfaceRefractiveIndexInside) + bytes(trace.surfaceRefractiveIndexOutside);
            residentBytes += bytes(trace.cellPointIndices) + bytes(trace.cellTypes) + bytes(trace.cellFaces);
            residentBytes += bytes(trace.barycentricFacePlanes) + bytes(trace.cellNeighborCells);
            residentBytes += bytes(trace.cellNeighborLocalFaces) + bytes(trace.cellFaceBoundaries);
            residentBytes += bytes(trace.cellVolumes) + bytes(trace.lumpedMaterialVertexVolumes);
            residentBytes += bytes(trace.cellVolumePrefix) + bytes(trace.sourceStrengthPrefix);
            residentBytes += bytes(trace.cellCenters) + bytes(trace.samplePoints);
            residentBytes += bytes(domain.boundaryTargetDomains) + bytes(domain.boundaryTargetCells)
                             + bytes(domain.boundaryTargetFaces);
            residentBytes += bytes(domain.boundaryReflectivities) + bytes(domain.boundarySourceRefractiveIndices)
                             + bytes(domain.boundaryTargetRefractiveIndices);
            auto const boundaryFaceCount = static_cast<std::uint64_t>(
                std::ranges::count_if(trace.cellFaceBoundaries, [](int const boundary) { return boundary > 0; }));
            double const sourceStrength = trace.sourceStrengthPrefix.empty() ? 0.0 : trace.sourceStrengthPrefix.back();
            result.push_back(
                {domain.id,
                 trace.numberOfCells,
                 residentBytes,
                 boundaryFaceCount,
                 sourceStrength,
                 domain.requestedRays});
        }
        return result;
    }
} // namespace hase::core
