#pragma once

#include <data/AseDomainGraph.hpp>

#include <algorithm>
#include <cmath>
#include <compare>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <ranges>
#include <span>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

namespace hase::core
{
    using WorkerId = std::uint32_t;
    using NodeId = std::uint32_t;

    /** @brief Stable identity of one indivisible domain-local statistical batch. */
    struct DomainRayPopulationId
    {
        data::DomainId domain{};
        data::RayPopulationId rayPopulationId{};
        std::uint32_t batchId{};

        constexpr auto operator<=>(DomainRayPopulationId const&) const = default;
    };

    /** @brief Domain-local primary launch metadata retained by boundary histories. */
    struct DomainWorkItem
    {
        DomainRayPopulationId id;
        std::uint64_t rayCount{};
        std::uint32_t rngSeed{};
    };

    /** @brief Location, capacity, and memory constraints of one device-owning worker. */
    struct WorkerDescriptor
    {
        WorkerId id{};
        NodeId node{};
        std::uint32_t localWorker{};
        double relativeCapacity{1.0};
        std::uint64_t deviceMemoryBytes{std::numeric_limits<std::uint64_t>::max()};
    };

    struct DomainCost
    {
        data::DomainId id{};
        std::uint64_t cellCount{};
        std::uint64_t residentBytes{};
        std::uint64_t boundaryFaceCount{};
        double sourceStrength{};
        std::optional<std::uint64_t> requestedRays;
    };

    struct DomainQuota
    {
        data::DomainId id{};
        std::uint64_t rayCount{};
        double sourceStrength{};
        double importanceWeight{};
    };

    struct DomainAssignment
    {
        DomainRayPopulationId id;
        WorkerId worker{};
        std::uint64_t rayCount{};
        double estimatedWork{};
    };

    /** @brief Static worker ownership for all domain/batch pairs. */
    class DomainSchedule
    {
    public:
        DomainSchedule() = default;

        explicit DomainSchedule(std::vector<DomainAssignment> assignments) : m_assignments(std::move(assignments))
        {
        }

        [[nodiscard]] std::vector<DomainAssignment> const& assignments() const
        {
            return m_assignments;
        }

        [[nodiscard]] WorkerId owner(DomainRayPopulationId const id) const
        {
            auto const found = std::ranges::find(m_assignments, id, &DomainAssignment::id);
            if(found == m_assignments.end())
                throw std::out_of_range("domain batch is absent from the schedule");
            return found->worker;
        }

        [[nodiscard]] std::vector<data::DomainId> requiredDomains(WorkerId const worker) const
        {
            std::vector<data::DomainId> result;
            for(auto const& assignment : m_assignments)
                if(assignment.worker == worker && std::ranges::find(result, assignment.id.domain) == result.end())
                    result.push_back(assignment.id.domain);
            std::ranges::sort(result);
            return result;
        }

    private:
        std::vector<DomainAssignment> m_assignments;
    };

    /**
     * @brief Resolve exact final domain quotas and their unbiased importance weights.
     * @throws std::invalid_argument for impossible or biased allocations.
     */
    [[nodiscard]] inline std::vector<DomainQuota> allocateDomainRays(
        std::vector<DomainCost> const& domains,
        std::uint64_t const totalRayCount)
    {
        if(domains.empty() || totalRayCount == 0u)
            throw std::invalid_argument("domain ray allocation requires domains and a positive total");

        std::vector<DomainQuota> result;
        result.reserve(domains.size());
        std::vector<double> allocationStrength(domains.size(), 0.0);
        std::uint64_t reserved = 0u;
        double automaticStrength = 0.0;
        std::size_t automaticSources = 0u;
        std::size_t automaticDomains = 0u;
        for(std::size_t index = 0u; index < domains.size(); ++index)
        {
            auto const& domain = domains[index];
            if(!std::isfinite(domain.sourceStrength) || domain.sourceStrength < 0.0)
                throw std::invalid_argument("domain source strength must be finite and non-negative");
            if(domain.requestedRays)
            {
                if(*domain.requestedRays == 0u)
                    throw std::invalid_argument("manual domain ASE rays must be positive");
                if(*domain.requestedRays > totalRayCount - std::min(totalRayCount, reserved))
                    throw std::invalid_argument("manual domain ASE rays exceed the global total");
                reserved += *domain.requestedRays;
            }
            else
            {
                ++automaticDomains;
                if(domain.sourceStrength > 0.0)
                {
                    allocationStrength[index] = domain.sourceStrength;
                    automaticStrength += domain.sourceStrength;
                    ++automaticSources;
                }
            }
            result.push_back({domain.id, domain.requestedRays.value_or(0u), domain.sourceStrength, 0.0});
        }
        if(reserved > totalRayCount)
            throw std::invalid_argument("manual domain ASE rays exceed the global total");
        auto const remainder = totalRayCount - reserved;
        if(remainder != 0u && automaticDomains == 0u)
            throw std::invalid_argument("no automatic domain can receive the remaining ASE rays");
        if(remainder != 0u && automaticSources == 0u)
        {
            for(std::size_t index = 0u; index < domains.size(); ++index)
                if(!domains[index].requestedRays)
                {
                    allocationStrength[index]
                        = static_cast<double>(std::max<std::uint64_t>(1u, domains[index].cellCount));
                    automaticStrength += allocationStrength[index];
                    ++automaticSources;
                }
        }
        if(remainder != 0u)
        {
            if(remainder < automaticSources)
                throw std::invalid_argument("global ASE ray count cannot sample every source domain");

            struct Fraction
            {
                std::size_t index;
                double remainder;
            };

            std::vector<Fraction> fractions;
            std::uint64_t assigned = 0u;
            for(std::size_t index = 0u; index < domains.size(); ++index)
            {
                if(domains[index].requestedRays || allocationStrength[index] == 0.0)
                    continue;
                double const exact = static_cast<double>(remainder) * allocationStrength[index] / automaticStrength;
                auto const count = std::max<std::uint64_t>(1u, static_cast<std::uint64_t>(std::floor(exact)));
                result[index].rayCount = count;
                assigned += count;
                fractions.push_back({index, exact - std::floor(exact)});
            }
            if(assigned > remainder)
            {
                std::ranges::sort(
                    fractions,
                    [](auto const& left, auto const& right)
                    {
                        if(left.remainder != right.remainder)
                            return left.remainder < right.remainder;
                        return left.index > right.index;
                    });
                for(auto const& fraction : fractions)
                    if(assigned > remainder && result[fraction.index].rayCount > 1u)
                    {
                        auto const excess = std::min(assigned - remainder, result[fraction.index].rayCount - 1u);
                        result[fraction.index].rayCount -= excess;
                        assigned -= excess;
                    }
            }
            else
            {
                std::ranges::sort(
                    fractions,
                    [](auto const& left, auto const& right)
                    {
                        if(left.remainder != right.remainder)
                            return left.remainder > right.remainder;
                        return left.index < right.index;
                    });
                for(std::uint64_t extra = remainder - assigned; extra > 0u; --extra)
                    ++result[fractions[(remainder - assigned - extra) % fractions.size()].index].rayCount;
            }
        }

        for(auto& quota : result)
        {
            if(quota.sourceStrength > 0.0 && quota.rayCount == 0u)
                throw std::invalid_argument("positive-source domain received no ASE rays");
            if(quota.rayCount > 0u)
                quota.importanceWeight
                    = quota.sourceStrength * static_cast<double>(totalRayCount) / static_cast<double>(quota.rayCount);
        }
        return result;
    }

    /** @brief Allocate one adaptive launch from the still-unlaunched part of each final quota. */
    [[nodiscard]] inline std::vector<std::uint32_t> allocateDomainLaunchCounts(
        std::span<DomainQuota const> const finalQuotas,
        std::span<std::uint64_t const> const completedCounts,
        std::uint64_t const launchRayCount)
    {
        if(finalQuotas.size() != completedCounts.size())
            throw std::invalid_argument("completed domain counts do not match the final quotas");

        std::uint64_t remainingTotal = 0u;
        std::vector<std::uint64_t> remaining(finalQuotas.size(), 0u);
        for(std::size_t domain = 0u; domain < finalQuotas.size(); ++domain)
        {
            if(completedCounts[domain] > finalQuotas[domain].rayCount)
                throw std::invalid_argument("completed domain count exceeds its final quota");
            remaining[domain] = finalQuotas[domain].rayCount - completedCounts[domain];
            remainingTotal += remaining[domain];
        }
        if(launchRayCount > remainingTotal || launchRayCount > std::numeric_limits<std::uint32_t>::max())
            throw std::invalid_argument("adaptive domain launch exceeds the remaining quota");

        struct Fraction
        {
            std::size_t domain;
            std::uint64_t remainder;
        };

        std::vector<Fraction> fractions;
        fractions.reserve(finalQuotas.size());
        std::vector<std::uint32_t> result(finalQuotas.size(), 0u);
        std::uint64_t assigned = 0u;
        for(std::size_t domain = 0u; domain < finalQuotas.size(); ++domain)
        {
            auto const product = launchRayCount * remaining[domain];
            auto const count = remainingTotal == 0u ? 0u : product / remainingTotal;
            result[domain] = static_cast<std::uint32_t>(count);
            assigned += count;
            fractions.push_back({domain, remainingTotal == 0u ? 0u : product % remainingTotal});
        }
        std::ranges::sort(
            fractions,
            [](auto const& left, auto const& right)
            {
                return left.remainder != right.remainder ? left.remainder > right.remainder
                                                         : left.domain < right.domain;
            });
        for(std::uint64_t extra = launchRayCount - assigned; extra > 0u; --extra)
            ++result[fractions[launchRayCount - assigned - extra].domain];
        return result;
    }

    /** @brief Allocate a surviving boundary population among neighbor routes from their weights. */
    [[nodiscard]] inline std::vector<std::uint32_t> allocateBoundaryRoutePopulations(
        std::span<double const> const routeWeights,
        std::uint32_t const survivingRayCount,
        std::span<std::uint32_t const> const routeCandidateCounts = {})
    {
        if(!routeCandidateCounts.empty() && routeCandidateCounts.size() != routeWeights.size())
            throw std::invalid_argument("boundary route weights and candidate counts must have matching sizes");
        std::vector<std::uint32_t> result(routeWeights.size(), 0u);
        double totalWeight = 0.0;
        std::size_t nonEmptyRoutes = 0u;
        for(double const weight : routeWeights)
        {
            if(!std::isfinite(weight) || weight < 0.0)
                throw std::invalid_argument("boundary route weights must be finite and non-negative");
            totalWeight += weight;
            nonEmptyRoutes += weight > 0.0 ? 1u : 0u;
        }
        if(survivingRayCount == 0u || totalWeight == 0.0)
            return result;
        if(nonEmptyRoutes > survivingRayCount)
            throw std::invalid_argument("surviving boundary population cannot represent every non-empty route");

        struct Fraction
        {
            std::size_t route;
            double remainder;
        };

        std::vector<Fraction> fractions;
        fractions.reserve(nonEmptyRoutes);
        std::uint32_t assigned = 0u;
        for(std::size_t route = 0u; route < routeWeights.size(); ++route)
        {
            if(routeWeights[route] == 0.0)
                continue;
            double const exact = static_cast<double>(survivingRayCount) * routeWeights[route] / totalWeight;
            result[route] = std::max<std::uint32_t>(1u, static_cast<std::uint32_t>(std::floor(exact)));
            assigned += result[route];
            fractions.push_back({route, exact - std::floor(exact)});
        }
        if(assigned > survivingRayCount)
        {
            std::ranges::sort(
                fractions,
                [](auto const& left, auto const& right)
                {
                    return left.remainder != right.remainder ? left.remainder < right.remainder
                                                             : left.route > right.route;
                });
            for(auto const& fraction : fractions)
                if(assigned > survivingRayCount && result[fraction.route] > 1u)
                {
                    auto const excess = std::min(assigned - survivingRayCount, result[fraction.route] - 1u);
                    result[fraction.route] -= excess;
                    assigned -= excess;
                }
        }
        else
        {
            std::ranges::sort(
                fractions,
                [](auto const& left, auto const& right)
                {
                    return left.remainder != right.remainder ? left.remainder > right.remainder
                                                             : left.route < right.route;
                });
            for(std::uint32_t extra = survivingRayCount - assigned; extra > 0u; --extra)
                ++result[fractions[(survivingRayCount - assigned - extra) % fractions.size()].route];
        }
        if(!routeCandidateCounts.empty())
        {
            std::uint32_t overflow = 0u;
            for(std::size_t route = 0u; route < result.size(); ++route)
                if(result[route] > routeCandidateCounts[route])
                {
                    overflow += result[route] - routeCandidateCounts[route];
                    result[route] = routeCandidateCounts[route];
                }
            using Choice = std::pair<double, std::size_t>;
            std::priority_queue<Choice> choices;
            for(std::size_t route = 0u; route < result.size(); ++route)
                if(routeWeights[route] > 0.0 && result[route] < routeCandidateCounts[route])
                    choices.emplace(routeWeights[route] / static_cast<double>(result[route] + 1u), route);
            while(overflow > 0u)
            {
                if(choices.empty())
                    throw std::invalid_argument("surviving boundary population exceeds positive route candidates");
                auto const [priority, route] = choices.top();
                static_cast<void>(priority);
                choices.pop();
                ++result[route];
                --overflow;
                if(result[route] < routeCandidateCounts[route])
                    choices.emplace(routeWeights[route] / static_cast<double>(result[route] + 1u), route);
            }
        }
        if(std::accumulate(result.begin(), result.end(), std::uint64_t{0u}) != survivingRayCount)
            throw std::runtime_error("boundary route population allocation did not preserve the surviving count");
        return result;
    }

    /** @brief Preserve the fixed passive-domain relaunch slots used by the SRM path. */
    [[nodiscard]] inline std::vector<std::uint32_t> reservePassiveDomainPopulationSlots(
        std::span<std::uint32_t const> const sourceCounts)
    {
        std::vector<std::uint32_t> result(sourceCounts.begin(), sourceCounts.end());
        std::uint64_t total = 0u;
        for(auto const count : sourceCounts)
            total += count;
        if(total < result.size())
            return result;

        for(std::size_t domain = 0u; domain < result.size(); ++domain)
        {
            if(result[domain] != 0u)
                continue;
            auto const donor = std::ranges::max_element(result);
            if(donor == result.end() || *donor <= 1u)
                throw std::runtime_error("cannot reserve a passive domain population slot");
            --*donor;
            result[domain] = 1u;
        }
        return result;
    }

    /**
     * @brief Build an exact worker-by-domain population matrix with fixed row and column totals.
     *
     * Worker totals preserve the histories already assigned to each worker, while domain totals
     * preserve the requested population after boundary resampling. Largest remainders make the
     * integer allocation deterministic without changing either marginal total.
     */
    [[nodiscard]] inline std::vector<std::vector<std::uint32_t>> distributeDomainPopulations(
        std::span<std::uint32_t const> const domainCounts,
        std::span<std::uint32_t const> const workerCounts)
    {
        if(domainCounts.empty() || workerCounts.empty())
            throw std::invalid_argument("domain population distribution requires domains and workers");
        auto sum = [](std::span<std::uint32_t const> const values)
        {
            std::uint64_t result = 0u;
            for(auto const value : values)
                result += value;
            return result;
        };
        auto const total = sum(domainCounts);
        if(total != sum(workerCounts))
            throw std::invalid_argument("domain and worker population totals must match");

        std::vector<std::vector<std::uint32_t>> result(
            workerCounts.size(),
            std::vector<std::uint32_t>(domainCounts.size(), 0u));
        std::vector<std::uint64_t> remainingWorkers(workerCounts.begin(), workerCounts.end());
        std::uint64_t remainingTotal = total;
        for(std::size_t domain = 0u; domain < domainCounts.size(); ++domain)
        {
            auto const domainCount = static_cast<std::uint64_t>(domainCounts[domain]);
            if(domain + 1u == domainCounts.size())
            {
                for(std::size_t worker = 0u; worker < workerCounts.size(); ++worker)
                    result[worker][domain] = static_cast<std::uint32_t>(remainingWorkers[worker]);
                break;
            }
            if(domainCount == 0u)
                continue;

            struct Fraction
            {
                std::size_t worker;
                std::uint64_t remainder;
            };

            std::vector<Fraction> fractions;
            fractions.reserve(workerCounts.size());
            std::uint64_t assigned = 0u;
            for(std::size_t worker = 0u; worker < workerCounts.size(); ++worker)
            {
                auto const product = domainCount * remainingWorkers[worker];
                auto const count = remainingTotal == 0u ? 0u : product / remainingTotal;
                result[worker][domain] = static_cast<std::uint32_t>(count);
                assigned += count;
                fractions.push_back({worker, remainingTotal == 0u ? 0u : product % remainingTotal});
            }
            std::ranges::sort(
                fractions,
                [](auto const& left, auto const& right)
                {
                    return left.remainder != right.remainder ? left.remainder > right.remainder
                                                             : left.worker < right.worker;
                });
            std::uint64_t extras = domainCount - assigned;
            for(auto const& fraction : fractions)
            {
                if(extras == 0u)
                    break;
                if(result[fraction.worker][domain] < remainingWorkers[fraction.worker])
                {
                    ++result[fraction.worker][domain];
                    --extras;
                }
            }
            if(extras != 0u)
                throw std::runtime_error("domain population rounding exhausted worker capacity");
            for(std::size_t worker = 0u; worker < workerCounts.size(); ++worker)
                remainingWorkers[worker] -= result[worker][domain];
            remainingTotal -= domainCount;
        }
        return result;
    }

    /** @brief Build the deterministic locality-aware static domain/batch schedule. */
    [[nodiscard]] DomainSchedule makeDomainSchedule(
        std::vector<WorkerDescriptor> const& workers,
        std::vector<DomainCost> const& domains,
        std::vector<DomainQuota> const& quotas,
        std::uint32_t numIndependentRayPopulations,
        std::span<data::AseDomainInterface const> interfaces = {});

    /** @brief Schedule concrete execution batches without changing statistical population identity. */
    [[nodiscard]] DomainSchedule makeDomainSchedule(
        std::vector<WorkerDescriptor> const& workers,
        std::vector<DomainCost> const& domains,
        std::span<DomainWorkItem const> work,
        std::span<data::AseDomainInterface const> interfaces = {});

    /** @brief Compute immutable scheduling statistics once from the prepared domain graph. */
    [[nodiscard]] std::vector<DomainCost> makeDomainCosts(data::AseDomainGraph const& graph);
} // namespace hase::core
