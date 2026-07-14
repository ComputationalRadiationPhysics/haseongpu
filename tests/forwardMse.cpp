#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <core/calcForwardPhiAse.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <ranges>
#include <string>
#include <type_traits>

TEST_CASE("forward PhiASE RSE includes zero-score histories", "[forward][rse]")
{
    // Four globally launched histories with cell scores [1, 3, 0, 0].
    // The forward cell estimator scales the per-history mean by totalVolume / cellVolume.
    double const sum = 4.0;
    double const sumSquares = 10.0;
    unsigned const rayCount = 4u;
    double const totalVolume = 8.0;
    double const cellVolume = 4.0;

    double const expectedRelativeStandardError = std::sqrt((rayCount * sumSquares / (sum * sum) - 1.0) / rayCount);
    double const expectedStandardError = expectedRelativeStandardError * (sum * totalVolume / (rayCount * cellVolume));

    CHECK(
        hase::core::calcForwardRelativeStandardError(sum, sumSquares, rayCount)
        == Catch::Approx(expectedRelativeStandardError));
    CHECK(
        hase::core::calcForwardStandardError(sum, sumSquares, rayCount, totalVolume, cellVolume)
        == Catch::Approx(expectedStandardError));
}

TEST_CASE("forward PhiASE RSE handles invalid and zero-score estimates", "[forward][rse]")
{
    CHECK(hase::core::calcForwardRelativeStandardError(1.0, 1.0, 1u) == std::numeric_limits<double>::max());
    CHECK(std::isnan(hase::core::calcForwardRelativeStandardError(0.0, 0.0, 2u)));
    CHECK(
        hase::core::calcForwardRelativeStandardError(std::numeric_limits<double>::infinity(), 1.0, 2u)
        == std::numeric_limits<double>::max());
    CHECK(hase::core::calcForwardStandardError(1.0, 1.0, 1u, 1.0, 1.0) == std::numeric_limits<double>::max());
    CHECK(hase::core::calcForwardStandardError(1.0, 1.0, 2u, 0.0, 1.0) == 0.0);
    CHECK(hase::core::calcForwardStandardError(0.0, 0.0, 2u, 1.0, 1.0) == 0.0);
    CHECK(hase::core::calcForwardStandardError(1.0, 1.0, 2u, 1.0, 0.0) == std::numeric_limits<double>::max());
    CHECK(
        hase::core::calcForwardStandardError(std::numeric_limits<double>::infinity(), 1.0, 2u, 1.0, 1.0)
        == std::numeric_limits<double>::max());
}

TEST_CASE("forward PhiASE beta-volume contribution uses double precision", "[forward][rse]")
{
    hase::core::BetaVolumeContribution contribution;
    auto const value = contribution(alpaka::Simd<double, 1u>{0.25}, alpaka::Simd<float, 1u>{0.5f});
    STATIC_REQUIRE(std::is_same_v<alpaka::trait::GetValueType_t<std::remove_cvref_t<decltype(value)>>, double>);
    CHECK(value[0] == Catch::Approx(0.125));
}

TEST_CASE("forward spectrum stratification balances discrete bins", "[forward][sampling]")
{
    constexpr unsigned spectrumSize = 7u;
    constexpr unsigned rayCount = 25u;
    std::array<unsigned, spectrumSize> visits{};
    for(unsigned ray = 0u; ray < rayCount; ++ray)
    {
        ++visits.at(hase::kernels::forward::stratifiedSpectrumIndex(spectrumSize, ray, rayCount, 3u));
    }

    auto const [minimum, maximum] = std::ranges::minmax_element(visits);
    CHECK(*maximum - *minimum <= 1u);
    CHECK(std::accumulate(visits.cbegin(), visits.cend(), 0u) == rayCount);
}

TEST_CASE("forward source stratification places one shifted point in each CDF interval", "[forward][sampling]")
{
    constexpr unsigned rayCount = 10u;
    constexpr double shift = 0.25;
    for(unsigned ray = 0u; ray < rayCount; ++ray)
    {
        double const target = hase::kernels::forward::stratifiedUnitInterval(ray, rayCount, shift);
        CHECK(target > static_cast<double>(ray) / rayCount);
        CHECK(target < static_cast<double>(ray + 1u) / rayCount);
    }
}

TEST_CASE("forward SRM environment controls are strict positive overrides", "[forward][srm]")
{
    auto const restore = [](char const* name, char const* value)
    {
        if(value == nullptr)
            unsetenv(name);
        else
            setenv(name, value, 1);
    };
    char const* oldMaxIterations = std::getenv("HASE_SRM_MAX_ITERATIONS");
    char const* oldDivergenceStreak = std::getenv("HASE_SRM_DIVERGENCE_STREAK");
    std::string const savedMaxIterations = oldMaxIterations == nullptr ? "" : oldMaxIterations;
    std::string const savedDivergenceStreak = oldDivergenceStreak == nullptr ? "" : oldDivergenceStreak;

    unsetenv("HASE_SRM_MAX_ITERATIONS");
    unsetenv("HASE_SRM_DIVERGENCE_STREAK");
    hase::core::ExperimentParameters experiment{};
    experiment.reflectionMaxIterations = 8u;
    auto const defaults = hase::core::resolveSrmControls(experiment);
    CHECK(defaults.maxIterations == 8u);
    CHECK(defaults.divergenceStreak == 3u);

    setenv("HASE_SRM_MAX_ITERATIONS", "11", 1);
    setenv("HASE_SRM_DIVERGENCE_STREAK", "4", 1);
    auto const overridden = hase::core::resolveSrmControls(experiment);
    CHECK(overridden.maxIterations == 11u);
    CHECK(overridden.divergenceStreak == 4u);

    setenv("HASE_SRM_DIVERGENCE_STREAK", "0", 1);
    CHECK_THROWS_AS(hase::core::resolveSrmControls(experiment), std::runtime_error);

    restore("HASE_SRM_MAX_ITERATIONS", oldMaxIterations == nullptr ? nullptr : savedMaxIterations.c_str());
    restore("HASE_SRM_DIVERGENCE_STREAK", oldDivergenceStreak == nullptr ? nullptr : savedDivergenceStreak.c_str());
}
