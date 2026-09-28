/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <core/reflectionTail.hpp>

#include <cmath>
#include <string>
#include <vector>

TEST_CASE("boundary gamma fit recovers a geometric reflected population", "[forward][reflection-tail]")
{
    std::vector<double> residualFractions;
    for(unsigned pass = 0u; pass < 24u; ++pass)
        residualFractions.push_back(std::pow(0.9, static_cast<double>(pass)));

    auto const fit = hase::core::fitBoundaryGamma(residualFractions);

    REQUIRE(fit.valid);
    CHECK(fit.gamma == Catch::Approx(0.9).epsilon(1.0e-12));
    CHECK(fit.standardError == Catch::Approx(0.0).margin(1.0e-12));
}

TEST_CASE("stationary subcritical reflection tail is completed geometrically", "[forward][reflection-tail]")
{
    std::vector<double> residualFractions;
    for(unsigned pass = 0u; pass < 24u; ++pass)
        residualFractions.push_back(std::pow(0.9, static_cast<double>(pass)));

    auto const estimate = hase::core::estimateBoundaryTail(residualFractions);

    CHECK_FALSE(estimate.divergent);
    REQUIRE(estimate.applicable);
    CHECK(estimate.tailFactor == Catch::Approx(9.0).epsilon(1.0e-12));
    CHECK(estimate.tailClosure == Catch::Approx(1.0).epsilon(1.0e-12));
}

TEST_CASE("confidently growing reflected population is divergent", "[forward][reflection-tail]")
{
    std::vector<double> residualFractions;
    for(unsigned pass = 0u; pass < 24u; ++pass)
        residualFractions.push_back(std::pow(1.02, static_cast<double>(pass)));

    auto const estimate = hase::core::estimateBoundaryTail(residualFractions);

    CHECK(estimate.divergent);
    CHECK_FALSE(estimate.applicable);
    CHECK(estimate.tailFactor == 0.0);
}

TEST_CASE("spectrally hardening reflected population refuses a stationary tail", "[forward][reflection-tail]")
{
    std::vector<double> residualFractions{1.0};
    for(unsigned pass = 1u; pass < 30u; ++pass)
    {
        double const gamma = 0.8 + 0.006 * static_cast<double>(pass);
        residualFractions.push_back(residualFractions.back() * gamma);
    }

    auto const estimate = hase::core::estimateBoundaryTail(residualFractions);

    CHECK_FALSE(estimate.divergent);
    CHECK_FALSE(estimate.applicable);
}

TEST_CASE("material integration accepts supported reflected ASE statuses", "[forward][reflection-tail]")
{
    hase::data::PhiAseResult result;

    result.boundaryStatus = hase::data::BoundaryStatus::disabled;
    CHECK_NOTHROW(hase::core::requireUsableBoundaryAseForIntegration(result, 0u));

    result.boundaryStatus = hase::data::BoundaryStatus::converged;
    CHECK_NOTHROW(hase::core::requireUsableBoundaryAseForIntegration(result, 0u));
}

TEST_CASE("material integration rejects an unresolved reflected ASE tally", "[forward][reflection-tail]")
{
    hase::data::PhiAseResult result;
    result.boundaryStatus = hase::data::BoundaryStatus::maxPasses;
    result.boundaryPasses = 162u;
    result.boundaryMaxPasses = 162u;
    result.boundaryRemainingFraction = 0.016;
    result.boundaryGamma = 0.9988;
    result.boundaryGammaStandardError = 0.0007;
    result.boundaryTailFactor = 818.0;
    result.boundaryTailClosure = 1.085;

    try
    {
        hase::core::requireUsableBoundaryAseForIntegration(result, 3u);
        FAIL("unresolved reflected ASE should stop material integration");
    }
    catch(std::runtime_error const& error)
    {
        std::string const message = error.what();
        CHECK(message.find("material step 4") != std::string::npos);
        CHECK(message.find("boundaryStatus=maxPasses") != std::string::npos);
        CHECK(message.find("boundaryPasses=162/162") != std::string::npos);
        CHECK(message.find("partial PhiASE tally was not integrated") != std::string::npos);
    }

    result.boundaryStatus = hase::data::BoundaryStatus::diverged;
    CHECK_THROWS_AS(hase::core::requireUsableBoundaryAseForIntegration(result, 3u), std::runtime_error);
}
