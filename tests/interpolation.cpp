/**
 * Copyright 2026 Tim Hanel
 **/

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <interpolation.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

std::vector<double> make_sorted_unique_x(std::size_t n, uint32_t seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);

    std::vector<double> x(n);
    for(auto& v : x)
        v = dist(rng);

    std::sort(x.begin(), x.end());

    // enforce strictly increasing (avoid dx==0)
    for(std::size_t i = 1; i < x.size(); ++i)
    {
        if(x[i] <= x[i - 1])
            x[i] = std::nextafter(x[i - 1], std::numeric_limits<double>::infinity());
    }
    return x;
}

std::vector<double> make_linear_y(std::vector<double> const& x, double a, double b)
{
    std::vector<double> y(x.size());
    for(std::size_t i = 0; i < x.size(); ++i)
        y[i] = a * x[i] + b;
    return y;
}

std::vector<double> reconstruct_interpolated_x(std::vector<double> const& x, unsigned nInterpolations)
{
    if(nInterpolations == 0)
        return {};
    if(nInterpolations == 1)
        return {x.front()};
    double const x_min = x.front();
    double const x_max = x.back();
    double const x_range = x_max - x_min;

    std::vector<double> xi(nInterpolations);
    for(unsigned i = 0; i < nInterpolations; ++i)
    {
        xi[i] = x_min + (double(i) * x_range) / double(nInterpolations - 1);
    }
    return xi;
}

TEST_CASE("interpolateLinear: exact for linear functions (fixed seeds)", "[interpolateLinear]")
{
    std::vector<std::size_t> const sizes = {100, 1000, 10000};
    std::vector<unsigned> const nInterps = {100, 1000, 10000, 100000};
    std::mt19937 rng(0); // fixed seed rng
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);
    // choose fixed m,b (fixed-seed generated over interval [-1000,1000])
    double const m = dist(rng);
    double const b = dist(rng);

    for(std::size_t n : sizes)
    {
        auto x = make_sorted_unique_x(n, /*seed=*/0xC0'FFEEu + static_cast<uint32_t>(n));
        auto y = make_linear_y(x, m, b);

        for(unsigned ni : nInterps)
        {
            auto y_interp = interpolateLinear(y, x, ni);
            REQUIRE(y_interp.size() == ni);

            auto x_interp = reconstruct_interpolated_x(x, ni);

            double max_abs_err = 0.0;
            unsigned worst_j = 0;

            for(unsigned j = 0; j < ni; ++j)
            {
                double const expected = m * x_interp[j] + b;
                double const err = std::abs(y_interp[j] - expected);

                if(err > max_abs_err)
                {
                    max_abs_err = err;
                    worst_j = j;
                }
            }
            CAPTURE(n, ni, worst_j, max_abs_err);
            // tolerance: depends on your formula; for correct interpolation it should be tiny
            REQUIRE(max_abs_err < 1e-8);
        }
    }
}

TEST_CASE("interpolateLinear: no overshoot between adjacent samples (fixed seeds)", "[interpolateLinear]")
{
    std::size_t const n = 1000;
    unsigned const ni = 100000;

    auto x = make_sorted_unique_x(n, /*seed=*/0xBAD'C0DEu);

    // fixed-seed y (random-looking but deterministic)
    std::mt19937 rng(0xDEAD'BEEFu);
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);
    std::vector<double> y(n);
    for(auto& v : y)
        v = dist(rng);

    auto y_interp = interpolateLinear(y, x, ni);
    auto x_interp = reconstruct_interpolated_x(x, ni);

    // Walk segments i and interpolated points j
    std::size_t i = 0;
    for(unsigned j = 0; j < ni; ++j)
    {
        while(i + 1 < x.size() && x_interp[j] >= x[i + 1])
        {
            ++i;
        }
        if(i + 1 >= x.size())
            break;

        double const lo = std::min(y[i], y[i + 1]);
        double const hi = std::max(y[i], y[i + 1]);

        // allow tiny epsilon for floating ops
        REQUIRE(y_interp[j] >= lo - 1e-9);
        REQUIRE(y_interp[j] <= hi + 1e-9);
    }
}
