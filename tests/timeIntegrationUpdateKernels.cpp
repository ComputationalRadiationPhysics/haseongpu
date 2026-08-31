#include <alpaka/alpaka.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <kernels/timeIntegrationUpdateKernels.hpp>

#include <cmath>
#include <cstdint>
#include <vector>

TEST_CASE("frozen-source derivative changes only with stage decay", "[time-integration][simd]")
{
    std::vector<std::uint8_t> const materialActive{1u, 0u};
    std::vector<double> const fluorescenceLifetimes{2.0, 0.0};
    hase::kernels::ComposeFrozenSourcesDerivative const compose{materialActive, fluorescenceLifetimes};

    alpaka::Simd<double, 4u> const betaVolume{0.2, 0.7, 0.4, 0.9};
    alpaka::Simd<double, 4u> const dndtPump{0.3, 0.8, 0.5, 1.0};
    alpaka::Simd<double, 4u> const dndtAse{0.1, 0.2, 0.4, 0.1};
    alpaka::Simd<unsigned, 4u> const materialIds{0u, 1u, 0u, 1u};

    auto const result = compose(betaVolume, dndtPump, dndtAse, materialIds);
    CHECK(result[0u] == Catch::Approx(0.3 - 0.1 - 0.2 / 2.0));
    CHECK(result[1u] == 0.0);
    CHECK(result[2u] == Catch::Approx(0.5 - 0.4 - 0.4 / 2.0));
    CHECK(result[3u] == 0.0);
}

TEST_CASE("frozen-source RK4 matches the scalar Julia update", "[time-integration][simd]")
{
    double const timeStep = 0.5;
    std::vector<std::uint8_t> const materialActive{1u};
    std::vector<double> const fluorescenceLifetimes{2.0};
    hase::kernels::ComposeFrozenSourcesDerivative const compose{materialActive, fluorescenceLifetimes};
    hase::kernels::CombineRungeKutta4 const combine{timeStep};

    alpaka::Simd<double, 4u> const beta{0.4, 0.4, 0.4, 0.4};
    alpaka::Simd<double, 4u> const pump{0.3, 0.3, 0.0, 0.0};
    alpaka::Simd<double, 4u> const ase{0.1, 0.0, 0.1, 0.0};
    alpaka::Simd<unsigned, 4u> const materialIds{0u, 0u, 0u, 0u};

    auto const k1 = compose(beta, pump, ase, materialIds);
    auto const k2 = compose(beta + 0.5 * timeStep * k1, pump, ase, materialIds);
    auto const k3 = compose(beta + 0.5 * timeStep * k2, pump, ase, materialIds);
    auto const k4 = compose(beta + timeStep * k3, pump, ase, materialIds);
    auto const result = combine(beta, k1, k2, k3, k4);

    auto const juliaStep = [timeStep](double y, double pumpRate, double aseRate)
    {
        double const source = pumpRate - aseRate;
        double const k1 = source - y / 2.0;
        double const k2 = source - (y + 0.5 * timeStep * k1) / 2.0;
        double const k3 = source - (y + 0.5 * timeStep * k2) / 2.0;
        double const k4 = source - (y + timeStep * k3) / 2.0;
        return y + (timeStep / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    };
    for(std::size_t lane = 0u; lane < 4u; ++lane)
        CHECK(result[lane] == Catch::Approx(juliaStep(beta[lane], pump[lane], ase[lane])));
}

TEST_CASE("exponential Euler update supports mixed-material SIMD packs", "[time-integration][simd]")
{
    std::vector<std::uint8_t> const materialActive{1u, 0u};
    std::vector<double> const fluorescenceLifetimes{2.0, 0.0};
    hase::kernels::ExponentialEulerUpdate const update{0.5, materialActive, fluorescenceLifetimes};

    alpaka::Simd<double, 4u> const betaVolume{0.2, 0.7, 0.4, 0.9};
    alpaka::Simd<double, 4u> const dndtPump{0.3, 0.8, 0.5, 1.0};
    alpaka::Simd<double, 4u> const dndtAse{0.1, 0.2, 0.4, 0.1};
    alpaka::Simd<unsigned, 4u> const materialIds{0u, 1u, 0u, 1u};

    auto const result = update(betaVolume, dndtPump, dndtAse, materialIds);
    double const decay = std::exp(-0.5 / 2.0);
    CHECK(result[0u] == Catch::Approx(2.0 * (0.3 - 0.1) * (1.0 - decay) + 0.2 * decay));
    CHECK(result[1u] == 0.0);
    CHECK(result[2u] == Catch::Approx(2.0 * (0.5 - 0.4) * (1.0 - decay) + 0.4 * decay));
    CHECK(result[3u] == 0.0);
}
