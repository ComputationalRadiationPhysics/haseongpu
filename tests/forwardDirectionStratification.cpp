#include <catch2/catch_test_macros.hpp>
#include <kernels/forward/rayBundling.hpp>

TEST_CASE("forward direction strata cover all octants", "[forward][direction]")
{
    using hase::core::Point;
    using hase::kernels::forward::directionStratum;
    CHECK(directionStratum(Point{1.0, 1.0, -1.0}) == 0u);
    CHECK(directionStratum(Point{-1.0, 1.0, -1.0}) == 1u);
    CHECK(directionStratum(Point{-1.0, -1.0, -1.0}) == 2u);
    CHECK(directionStratum(Point{1.0, -1.0, -1.0}) == 3u);
    CHECK(directionStratum(Point{1.0, 1.0, 1.0}) == 4u);
    CHECK(directionStratum(Point{-1.0, 1.0, 1.0}) == 5u);
    CHECK(directionStratum(Point{-1.0, -1.0, 1.0}) == 6u);
    CHECK(directionStratum(Point{1.0, -1.0, 1.0}) == 7u);
}
