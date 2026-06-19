#include <catch2/catch_test_macros.hpp>
#include <hase/version.hpp>

#include <string>

TEST_CASE("HASEonGPU version macros expose string and comparable numeric version", "[version]")
{
    static_assert(HASEONGPU_VERSION_MAJOR == 2);
    static_assert(HASEONGPU_VERSION_MINOR == 0);
    static_assert(HASEONGPU_VERSION_PATCH == 0);
    static_assert(HASEONGPU_VERSION == HASEONGPU_VERSION_NUMBER(2, 0, 0));

    CHECK(std::string{HASEONGPU_VERSION_STRING} == "2.0.0");
}
