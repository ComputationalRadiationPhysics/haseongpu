/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#include <catch2/catch_test_macros.hpp>
#include <core/logging.hpp>

TEST_CASE("default logging keeps operational progress enabled")
{
    CHECK((hase::core::verbosity & V_ERROR) != 0u);
    CHECK((hase::core::verbosity & V_WARNING) != 0u);
    CHECK((hase::core::verbosity & V_INFO) != 0u);
    CHECK((hase::core::verbosity & V_STAT) != 0u);
    CHECK((hase::core::verbosity & V_PROGRESS) != 0u);
}

TEST_CASE("debug logging follows its build option")
{
#ifdef HASE_ENABLE_DEBUG_LOGGING
    CHECK((hase::core::verbosity & V_DEBUG) != 0u);
#else
    CHECK((hase::core::verbosity & V_DEBUG) == 0u);
#endif
}
