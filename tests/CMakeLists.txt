cmake_minimum_required(VERSION 3.24)

include(FetchContent)

# Catch2 v3
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.6.0
)
FetchContent_MakeAvailable(Catch2)

# Collect all test translation units in this folder (cpp + cu)
file(GLOB TEST_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cu"
)

add_executable(hase_tests ${TEST_SOURCES})

# Inherit all include dirs, compile options, defines, and link libs via hase_core
target_link_libraries(hase_tests PRIVATE
        hase_core
        Catch2::Catch2WithMain
)

# If your tests include CUDA code and need separable compilation too:
set_target_properties(hase_tests PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${HASE_CUDA_ARCHITECTURES}"
)

include(CTest)
include(Catch)

catch_discover_tests(hase_tests)