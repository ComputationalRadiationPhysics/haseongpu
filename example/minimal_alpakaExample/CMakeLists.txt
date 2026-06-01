cmake_minimum_required(VERSION 3.24)
project(MinimalExample)
# install alpaka
find_package(alpaka CONFIG)
if(NOT alpaka_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        alpaka
        GIT_REPOSITORY https://github.com/alpaka-group/alpaka3.git
        GIT_TAG 245c74bf96d751cc52ec3f083934cf58d2b7313b
    )
    FetchContent_MakeAvailable(alpaka)
endif()
# --- main executable: main.cu + core lib ---
add_executable(minimalExample "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
target_link_libraries(minimalExample PRIVATE alpaka::alpaka)
alpaka_finalize(minimalExample)
