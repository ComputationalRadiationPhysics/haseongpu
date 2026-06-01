cmake_minimum_required(VERSION 3.24)
project(HASEonGPU)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
include(CMakePrintHelpers)

# dependencies
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30")
    cmake_policy(SET CMP0167 NEW)
endif()

find_package(Threads REQUIRED)
find_package(CUDAToolkit REQUIRED COMPONENTS cudart)

# build options
option(HASE_BUILD_PhiAse "Build the main calcPhiAse - HaseOnGPU Binary" ON)
option(
    HASE_SELECT_BACKEND_ALPAKA
    "Use the cmake options provided by alpaka in order to select backend and executor manually"
    OFF
)
option(HASE_CUDA_FLUSHTOZERO "Set flush to zero for GPU" OFF)
option(HASE_CUDA_FASTMATH "Enable fast-math" ON)
option(HASE_CUDA_SHOW_REGISTER "Show kernel registers and create PTX" OFF)
option(
    HASE_CUDA_KEEP_FILES
    "Keep all intermediate files that are generated during internal compilation steps (folder: nvcc_tmp)"
    OFF
)
option(
    HASE_CUDA_SHOW_CODELINES
    "Show kernel lines in cuda-gdb and cuda-memcheck"
    OFF
)
set(DISABLE_MPI
    "AUTO"
    CACHE STRING
    "MPI support selection: AUTO -> MPI detection, ON -> MPI disabled, or OFF -> MPI enabled"
)

set_property(CACHE DISABLE_MPI PROPERTY STRINGS AUTO ON OFF)

if(NOT DISABLE_MPI MATCHES "^(AUTO|ON|OFF)$")
    message(
        FATAL_ERROR
        "Invalid value for DISABLE_MPI='${DISABLE_MPI}'. "
        "Allowed values are: AUTO, ON, OFF."
    )
endif()
option(HASE_TESTING "Build unit tests" OFF)
option(HASE_ENABLE_PYTHON "Build Python bindings" OFF)

set(HASE_CUDA_ARCHITECTURES native CACHE STRING "CUDA Architectures")

# search for backends (dependent) -- important that this runs before making alpaka available
if(NOT HASE_SELECT_BACKEND_ALPAKA)
    include(cmake/searchBackends.cmake)
endif()

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

# CUDA compiler Flags
set(HASE_NVCC_FLAGS "")

if(HASE_CUDA_SHOW_CODELINES)
    set(HASE_CUDA_KEEP_FILES ON)
    list(APPEND HASE_NVCC_FLAGS --source-in-ptx -Xcompiler -rdynamic -lineinfo)
endif()

if(HASE_CUDA_FLUSHTOZERO)
    list(APPEND HASE_NVCC_FLAGS --ftz=true)
else()
    list(APPEND HASE_NVCC_FLAGS --ftz=false)
endif()

if(HASE_CUDA_FASTMATH)
    list(APPEND HASE_NVCC_FLAGS --use_fast_math)
endif()

if(HASE_CUDA_KEEP_FILES)
    set(NVCC_TMP_DIR "${PROJECT_BINARY_DIR}/nvcc_tmp")
    file(MAKE_DIRECTORY "${NVCC_TMP_DIR}")
    list(APPEND HASE_NVCC_FLAGS --keep --keep-dir "${NVCC_TMP_DIR}")
endif()

cmake_print_variables(HASE_NVCC_FLAGS)

# C++ compiler Flags
set(HASE_GCC_FLAGS -Wall -Wextra)
set(HASE_CLANG_FLAGS -Wall -Wextra)

# warning for Thrust-related Bug
function(warn_path_cuda PATH_VARIABLE)
    string(TOLOWER "$ENV{${PATH_VARIABLE}}" CPLUS_TOLOWER)
    string(FIND "${CPLUS_TOLOWER}" "cuda" CPLUS_TOLOWER_FOUND)
    if(NOT "${CPLUS_TOLOWER_FOUND}" MATCHES "-1")
        message(
            WARNING
            "Your ${PATH_VARIABLE} seems to contain CUDA includes."
            "This might result in the following error:\n"
            "'error: kernel launches from templates are not allowed in system files'\n"
            "To fix the problem, remove the CUDA includes from the path."
            "(A better place for these includes might be your CPATH)\n"
            "source: https://code.google.com/p/thrust/issues/detail?id=359#c5 \n"
            "issue:  https://github.com/ComputationalRadiationPhysics/haseongpu/issues/26"
            "CPATH:  https://gcc.gnu.org/onlinedocs/cpp/Environment-Variables.html \n"
        )
    endif()
endfunction()

warn_path_cuda("CPLUS_INCLUDE_PATH")
warn_path_cuda("C_INCLUDE_PATH")

# build type (debug, release)
option(HASE_RELEASE "Build release version, disables all runtime asserts" ON)
if(HASE_RELEASE)
    message(STATUS "Release version")
    set(CMAKE_BUILD_TYPE Release)
endif()

# core library with all sources
add_library(hase_core INTERFACE)
add_library(hase::core ALIAS hase_core)
#target_sources(hase_core PRIVATE ${HASE_CORE_SOURCES})
target_include_directories(
    hase_core
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

set_target_properties(
    hase_core
    PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${HASE_CUDA_ARCHITECTURES}"
)

if(DISABLE_MPI STREQUAL "AUTO")
    find_package(MPI QUIET)
    if(NOT MPI_FOUND)
        message(
            STATUS
            "MPI is not found on your system - compilation will proceed without MPI Support."
        )
    endif()
endif()
if(DISABLE_MPI STREQUAL "OFF" OR MPI_FOUND)
    target_compile_definitions(hase_core INTERFACE MPI_FOUND)
    target_link_libraries(hase_core INTERFACE MPI::MPI_CXX)
else()
    message(STATUS "MPI disabled.")
    target_compile_definitions(hase_core INTERFACE DISABLE_MPI)
endif()

target_link_libraries(hase_core INTERFACE alpaka::alpaka CUDA::cudart)

file(
    GLOB_RECURSE HASE_RUNTIME_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu"
)

list(
    REMOVE_ITEM
    HASE_RUNTIME_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
)
add_library(hase STATIC ${HASE_RUNTIME_SOURCES})
add_library(hase::hase ALIAS hase)
set_target_properties(hase PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(hase PUBLIC hase::core)
alpaka_finalize(hase)

if(HASE_BUILD_PhiAse)
    # --- main executable: main.cpp ---
    add_executable(calcPhiASE "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
    set_target_properties(calcPhiASE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(calcPhiASE PRIVATE hase::hase)
    alpaka_finalize(calcPhiASE)
endif()

# enable Testing
if(HASE_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
if(HASE_ENABLE_PYTHON)
    set_target_properties(hase PROPERTIES POSITION_INDEPENDENT_CODE ON)
    find_package(pybind11 REQUIRED)
    pybind11_add_module(HASEonGPU HASEonGPU_Bindings/module.cpp)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
    set_target_properties(
        HASEonGPU
        PROPERTIES BUILD_RPATH_USE_ORIGIN ON INSTALL_RPATH_USE_LINK_PATH ON
    )
    target_link_libraries(HASEonGPU PRIVATE hase)
    install(TARGETS HASEonGPU LIBRARY DESTINATION HASEonGPU_Bindings)
    install(TARGETS calcPhiASE RUNTIME DESTINATION HASEonGPU_Bindings)
    install(FILES HASEonGPU.py DESTINATION .)
    install(
        DIRECTORY HASEonGPU_Bindings
        DESTINATION .
        FILES_MATCHING
        PATTERN "__init__.py"
    )
endif()
