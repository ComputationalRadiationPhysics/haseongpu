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

# build options
option(HASE_BUILD_PhiAse "Build the main calcPhiAse - HaseOnGPU Binary" ON)
option(
    HASE_BUILD_RELEASE
    "Build release version; forces Release build type and release optimization options."
    ON
)
option(
    HASE_SELECT_BACKEND_ALPAKA
    "Use the cmake options provided by alpaka in order to select backend and executor manually"
    OFF
)
option(
    HASE_USE_SYSTEM_ALPAKA
    "Use an existing alpaka package from alpaka_DIR/CMAKE_PREFIX_PATH instead of fetching the pinned version"
    OFF
)
option(HASE_CUDA_FLUSHTOZERO "Set flush to zero for GPU" OFF)
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
option(
    HASE_ENABLE_PYTHON
    "Build Python bindings - requires python dependency."
    ON
)
option(HASE_BENCHMARK "Enable scoped PhiASE benchmark CSV output" OFF)

set(HASE_CUDA_ARCHITECTURES native CACHE STRING "CUDA Architectures")

if(HASE_BUILD_RELEASE)
    message(
        STATUS
        "HASE_BUILD_RELEASE=ON: forcing CMAKE_BUILD_TYPE=Release and enabling release optimization options. "
        "Disable HASE_BUILD_RELEASE to keep user-defined CMAKE_BUILD_TYPE or optimization settings."
    )
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# search for backends (dependent) -- important that this runs before making alpaka available
if(NOT HASE_SELECT_BACKEND_ALPAKA)
    include(cmake/searchBackends.cmake)
endif()

if(NOT alpaka_DEP_CUDA)
    set(alpaka_EXEC_GpuCuda OFF CACHE BOOL "Enable/Disable CUDA executor" FORCE)
endif()
if(NOT alpaka_DEP_HIP)
    set(alpaka_EXEC_GpuHip OFF CACHE BOOL "Enable/Disable HIP executor" FORCE)
endif()
if(NOT alpaka_DEP_ONEAPI)
    set(alpaka_EXEC_OneApi
        OFF
        CACHE BOOL
        "Enable/Disable Intel OneAPI SYCL executor"
        FORCE
    )
endif()

if(alpaka_DEP_HIP)
    set(alpaka_RELOCATABLE_DEVICE_CODE
        ON
        CACHE STRING
        "Enable relocatable device code so HIP kernels can call device functions from other translation units"
        FORCE
    )
endif()

include(cmake/findAlpaka.cmake)
# CUDA compiler Flags
set(HASE_NVCC_FLAGS "")
set(HASE_CLANG_CUDA_FLAGS "")

if(HASE_CUDA_SHOW_CODELINES)
    set(HASE_CUDA_KEEP_FILES ON)
    list(APPEND HASE_NVCC_FLAGS --source-in-ptx -Xcompiler -rdynamic -lineinfo)
endif()

if(HASE_CUDA_FLUSHTOZERO)
    list(APPEND HASE_NVCC_FLAGS --ftz=true)
    list(APPEND HASE_CLANG_CUDA_FLAGS -fcuda-flush-denormals-to-zero)
else()
    list(APPEND HASE_NVCC_FLAGS --ftz=false)
endif()

if(HASE_BUILD_RELEASE)
    list(APPEND HASE_NVCC_FLAGS --use_fast_math)
    list(APPEND HASE_CLANG_CUDA_FLAGS -ffast-math -ffp-contract=fast)
endif()
if(HASE_BUILD_RELEASE AND alpaka_DEP_CUDA)
    list(
        APPEND
        HASE_NVCC_FLAGS
        -Xcompiler=-march=native
        -Xcompiler=-mtune=native
    )
    list(
        APPEND
        HASE_CLANG_CUDA_FLAGS
        -Xarch_host=-march=native
        -Xarch_host=-mtune=native
    )
endif()
if(HASE_CUDA_KEEP_FILES)
    set(NVCC_TMP_DIR "${PROJECT_BINARY_DIR}/nvcc_tmp")
    file(MAKE_DIRECTORY "${NVCC_TMP_DIR}")
    list(APPEND HASE_NVCC_FLAGS --keep --keep-dir "${NVCC_TMP_DIR}")
    list(APPEND HASE_CLANG_CUDA_FLAGS -save-temps)
endif()

cmake_print_variables(HASE_NVCC_FLAGS)
cmake_print_variables(HASE_CLANG_CUDA_FLAGS)

# C++ compiler Flags
set(HASE_GCC_FLAGS -Wall -Wextra)
set(HASE_CLANG_FLAGS -Wall -Wextra)

if(HASE_BUILD_RELEASE)
    list(APPEND HASE_GCC_FLAGS -march=native -mtune=native)
    list(APPEND HASE_CLANG_FLAGS -march=native -mtune=native)
endif()

# core library with all sources
add_library(hase_core INTERFACE)
add_library(hase::core ALIAS hase_core)

target_compile_options(
    hase_core
    INTERFACE
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:${HASE_GCC_FLAGS}>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Clang>>:${HASE_CLANG_FLAGS}>
        $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:AppleClang>>:${HASE_CLANG_FLAGS}>
)

if(HASE_BUILD_RELEASE AND alpaka_DEP_HIP)
    target_compile_options(
        hase_core
        INTERFACE
            $<$<COMPILE_LANGUAGE:HIP>:-Xarch_host=-march=native>
            $<$<COMPILE_LANGUAGE:HIP>:-Xarch_host=-mtune=native>
    )
endif()
target_compile_options(
    hase_core
    INTERFACE
        $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:${HASE_NVCC_FLAGS}>
        $<$<COMPILE_LANG_AND_ID:CUDA,Clang>:${HASE_CLANG_CUDA_FLAGS}>
)
#target_sources(hase_core PRIVATE ${HASE_CORE_SOURCES})
target_include_directories(
    hase_core
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Benchmark>
)

set_target_properties(
    hase_core
    PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${HASE_CUDA_ARCHITECTURES}"
)

set(HASE_MPI_PROBE_ALLOWED ON)
if(DISABLE_MPI STREQUAL "AUTO")
    find_program(HASE_MPI_CXX_WRAPPER NAMES mpicxx mpic++ mpiCC)
    if(HASE_MPI_CXX_WRAPPER)
        execute_process(
            COMMAND "${HASE_MPI_CXX_WRAPPER}" --showme:compile
            RESULT_VARIABLE HASE_MPI_SHOWME_RESULT
            OUTPUT_VARIABLE HASE_MPI_COMPILE_FLAGS
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(NOT HASE_MPI_SHOWME_RESULT EQUAL 0)
            execute_process(
                COMMAND "${HASE_MPI_CXX_WRAPPER}" -show
                RESULT_VARIABLE HASE_MPI_SHOWME_RESULT
                OUTPUT_VARIABLE HASE_MPI_COMPILE_FLAGS
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
        endif()
        if(HASE_MPI_SHOWME_RESULT EQUAL 0)
            separate_arguments(
                HASE_MPI_COMPILE_ARGS
                UNIX_COMMAND
                "${HASE_MPI_COMPILE_FLAGS}"
            )
            foreach(HASE_MPI_COMPILE_ARG IN LISTS HASE_MPI_COMPILE_ARGS)
                if(HASE_MPI_COMPILE_ARG MATCHES "^-I(.+)")
                    set(HASE_MPI_INCLUDE_DIR "${CMAKE_MATCH_1}")
                    if(NOT EXISTS "${HASE_MPI_INCLUDE_DIR}")
                        message(
                            WARNING
                            "MPI wrapper '${HASE_MPI_CXX_WRAPPER}' references missing include directory "
                            "'${HASE_MPI_INCLUDE_DIR}'. MPI support will be disabled."
                        )
                        set(HASE_MPI_PROBE_ALLOWED OFF)
                    endif()
                endif()
            endforeach()
        endif()
    endif()
    if(HASE_MPI_PROBE_ALLOWED)
        find_package(MPI COMPONENTS CXX QUIET)
    endif()
    if(NOT MPI_FOUND)
        message(
            STATUS
            "MPI is not found on your system - compilation will proceed without MPI Support."
        )
    endif()
elseif(DISABLE_MPI STREQUAL "OFF")
    find_package(MPI COMPONENTS CXX REQUIRED)
endif()
if(MPI_FOUND)
    target_compile_definitions(hase_core INTERFACE MPI_FOUND)
    target_link_libraries(hase_core INTERFACE MPI::MPI_CXX)
else()
    message(STATUS "MPI disabled.")
    target_compile_definitions(hase_core INTERFACE DISABLE_MPI)
endif()

target_link_libraries(hase_core INTERFACE alpaka::alpaka)
if(TARGET CUDA::cudart)
    target_link_libraries(hase_core INTERFACE CUDA::cudart)
endif()

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
    "${CMAKE_CURRENT_SOURCE_DIR}/src/alpakaUtils/backendNamesLibrary.cpp"
)
add_library(hase STATIC ${HASE_RUNTIME_SOURCES})
add_library(hase::hase ALIAS hase)
set_target_properties(hase PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(hase PUBLIC hase::core)
if(HASE_BENCHMARK)
    target_sources(
        hase
        PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Benchmark/benchmark.cpp"
    )
    target_compile_definitions(hase_core INTERFACE HASE_ENABLE_BENCHMARK)
endif()
alpaka_finalize(hase)

include(cmake/AlpakaBackendNames.cmake)
# @TODO think about not exposing this at all for the alpaka release
if(HASE_BUILD_PhiAse)
    add_executable(calcPhiASE "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
    set_target_properties(
        calcPhiASE
        PROPERTIES
            BUILD_RPATH_USE_ORIGIN ON
            CUDA_SEPARABLE_COMPILATION ON
            INSTALL_RPATH_USE_LINK_PATH ON
    )
    target_link_libraries(calcPhiASE PRIVATE hase::hase)
    alpaka_finalize(calcPhiASE)
endif()

# Testing
if(HASE_TESTING)
    enable_testing()
    add_subdirectory(tests)
    # Python
endif()
if(HASE_ENABLE_PYTHON)
    include(cmake/pythonBindings.cmake)
endif()
