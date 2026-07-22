cmake_minimum_required(VERSION 3.24)

# Keep the public C++ version header as the single source of truth and read it
# before project() so CMake package metadata cannot drift from the installed API.
file(
    STRINGS
    "${CMAKE_CURRENT_LIST_DIR}/include/hase/version.hpp"
    HASEONGPU_VERSION_MAJOR_HPP
    REGEX "#define HASEONGPU_VERSION_MAJOR "
)
file(
    STRINGS
    "${CMAKE_CURRENT_LIST_DIR}/include/hase/version.hpp"
    HASEONGPU_VERSION_MINOR_HPP
    REGEX "#define HASEONGPU_VERSION_MINOR "
)
file(
    STRINGS
    "${CMAKE_CURRENT_LIST_DIR}/include/hase/version.hpp"
    HASEONGPU_VERSION_PATCH_HPP
    REGEX "#define HASEONGPU_VERSION_PATCH "
)

string(
    REGEX MATCH
    "([0-9]+)"
    HASEONGPU_VERSION_MAJOR
    "${HASEONGPU_VERSION_MAJOR_HPP}"
)
string(
    REGEX MATCH
    "([0-9]+)"
    HASEONGPU_VERSION_MINOR
    "${HASEONGPU_VERSION_MINOR_HPP}"
)
string(
    REGEX MATCH
    "([0-9]+)"
    HASEONGPU_VERSION_PATCH
    "${HASEONGPU_VERSION_PATCH_HPP}"
)

project(
    HASEonGPU
    VERSION
        ${HASEONGPU_VERSION_MAJOR}.${HASEONGPU_VERSION_MINOR}.${HASEONGPU_VERSION_PATCH}
)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
include(CMakePrintHelpers)
include(GNUInstallDirs)

# dependencies
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30")
    cmake_policy(SET CMP0167 NEW)
endif()

find_package(Threads REQUIRED)

# A thin frontend must not mistake scikit-build's temporary compiler/prefix
# cache entries for deliberate changes to an existing native runtime. Forward
# only definitions the caller explicitly supplied through CMAKE_ARGS (the same
# interface used by hase-configure).
set(HASE_INITIAL_CACHE_VARIABLES "")
if(DEFINED ENV{CMAKE_ARGS} AND NOT "$ENV{CMAKE_ARGS}" STREQUAL "")
    separate_arguments(
        HASE_EXPLICIT_CMAKE_ARGUMENTS
        NATIVE_COMMAND
        "$ENV{CMAKE_ARGS}"
    )
    foreach(HASE_EXPLICIT_CMAKE_ARGUMENT IN LISTS HASE_EXPLICIT_CMAKE_ARGUMENTS)
        if(HASE_EXPLICIT_CMAKE_ARGUMENT MATCHES "^-D([^:=]+)(:[^=]+)?=")
            list(APPEND HASE_INITIAL_CACHE_VARIABLES "${CMAKE_MATCH_1}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES HASE_INITIAL_CACHE_VARIABLES)
endif()

# build options
option(
    HASE_BUILD_RELEASE
    "Build release version; forces Release build type and release optimization options."
    ON
)
option(
    HASE_NATIVE_OPTIMIZATIONS
    "Enable -march=native/-mtune=native for host-tuned source builds. Disable for redistributable binaries."
    OFF
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
option(
    HASE_ENABLE_PYTHON
    "Generate/install the HASEonGPU Python frontend and runtime metadata"
    ON
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
option(HASE_BENCHMARK "Enable scoped PhiASE benchmark CSV output" OFF)
option(
    HASE_FORWARD_LOGGING
    "Forward captured calcPhiASE stdout and stderr through the Python openPMD launcher."
    OFF
)
set(HASE_RUNTIME_DIR
    "${CMAKE_SOURCE_DIR}/build"
    CACHE PATH
    "Shared HASE native runtime build used by the installed Python frontend"
)
if(NOT HASE_RUNTIME_DIR)
    set(HASE_RUNTIME_DIR
        "${CMAKE_SOURCE_DIR}/build"
        CACHE PATH
        "Shared HASE native runtime build"
        FORCE
    )
endif()
option(HASE_BUILD_RUNTIME "Build the native HASE runtime in this CMake tree" ON)
if(HASE_BUILD_RUNTIME)
    set(HASE_RUNTIME_DIR
        "${CMAKE_BINARY_DIR}"
        CACHE PATH
        "Native HASE runtime build"
        FORCE
    )
endif()

if(NOT HASE_BUILD_RUNTIME)
    if(NOT HASE_ENABLE_PYTHON)
        message(
            FATAL_ERROR
            "HASE_BUILD_RUNTIME=OFF requires HASE_ENABLE_PYTHON=ON because this mode only exports the Python frontend."
        )
    endif()
    get_filename_component(
        HASE_RUNTIME_DIR
        "${HASE_RUNTIME_DIR}"
        ABSOLUTE
        BASE_DIR "${CMAKE_SOURCE_DIR}"
    )
    set(HASE_RUNTIME_DIR
        "${HASE_RUNTIME_DIR}"
        CACHE PATH
        "External HASE native runtime root recorded for the installed Python frontend"
        FORCE
    )
    if(HASE_RUNTIME_DIR STREQUAL CMAKE_BINARY_DIR)
        message(
            FATAL_ERROR
            "HASE_BUILD_RUNTIME=OFF requires HASE_RUNTIME_DIR to differ from the Python frontend build directory."
        )
    endif()

    set(HASE_RUNTIME_FORWARDED_VARIABLES "")
    foreach(HASE_CACHE_VARIABLE IN LISTS HASE_INITIAL_CACHE_VARIABLES)
        if(
            HASE_CACHE_VARIABLE MATCHES "^HASE_"
            OR HASE_CACHE_VARIABLE MATCHES "^alpaka_"
            OR HASE_CACHE_VARIABLE
                MATCHES
                "^CMAKE_(BUILD_TYPE|INTERPROCEDURAL_OPTIMIZATION|BUILD_RPATH_USE_ORIGIN|INSTALL_RPATH|INSTALL_RPATH_USE_LINK_PATH|PREFIX_PATH|TOOLCHAIN_FILE)$"
            OR HASE_CACHE_VARIABLE
                MATCHES
                "^CMAKE_(C|CXX|CUDA|HIP)_(COMPILER|COMPILER_LAUNCHER|FLAGS|FLAGS_[A-Z]+|ARCHITECTURES)$"
            OR HASE_CACHE_VARIABLE MATCHES "^(ADIOS2|HDF5|CUDAToolkit)_DIR$"
            OR HASE_CACHE_VARIABLE MATCHES "^(openPMD|CUDAToolkit)_ROOT$"
            OR HASE_CACHE_VARIABLE MATCHES "^(BUILD_TESTING|DISABLE_MPI)$"
        )
            list(
                APPEND
                HASE_RUNTIME_FORWARDED_VARIABLES
                "${HASE_CACHE_VARIABLE}"
            )
        endif()
    endforeach()
    list(
        REMOVE_ITEM
        HASE_RUNTIME_FORWARDED_VARIABLES
        HASE_BUILD_RUNTIME
        HASE_RUNTIME_DIR
        HASE_FRONTEND_EXPLICIT_CACHE_VARIABLES
    )
    list(REMOVE_DUPLICATES HASE_RUNTIME_FORWARDED_VARIABLES)
    list(SORT HASE_RUNTIME_FORWARDED_VARIABLES)

    set(HASE_RUNTIME_NEEDS_CONFIGURE FALSE)
    if(EXISTS "${HASE_RUNTIME_DIR}/CMakeCache.txt")
        if(HASE_RUNTIME_FORWARDED_VARIABLES)
            load_cache(
                "${HASE_RUNTIME_DIR}"
                READ_WITH_PREFIX HASE_RUNTIME_CACHE_
                ${HASE_RUNTIME_FORWARDED_VARIABLES}
            )
        endif()
        foreach(HASE_CACHE_VARIABLE IN LISTS HASE_RUNTIME_FORWARDED_VARIABLES)
            set(HASE_RUNTIME_CACHE_VARIABLE
                "HASE_RUNTIME_CACHE_${HASE_CACHE_VARIABLE}"
            )
            set(HASE_RUNTIME_OPTION_CHANGED FALSE)
            if(DEFINED ${HASE_RUNTIME_CACHE_VARIABLE})
                if(
                    NOT "${${HASE_RUNTIME_CACHE_VARIABLE}}"
                        STREQUAL
                        "${${HASE_CACHE_VARIABLE}}"
                )
                    set(HASE_RUNTIME_OPTION_CHANGED TRUE)
                endif()
            elseif(NOT "${${HASE_CACHE_VARIABLE}}" STREQUAL "")
                set(HASE_RUNTIME_OPTION_CHANGED TRUE)
            endif()
            if(HASE_RUNTIME_OPTION_CHANGED)
                message(
                    STATUS
                    "Shared runtime option changed: ${HASE_CACHE_VARIABLE}='${${HASE_CACHE_VARIABLE}}'"
                )
                set(HASE_RUNTIME_NEEDS_CONFIGURE TRUE)
                break()
            endif()
        endforeach()
    else()
        set(HASE_RUNTIME_NEEDS_CONFIGURE TRUE)
    endif()

    if(HASE_RUNTIME_NEEDS_CONFIGURE)
        set(HASE_RUNTIME_INITIAL_CACHE
            "${CMAKE_BINARY_DIR}/hase-runtime-initial-cache.cmake"
        )
        file(
            WRITE
            "${HASE_RUNTIME_INITIAL_CACHE}"
            "# Generated by the HASE Python frontend configuration.\n"
        )
        foreach(HASE_CACHE_VARIABLE IN LISTS HASE_RUNTIME_FORWARDED_VARIABLES)
            get_property(
                HASE_CACHE_TYPE
                CACHE "${HASE_CACHE_VARIABLE}"
                PROPERTY TYPE
            )
            if(
                NOT HASE_CACHE_TYPE
                OR HASE_CACHE_TYPE STREQUAL "UNINITIALIZED"
                OR HASE_CACHE_TYPE STREQUAL "INTERNAL"
            )
                set(HASE_CACHE_TYPE STRING)
            endif()
            file(
                APPEND
                "${HASE_RUNTIME_INITIAL_CACHE}"
                "set(${HASE_CACHE_VARIABLE} [==[${${HASE_CACHE_VARIABLE}}]==] CACHE ${HASE_CACHE_TYPE} \"Forwarded by the HASE Python frontend\" FORCE)\n"
            )
        endforeach()
        file(
            APPEND
            "${HASE_RUNTIME_INITIAL_CACHE}"
            "set(HASE_BUILD_RUNTIME ON CACHE BOOL \"Build the shared native HASE runtime\" FORCE)\n"
            "set(HASE_ENABLE_PYTHON ON CACHE BOOL \"Build matching Python runtime helpers\" FORCE)\n"
        )

        set(HASE_RUNTIME_CONFIGURE_COMMAND
            "${CMAKE_COMMAND}"
            -S
            "${CMAKE_SOURCE_DIR}"
            -B
            "${HASE_RUNTIME_DIR}"
            -C
            "${HASE_RUNTIME_INITIAL_CACHE}"
        )
        if(NOT EXISTS "${HASE_RUNTIME_DIR}/CMakeCache.txt")
            list(APPEND HASE_RUNTIME_CONFIGURE_COMMAND -G "${CMAKE_GENERATOR}")
            if(CMAKE_GENERATOR_PLATFORM)
                list(
                    APPEND
                    HASE_RUNTIME_CONFIGURE_COMMAND
                    -A
                    "${CMAKE_GENERATOR_PLATFORM}"
                )
            endif()
            if(CMAKE_GENERATOR_TOOLSET)
                list(
                    APPEND
                    HASE_RUNTIME_CONFIGURE_COMMAND
                    -T
                    "${CMAKE_GENERATOR_TOOLSET}"
                )
            endif()
        endif()
        message(
            STATUS
            "Configuring shared native HASE runtime: ${HASE_RUNTIME_DIR}"
        )
        execute_process(
            COMMAND ${HASE_RUNTIME_CONFIGURE_COMMAND}
            RESULT_VARIABLE HASE_RUNTIME_CONFIGURE_RESULT
            COMMAND_ECHO STDOUT
        )
        if(NOT HASE_RUNTIME_CONFIGURE_RESULT EQUAL 0)
            message(
                FATAL_ERROR
                "Failed to configure the shared native HASE runtime."
            )
        endif()
    endif()

    message(STATUS "Building shared native HASE runtime: ${HASE_RUNTIME_DIR}")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" --build "${HASE_RUNTIME_DIR}"
        RESULT_VARIABLE HASE_RUNTIME_BUILD_RESULT
        COMMAND_ECHO STDOUT
    )
    if(NOT HASE_RUNTIME_BUILD_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to build the shared native HASE runtime.")
    endif()
    if(
        NOT EXISTS "${HASE_RUNTIME_DIR}/calcPhiASE"
        AND NOT EXISTS "${HASE_RUNTIME_DIR}/python/pyInclude/_runtime/calcPhiASE"
    )
        message(
            FATAL_ERROR
            "The shared runtime '${HASE_RUNTIME_DIR}' did not produce calcPhiASE."
        )
    endif()

    message(
        STATUS
        "Configuring Python frontend for runtime: ${HASE_RUNTIME_DIR}"
    )
    include(cmake/pythonRuntime.cmake)
    return()
endif()

set(HASE_CUDA_ARCHITECTURES native CACHE STRING "CUDA Architectures")
if(
    DEFINED CMAKE_CUDA_ARCHITECTURES
    AND HASE_CUDA_ARCHITECTURES STREQUAL "native"
)
    set(HASE_CUDA_ARCHITECTURES
        "${CMAKE_CUDA_ARCHITECTURES}"
        CACHE STRING
        "CUDA Architectures"
        FORCE
    )
endif()

if(HASE_BUILD_RELEASE)
    message(
        STATUS
        "HASE_BUILD_RELEASE=ON: forcing CMAKE_BUILD_TYPE=Release and enabling release optimization options. "
        "Disable HASE_BUILD_RELEASE to keep user-defined CMAKE_BUILD_TYPE or optimization settings."
    )
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(HASE_NATIVE_OPTIMIZATIONS)
    message(
        STATUS
        "HASE_NATIVE_OPTIMIZATIONS=ON: enabling host-specific CPU tuning. "
        "Disable this option for redistributable binaries."
    )
endif()

# search for backends (dependent) -- important that this runs before making alpaka available
if(NOT HASE_SELECT_BACKEND_ALPAKA)
    include(cmake/searchBackends.cmake)
endif()

if(alpaka_DEP_CUDA AND HASE_CUDA_ARCHITECTURES STREQUAL "native")
    # Probe nvidia-smi to check whether a GPU is visible and "native" can work.
    execute_process(
        COMMAND nvidia-smi -L
        RESULT_VARIABLE HASE_NVIDIA_SMI_RESULT
        OUTPUT_QUIET
        ERROR_QUIET
    )

    if(NOT HASE_NVIDIA_SMI_RESULT EQUAL 0)
        # If no GPU is visible, fall back to a static architecture list.
        if(NOT DEFINED CUDAToolkit_VERSION)
            find_package(CUDAToolkit QUIET)
        endif()

        if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL "13.0")
            # CUDA 13 removed Maxwell/Pascal/Volta support, so no 50/60/70.
            set(_hase_default_cuda_archs "75;80;90")
        else()
            set(_hase_default_cuda_archs "60;70;80;90")
        endif()

        message(
            STATUS
            "No visible NVIDIA GPU found; falling back to CUDA architectures "
            "${_hase_default_cuda_archs}."
        )

        set(HASE_CUDA_ARCHITECTURES
            "${_hase_default_cuda_archs}"
            CACHE STRING
            "CUDA Architectures"
            FORCE
        )

        unset(_hase_default_cuda_archs)
    endif()

    unset(HASE_NVIDIA_SMI_RESULT)
endif()

if(alpaka_DEP_CUDA)
    set(CMAKE_CUDA_ARCHITECTURES
        "${HASE_CUDA_ARCHITECTURES}"
        CACHE STRING
        "CUDA Architectures"
        FORCE
    )
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

if(alpaka_DEP_CUDA AND NOT TARGET alpaka::cuda)
    if(NVCC_EXECUTABLE)
        set(HASE_NVCC_EXECUTABLE "${NVCC_EXECUTABLE}")
    else()
        find_program(
            HASE_NVCC_EXECUTABLE
            NAMES nvcc
            HINTS
                "${CMAKE_CUDA_COMPILER}"
                "${CUDAToolkit_ROOT}/bin"
                "$ENV{CUDAToolkit_ROOT}/bin"
                "$ENV{CUDA_HOME}/bin"
                "$ENV{CUDA_PATH}/bin"
                "$ENV{CUDA_ROOT}/bin"
            DOC "Path to NVIDIA nvcc compiler"
        )
    endif()

    if(HASE_NVCC_EXECUTABLE)
        message(
            WARNING
            "[alpaka CUDA] nvcc exists, but CMake did not accept CUDA as a language.\n"
            "[alpaka CUDA] This usually means a configuration/compatibility problem,\n"
            "[alpaka CUDA] e.g. CUDA does not support the host compiler.\n"
            "[alpaka CUDA] nvcc found at: ${HASE_NVCC_EXECUTABLE}\n"
            "[alpaka CUDA] Check ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log"
        )
    endif()
endif()
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
if(HASE_BUILD_RELEASE AND HASE_NATIVE_OPTIMIZATIONS AND alpaka_DEP_CUDA)
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

if(HASE_BUILD_RELEASE AND HASE_NATIVE_OPTIMIZATIONS)
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

if(HASE_BUILD_RELEASE AND HASE_NATIVE_OPTIMIZATIONS AND alpaka_DEP_HIP)
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

set(HASE_MPI_ENABLED OFF)
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
                            "MPI wrapper ${HASE_MPI_CXX_WRAPPER} references missing include directory "
                            "${HASE_MPI_INCLUDE_DIR}. MPI support will be disabled."
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
    if(MPI_FOUND AND HASE_MPI_PROBE_ALLOWED)
        set(HASE_MPI_ENABLED ON)
    else()
        message(
            STATUS
            "MPI is not found on your system - compilation will proceed without MPI Support."
        )
    endif()
elseif(DISABLE_MPI STREQUAL "OFF")
    find_package(MPI COMPONENTS CXX REQUIRED)
    set(HASE_MPI_ENABLED ON)
endif()
if(HASE_MPI_ENABLED)
    target_compile_definitions(hase_core INTERFACE MPI_FOUND)
    target_link_libraries(hase_core INTERFACE MPI::MPI_CXX)
else()
    message(STATUS "MPI disabled.")
    target_compile_definitions(hase_core INTERFACE DISABLE_MPI)
endif()

set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
set(CMAKE_INSTALL_RPATH "$ORIGIN")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set(HASE_BUILD_TREE_RPATH
    "$ORIGIN/lib;$ORIGIN/_deps/adios2-build/lib;$ORIGIN/../../lib;$ORIGIN/../../_deps/adios2-build/lib"
)
set(HASE_INSTALL_RPATH "$ORIGIN;$ORIGIN/../lib")
set(HASE_INSTALL_LIB_RPATH "$ORIGIN")
set(HASE_OPENPMD_RUNTIME_RPATH
    ""
    CACHE STRING
    "Additional runtime search paths for a non-vendored openPMD-api provider"
)

include(cmake/findOpenPmd.cmake)
set(HASE_OPENPMD_PROVIDER_RPATH ${HASE_OPENPMD_RUNTIME_RPATH})
function(hase_append_openpmd_rpath_from_target target)
    if(NOT TARGET ${target})
        return()
    endif()

    get_target_property(HASE_OPENPMD_PROVIDER_TARGET_TYPE ${target} TYPE)
    if(HASE_OPENPMD_PROVIDER_TARGET_TYPE STREQUAL "INTERFACE_LIBRARY")
        return()
    endif()

    get_target_property(HASE_OPENPMD_PROVIDER_IMPORTED ${target} IMPORTED)
    if(HASE_OPENPMD_PROVIDER_IMPORTED)
        foreach(
            HASE_OPENPMD_PROVIDER_LOCATION_PROPERTY
            IN
            ITEMS
                IMPORTED_LOCATION
                IMPORTED_LOCATION_RELEASE
                IMPORTED_LOCATION_RELWITHDEBINFO
                IMPORTED_LOCATION_MINSIZEREL
                IMPORTED_LOCATION_DEBUG
        )
            get_target_property(
                HASE_OPENPMD_PROVIDER_LOCATION
                ${target}
                ${HASE_OPENPMD_PROVIDER_LOCATION_PROPERTY}
            )
            if(HASE_OPENPMD_PROVIDER_LOCATION)
                get_filename_component(
                    HASE_OPENPMD_PROVIDER_LOCATION_DIR
                    "${HASE_OPENPMD_PROVIDER_LOCATION}"
                    DIRECTORY
                )
                list(
                    APPEND
                    HASE_OPENPMD_PROVIDER_RPATH
                    "${HASE_OPENPMD_PROVIDER_LOCATION_DIR}"
                )
            endif()
        endforeach()
    else()
        list(APPEND HASE_OPENPMD_PROVIDER_RPATH "$<TARGET_FILE_DIR:${target}>")
    endif()
    set(HASE_OPENPMD_PROVIDER_RPATH
        "${HASE_OPENPMD_PROVIDER_RPATH}"
        PARENT_SCOPE
    )
endfunction()

hase_append_openpmd_rpath_from_target(openPMD::openPMD)
hase_append_openpmd_rpath_from_target(openPMD)
if(HASE_OPENPMD_PROVIDER_RPATH)
    list(REMOVE_DUPLICATES HASE_OPENPMD_PROVIDER_RPATH)
    list(APPEND HASE_INSTALL_RPATH ${HASE_OPENPMD_PROVIDER_RPATH})
    list(APPEND HASE_BUILD_TREE_RPATH ${HASE_OPENPMD_PROVIDER_RPATH})
endif()
foreach(HASE_BUILD_RPATH_TARGET IN ITEMS openPMD hdf5-shared)
    if(TARGET ${HASE_BUILD_RPATH_TARGET})
        list(
            APPEND
            HASE_BUILD_TREE_RPATH
            "$<TARGET_FILE_DIR:${HASE_BUILD_RPATH_TARGET}>"
        )
    endif()
endforeach()
target_link_libraries(hase_core INTERFACE openPMD::openPMD)
target_compile_definitions(
    hase_core
    INTERFACE
        HASE_ENABLE_OPENPMD
        HASE_OPENPMD_FILE_EXTENSION=\"${HASE_OPENPMD_FILE_EXTENSION}\"
)
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
    "${CMAKE_CURRENT_SOURCE_DIR}/src/openpmd_main.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/openpmd/backendProbe.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/alpakaUtils/backendNamesLibrary.cpp"
)
add_library(hase STATIC ${HASE_RUNTIME_SOURCES})
add_library(hase::hase ALIAS hase)
set_target_properties(
    hase
    PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${HASE_CUDA_ARCHITECTURES}"
)
target_link_libraries(hase PUBLIC hase::core)
if(HASE_BENCHMARK)
    target_sources(
        hase
        PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Benchmark/benchmark.cpp"
    )
    target_compile_definitions(hase_core INTERFACE HASE_ENABLE_BENCHMARK)
endif()
if(HASE_FORWARD_LOGGING)
    target_compile_definitions(hase_core INTERFACE HASE_FORWARD_LOGGING)
endif()
alpaka_finalize(hase)

include(cmake/AlpakaBackendNames.cmake)
include(cmake/OpenPmdBackendProbe.cmake)
add_executable(calcPhiASE "${CMAKE_CURRENT_SOURCE_DIR}/src/openpmd_main.cpp")
set_target_properties(
    calcPhiASE
    PROPERTIES
        BUILD_RPATH_USE_ORIGIN ON
        BUILD_RPATH "${HASE_BUILD_TREE_RPATH}"
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${HASE_CUDA_ARCHITECTURES}"
        INSTALL_RPATH "${HASE_INSTALL_RPATH}"
        INSTALL_RPATH_USE_LINK_PATH ON
)
target_link_libraries(calcPhiASE PRIVATE hase::hase)
foreach(HASE_ADIOS2_RUNTIME_MODULE IN ITEMS cmselect cmsockets)
    if(TARGET ${HASE_ADIOS2_RUNTIME_MODULE})
        add_dependencies(calcPhiASE ${HASE_ADIOS2_RUNTIME_MODULE})
    endif()
endforeach()
alpaka_finalize(calcPhiASE)

if(HASE_ENABLE_PYTHON)
    include(cmake/pythonRuntime.cmake)
endif()

# Testing
if(HASE_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
