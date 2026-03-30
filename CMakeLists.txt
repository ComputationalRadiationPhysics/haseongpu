cmake_minimum_required(VERSION 3.24)
project(HASEonGPU LANGUAGES CUDA CXX C)
set(CMAKE_CXX_STANDARD 20)
include(CMakePrintHelpers)

# dependencies
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30")
    cmake_policy(SET CMP0167 NEW)
endif()

find_package(Threads REQUIRED)
find_package(CUDAToolkit REQUIRED COMPONENTS cudart)

# build options
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
option(
    DISABLE_MPI
    "toggle that allows compilation of haseongpu without requiring the MPI dependency"
    OFF
)
option(HASE_TESTING "Build unit tests" OFF)
option(HASE_ENABLE_PYTHON "Build Python bindings" OFF)

set(HASE_CUDA_ARCHITECTURES native CACHE STRING "CUDA Architectures")

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
add_library(hase_core STATIC)

file(
    GLOB_RECURSE HASE_CORE_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu"
)
set_target_properties(hase_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

file(
    GLOB_RECURSE HASE_CORE_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu"
)
list(REMOVE_ITEM HASE_CORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc")
target_sources(hase_core PRIVATE ${HASE_CORE_SOURCES})

target_include_directories(
    hase_core
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

set_target_properties(
    hase_core
    PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${HASE_CUDA_ARCHITECTURES}"
)

target_compile_options(
    hase_core
    PRIVATE
        "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:${HASE_NVCC_FLAGS}>"
        "$<$<COMPILE_LANG_AND_ID:CXX,GNU>:${HASE_GCC_FLAGS}>"
        "$<$<COMPILE_LANG_AND_ID:CXX,Clang>:${HASE_CLANG_FLAGS}>"
)

target_link_libraries(
    hase_core
    PUBLIC
        Threads::Threads
        CUDA::cudart
)



# enable/disable mpi dependency based on cmake option
if(NOT DISABLE_MPI)
    find_package(MPI QUIET)
    if(MPI_FOUND)
        target_compile_definitions(hase_core PUBLIC MPI_FOUND)
        target_link_libraries(hase_core PUBLIC MPI::MPI_CXX)
    else()
        message(STATUS "MPI is not found on your system. You can compile HASEonGPU without MPI by using DISABLE_MPI=ON")
    endif()
    option(
        USE_GRAYBAT
        "Use Graph Approach for Highly Generic Communication Schemes Based on Adaptive Topologies (Graybat)"
        OFF
    )
    if(USE_GRAYBAT)
        include(FetchContent)

        FetchContent_Declare(
            graybat
            GIT_REPOSITORY
                https://github.com/ComputationalRadiationPhysics/graybat.git
            GIT_TAG dev-c++11
            GIT_SHALLOW TRUE
        )

        # Fetch sources, but do NOT add_subdirectory() / MakeAvailable()
        FetchContent_GetProperties(graybat)
        if(NOT graybat_POPULATED)
            FetchContent_Populate(graybat)
        endif()

        # Graybat ships graybatConfig.cmake in its repo root, so we can use it directly
        find_package(
            graybat
            1.0.0
            REQUIRED
            CONFIG
            PATHS "${graybat_SOURCE_DIR}"
        )

        # Example integration with your target
        target_include_directories(
            hase_core
            SYSTEM
            PRIVATE ${graybat_INCLUDE_DIRS}
        )
        target_sources(hase_core PRIVATE ${graybat_GENERATED_FILES})

        target_compile_definitions(hase_core PRIVATE ${graybat_DEFINITIONS})

        target_link_libraries(hase_core PRIVATE ${graybat_LIBRARIES})
        target_compile_definitions(hase_core PUBLIC USE_GRAYBAT)
    endif()
else()
    message(STATUS "MPI explicitly disabled via DISABLE_MPI=ON")
    target_compile_definitions(hase_core PUBLIC DISABLE_MPI)
endif()

# enable Testing
if(HASE_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
# (keep your CUDA props / flags / link libs on hase_core as discussed)

# --- main executable: main.cu + core lib ---
add_executable(calcPhiASE "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc")

target_link_libraries(calcPhiASE PRIVATE hase_core)

if(HASE_ENABLE_PYTHON)
    set_target_properties(hase_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
    find_package(pybind11 REQUIRED)
    pybind11_add_module(HASEonGPU HASEonGPU_Bindings/module.cpp)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
    set_target_properties(
        HASEonGPU
        PROPERTIES BUILD_RPATH_USE_ORIGIN ON INSTALL_RPATH_USE_LINK_PATH ON
    )
    target_link_libraries(HASEonGPU PRIVATE hase_core)
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
