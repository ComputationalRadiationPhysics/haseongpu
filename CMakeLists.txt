cmake_minimum_required(VERSION 3.24)
project(HASEonGPU LANGUAGES CUDA CXX C)

include(CMakePrintHelpers)

# dependencies
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30")
  cmake_policy(SET CMP0167 NEW)
endif()

find_package(Boost 1.70.0 REQUIRED COMPONENTS program_options)
find_package(Threads REQUIRED)
find_package(CUDAToolkit REQUIRED COMPONENTS cudart)

# build options
option(HASE_CUDA_FLUSHTOZERO "Set flush to zero for GPU" OFF)
option(HASE_CUDA_FASTMATH "Enable fast-math" ON)
option(HASE_CUDA_SHOW_REGISTER "Show kernel registers and create PTX" OFF)
option(HASE_CUDA_KEEP_FILES
        "Keep all intermediate files that are generated during internal compilation steps (folder: nvcc_tmp)"
        OFF
)
option(HASE_CUDA_SHOW_CODELINES "Show kernel lines in cuda-gdb and cuda-memcheck" OFF)
option(DISABLE_MPI "toggle that allows compilation of haseongpu without requiring the MPI dependency" OFF)
option(HASE_TESTING "Build unit tests" OFF)

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
    message(WARNING
      "Your ${PATH_VARIABLE} seems to contain CUDA includes."
      "This might result in the following error:\n"
      "'error: kernel launches from templates are not allowed in system files'\n"
      "To fix the problem, remove the CUDA includes from the path."
      "(A better place for these includes might be your CPATH)\n"
      "source: https://code.google.com/p/thrust/issues/detail?id=359#c5 \n"
      "issue:  https://github.com/ComputationalRadiationPhysics/haseongpu/issues/26"
      "CPATH:  https://gcc.gnu.org/onlinedocs/cpp/Environment-Variables.html \n")
  endif()
endfunction()

warn_path_cuda("CPLUS_INCLUDE_PATH")
warn_path_cuda("C_INCLUDE_PATH")

# build type (debug, release)
option(HASE_RELEASE "Build release version, disables all runtime asserts" OFF)
if(HASE_RELEASE)
  message(STATUS "Release version")
  set(CMAKE_BUILD_TYPE Release)
else()
  message(STATUS "Debug version")
  set(CMAKE_BUILD_TYPE Debug)
endif()


# core library with all sources
add_library(hase_core STATIC)

file(GLOB_RECURSE HASE_CORE_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu"
)
list(REMOVE_ITEM HASE_CORE_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cu"
)
target_sources(hase_core PRIVATE ${HASE_CORE_SOURCES})

target_include_directories(hase_core PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

set_target_properties(hase_core PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "${HASE_CUDA_ARCHITECTURES}"
)

target_compile_options(hase_core PRIVATE
        "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:${HASE_NVCC_FLAGS}>"
        "$<$<COMPILE_LANG_AND_ID:CXX,GNU>:${HASE_GCC_FLAGS}>"
        "$<$<COMPILE_LANG_AND_ID:CXX,Clang>:${HASE_CLANG_FLAGS}>"
)

target_link_libraries(hase_core PUBLIC
        Boost::program_options
        Threads::Threads
        CUDA::cudart
)

# enable/disable mpi dependency based on cmake option
if(NOT DISABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS C)
  target_link_libraries(hase_core PUBLIC MPI::MPI_C)
  target_compile_definitions(hase_core PUBLIC OMPI_SKIP_MPICXX)
else()
  target_compile_definitions(hase_core PUBLIC DISABLE_MPI)
endif()
# enable Testing
if(HASE_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()
# (keep your CUDA props / flags / link libs on hase_core as discussed)


# --- main executable: main.cu + core lib ---
add_executable(calcPhiASE "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cu")

target_link_libraries(calcPhiASE PRIVATE hase_core)