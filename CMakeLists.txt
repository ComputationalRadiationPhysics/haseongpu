cmake_minimum_required(VERSION 3.24)
project(HASEonGPU LANGUAGES CUDA CXX C)

include(CMakePrintHelpers)

# dependencies
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30")
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost 1.70.0 REQUIRED COMPONENTS program_options)
find_package(MPI REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)
find_package(CUDAToolkit REQUIRED COMPONENTS cudart)


# build options
set(HASE_CUDA_FLUSHTOZERO OFF CACHE BOOL "Set flush to zero for GPU")
set(HASE_CUDA_FASTMATH ON CACHE BOOL "Enable fast-math" )
set(HASE_CUDA_SHOW_REGISTER OFF CACHE BOOL "Show kernel registers and create PTX")
set(HASE_CUDA_KEEP_FILES OFF CACHE BOOL "Keep all intermediate files that are generated during internal compilation steps (folder: nvcc_tmp)")
set(HASE_CUDA_SHOW_CODELINES OFF CACHE BOOL "Show kernel lines in cuda-gdb and cuda-memcheck")
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
  make_directory("${NVCC_TMP_DIR}")
  list(APPEND HASE_NVCC_FLAGS --keep-dir "${NVCC_TMP_DIR}")
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

# calcPhiASE executable target
add_executable(calcPhiASE
  src/calc_phi_ase.cu
  src/calc_phi_ase_mpi.cc
  src/calc_phi_ase_threaded.cu
  src/calc_sample_gain_sum.cu
  src/cuda_utils.cu
  src/for_loops_clad.cu
  src/geometry.cu
  src/importance_sampling.cu
  src/interpolation.cu
  src/logging.cu
  src/main.cu
  src/map_rays_to_prisms.cu
  src/mesh.cu
  src/mt19937ar.cu
  src/parser.cu
  src/progressbar.cu
  src/propagate_ray.cu
  src/ray_histogram.cu
  src/reflection.cu
  src/write_matlab_output.cu
  src/write_to_file.cu
  src/write_to_vtk.cu
)

target_include_directories(calcPhiASE PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
  
set_target_properties(calcPhiASE PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON
  CUDA_ARCHITECTURES "${HASE_CUDA_ARCHITECTURES}"
)

target_compile_options(calcPhiASE PRIVATE
  "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:${HASE_NVCC_FLAGS}>"
  "$<$<COMPILE_LANG_AND_ID:CXX,GNU>:${HASE_GCC_FLAGS}>"
  "$<$<COMPILE_LANG_AND_ID:CXX,Clang>:${HASE_CLANG_FLAGS}>"
)

target_link_libraries(calcPhiASE PRIVATE
  Boost::program_options
  MPI::MPI_C
  Threads::Threads
  CUDA::cudart
)


# examples target
add_custom_target(
  examples
  ALL
  COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/example/matlab_example/bin
  COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/example/c_example/bin
  COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/example/c_example/output
  COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/calcPhiASE ${CMAKE_CURRENT_SOURCE_DIR}/example/matlab_example/bin
  COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/calcPhiASE ${CMAKE_CURRENT_SOURCE_DIR}/example/c_example/bin
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/calcPhiASE.m ${CMAKE_CURRENT_SOURCE_DIR}/example/matlab_example
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/calcPhiASE.m ${CMAKE_CURRENT_SOURCE_DIR}/.
  DEPENDS calcPhiASE
)
