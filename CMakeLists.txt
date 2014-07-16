project(HASEonGPU)
cmake_minimum_required(VERSION 2.8.10.1)

# helper for libs and packages
set(CMAKE_PREFIX_PATH "/usr/lib/x86_64-linux-gnu/"
  "$ENV{MPI_ROOT}" "$ENV{CUDA_ROOT}" "$ENV{BOOST_ROOT}")

###############################################################################
# CUDA
###############################################################################
find_package(CUDA 5.0 REQUIRED)

if(CUDA_VERSION VERSION_LESS 5.5)
  message(STATUS "CUDA Toolkit < 5.5 detected. We strongly recommend to still "
    "use CUDA 5.5+ drivers (319.82 or higher)!")
endif(CUDA_VERSION VERSION_LESS 5.5)

set(CUDA_ARCH sm_20 CACHE STRING "Set GPU architecture")
string(COMPARE EQUAL ${CUDA_ARCH} "sm_10" IS_CUDA_ARCH_UNSUPPORTED)
string(COMPARE EQUAL ${CUDA_ARCH} "sm_11" IS_CUDA_ARCH_UNSUPPORTED)
string(COMPARE EQUAL ${CUDA_ARCH} "sm_12" IS_CUDA_ARCH_UNSUPPORTED)
string(COMPARE EQUAL ${CUDA_ARCH} "sm_13" IS_CUDA_ARCH_UNSUPPORTED)

if(IS_CUDA_ARCH_UNSUPPORTED)
  message(FATAL_ERROR "Unsupported CUDA architecture ${CUDA_ARCH} specified. "
    "SM 2.0 or higher is required.")
endif(IS_CUDA_ARCH_UNSUPPORTED)

set(CUDA_FTZ "--ftz=false" CACHE STRING "Set flush to zero for GPU")

set(CUDA_MATH --use_fast_math CACHE STRING "Enable fast-math" )
option(CUDA_SHOW_REGISTER "Show kernel registers and create PTX" OFF)
option(CUDA_KEEP_FILES "Keep all intermediate files that are generated during internal compilation steps (folder: nvcc_tmp)" OFF)
option(CUDA_SHOW_CODELINES "Show kernel lines in cuda-gdb and cuda-memcheck" OFF)

if(CUDA_SHOW_CODELINES)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" --source-in-ptx -Xcompiler -rdynamic -lineinfo)
  set(CUDA_KEEP_FILES ON CACHE BOOL "activate keep files" FORCE)
endif(CUDA_SHOW_CODELINES)

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${nvcc_flags} -arch=${CUDA_ARCH} ${CUDA_MATH} ${CUDA_FTZ})
if(CUDA_SHOW_REGISTER)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" -Xptxas=-v)
endif(CUDA_SHOW_REGISTER)

if(CUDA_KEEP_FILES)
  make_directory("${PROJECT_BINARY_DIR}/nvcc_tmp")
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" --keep --keep-dir "${PROJECT_BINARY_DIR}/nvcc_tmp")
endif(CUDA_KEEP_FILES)


###############################################################################
# Boost
###############################################################################
find_package(Boost 1.48.0 REQUIRED COMPONENTS program_options filesystem)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})

# nvcc + boost 1.55 work around
#if(Boost_VERSION EQUAL 105500)
#  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}; \"-DBOOST_NOINLINE=__attribute__((noinline))\" ")
#endif(Boost_VERSION EQUAL 105500)


################################################################################
# Find MPI
################################################################################
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
set(LIBS ${LIBS} ${MPI_C_LIBRARIES})

# bullxmpi fails if it can not find its c++ counter part
if(MPI_CXX_FOUND)
  set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
endif(MPI_CXX_FOUND)


################################################################################
# Find PThreads
################################################################################
find_package(Threads REQUIRED)
set(LIBS ${LIBS} ${CMAKE_THREAD_LIBS_INIT})


################################################################################
# Compiler Flags
################################################################################
# GNU
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Waddress")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Warray-bounds")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wchar-subscripts")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcomment")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wformat")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmain")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmissing-braces")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wparentheses")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsequence-point")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-compare")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-aliasing")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-overflow=1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wswitch")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wtrigraphs")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wuninitialized")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-function")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-label")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-value")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wvolatile-register-var")
  # new warning in gcc 4.8 (flag ignored in previous version)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
  # ICC
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_NO_VARIADIC_TEMPLATES")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_NO_CXX11_VARIADIC_TEMPLATES")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_NO_FENV_H")
  # PGI
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Minform=inform")
endif()


################################################################################
# Build type (debug, release)
################################################################################
option(HASE_RELEASE "Build release version, disables all runtime asserts" OFF)
if(HASE_RELEASE)
  message(STATUS "Release version")

  set(CMAKE_BUILD_TYPE Release)
  add_definitions(-DNDEBUG)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" "-Xcompiler=-pthread")
else(HASE_RELEASE)
  message(STATUS "Debug version")

  set(CMAKE_BUILD_TYPE Debug)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" -g "-Xcompiler=-g,-pthread")
endif(HASE_RELEASE)


###############################################################################
# Executables
###############################################################################

set(HASEonGPU_NAME "calcPhiASE")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB SRCFILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c*")

set(CUDA_SEPARABLE_COMPILATION ON)
cuda_add_executable(${HASEonGPU_NAME} ${SRCFILES})

target_link_libraries(${HASEonGPU_NAME} ${LIBS} ${CUDA_CUDART_LIBRARY})



###############################################################################
# more targets (make examples, make new, make doc)
###############################################################################

add_custom_target(
  examples
  ALL
  COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/example/matlab_example/bin
  COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/example/c_example/bin
  COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/example/c_example/output
  COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/${HASEonGPU_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/example/matlab_example/bin
  COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/${HASEonGPU_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/example/c_example/bin
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/${HASEonGPU_NAME}.m ${CMAKE_CURRENT_SOURCE_DIR}/example/matlab_example
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/src/${HASEonGPU_NAME}.m ${CMAKE_CURRENT_SOURCE_DIR}/.
  )

add_custom_target(
  new
  COMMAND make clean
  COMMAND make
  )

# for later...
#
#find_package(Doxygen)
#if(DOXYGEN_FOUND)
#  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxygen.conf ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
#  add_custom_target(doc
#    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#    COMMENT "Generating API documentation with Doxygen" VERBATIM
#    )
#endif(DOXYGEN_FOUND)
