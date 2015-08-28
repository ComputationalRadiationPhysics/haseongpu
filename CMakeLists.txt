project(HASEonGPU)
cmake_minimum_required(VERSION 3.0.1)

# helper for libs and packages
set(CMAKE_PREFIX_PATH "/usr/lib/x86_64-linux-gnu/"
  "$ENV{MPI_ROOT}" "$ENV{CUDA_ROOT}" "$ENV{BOOST_ROOT}")

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8.2)
  message(FATAL_ERROR "GCC version must be at least 4.8.2")
endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8.2)


###############################################################################
# CUDA
###############################################################################
find_package(CUDA 5.5 REQUIRED)

set(CUDA_ARCH sm_20 CACHE STRING "Set GPU architecture")
# Separate C++ from CUDA file compilation (So we can use our own compiler)
set(CUDA_SEPARABLE_COMPILATION ON)
# Do not propagate CXX flags to the nvcc host compiler (especially -std=c++11)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)

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



################################################################################
# Find MPI
################################################################################
find_package(MPI)
include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
set(LIBS ${LIBS} ${MPI_C_LIBRARIES})

# bullxmpi fails if it can not find its c++ counter part
if(MPI_CXX_FOUND)
  set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
endif(MPI_CXX_FOUND)


###############################################################################
# Boost
###############################################################################
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
## The two paths might need to be set on windows
## the following example paths are for a default installation of the pre-compiled
## libraries from http://boost.teeks99.com/ (does NOT include Boost.MPI)
# set(BOOST_LIBRARYDIR "C:/local/boost_1_58_0/lib64-msvc-12.0")
# set(BOOST_ROOT "C:/local/boost_1_58_0")
find_package(Boost 1.50.0 REQUIRED COMPONENTS program_options system filesystem)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})

find_package(Boost 1.50.0 OPTIONAL_COMPONENTS mpi)
if(Boost_MPI_FOUND)
    # The define can be used in the code to disable MPI features
    add_definitions(-DBOOST_MPI_FOUND)
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${Boost_LIBRARIES})
endif()

if(Boost_VERSION EQUAL 105500)
  message(FATAL_ERROR "Boost 1.55 is not supported. Use a lower or a higher version instead!")
endif(Boost_VERSION EQUAL 105500)



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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnarrowing")
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
  
  # Visual Studio (Microsoft)
 elseif(MSVC)
    # This loop sets all linking to static (MT)
    foreach(flag_var
    CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        if(${flag_var} MATCHES "/MD")
            string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
        endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
endif()


###############################################################################
# Warning for Thrust-related Bug
###############################################################################
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
  endif(NOT "${CPLUS_TOLOWER_FOUND}" MATCHES "-1")
endfunction(warn_path_cuda)

warn_path_cuda("CPLUS_INCLUDE_PATH")
warn_path_cuda("C_INCLUDE_PATH")


###############################################################################
# Executables
###############################################################################
set(HASEonGPU_NAME "calcPhiASE")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
            ${CMAKE_CURRENT_SOURCE_DIR}/include/graybat)


cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB SRCFILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c*")

cuda_add_executable(${HASEonGPU_NAME} ${SRCFILES})

target_link_libraries(${HASEonGPU_NAME} ${LIBS} ${CUDA_CUDART_LIBRARY})


###############################################################################
# more targets (make examples, make new, make doc)
###############################################################################
if(CMAKE_COMPILER_IS_GNUCXX) # Only useful on Linux
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
  DEPENDS ${HASEonGPU_NAME}
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
endif(CMAKE_COMPILER_IS_GNUCXX)
