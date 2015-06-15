project(genericCommunicator)
cmake_minimum_required(VERSION 2.8.5)


###############################################################################
# CMAKE module path in utils/cmake/modules/
###############################################################################
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/utils/cmake/modules/" ${CMAKE_MODULE_PATH})

###############################################################################
# METIS
###############################################################################
find_package(METIS 5.1.0)
include_directories(SYSTEM ${METIS_INCLUDE_DIRS})
set(LIBS ${LIBS} ${METIS_LIBRARIES})

###############################################################################
# Boost
###############################################################################
find_package(Boost 1.58.0 COMPONENTS mpi unit_test_framework serialization REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})

################################################################################
# MPI
################################################################################
find_package(MPI REQUIRED)
  include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
set(LIBS ${LIBS} ${MPI_C_LIBRARIES})

# bullxmpi fails if it can not find its c++ counter part
if(MPI_CXX_FOUND)
  set(LIBS ${LIBS} ${MPI_CXX_LIBRARIES})
endif(MPI_CXX_FOUND)

###############################################################################
# Compiler Flags
###############################################################################
string(FIND ${CMAKE_CXX_COMPILER} "clang++" CXX_COMPILER_IS_CLANG++)
string(FIND ${CMAKE_CXX_COMPILER} "g++" CXX_COMPILER_IS_G++)

if(${CXX_COMPILER_IS_G++})
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(${CXX_COMPILER_IS_G++})

if(${CXX_COMPILER_IS_CLANG++})
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif(${CXX_COMPILER_IS_CLANG++})

###############################################################################
# Build executables
###############################################################################
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
		   ${CMAKE_CURRENT_SOURCE_DIR}/include/dout/include )
#add_executable(Main EXCLUDE_FROM_ALL src/main.cc)
add_executable(GoL EXCLUDE_FROM_ALL example/gol.cc )
add_executable(ZMQ1 EXCLUDE_FROM_ALL example/zmq/request_reply.cc )
add_executable(Test EXCLUDE_FROM_ALL utils/unittests/CageUT.cc )



# add the libraries (Boost...) to this executable
target_link_libraries(GoL ${LIBS})
target_link_libraries(ZMQ1 ${LIBS})
target_link_libraries(Test ${LIBS})


###############################################################################
# Custom targets
###############################################################################

# doc
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.conf ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc 
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)
