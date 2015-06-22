GrayBat
=======

<b>Gr</b>aph <b>A</b>pproach  for Highl<b>y</b>  Generic Communication
Schemes <b>B</b>ased on <b>A</b>daptive <b>T</b>opologies

##Description##

**GrayBat** is a C++ library that presents a graph-based communication
approach, which enables a mapping of algorithms to communication
patterns and further a mapping of these communication patterns to
varying hardware architectures. Therefore, a flexible and configurable
communication approach for parallel and distributed
applications. These mappings are established as an intermediate layer
between an application and communication libraries and are dynamically
adptable during run-time.

##Documentation##

Have a look at the documentation that is available [here](https://ComputationalRadiationPhysics.github.io/graybat)

##Referencing##

GrayBat is a scientific project. If you **present and/or publish** scientific
results that used GrayBat, you should set this as a **reference**.

##Software License##


GrayBat  is licensed under the <b>GPLv3+</b>. Please refer to our [LICENSE.md](LICENSE.md)


##Dependencies##

 * cmake 3.0.2
 * Boost 1.57.0
 * OpenMPI 1.8.0 (mpi communication policy)
 * g++ 4.9.2
 * metis 5.1 (graph partitioning)

##Compiling##

 * Clone the repository: `git clone https://github.com/computationalradiationphysics/graybat.git`
 * Change directory: `cd graybat`
 * Init git submodules: `git submodule init && git submodule update`
 * Create the build directory: `mkdir -p build`
 * Change to build directory: `cd build`
 * Set compiler: `export CXX=[g++,clang++]`
 * Create Makefile `cmake ..`
 * Build project : `make [target]`

##Tested Compilers##

 * clang 3.5
 * g++ 4.9.2

### Current Compilation Status:

| *branch* | *state* | *description* |
| -------- | --------| ------------- |
| **master** | [![Build Status](http://haseongpu.mooo.com/api/badge/github.com/erikzenker/GrayBat/status.svg?branch=master)](http://haseongpu.mooo.com/github.com/erikzenker/GrayBat) |  stable releases |
| **dev**  | [![Build Status](http://haseongpu.mooo.com/api/badge/github.com/erikzenker/GrayBat/status.svg?branch=dev)](http://haseongpu.mooo.com/github.com/erikzenker/GrayBat) |development branch |

##Example Targets##

 * **GoL**: Game of Life simulation

 * **doc**: Build documentation in doc/

 * **clean**: Cleanup build directory

##Literature##
 * Talk by Erik Zenker of his diploma defence [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.16306.svg)](http://dx.doi.org/10.5281/zenodo.16306)


##Authors##

 * Erik Zenker (erikzenker@posteo.de)

