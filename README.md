HASEonGPU 
=======================================

<b>H</b>igh performance <b>A</b>mplified <b>S</b>pontaneous <b>E</b>mission <b>on GPU</b>

Description
-----------


Referencing
-----------

   HASEonGPU is a *scientific project*. If you **present and/or publish** scientific
   results that used HASEonGPU, you should set this as a **reference**.

Software License
----------------

   *HASEonGPU* is licensed under the **GPLv3+**.
   Please refer to our [LICENSE.md](LICENSE.md)


Dependencies
------------

   + Software:
     + make
     + cmake >= 3.0.1
     + gcc >= 4.8.2
     + cuda >= 5.5
     + boost >= 1.50 ( **not** 1.55, due to a regression bug )

   + Optional:
     + octave / matlab
     + paraview
     + OpenMPI 1.8 or compatible
   
   + Hardware:
     + Nvidia device >= Compute capability 2.0 (at least fermi generation)


Compiling
---------

   + clone the repository: `git clone https://github.com/computationalradiationphysics/haseongpu.git`
   + create the build directory: `mkdir haseongpu/build`
   + go to build directory: `cd haseongpu/build`
   + create Makefile `cmake ..`
   + build project : `make`

### Current Compilation Status:

| *branch* | *state* | *description* |
| -------- | --------| ------------- |
| **master** | [![Build Status](http://haseongpu.mooo.com/api/badge/github.com/ComputationalRadiationPhysics/haseongpu/status.svg?branch=master)](http://haseongpu.mooo.com/github.com/ComputationalRadiationPhysics/haseongpu) | our stable new releases |
| **dev**  | [![Build Status](http://haseongpu.mooo.com/api/badge/github.com/ComputationalRadiationPhysics/haseongpu/status.svg?branch=dev)](http://haseongpu.mooo.com/github.com/ComputationalRadiationPhysics/haseongpu) | our development branch |
 

Usage
-----

   + MATLAB compatible interface
   + C-Application interface


### Quick MATLAB laser pump example

  A small example for Phi ASE calculation
  with a pumped crystal. The simulation
  can be started by the following:

  1. follow the compile instructions above
  2. change path `cd example/matlab_example/`
  3. run : `matlab laserPumpCladdingExample`
  4. watch progress
  5. take a look at the results (`*.vtk`) with paraview


### Quick C-Application laser pump example 

  1. follow the compile instructions above
  2. change path `cd example/c_example/`
  3. run : `./bin/calcPhiASE --input-path=./input/cylindrical --min-rays=10000`
  4. watch progress
  5. take a look at the results in the output directory

### MATLAB compatible interface

   + Add calcPhiASE.m to your matlab path
   + Call calcPhiASE from your matlab script 
     like following:

```matlab
[phiASE, MSE, nRays] = calcPhiASE(
                                    points,  
                                    trianglePointIndices,           
                                    betaCells,  
                                    betaVolume,  
                                    claddingCellTypes,  
                                    claddingNumber,  
                                    claddingAbsorption,  
                                    useReflections,  
                                    refractiveIndices,  
                                    reflectivities,  
                                    triangleNormalsX,  
                                    triangleNormalsY,  
                                    triangleNeighbors,  
                                    triangleSurfaces,  
                                    triangleCenterX,  
                                    triangleCenterY,  
                                    triangleNormalPoint,  
                                    forbiddenEdge,  
                                    minRaysPerSample,  
                                    maxRaysPerSample,  
                                    mseThreshold,  
                                    repetitions,  
                                    adaptiveSteps,
                                    nTot,  
                                    thickness,  
                                    laserParameter,  
                                    crystal,  
                                    numberOfLevels,  
                                    deviceMode,
                                    parallelMode,  
                                    maxGPUs,  
                                    nPerNode
                                );  
```
  + The returned values are represented as
    two-dimensional matrices in which columns are
    slice indices(levels) and rows are point
    indices. The value for the ith point and jth 
    slice can then be optained by MATLAB with:
    	  
	  value = values(i,j);


### Input argument description

In the following all arguments of the MATLAB call are described.
You will find on each point a head with datatype (an array when
in brackets []), the size of the array and to which set of
numbers the array belongs.

   + __points__ [float], in {0, ..., numberOfPoints}, size = numberOfPoints  
     The coordinates of the triangle vertices. All x coordinates followed by all
     of the y coordinates of the triangle vertices  
     structure: [x_1, x_2, ... x_n, y_1, y_2, ... y_n] (n == numberOfPoints)

   + __trianglePointIndices__ [int] in {0, ..., numberOfPoints}, size = numberOfTriangles * 3  
     Contains the indices to access the "points" datastructure 
     (each triangle has 3 points as vertices). Each entry is an
     index from 0 to numberOfPoints, corresponding to the positions 
     of a vertex in "points".
     Structure as follows:  
     [ triangle1A, triangle2A, ... triangleNA, triangle1B, triangle2B, ... triangleNB, triangle1C, ... ]  
     i.e. for triangles with vertices A,B,C there are all the indices
     of the A-vertices, followed by all the B and C vertices.

   + __betaCells__ [float]  
     Stimulus in the sample points.

   + __betaVolume__ [float], size = numberOfTriangles * numberOfLevels - 1  
     Stimulus in the volume (prisms).
     Beta values for all prisms ordered accordingly to the prismIDs:
     prismID = triangleID + layer * numberOfTriangles.
     Therefore, all betaValues for a layer are grouped together

   + __claddingCellTypes__ [int], size = numberOfTriangles  
     Sets cladding index for triangles {0,1,2,...}

   + __claddingNumber__ unsigned, size = 1  
     Set which cladding to use

   + __claddingAbsorption__ float, size = 1  
     Absorption coefficient of cladding

   + __useReflections__ bool, size = 1  
     Switch to activate reflections.

   + __refractiveIndices__ [float], size = 4  
     Describes the refractive indices of the active
     gain medium top and bottom planes.
     It is structured as follows:
     {bottomInside, bottomOutside, topInside, topOutside}
     bottomInside = topInside (because it is the same medium)

   + __reflectivities__ [float], in {0, ...,1}, size = 2 * numberOfTriangles  
     Defines the reflectivities of prism planes.
     First the reflectivities of bottom plane and then the reflectivities
     of top plane. Both it ordered by TriangleID.

   + __triangleNormalsX__ [float], size = numberOfTriangles * 3  
     The x coordinate of the normal vectors for each triangle edge.
     It is ordered as follows:  
     [ triangle1_edge0, triangle2_edge0, ... triangleN_edge0, triangle1_edge1, triangle2_edge1, ... ]  
     i.e. all first edges of each triangle, followed by all second edges of each triangle and so on.

   + __triangleNormalsY__ [float], size = numberOfTriangles * 3  
     The y coordinate of the normal vectors for each triangle edge.
     It is ordered as follows:  
     [ triangle1_edge0, triangle2_edge0, ... triangleN_edge0, triangle1_edge1, triangle2_edge1, ... ]  
     i.e. all first edges of each triangle, followed by all second edges of each triangle and so on.

   + __triangleNeighbors__ [int], in {-1,0,1,2,4}, size = 5  
     Describes the neighnor relation of triangles to each other.
     Each entry corresponds to a triangleID (see "triangles") which
     is adjacent to the current triangle and edge.
     Structure is similar to "forbidden":  
     [ triangle1_edge0, triangle2_edge0, ... triangleN_edge0, triangle1_edge1, triangle2_edge1, ... ]

   + __triangleSurfaces__ [float], size = numberOfTriangles  
     The sizes of the surfaces of each triangle, ordered by the triangleID.

   + __triangleCenterX__ [float], size = numberOfTriangles  
     The x coordinates of the center points for each triangle
     ordered by TriangleID.

   + __triangleCenterY__ [float], size = numberOfTriangles  
     The y coordinates of the center points for each triangle
     ordered by TriangleID.

   + __triangleNormalPoint [unsigned]__, in {0, ...,  numberOfPoints}, size = numberOfTriangles * 3  
     Contains indices to the point where the
     triangleNormalVectors start. For each Triangle 3 points (3 edges)
     are stored in this list. Indices point to locations in "points" 
     (i.e. normal vectors start at triangle vertices!)l
     Structure is VERY similar to triangles: 
     [ triangle1_p0, triangle2_p0, ... triangleN_p0, triangle1_p1, triangle2_p1, ... ]

   + __forbiddenEdge__ [int], in {-1,0,1,2,4}, size = 5  
     Describes the relation of edge indices of adjacent triangles
     -1 means, there is no adjacent triangle to that edge
     0,1,2 describes the index of the edge as seen from the ADJACENT triangle
     Order of data is similar to normalVec:
     [ triangle1_edge0, triangle2_edge0, ... triangleN_edge0, triangle1_edge1, triangle2_edge1, ... ]
     i.e. all first edges of each triangle, followed by all second edges of each triangle and so on.

   + __minRaysPerSample__ unsigned, size = 1  
     Minimal number of rays for adaptive sampling

   + __maxRaysPerSample__ unsigned, size = 1  
     Maximal number of rays for adaptive sampling
 
   + __mseThreshold__ float, size = 1  
      Sets the maximal MSE of the ASE value.
      If a sample-point does not reach this MSE-threshold, the number 
      of rays per sample-point will be increased upto maxRaysPerSample or
      resampled with repetitive sampling.

   + __repetitions__ unsigned, size = 1  
     Sets the number of maximal repetitions when the
     mseThreshold was not reached.

   + __adaptiveSteps__ unsigned, size = 1  
     Sets the number of adaptive steps. The range between min-rays and max-rays will be
     split into that many parts. Setting it to 1 will result in no adaptive
     steps and only minRaysPerSample is used.

   + __nTot__ float, size = 1  
     Doping of the active gain medium

   + __thickness__ float, size = 1  
     Thickness of one prism level of the mesh.

   + __laserParameter__ [float]   
     Is a structure for the laser parameters (intensities sigma, wavelength lambda)
     s_ems corresponds to l_ems and s_abs to l_abs
     struct(s_abs, VALUES, s_ems, VALUES, l_abs, VALUES, l_ems, VALUES)

   + __crystal__ [float]  
     Is a structure for the crystal parameters 
     crystal.tfluo describes the crystalFluorescence of the active gain medium.

   + __numberOfLevels__ unsigned, size = 1  
     Total number of levels of the mesh. Thus the total thickness
     of the mesh is thickness * numberOfLevels!

   + __parallelMode__  
    + 'mpi'      use mpi to distribute workload (use nPerNode)
    + 'threaded' use pthreads to distribute workload locally
    + 'graybat'  use experimental graybat api (internally mpi)

   + __deviceMode__
    + 'cpu'      use cpu algorithm (does not have all features)
    + 'gpu'      use gpu algorithm

   + __maxGpus__ unsigned, size = 1  
     Maximal number of GPUs for threaded case

   + __nPerNode__  
     Number of devices per mpi-node


Synopsis
--------

   + Command:
   
        ./bin/calcPhiASE [OPTIONS] 

   + Options:
   
    --input-path 
      Path to the experiment location.
      This folder contains several .txt files usually
      generated by an matlab script. The content of this
      .txt files contains all experiment data you need
      to run one experiment.
        
    --output-path
      Path to a writable location. Is used to write
      input and output for matlab script.

    --parallel-mode=[graybat|threaded|mpi]
      Defines the method of parallelization to start the
      simulation with. Mode "threaded" uses pthreads on a single
      node. Mode "mpi" is a parallel mpi implementation for clusters. Mode
      "graybat" is similar to "mpi", but uses the communication framework
      [GrayBat](https://github.com/ComputationalRadiationPhysics/graybat)
      instead of plain MPI. Note, that the last 2 parameters
      are currently only available when using `--device-mode=gpu`

    --device-mode=[cpu|gpu]  
      Defines on which hardware the simulation will run.
      Mode "cpu" is the original
      algorithm based on single core cpu.
      Mode "gpu" uses nVIDIA CUDA GPUs, that can be parallelized either
      with Pthreads or MPI.

    --min-rays=  
      Sets the minimum number of rays per sample-point in the
      crystal structure.

    --max-rays=  
      Sets the maximal number of rays per sample-point. The number
      of rays per sample-point will vary between minimum and
      maximum number of rays in dependance of a MSE-Threshold.
      (see --mse-threshold)

    --ngpus=  
      Set the number of gpus to use. "mpi" parallel-mode should set this
      to 1 and "threaded" to the maximum number
      of GPUs on the node. If you don't set it, it will
      be set to the maximum automatically.

    --min-sample-i=  
      Index of the first sample point (normally 0).

    --max-sample-i=  
      Index of the last sample point (numberOfSample - 1).
          
    --verbosity=  
      Add the following for different verbosity levels:
      0  : quiet
      1  : error
      2  : warning
      4  : info
      8  : statistics
      16 : debug
      32 : progress-bar
      Levels the verbosity level is interpreted as a bitmask and 
      can be composed of different levels.

    --reflection  
      Use reflection on upper and lower plane of gain
      medium. Maximal number of reflections will be
      calculated

    --mse-threshold=  
      Algorithm tries to stay under this threshold
      by adaptive and repetitive sampling.

    --repetitions=  
      Number of repetitions, that will be done
      when mse-threshold was not met.

    --adaptive-steps= 
      Number of adaptive steps. The range between min-rays and max-rays will be
      split into that many parts. Setting --adaptive=1 will result in no adaptive
      steps and only min-rays is used.

    --spectral-resolution= 
      Resolution of absorption and emission spectrum to which the
      input spectrum will be interpolated linear.  Interpolation is
      used to distribute spectrum values equidistant over the
      wavelength.  Omitting this option or setting a to small
      resolutionwill set the lambda resolution to the maximum number
      of absorption or emission values.

    --config=
      Specify the location of an optional configuration file. This file holds
      simple key=value pairs that have the same names as the command-line
      parameters. If values are defined in the file and on the command line,
      values from the file are overruled. It is not possible to specify
      --help or --config= in the config file itself.
      For an example configuration, see the C example (`calcPhiASE.cfg`).

C-Application Templates
-----------------------

+ 4 GPUs, 10K to 100K Rays, 4 Repetitions   

        ./bin/calcPhiASE --input-path=/input/  
	       	 --output-path=/tmp/  
    		 --parallel-mode=threaded  
    		 --min-rays=10000  
    		 --max-rays=100000  
    		 --reflection  
    		 --repetitions=4  
    		 --adaptive-steps=4  
    		 --ngpus=4  
    		 --min-sample-i=0  
    		 --max-sample-i=1234  
    		 --mse-threshold=0.05  

+ MPI with 4 GPUs per node   

        mpiexec -npernode 4 ./bin/calcPhiASE  --input-path=/input/  
		     --output-path=/tmp/  
    		 --parallel-mode=mpi  
    		 --min-rays=10000  
    		 --max-rays=100000  
    		 --reflection  
    		 --repetitions=4  
    		 --adaptive-steps=4  
    	   --ngpus=1  
    		 --min-sample-i=0  
    		 --max-sample-i=1234  
    		 --mse-threshold=0.05  
		 
Active Team
----------

### Scientific Supervision

  + Dr. Michael Bussmann
  + Dr. Daniel Albach

### Maintainers and core developers

  + Erik Zenker
  + Carlchristian Eckert      

### Participants, Former Members and Thanks

  + Marius Melzer
  + Frank Liebold
  + Daniel Höhne


File Descriptions
-----------------

 - `CMakeLists.txt` cmake file to generate a Makefile
 - `src/` folder containing all the source code that is not a header
  - `src/map_rays_to_prisms.cu` CUDA code to generate a schedule of which ray will be launched from which prism of the gain medium
  - `src/calcPhiASE.m` MATLAB adapter script 
  - `src/logging.cu` creates nicely readable output based on log-levels
  - `src/importance_sampling.cu` CUDA parallelized importance sampling
  - `src/geometry.cu` basic 3D geometry calculations
  - `src/interpolation.cu` interpolation functions for wavelengths of polychromatic laser pulses
  - `src/reflection.cu` CUDA functions to calculate reflections inside the gain medium
  - `src/progressbar.cu` progressbar for the command line
  - `src/parser.cu` parsing of command line arguments and textfiles
  - `src/for_loops_clad.cu` old CPU code for ASE calculation
  - `src/calc_phi_ase_mpi.cc` MPI workload distribution. Code for Master and Slaves
  - `src/write_to_file.cu` writing formatted data to a file
  - `src/ray_histogram.cu` print a histogram of the adaptive ray count to command line
  - `src/calc_phi_ase_threaded.cu` pthreads workload distribution
  - `src/mt19937ar.cu` CPU code for Mersenne Twister PRNG used by for_loops_clad.cu
  - `src/write_to_vtk.cu` generate VTK-files from simulation results (deprecated)
  - `src/propagate_ray.cu` CUDA code to propagate a single ray through the prism mesh structure
  - `src/mesh.cu` class that holds the information and all parameters about the gain medium mesh
  - `src/cuda_utils.cu` utility functions (getFreeDevice)
  - `src/calc_sample_gain_sum.cu` CUDA code to calculate all the rays for a single sample point
  - `src/calc_phi_ase.cu` CUDA code to calculate ASE for all the sample points 
  - `src/main.cu` main entry file
  - `src/write_matlab_output.cu` generate MATLAB-readable matrixes from the simulation data
 - `example/` folder that contains executable examples
  - `example/c_example/` folder that contains the commandline-example
   - `example/c_example/input/` folder that contains input for 2 different experiments
    - `example/c_example/input/cylindrical/` folder that contains data for the cylindrical gain medium. For details on the files, see detailed information above (Input argument description)
    - `example/c_example/input/cuboid/` example input with a cuboid gain medium. contents similar to cylindrical example.
   - `example/c_example/output/` folder to gather the output
  - `example/matlab_example/` folder that contains the input data for the matlab example
   - `example/matlab_example/lambda_e.txt` emission wavelengths
   - `example/matlab_example/sigma_e.txt` emission crosssection
   - `example/matlab_example/pt.mat` sampling points and delaynay-triangles of the gain medium
   - `example/matlab_example/set_variables.m` generate information about the mesh
   - `example/matlab_example/vtk_wedge.m` generate a VTK file from the mesh
   - `example/matlab_example/laserPumpCladdingExample.m` experimental setup. Run this file to see the progress of a whole experiment
   - `example/matlab_example/sve.mat` 
   - `example/matlab_example/sigma_a.txt` absorption crosssection
   - `example/matlab_example/gain.m` calculate gain distribution inside the gain medium
   - `example/matlab_example/beta_int3.m` utility function to calculate gain distribution
   - `example/matlab_example/extract_gain_map.m` calculate the gain for the sample point used in the actual measurement
   - `example/matlab_example/beta_int.m` utility function to calculate gain distribution 
   - `example/matlab_example/lambda_a.txt` absorption wavelengths
 - `include/` folder containing all the header source code
  - `include/calc_phi_ase_mpi.hpp`  header for calc_phi_ase_mpi.cu
  - `include/mesh.hpp` header for mesh.cu
  - `include/importance_sampling.hpp` header for importance_sampling.cu
  - `include/ray_histogram.hpp` header for ray_histogram.cu
  - `include/for_loops_clad.hpp` header for for_loops_clad.cu
  - `include/mt19937ar.hpp` header for mt19937ar.cu
  - `include/calc_phi_ase.hpp` header for calc_phi_ase.cu
  - `include/write_matlab_output.hpp` header for write_matlab_output.cu
  - `include/cuda_utils.hpp` header for cuda_utils.cu
  - `include/logging.hpp` header for logging.cu
  - `include/cudachecks.hpp` Macros to check the success state of CUDA calls
  - `include/reflection.hpp` header for reflection.cu
  - `include/parser.hpp` header for parser.cu
  - `include/map_rays_to_prisms.hpp` header for map_rays_to_prisms.cu
  - `include/calc_phi_ase_threaded.hpp` header for calc_phi_ase_threaded.cu
  - `include/thrust_device_vector_nowarn.hpp` wrapper to switch off compiler warning that is produced by 3rd party library (CUDA Thrust)
  - `include/propagate_ray.hpp` header for propagate_ray.cu
  - `include/thrust_host_vector_nowarn.hpp` wrapper to switch off compiler warning that is produced by 3rd party library (CUDA Thrust)
  - `include/calc_sample_gain_sum.hpp` header for calc_sample_gain_sum.cu
  - `include/interpolation.hpp` header for interpolation.cu
  - `include/version.hpp` version information for HASEonGPU
  - `include/geometry.hpp` header for geometry.cu
  - `include/write_to_file.hpp` header for write_to_file.cu
  - `include/types.hpp` type definitions for HASEonGPU
  - `include/progressbar.hpp` header for progressbar.cu
  - `include/nan_fix.hpp` wrapper to allow usage of `isnan()` in a template
  - `include/write_to_vtk.hpp` header for write_to_vtk.cu
 - `LICENSE.md`  additional licensing information
 - `README.md` this README file
 - `REFERENCE.md`  Referencing information
 - `COPYING` Full License information
 - `utils/` folder that contains utility files
  - `utils/cmake/`
   - `utils/cmake/modules/` 3rd Party CMAKE module that was modified to circumvent a bug where the NVCC linker would crash on unknown arguments

