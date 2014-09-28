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
     + cmake 2.8.10
     + gcc 4.4.1
     + cuda 5.5

   + Optional:
     + octave / matlab
     + paraview
   
   + Hardware:
     + Nvidia device >= Compute capability 2.0 (at least fermi generation)


Compiling
---------

   + clone the repository: `git clone https://github.com/computationalradiationphysics/haseongpu.git`
   + create the build directory: `mkdir haseongpu/build`
   + go to build directory: `cd haseongpu/build`
   + create Makefile `cmake ..`
   + build project : `make`
 

Usage
-----

   + MATLAB compatible interface
   + C-Application interface


### Quick MATLAB laser pump example

  A small example for Phi ASE calculation
  with a pumped crystal. The simulation
  can be started by the following:

  1. follow the compile instructions above
  2. change path "cd example/matlab_example/"
  3. run : `matlab laserPumpCladdingExample`
  4. watch progress
  5. take a look at the results (*.vtk) with paraview 


### Quick C-Application laser pump example 

  1. follow the compile instructions above
  2. change path "cd example/c_example/"
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

    --parallel-mode=[|threaded|mpi]  
      Defines the method of parallelization to start the
      simulation with. Mode "threaded" uses pthreads on a single
      node. Mode "mpi" is a parallel mpi
      implementation for clusters. Note, that this parameter
      is currently only available when using `--device-mode=gpu`

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

    --spectral-resolution= 
      Resolution of absorption and emission spectrum to which the
      input spectrum will be interpolated linear.  Interpolation is
      used to distribute spectrum values equidistant over the
      wavelength.  Omitting this option or setting a to small
      resolutionwill set the lambda resolution to the maximum number
      of absorption or emission values.

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

