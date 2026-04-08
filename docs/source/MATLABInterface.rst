MATLAB Interface
================

This page documents the MATLAB-compatible interface of HASEonGPU.

The MATLAB interface is mainly intended for compatibility with existing MATLAB
or Octave workflows. If this is your first time using HASEonGPU, the
:doc:`Python Interface <pythonInterface>` is usually the better starting point.

For general setup and dependency information, see
:doc:`Getting Started <gettingStarted>`.

Compilation
-----------

The MATLAB interface requires a compiled ``calcPhiASE`` binary, see :doc:`compilation`.

Installation
------------

No separate installation step is required.

Add ``calcPhiASE.m`` to your MATLAB path, or run your MATLAB script from a
location where the wrapper can be found.

Usage
-----

Call ``calcPhiASE`` from a MATLAB script:

.. code-block:: matlab

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

The returned values are represented as two-dimensional matrices, where rows
correspond to point indices and columns correspond to slice indices.

Example
-------

A complete MATLAB example using laser pumping and cladding is provided under:

.. code-block:: text

   example/matlab_example/laserPumpCladdingExample.m

Typical workflow:

#. Compile HASEonGPU if needed
#. Change into ``example/matlab_example/``
#. Run ``matlab laserPumpCladdingExample``
#. Inspect the generated ``.vtk`` files with ParaView

Input Argument Description
--------------------------

The following section describes all arguments of the MATLAB call.

``points``
^^^^^^^^^^

* Element type: ``float``
* Value range: ``{0, ..., max(maxPointX, maxPointY)}``
* Size: ``2 * numberOfPoints``
* Description:
  Coordinates of the 2D mesh points. All x coordinates are followed by all y
  coordinates.

* Structure::

     [x_1, x_2, ... x_n, y_1, y_2, ... y_n]

``trianglePointIndices``
^^^^^^^^^^^^^^^^^^^^^^^^

* Element type: ``int``
* Value range: ``{0, ..., numberOfPoints}``
* Size: ``3 * numberOfTriangles``
* Description:
  Point indices defining the triangles of the mesh. Each triangle is described
  by three point indices.

* Structure::

     [triangle1_p1, triangle2_p1, ... triangleN_p1,
      triangle1_p2, triangle2_p2, ... triangleN_p2,
      triangle1_p3, triangle2_p3, ... triangleN_p3]

``betaCells``
^^^^^^^^^^^^^

* Element type: ``float``
* Size: ``numberOfPoints``
* Description:
  Stimulus values at the sample points.

``betaVolume``
^^^^^^^^^^^^^^

* Element type: ``float``
* Size: ``numberOfTriangles * (numberOfLevels - 1)``
* Description:
  Stimulus values in the prism volume.

  Prism values are ordered according to the prism ID:

  ``prismID = triangleID + layer * numberOfTriangles``

  Therefore, all prism values of one layer are grouped together.

``claddingCellTypes``
^^^^^^^^^^^^^^^^^^^^^

* Element type: ``int``
* Size: ``numberOfTriangles``
* Description:
  Cladding type index for each triangle.

``claddingNumber``
^^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Size: ``1``
* Description:
  Selects which cladding to use.

``claddingAbsorption``
^^^^^^^^^^^^^^^^^^^^^^

* Type: ``float``
* Size: ``1``
* Description:
  Absorption coefficient of the cladding.

``useReflections``
^^^^^^^^^^^^^^^^^^

* Type: ``bool``
* Size: ``1``
* Description:
  Enables or disables reflections at the top and bottom planes.

``refractiveIndices``
^^^^^^^^^^^^^^^^^^^^^

* Element type: ``float``
* Size: ``4``
* Description:
  Refractive indices of the active gain medium and its surrounding media.

* Structure::

     [bottomInside, bottomOutside, topInside, topOutside]

``reflectivities``
^^^^^^^^^^^^^^^^^^

* Element type: ``float``
* Value range: ``{0, ..., 1}``
* Size: ``2 * numberOfTriangles``
* Description:
  Reflectivities of the prism planes. First the reflectivities of the bottom
  plane, then the reflectivities of the top plane, both ordered by triangle
  index.

``triangleNormalsX``
^^^^^^^^^^^^^^^^^^^^

* Element type: ``float``
* Size: ``3 * numberOfTriangles``
* Description:
  x-components of the normal vectors for each triangle edge.

* Structure::

     [triangle0_edge0, triangle1_edge0, ..., triangleN_edge0,
      triangle0_edge1, triangle1_edge1, ..., triangleN_edge1,
      triangle0_edge2, triangle1_edge2, ..., triangleN_edge2]

``triangleNormalsY``
^^^^^^^^^^^^^^^^^^^^

* Element type: ``float``
* Size: ``3 * numberOfTriangles``
* Description:
  y-components of the normal vectors for each triangle edge.

* Structure::

     [triangle0_edge0, triangle1_edge0, ..., triangleN_edge0,
      triangle0_edge1, triangle1_edge1, ..., triangleN_edge1,
      triangle0_edge2, triangle1_edge2, ..., triangleN_edge2]

``triangleNeighbors``
^^^^^^^^^^^^^^^^^^^^^

* Element type: ``int``
* Value range: ``{-1, 0, ..., numberOfTriangles - 1}``
* Size: ``3 * numberOfTriangles``
* Description:
  Neighbor triangle indices for each triangle edge.
  ``-1`` means that no neighboring triangle exists at that edge.

``triangleSurfaces``
^^^^^^^^^^^^^^^^^^^^

* Element type: ``float``
* Size: ``numberOfTriangles``
* Description:
  Surface area of each triangle.

``triangleCenterX``
^^^^^^^^^^^^^^^^^^^

* Element type: ``float``
* Size: ``numberOfTriangles``
* Description:
  x-coordinates of the triangle center points.

``triangleCenterY``
^^^^^^^^^^^^^^^^^^^

* Element type: ``float``
* Size: ``numberOfTriangles``
* Description:
  y-coordinates of the triangle center points.

``triangleNormalPoint``
^^^^^^^^^^^^^^^^^^^^^^^

* Element type: ``unsigned``
* Value range: ``{0, ..., numberOfPoints}``
* Size: ``3 * numberOfTriangles``
* Description:
  Point indices where the corresponding triangle normal vectors start.
  For each triangle, one point index is stored for each of its three edges.

``forbiddenEdge``
^^^^^^^^^^^^^^^^^

* Element type: ``int``
* Value range: ``{-1, 0, 1, 2}``
* Size: ``3 * numberOfTriangles``
* Description:
  Edge indices of the adjacent triangles.
  ``-1`` means there is no adjacent triangle at that edge.
  ``0``, ``1`` and ``2`` denote the corresponding edge index in the adjacent
  triangle.

``minRaysPerSample``
^^^^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Size: ``1``
* Description:
  Minimum number of rays used for adaptive sampling.

``maxRaysPerSample``
^^^^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Size: ``1``
* Description:
  Maximum number of rays used for adaptive sampling.

``mseThreshold``
^^^^^^^^^^^^^^^^

* Type: ``float``
* Size: ``1``
* Description:
  Sets the maximal MSE of the ASE value. If a sample point does not reach this
  threshold, the number of rays per sample point is increased up to
  ``maxRaysPerSample`` or resampled with repetitive sampling.

``repetitions``
^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Size: ``1``
* Description:
  Maximum number of repetitions if the MSE threshold is not reached.

``adaptiveSteps``
^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Size: ``1``
* Description:
  Number of adaptive steps. The range between minimum and maximum ray count is
  split into that many parts. Setting it to ``1`` disables adaptive stepping
  and uses only ``minRaysPerSample``.

``nTot``
^^^^^^^^

* Type: ``float``
* Size: ``1``
* Description:
  Doping of the active gain medium.

``thickness``
^^^^^^^^^^^^^

* Type: ``float``
* Size: ``1``
* Description:
  Thickness of one prism layer of the mesh.

``laserParameter``
^^^^^^^^^^^^^^^^^^

* Type: ``struct``
* Description:
  Structure containing the laser material spectra.
  It provides the absorption and emission spectrum values together with the
  corresponding wavelength values.

* Required fields:

  * ``l_abs``
  * ``l_ems``
  * ``s_abs``
  * ``s_ems``
  * ``l_res``

* Field layouts:

  * ``l_abs``
    * Container type: vector
    * Element type: ``double``
    * Supported layout: flat
    * Description:
      Wavelength values of the absorption spectrum in ``nm``.
      Accessed in MATLAB as ``laserParameter.l_abs``.

  * ``l_ems``
    * Container type: vector
    * Element type: ``double``
    * Supported layout: flat
    * Description:
      Wavelength values of the emission spectrum in ``nm``.
      Accessed in MATLAB as ``laserParameter.l_ems``.

  * ``s_abs``
    * Container type: vector
    * Element type: ``double``
    * Supported layout: flat
    * Description:
      Values of the absorption spectrum in ``cm^2``, corresponding to
      ``l_abs``.
      Accessed in MATLAB as ``laserParameter.s_abs``.

  * ``s_ems``
    * Container type: vector
    * Element type: ``double``
    * Supported layout: flat
    * Description:
      Values of the emission spectrum in ``cm^2``, corresponding to
      ``l_ems``.
      Accessed in MATLAB as ``laserParameter.s_ems``.

  * ``l_res``
    * Type: scalar
    * Element type: ``double``
    * Description:
      Resolution used for linear interpolation of the spectrum.
      Accessed in MATLAB as ``laserParameter.l_res``.


``crystal``
^^^^^^^^^^^

* Type: ``struct``
* Description:
  Structure containing the crystal parameters.

* Required fields:

  * ``tfluo``

* Field layouts:

  * ``tfluo``
    * Type: scalar
    * Element type: ``double``
    * Description:
      Fluorescence lifetime of the active medium in seconds.
      Accessed in MATLAB as ``crystal.tfluo``.

``numberOfLevels``
^^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Size: ``1``
* Description:
  Total number of mesh levels in z-direction.

``deviceMode``
^^^^^^^^^^^^^^

* Type: ``str``
* Allowed values:

  * ``"cpu"``
  * ``"gpu"``

* Description:
  Selects whether the CPU or GPU implementation is used.
  The CPU implementation is mostly deprecated.

``parallelMode``
^^^^^^^^^^^^^^^^
* Type: ``str``
* Allowed values:

  * ``"threaded"``: threads distribute samples across GPUs/devices within one node
  * ``"mpi"``: distributes samples across nodes


* Description:
  Selects the parallelization mode.

``maxGPUs``
^^^^^^^^^^^

* Type: ``unsigned``
* Description:
  Maximum number of GPUs used by a single process.

  In ``threaded`` mode, this usually corresponds to the total number of GPUs
  used on the local node.

  In ``mpi`` mode, this corresponds to the number of GPUs used per MPI process
  rank. Therefore, if ``nPerNode = 1 (default)`` , it corresponds to the number of
  GPUs used per node.

``nPerNode``
^^^^^^^^^^^^

* Type: ``unsigned``
* Description:
  Number of MPI process ranks launched per node (default = 1).
