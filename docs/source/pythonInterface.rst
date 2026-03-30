Python Interface
================

This page documents the Python interface of HASEonGPU.

The Python interface is the recommended starting point for new users. It
provides a library call that can be integrated into custom workflows,
and the provided examples demonstrate this usage.

For general setup and dependency information, see
:doc:`Getting Started <gettingStarted>`.

Installation
------------

Install the haseongpu Python library from the repository root:

.. code-block:: bash

   pip install -e .

Usage
-----

The main Python entry point is:

``calcPhiASE(...)``

For single-node runs, the Python binding is used directly without disk I/O.

If ``parallelMode`` is set to ``"mpi"`` or ``"graybat"``, the Python interface
writes input data to a temporary directory, launches the external
``calcPhiASE`` executable via ``mpiexec``, and reads the results back into
Python.

In typical workflows, ``calcPhiASE(...)`` is used as part of a larger simulation
pipeline. This usually includes preparing the mesh and material parameters,
computing the gain distribution from laser pumping, and evaluating the ASE flux
and related quantities over one or multiple time steps.

Such a workflow is demonstrated by the accompanying example:

.. code-block:: text

   example/python_example/laserPumpCladdingExample.py



Example
-------

A complete Python example using laser pumping and cladding is provided under:

.. code-block:: text

   example/python_example/laserPumpCladdingExample.py

Minimal example (only using the library function call):

.. code-block:: python

    from HASEonGPU import calcPhiASE

    phi_ASE, mse_values, n_rays = calcPhiASE(
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
    )

Input Conventions
-----------------

The Python interface supports two array layout conventions most array-like
arguments (see **input argument description** below for details):

* **flat layout**
* **matrix layout**

The selected layout must be used consistently across all array-like inputs of a
single ``calcPhiASE(...)`` call.

The functions return type depends on the input-layout.

If flat layout is used for the input arrays, array-like outputs are returned in
flat layout.

If matrix layout is used for the input arrays, array-like outputs are returned
in matrix layout.

Matrix array layout:

.. code-block:: text

   (N_points, 2) -> [[x0, y0], [x1, y1], ...]

Flat array layout:

.. code-block:: text

   (2, N_points) -> [x0, x1, ..., xn, y0, y1, ..., yn]

Input Argument Description
--------------------------

The following section describes all arguments of the Python call.

``points``
^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``float``
* Value range: ``{0, ..., max(maxPointX, maxPointY)}``
* Description:
  Coordinates of the 2D mesh points.

* Supported layouts:

  * flat
    * Size: ``2 * numberOfPoints``
    * Structure::

         [x_1, x_2, ... x_n, y_1, y_2, ... y_n]

  * matrix
    * Size: ``(numberOfPoints, 2)``
    * Structure::

         [[x0, y0], [x1, y1], ...]

``trianglePointIndices``
^^^^^^^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``int``
* Value range: ``{0, ..., numberOfPoints}``
* Description:
  Contains the point indices defining the triangles of the mesh.
  Each triangle is described by three point indices.

* Supported layouts:

  * flat
    * Size: ``3 * numberOfTriangles``
    * Structure::

         [triangle1_p1, triangle2_p1, ... triangleN_p1,
          triangle1_p2, triangle2_p2, ... triangleN_p2,
          triangle1_p3, triangle2_p3, ... triangleN_p3]

  * matrix
    * Size: ``(numberOfTriangles, 3)``
    * Structure::

         [[p1, p2, p3],
          [p1, p2, p3],
          ...]

``betaCells``
^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``float``
* Description:
  Stimulus values at the sample points.

* Supported layout:

  * flat
    * Size: ``numberOfPoints``
    * Structure::

         [beta_1, beta_2, ..., beta_n]

``betaVolume``
^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``float``
* Description:
  Stimulus values in the prism volume.

* Supported layouts:

  * flat
    * Size: ``numberOfTriangles * (numberOfLevels - 1)``
    * Structure::

         [layer0_triangle0, layer0_triangle1, ..., layer0_triangleN,
          layer1_triangle0, layer1_triangle1, ..., layer1_triangleN,
          ...]

  * matrix
    * Size: ``(numberOfTriangles, numberOfLevels - 1)``
    * Structure::

         [[triangle0_layer0, triangle0_layer1, ..., triangle0_layerM],
          [triangle1_layer0, triangle1_layer1, ..., triangle1_layerM],
          ...]

``claddingCellTypes``
^^^^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``int``
* Description:
  Cladding type index for each triangle.

* Supported layout:

  * flat
    * Size: ``numberOfTriangles``
    * Structure::

         [cladding_0, cladding_1, ..., cladding_N]

``claddingNumber``
^^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Description:
  Selects which cladding to use.

``claddingAbsorption``
^^^^^^^^^^^^^^^^^^^^^^

* Type: ``float``
* Description:
  Absorption coefficient of the cladding.

``useReflections``
^^^^^^^^^^^^^^^^^^

* Type: ``bool``
* Description:
  Enables or disables reflections at the top and bottom planes.

``refractiveIndices``
^^^^^^^^^^^^^^^^^^^^^

* Type: ``[float]``
* Accepted input:

  * NumPy array
  * Python list

* Description:
  Refractive indices of the active gain medium at the top and bottom planes.

* Supported layout:

  * flat
    * Size: ``4``
    * Structure::

         [bottomInside, bottomOutside, topInside, topOutside]

``reflectivities``
^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``float``
* Value range: ``{0, ..., 1}``
* Description:
  Reflectivities of the prism planes. First the reflectivities of the bottom
  plane, then the reflectivities of the top plane, both ordered by triangle
  index.

* Supported layout:

  * flat
    * Size: ``2 * numberOfTriangles``
    * Structure::

         [bottom_triangle0, bottom_triangle1, ..., bottom_triangleN,
          top_triangle0, top_triangle1, ..., top_triangleN]

``triangleNormalsX``
^^^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``float``
* Description:
  x-components of the normal vectors for each triangle edge.

* Supported layouts:

  * flat
    * Size: ``3 * numberOfTriangles``
    * Structure::

         [triangle0_edge0, triangle1_edge0, ..., triangleN_edge0,
          triangle0_edge1, triangle1_edge1, ..., triangleN_edge1,
          triangle0_edge2, triangle1_edge2, ..., triangleN_edge2]

  * matrix
    * Size: ``(numberOfTriangles, 3)``
    * Structure::

         [[edge0_x, edge1_x, edge2_x],
          [edge0_x, edge1_x, edge2_x],
          ...]

``triangleNormalsY``
^^^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``float``
* Description:
  y-components of the normal vectors for each triangle edge.

* Supported layouts:

  * flat
    * Size: ``3 * numberOfTriangles``
    * Structure::

         [triangle0_edge0, triangle1_edge0, ..., triangleN_edge0,
          triangle0_edge1, triangle1_edge1, ..., triangleN_edge1,
          triangle0_edge2, triangle1_edge2, ..., triangleN_edge2]

  * matrix
    * Size: ``(numberOfTriangles, 3)``
    * Structure::

         [[edge0_y, edge1_y, edge2_y],
          [edge0_y, edge1_y, edge2_y],
          ...]

``triangleNeighbors``
^^^^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``int``
* Description:
  Neighbor triangle indices for each triangle edge.
  ``-1`` means that no neighboring triangle exists at that edge.

* Supported layouts:

  * flat
    * Size: ``3 * numberOfTriangles``
    * Structure::

         [triangle0_edge0, triangle1_edge0, ..., triangleN_edge0,
          triangle0_edge1, triangle1_edge1, ..., triangleN_edge1,
          triangle0_edge2, triangle1_edge2, ..., triangleN_edge2]

  * matrix
    * Size: ``(numberOfTriangles, 3)``
    * Structure::

         [[neighbor_e0, neighbor_e1, neighbor_e2],
          [neighbor_e0, neighbor_e1, neighbor_e2],
          ...]

``triangleSurfaces``
^^^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``float``
* Description:
  Surface area of each triangle.

* Supported layout:

  * flat
    * Size: ``numberOfTriangles``
    * Structure::

         [surface_0, surface_1, ..., surface_N]

``triangleCenterX``
^^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``float``
* Description:
  x-coordinates of the triangle center points.

* Supported layout:

  * flat
    * Size: ``numberOfTriangles``
    * Structure::

         [x_0, x_1, ..., x_N]

``triangleCenterY``
^^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``float``
* Description:
  y-coordinates of the triangle center points.

* Supported layout:

  * flat
    * Size: ``numberOfTriangles``
    * Structure::

         [y_0, y_1, ..., y_N]

``triangleNormalPoint``
^^^^^^^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``unsigned``
* Value range: ``{0, ..., numberOfPoints}``
* Description:
  Point indices where the corresponding triangle normal vectors start.
  For each triangle, one point index is stored for each of its three edges.

* Supported layouts:

  * flat
    * Size: ``3 * numberOfTriangles``
    * Structure::

         [triangle0_edge0, triangle1_edge0, ..., triangleN_edge0,
          triangle0_edge1, triangle1_edge1, ..., triangleN_edge1,
          triangle0_edge2, triangle1_edge2, ..., triangleN_edge2]

  * matrix
    * Size: ``(numberOfTriangles, 3)``
    * Structure::

         [[p_e0, p_e1, p_e2],
          [p_e0, p_e1, p_e2],
          ...]

``forbiddenEdge``
^^^^^^^^^^^^^^^^^

* Container type: ``list`` or ``ndarray``
* Element type: ``int``
* Value range: ``{-1, 0, 1, 2}``
* Description:
  Edge indices of the adjacent triangles.
  ``-1`` means there is no adjacent triangle at that edge.
  ``0``, ``1`` and ``2`` denote the corresponding edge index in the adjacent
  triangle.

* Supported layouts:

  * flat
    * Size: ``3 * numberOfTriangles``
    * Structure::

         [triangle0_edge0, triangle1_edge0, ..., triangleN_edge0,
          triangle0_edge1, triangle1_edge1, ..., triangleN_edge1,
          triangle0_edge2, triangle1_edge2, ..., triangleN_edge2]

  * matrix
    * Size: ``(numberOfTriangles, 3)``
    * Structure::

         [[edge_e0, edge_e1, edge_e2],
          [edge_e0, edge_e1, edge_e2],
          ...]

``minRaysPerSample``
^^^^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Description:
  Minimum number of rays used for adaptive sampling.

``maxRaysPerSample``
^^^^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Description:
  Maximum number of rays used for adaptive sampling.

``mseThreshold``
^^^^^^^^^^^^^^^^

* Type: ``float``
* Description:
  Sets the maximal MSE of the ASE value. If a sample point does not reach this
  threshold, the number of rays per sample point is increased up to
  ``maxRaysPerSample`` or resampled with repetitive sampling.

``repetitions``
^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Description:
  Maximum number of repetitions if the MSE threshold is not reached.

``adaptiveSteps``
^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
* Description:
  Sets the number of adaptive steps. The range between minimum and maximum ray
  count is split into that many parts. Setting it to ``1`` disables adaptive
  stepping and uses only ``minRaysPerSample``.

``nTot``
^^^^^^^^

* Type: ``float``
* Description:
  Doping of the active gain medium.

``thickness``
^^^^^^^^^^^^^

* Type: ``float``
* Description:
  Thickness of one prism layer of the mesh.

``laserParameter``
^^^^^^^^^^^^^^^^^^

* Type: ``dict``
* Description:
  Structure containing the laser parameters.

* Required fields:

  * ``l_abs``
  * ``l_ems``
  * ``s_abs``
  * ``s_ems``
  * ``l_res``

* Field layouts:

  * ``l_abs``
    * Container type: ``list`` or ``ndarray``
    * Element type: ``float``
    * Supported layout: flat

  * ``l_ems``
    * Container type: ``list`` or ``ndarray``
    * Element type: ``float``
    * Supported layout: flat

  * ``s_abs``
    * Container type: ``list`` or ``ndarray``
    * Element type: ``float``
    * Supported layout: flat

  * ``s_ems``
    * Container type: ``list`` or ``ndarray``
    * Element type: ``float``
    * Supported layout: flat

  * ``l_res``
    * Type: scalar

``crystal``
^^^^^^^^^^^

* Type: structure-like input
* Description:
  Structure containing the crystal parameters.

``numberOfLevels``
^^^^^^^^^^^^^^^^^^

* Type: ``unsigned``
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

  * ``"threaded"`` - threads distribute samples across GPUs/devices within one node
  * ``"mpi"`` - distributes samples across nodes;
  * ``"graybat"`` - similar to MPI but supports a more refined node topology (experimental)

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