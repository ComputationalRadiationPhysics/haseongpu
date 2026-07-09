Topology
========

``MeshTopology`` describes the spatial discretization used by HASEonGPU.  It is
the user-facing replacement for manually constructing point arrays, triangle
indices, derived triangle geometry, and z-level metadata.

The public import names are:

.. code-block:: python

   from HASEonGPU import GainMedium, Grid, MeshTopology, VolumeTopology

``GainMediumGeometry`` is an alias for ``MeshTopology``.

Grid
----

``Grid`` is the shortest path to a usable topology for rectangular media.

.. code-block:: python

   grid = Grid(
       xExtent=4.0,
       yExtent=4.0,
       zExtent=0.7,
       tileSizeX=0.25,
       tileSizeZ=0.7 / 9.0,
   )
   topology = MeshTopology.fromGrid(grid)

``tileSizeY`` defaults to ``tileSizeX``.  ``tileSizeZ`` defaults to
``tileSizeX`` if it is not supplied.

Useful ``Grid`` members:

* ``numberOfLevels``: number of z-levels generated from ``zExtent`` and
  ``tileSizeZ``.
* ``thickness``: z-spacing :math:`\Delta z` between levels.  This is the value
  passed to HASEonGPU as the layer thickness.
* ``constructPoints()``: returns the transverse point array with shape
  ``(numberOfPoints, 2)``.

MeshTopology Construction
-------------------------

From a grid:

.. code-block:: python

   topology = MeshTopology.fromGrid(grid)

From points:

.. code-block:: python

   points = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
   topology = MeshTopology.fromPoints(points, numberOfLevels=5).withThickness(0.1)

``fromPoints`` performs a Delaunay triangulation of the 2D points.

From a legacy planar file:

.. code-block:: python

   topology = MeshTopology.fromFile("medium.vtk")
   topology = MeshTopology.fromFile("mesh.msh", numberOfLevels=5, thickness=0.1)

Supported ``MeshTopology`` mesh formats are:

* ``vtk`` or ``legacy-vtk``: legacy ASCII VTK unstructured grids containing
  wedge cells are accepted.  The importer reads the wedge cells, extracts the
  transverse ``(x, y)`` points from the first z-level, reconstructs the base
  triangle connectivity from the lower three nodes of each wedge, and infers
  ``numberOfLevels`` and ``thickness`` from the point z-coordinates.
* ``msh`` or ``gmsh``: gmsh 2D triangle meshes are supported through the Python
  ``gmsh`` package.  The mesh must be planar and must contain triangle
  elements.  ``numberOfLevels`` and ``thickness`` are required.

Tet4 volume files, including closed 3D STL surfaces, are loaded as
``VolumeTopology``:

.. code-block:: python

   topology = VolumeTopology.fromFile("mesh.stl", meshSize=0.05)
   topology = VolumeTopology.fromFile("mesh.msh")
   topology = VolumeTopology.fromFile("mesh.vtk")

STL volume import expects a closed 3D surface suitable for Tet4 volume meshing
and emits a warning that HASEonGPU does not run a full tetrahedral mesh
validation pass.

``MeshTopology.fromFile(...)`` and ``VolumeTopology.fromFile(...)`` can override
auto-detection with ``format=``.  This is useful for temporary files or
non-standard extensions:

.. code-block:: python

   topology = MeshTopology.fromFile(
       "surface.mesh",
       format="gmsh",
       numberOfLevels=8,
       thickness=0.05,
   )
   volume = VolumeTopology.fromFile("volume.mesh", format="stl", meshSize=0.05)

The topology importer only creates geometry and layer metadata.  It does not
populate material arrays such as ``betaCells`` or ``reflectivities`` except for
metadata that can later be used by ``GainMedium``.  For VTK files that contain
both geometry and HASEonGPU field data, load a full gain medium instead:

.. code-block:: python

   medium = GainMedium.fromVtk("medium.vtk")
   medium = GainMedium.fromFile("medium.vtk")

``GainMedium.fromVtk(...)`` reads the same legacy wedge geometry as
``MeshTopology.fromFile(...)`` and then imports recognized HASEonGPU fields.
``betaCells`` is read from point data, ``betaVolume`` from cell data, and
``claddingCellTypes``, ``refractiveIndices``, ``reflectivities``, ``nTot``,
``crystalTFluo``, ``claddingNumber``, and ``claddingAbsorption`` from VTK
``FIELD`` arrays when present.  ``numberOfLevels`` and ``thickness`` may be
passed as overrides after loading.

The gmsh importer can also map physical groups whose names contain
``cladding`` to ``claddingCellTypes`` when the topology is used by
``GainMedium``.  The stored value is the gmsh physical tag for triangles in
matching two-dimensional physical groups; all other triangles keep the default
cladding type ``0``.

Shape and Size Queries
----------------------

``MeshTopology`` exposes the dimensions needed by other objects:

.. code-block:: python

   topology.numberOfPoints
   topology.numberOfTriangles
   topology.numberOfLevels(10)        # sets levels and returns topology
   topology.numberOfPrisms            # triangles * (levels - 1)
   topology.levelCoordinates()        # z coordinates

``numberOfLevels`` is both a construction parameter name and a setter method on
``MeshTopology``.  After construction, use it as:

.. code-block:: python

   topology.numberOfLevels(10).withThickness(0.05)

Use ``GainMedium`` property metadata for array shapes:

.. code-block:: python

   medium = topology.asGainMedium()
   medium.get("betaCells").expectedShape       # (numberOfPoints, levels)
   medium.get("betaVolume").expectedShape      # (numberOfTriangles, levels - 1)
   medium.get("reflectivities").expectedShape  # (numberOfTriangles, 2)

Index Utilities
---------------

These helpers convert physical coordinates to topology indices:

.. code-block:: python

   point_index = topology.pointIndexAt(x=0.0, y=0.0)
   level_index = topology.levelIndexAt(z=0.1)
   point_index, level_index = topology.betaCellIndexAt(0.0, 0.0, 0.1)
   flat_index = topology.betaCellIndexAt(0.0, 0.0, 0.1, flat=True)

``betaCellIndexAt`` is useful when initializing or inspecting one
:math:`\beta_i` entry in the ``betaCells`` array.

Conversion to GainMedium
------------------------

.. code-block:: python

   medium = topology.asGainMedium()

This is equivalent to:

.. code-block:: python

   medium = GainMedium(topology=topology)
