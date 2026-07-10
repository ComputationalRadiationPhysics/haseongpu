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

``Grid`` is the shortest path to a usable topology for three-dimensional
rectangular media.

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

From a legacy planar file:

.. code-block:: python

   topology = MeshTopology.fromFile("mesh.msh", numberOfLevels=5, thickness=0.1)

Supported ``MeshTopology`` mesh formats are:

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

``GainMedium.fromVtk(...)`` expects a Tet4 VTK unstructured grid, the same
format accepted by ``VolumeTopology.fromVtk(...)``. ``betaCells`` is read from
point or cell data, ``betaVolume`` from cell data, and ``claddingCellTypes``,
``refractiveIndices``, ``reflectivities``, ``nTot``, ``crystalTFluo``,
``claddingNumber``, and ``claddingAbsorption`` from VTK ``FIELD`` arrays when
present. ``numberOfLevels`` and ``thickness`` are metadata in the Tet4 file and
are not accepted as import overrides.

The gmsh importer can also map physical groups whose names contain
``cladding`` to ``claddingCellTypes`` when the topology is used by
``GainMedium``.  The stored value is the gmsh physical tag for triangles in
matching two-dimensional physical groups; all other triangles keep the default
cladding type ``0``.

Domains
-------

Domains are integer labels attached to topology entities.  Cell domains label
volume cells; surface domains label boundary faces.  The labels are independent
of mesh ordering, so they are a stable way to assign material regions, cladding
groups, boundary optics, or later solver options to named parts of a mesh.

For gmsh input, HASEonGPU keeps physical group names as domain metadata:
three-dimensional physical groups become cell-domain names, and
two-dimensional physical groups become surface-domain names.  The names can be
resolved later, so code can refer to ``"gain"`` or ``"crystal_exit"`` instead
of hard-coding a physical tag.  Domains can also be assigned directly in Python
when a file format does not contain the required group labels.

.. code-block:: python

   from HASEonGPU import GainMedium, SurfaceOptics, VolumeTopology

   topology = (
       VolumeTopology.fromFile("crystal.vtk")
       .withCellDomains({"where": "all", "domain": 1, "name": "gain"})
       .withSurfaceDomains(
           [
               {"where": "z_min", "domain": 10, "name": "entry"},
               {"where": "z_max", "domain": 11, "name": "exit"},
           ]
       )
   )

   medium = GainMedium(topology).withSurfaceOptics(
       {
           "entry": SurfaceOptics(reflectivity=0.0, n_inside=1.83, n_outside=1.0),
           "exit": SurfaceOptics(reflectivity=0.0, n_inside=1.83, n_outside=1.0),
       }
   )

``withCellDomains(...)`` accepts assignments for cell indices, ``where="all"``,
gmsh physical names, or gmsh physical tags.  ``withSurfaceDomains(...)`` accepts
face indices, exterior z-plane selectors such as ``where="z_min"`` and
``where="z_max"``, all exterior faces, gmsh physical names, or gmsh physical
tags.  Surface assignments reject internal faces by default; pass
``allowInternal=True`` only when the caller intentionally labels internal
interfaces.

``SurfaceOptics`` fills the backend arrays used for reflective boundaries:
explicit reflectivity, refractive index inside the surface, and refractive
index outside the surface.  A reflectivity of ``0.0`` still allows total
internal reflection when ``n_inside`` and ``n_outside`` make the incident angle
supercritical; this matches the legacy reflection model.

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
