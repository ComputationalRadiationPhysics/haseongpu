Volume topology
===============

``VolumeTopology`` is the geometry contract for current HASEonGPU transport.
It stores an explicit unstructured Tet4 mesh; ASE, pump, and time-dependent
state are cell-centered. It does not store excitation, material constants,
spectra, or boundary optics. Only VTK cell type ``10`` (Tet4) is supported.

Construction
------------

Load gmsh, VTK, or a closed STL surface:

.. code-block:: python

   from HASEonGPU import VolumeTopology

   topology = VolumeTopology.fromFile("crystal.msh")
   topology = VolumeTopology.fromFile("crystal.vtk")
   topology = VolumeTopology.fromFile("crystal.stl", meshSize=0.05)

Use ``format=`` when a filename has a non-standard extension. ``fromVtk``
expects an ASCII VTK unstructured grid. ``fromStl`` uses gmsh to tetrahedralize
a closed three-dimensional surface and warns that HASEonGPU does not perform a
complete mesh-validity proof.

Construct a topology directly when another mesher already provides arrays:

.. code-block:: python

   topology = VolumeTopology.fromTetrahedra(
       points,                 # (numberOfPoints, 3)
       cellPointIndices,       # (numberOfCells, 4)
       cellDomains=cell_ids,   # optional (numberOfCells,)
       faceBoundaries=faces,   # optional (numberOfCells, 4)
   )

Derived geometry
----------------

Construction derives and validates the arrays needed by transport:

``cellPointIndices``
   Four point indices per cell.

``facePointIndices``
   Three point indices for each of the four local faces.

``neighborCells`` and ``neighborLocalFaces``
   The adjacent cell and its matching local face, or a negative value at the
   exterior boundary.

``cellCenters`` and ``cellVolumes``
   Cell geometry used for state placement, source weighting, pump-rate
   normalization, and volume-weighted post-processing.

``faceCenters``, ``faceNormals``, and ``faceAreas``
   Oriented face geometry used by boundary optics and pump injection.

The main size queries are ``numberOfPoints``, ``numberOfCells``,
``numberOfFacesPerCell``, and ``numberOfSamplePoints``. ``samplePoints`` equals
the cell centers in an explicit volume topology.

Named domains
-------------

Domains are positive integer labels. Cell domains identify volume regions;
surface domains identify faces used by pump injection, relays, or optical
boundaries. gmsh physical names are retained and can be used instead of numeric
tags. A domain label is geometry metadata: it acquires physical meaning when a
``GainMedium``, pump injector, or relay refers to it.

.. code-block:: python

   topology = (
       topology
       .withCellDomains(where="all", domain=1, name="gain")
       .withSurfaceDomains([
           {"where": "z_min", "domain": 10, "name": "pump_input"},
           {"where": "z_max", "domain": 11, "name": "pump_output"},
       ])
   )

``withCellDomains`` accepts cell indices, ``where="all"``, and gmsh physical
names or tags. ``withSurfaceDomains`` additionally accepts face indices,
``z_min``, ``z_max``, and all exterior faces. It rejects internal faces unless
``allowInternal=True`` is explicit. Both methods return a copied topology, so
the input object remains unchanged.

Resolve names with ``cellDomainMap()`` or ``surfaceDomainMap()``:

.. code-block:: python

   entry_id = topology.surfaceDomainMap().resolve("pump_input")

Boundary optics belong to ``GainMedium`` because they are physical fields, not
connectivity. The topology supplies only the name that connects an optical
assignment to a set of faces:

.. code-block:: python

   from HASEonGPU import GainMedium, SurfaceOptics

   medium = GainMedium(topology).with_surface_optics({
       "pump_input": SurfaceOptics(
           reflectivity=0.0, n_inside=1.83, n_outside=1.0
       )
   })

See :doc:`gain_medium` for ``SurfaceOptics`` syntax and
:ref:`ase-surface-reflections` for the implemented boundary physics.

VTK state input
---------------

``VolumeTopology.fromVtk`` reads geometry only. Use ``GainMedium.fromVtk``
when a Tet4 VTK file also contains ``betaVolume`` and supported physical
fields:

.. code-block:: python

   topology = VolumeTopology.fromVtk("geometry.vtk")
   medium = GainMedium.fromVtk("prepared-state.vtk")

This distinction prevents a geometry loader from silently becoming a material
or initial-state loader.

Legacy extruded topology
------------------------

``Grid`` and ``MeshTopology`` remain available for compatibility with planar
triangle meshes extruded into wedge layers. They expose legacy queries such as
``numberOfTriangles``, ``numberOfLevels``, and ``numberOfPrisms``. New
forward-volume simulations should use ``VolumeTopology``; VTK Tet4 input is
intentionally rejected by ``MeshTopology.fromVtk``.
