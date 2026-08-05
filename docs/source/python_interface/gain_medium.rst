GainMedium
==========

``GainMedium`` attaches physical fields and evolving excitation to a topology.
For ``VolumeTopology``, all transport state has exactly one value per Tet4 cell.
Geometry remains in ``medium.topology`` and is not duplicated in the fields.
Spectra and solver controls are separate objects because they describe material
response and numerical execution rather than mesh state.

.. code-block:: python

   import numpy as np
   from HASEonGPU import GainMedium

   medium = GainMedium(topology).withPhysicalProperties(
       betaVolume=np.zeros(topology.numberOfCells),
       claddingCellTypes=np.zeros(topology.numberOfCells, dtype=np.uint32),
       nTot=2.776e20,
       crystalTFluo=9.41e-4,
       claddingNumber=1,
       claddingAbsorption=5.5,
   )

Built-in fields
---------------

``betaVolume``
   Dimensionless, cell-centered excited-state fraction :math:`\beta_j`. This is
   the authoritative state used by ASE, pump transport, and time integration.

``claddingCellTypes``
   Unsigned cell type used to select the cells governed by
   ``claddingAbsorption`` rather than the active-medium gain coefficient. A cell
   is treated as cladding when its value equals ``claddingNumber``; other cells
   use :math:`\beta_j`, ``nTot``, and the spectral cross sections.

``nTot``
   Total active-ion concentration :math:`N_{\mathrm{tot}}` in ``cm^-3``.

``crystalTFluo``
   Fluorescence lifetime :math:`\tau` in seconds.

``claddingNumber`` and ``claddingAbsorption``
   Selected cladding type and its absorption coefficient. The coefficient's
   inverse-length unit must match the geometry length unit.

``surfaceReflectivity``, ``surfaceRefractiveIndexInside``, and ``surfaceRefractiveIndexOutside``
   Arrays indexed by positive surface-domain id. Prefer
   ``withSurfaceOptics`` to construct these arrays by domain name.

``refractiveIndices`` and ``reflectivities``
   Compatibility fields for the planar top/bottom surface layout. Explicit
   Tet4 boundaries use the domain-indexed surface fields above.

The equations that consume these fields are collected in
:doc:`../theoryAndModel`; the field-to-symbol table there maps the Python names
to the ASE, pump, and population-rate equations.

Inspect and set fields
----------------------

``get(name)`` returns a property handle with dtype, expected shape, value, and
metadata:

.. code-block:: python

   beta = medium.get("betaVolume")
   print(beta.dtype, beta.expectedShape, beta.value)
   beta.value = np.zeros(beta.expectedShape)

   medium.set("nTot", 2.776e20)
   for item in medium.listProperties():
       print(item["name"], item["expectedShape"], item["isSet"])

Array input should use the reported primitive shape. Multi-dimensional arrays
are flattened in Fortran order at the backend boundary. A one-dimensional
array for a multi-dimensional legacy field is ambiguous; wrap an already
canonical flat array with ``backendFlat(values)``.

Surface optics
--------------

.. code-block:: python

   from HASEonGPU import SurfaceOptics

   medium.withSurfaceOptics({
       "input": SurfaceOptics(
           reflectivity=0.0, n_inside=1.83, n_outside=1.0
       ),
       "mirror": SurfaceOptics(
           reflectivity=0.98, n_inside=1.83, n_outside=1.0
       ),
   })

Names are resolved through ``topology.surfaceDomainMap()``. Index zero is
unused because physical domain ids are positive. A configured reflectivity of
zero still permits total internal reflection when ``n_inside`` and ``n_outside``
make the incident angle supercritical. See :ref:`ase-surface-reflections` for
the implemented model and its limits.

Custom openPMD fields
---------------------

``defineField`` writes application data next to HASE records. It does not make
the current backend consume that field. Use ``entity="cell"`` for an explicit
volume cell array and provide unit metadata:

.. code-block:: python

   medium.defineField(
       "temperature",
       entity="cell",
       values=np.full(topology.numberOfCells, 300.0),
       unit="K",
       unitDimension=(0, 0, 0, 0, 1, 0, 0),
   )

Other supported entity axes include ``point``, ``("cell", "local_vertex")``,
``("cell", "local_side")``, ``interface``, ``wavelength``, and ``domain``. Marking a
custom field ``backendRequired=True`` causes transport validation to reject a
backend that does not declare support instead of silently treating it as
analysis-only metadata.

VTK round trip
--------------

.. code-block:: python

   medium.toVtk("prepared-state.vtk")
   restored = GainMedium.fromVtk("prepared-state.vtk")

These methods use ASCII Tet4 VTK. ``GainMedium.fromVtk`` restores supported
geometry and field records; use ``VolumeTopology.fromVtk`` when only geometry
is wanted.
