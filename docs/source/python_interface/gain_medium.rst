GainMedium
==========

``GainMedium`` combines a ``MeshTopology`` with the material and state arrays
that the ASE calculation needs.

.. code-block:: python

   from HASEonGPU import GainMedium

   medium = GainMedium(topology=topology)

Required Fields
---------------

Assign the built-in transport fields directly on the medium.  For arrays, ask
the field for its primitive shape first so the code follows the topology:

.. code-block:: python

   medium.get("betaCells").value = np.zeros(medium.get("betaCells").expectedShape)
   medium.get("betaVolume").value = np.zeros(medium.get("betaVolume").expectedShape)
   medium.get("claddingCellTypes").value = np.zeros(
       medium.get("claddingCellTypes").expectedShape, dtype=np.uint32
   )
   medium.get("refractiveIndices").value = np.asarray([2.0, 1.0, 2.0, 1.0], dtype=np.float32)
   medium.get("reflectivities").value = np.zeros(
       medium.get("reflectivities").expectedShape, dtype=np.float32
   )
   medium.get("nTot").value = 2.776e20
   medium.get("crystalTFluo").value = 9.41e-4
   medium.get("claddingNumber").value = 1
   medium.get("claddingAbsorption").value = 5.5

The transport writer reads these named fields from the medium; examples do not
construct a separate adapter object for backend input.

Property Reference
------------------

``betaCells``
   Excited-state fraction :math:`\beta_i` at topology points and z-levels.
   Matrix shape: ``(numberOfPoints, numberOfLevels)``.

``betaVolume``
   Prism-centered excited-state fraction :math:`\beta_j` used by the ASE ray
   integration.  Matrix shape: ``(numberOfTriangles, numberOfLevels - 1)``.

``claddingCellTypes``
   Triangle-wise cladding type index.  Shape: ``(numberOfTriangles,)``.

``refractiveIndices``
   Four refractive indices:
   ``[bottomInside, bottomOutside, topInside, topOutside]``.

``reflectivities``
   Surface reflectivity per triangle.  Matrix shape:
   ``(numberOfTriangles, 2)`` where column 0 is bottom and column 1 is top.

``nTot``
   Total active-ion concentration :math:`N_{\mathrm{tot}}` in ``cm^-3``.

``crystalTFluo``
   Fluorescence lifetime :math:`\tau`.

``claddingNumber``
   Cladding type selected for cladding absorption handling.

``claddingAbsorption``
   Absorption coefficient of the selected cladding.

Shape and Metadata Utilities
----------------------------

``get(name)`` returns a property wrapper:

.. code-block:: python

   prop = medium.get("betaCells")
   prop.name
   prop.description
   prop.dtype
   prop.expectedShape
   prop.value
   prop.meta()

``listProperties()`` returns metadata for all known physical properties:

.. code-block:: python

   for prop in medium.listProperties():
       print(prop["name"], prop["expectedShape"], prop["isSet"])

``set(name, value)`` validates and stores one property:

.. code-block:: python

   medium.set("nTot", 2.776e20)
   medium.set("betaCells", np.zeros(medium.get("betaCells").expectedShape))

Arrays can be supplied either in matrix shape or flat Fortran order.  Stored
arrays are flattened internally for the HASEonGPU binding.

Indexing Helpers
----------------

``GainMedium`` forwards beta-cell coordinate lookups for :math:`\beta_i` to
its topology:

.. code-block:: python

   i, k = medium.betaCellIndexAt(x=0.0, y=0.0, z=0.0)
   beta = medium.get("betaCells").value.reshape(
       medium.get("betaCells").expectedShape,
       order="F",
   )
   beta[i, k] = 0.5
   medium.get("betaCells").value = beta

Convenience Dimensions
----------------------

.. code-block:: python

   medium.numberOfPoints
   medium.numberOfTriangles
   medium.numberOfPrisms
   medium.numberOfLevels

``emptyBetaCells(fill=0.0)`` creates a correctly shaped beta array
:math:`\beta_i`:

.. code-block:: python

   medium.get("betaCells").value = medium.emptyBetaCells(fill=0.0)

Custom Fields
-------------

Custom fields are primarily an openPMD extension point.  They let a Python
workflow write additional mesh records next to the HASEonGPU records so
downstream analysis tools, coupled codes, or future backends can read them.
The current ASE backend ignores custom records unless a backend explicitly opts
in.

``GainMedium.defineField(...)`` creates one additional openPMD mesh record.
Choose the entity from the location of the data:

* ``"point"`` for arrays shaped like ``(numberOfPoints, numberOfLevels)``
* ``"prism"`` for arrays shaped like ``(numberOfTriangles, numberOfLevels - 1)``
* ``"triangle"`` for arrays shaped like ``(numberOfTriangles,)``

Always provide unit metadata when the field has a physical meaning.  If omitted,
the transport writes ``unitSI=1.0`` and
``unitDimension=unitDimension.dimensionless``.  Unit dimensions follow the
standard seven-entry openPMD tuple.

.. code-block:: python

   import numpy as np
   from HASEonGPU import unitDimension

   temperature = np.full(medium.get("betaVolume").expectedShape, 300.0)

   medium.defineField(
       "temperature",
       entity="prism",
       values=temperature,
       unit="K",
       unitSI=1.0,
       unitDimension=(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
   )

Inheritance-based declarations use the same openPMD metadata on
``PrimitiveFieldSpec``:

.. code-block:: python

   import numpy as np
   from HASEonGPU import PrimitiveFieldSpec, PrismSchema

   class ThermalPrism(PrismSchema):
       temperature = PrimitiveFieldSpec(
           "temperature",
           "custom_temperature",
           np.float64,
           unit="K",
           unitSI=1.0,
           unitDimension=(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
           backendRequired=False,
       )

   medium.withPrimitiveSchema(ThermalPrism, temperature=temperature)

Predefined openPMD dimension tuples are available from ``unitDimension`` for
HASEonGPU fields and common dimensionless records.
