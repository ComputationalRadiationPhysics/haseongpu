GainMedium
==========

``GainMedium`` combines a ``MeshTopology`` with the material and state arrays
that the ASE calculation needs.

.. code-block:: python

   from HASEonGPU import GainMedium

   medium = GainMedium(topology=topology)

Physical Properties
-------------------

The most important pattern is to ask the property object for its expected
shape, then create an array with that shape:

.. code-block:: python

   beta_shape = medium.get("betaCells").expectedShape
   medium.get("betaCells").value = np.zeros(beta_shape)

You can set several properties at once:

.. code-block:: python

   medium.withPhysicalProperties(
       betaCells=np.zeros(medium.get("betaCells").expectedShape),
       betaVolume=np.zeros(medium.get("betaVolume").expectedShape),
       claddingCellTypes=np.zeros(medium.get("claddingCellTypes").expectedShape, dtype=np.uint32),
       refractiveIndices=[2.0, 1.0, 2.0, 1.0],
       reflectivities=np.zeros(medium.get("reflectivities").expectedShape),
       nTot=2.776e20,
       crystalTFluo=9.41e-4,
       claddingNumber=1,
       claddingAbsorption=5.5,
   )

``withPhysicalProperties`` returns the same ``GainMedium`` instance, so it can
be chained.

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
   ``(2, numberOfTriangles)`` where row 0 is bottom and row 1 is top.

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
