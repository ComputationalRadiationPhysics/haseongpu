SpectralDecomposition
=====================

``SpectralDecomposition`` stores absorption and emission spectra,
:math:`\sigma_a(\lambda)` and :math:`\sigma_e(\lambda)`.  It is an alias for
``CrossSectionData``.

.. code-block:: python

   from HASEonGPU import SpectralDecomposition

   spectra = SpectralDecomposition(
       wavelengthsAbsorption=[900.0, 910.0],
       crossSectionAbsorption=[1.1e-21, 1.2e-21],
       wavelengthsEmission=[1020.0, 1030.0],
       crossSectionEmission=[2.0e-20, 2.48e-20],
       resolution=2,
   )

The wavelength arrays :math:`\lambda` and matching cross-section arrays must
have the same length.  ``resolution`` is the spectral interpolation resolution
passed to the ASE calculation.

Constructors
------------

``SpectralDecomposition(...)``
   Builds spectral data from explicit arrays.

``SpectralDecomposition.monochromatic(...)``
   Creates a single-wavelength :math:`\lambda` data set:

.. code-block:: python

   spectra = SpectralDecomposition.monochromatic(
       wavelength=940e-9,
       crossSectionAbsorption=1.2e-21,
       crossSectionEmission=2.0e-20,
   )

``SpectralDecomposition.fromDirectory(path, resolution=1000)``
   Loads four text files from a directory:

* ``lambda_a.txt``
* ``sigma_a.txt``
* ``lambda_e.txt``
* ``sigma_e.txt``

Interpolation Utilities
-----------------------

.. code-block:: python

   sigma_abs = spectra.absorptionAt(940e-9)
   sigma_ems = spectra.emissionAt(1030.0)

The interpolation helpers accept wavelengths :math:`\lambda` in the same unit
as the stored data.  They also handle the common case where stored spectra are
in ``nm`` and a pump wavelength is supplied in ``m``.

Conversion Utilities
--------------------

``toDict()`` returns the backend-compatible laser property dictionary:

.. code-block:: python

   data = spectra.toDict()
   data["l_abs"]
   data["s_abs"]
   data["l_ems"]
   data["s_ems"]
   data["l_res"]

``toLaserProperties()`` wraps the same data in ``LaserProperties``:

.. code-block:: python

   laser = spectra.toLaserProperties()
   laser.maxSigmaA
   laser.maxSigmaE

``LaserProperties`` remains available for compatibility workflows, but normal
simulation code passes ``SpectralDecomposition`` directly to
``PumpProperties``, ``PhiASE``, or ``Simulation``.
