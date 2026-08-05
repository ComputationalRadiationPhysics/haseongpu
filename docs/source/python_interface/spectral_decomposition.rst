Cross-section data
==================

``CrossSectionData`` stores the wavelength-dependent absorption and emission
cross sections :math:`\sigma_a(\lambda)` and :math:`\sigma_e(\lambda)` used by
ASE and pump transport. ``SpectralDecomposition`` is a compatibility alias for
the same class; new examples use ``CrossSectionData``.

The cross sections are stored in ``cm^2``. Wavelength values retain the unit in
which they are supplied, so each wavelength table must be internally
consistent.

Load a measured spectrum
------------------------

``fromDirectory`` loads the four text files used by the
:doc:`laserPumpCladding tutorial <../laserPumpCladding>`:

.. code-block:: python

   from HASEonGPU import CrossSectionData

   ase_spectra = CrossSectionData.fromDirectory(
       "spectra",
       resolution=1000,
   )

The directory must contain:

* ``lambda_a.txt`` and ``sigma_a.txt`` for absorption;
* ``lambda_e.txt`` and ``sigma_e.txt`` for emission.

The wavelength and cross-section arrays in each pair must have equal lengths.
``resolution`` controls the numerical spectral table passed to the ASE backend;
it does not resample the input files or create additional measured information.

Explicit and monochromatic construction
---------------------------------------

Construct the same object directly from arrays when the data does not use the
four-file layout:

.. code-block:: python

   spectra = CrossSectionData(
       wavelengthsAbsorption=[930.0, 940.0, 950.0],
       crossSectionAbsorption=[6.7e-21, 7.8e-21, 8.1e-21],
       wavelengthsEmission=[930.0, 940.0, 950.0],
       crossSectionEmission=[1.5e-21, 1.9e-21, 2.4e-21],
       resolution=1000,
   )

A monochromatic physical pump needs one wavelength and the two material cross
sections at that wavelength:

.. code-block:: python

   wavelength = 940e-9
   pump_cross_sections = CrossSectionData.monochromatic(
       wavelength=wavelength,
       crossSectionAbsorption=spectra.absorptionAt(wavelength),
       crossSectionEmission=spectra.emissionAt(wavelength),
   )

``absorptionAt`` and ``emissionAt`` linearly interpolate sorted wavelength
samples. They also recognize the common case in which a material table is in
nanometers but the query is in meters. Outside that clear magnitude-based
conversion, pass the query in the same unit as the stored wavelengths.

Role in a simulation
--------------------

The full spectrum and the pump-specific data have different purposes:

.. code-block:: python

   phi_ase = PhiASE(spectralProperties=ase_spectra)
   pump = Pump(
       total_power=16_000.0,
       spectrum=PumpSpectrum.monochromatic(wavelength),
       cross_sections=pump_cross_sections,
   )
   simulation = Simulation(
       gain_medium=medium,
       phi_ase=phi_ase,
       time_integrator=FrozenPhiAseRungeKutta4(),
       time_step_size=2e-5,
       cross_sections=ase_spectra,
   )

``PhiASE`` samples the complete absorption/emission spectrum. ``Pump`` uses its
own spectrum to choose pump wavelengths and its own ``cross_sections`` to
evaluate material interaction at those wavelengths. Passing ``ase_spectra`` to
``Simulation`` supplies the spectral data used by ASE throughout the compiled
time loop.

Low-level compatibility
-----------------------

``toDict`` exposes the historical ``l_abs``, ``s_abs``, ``l_ems``, ``s_ems``,
and ``l_res`` names, and ``toLaserProperties`` wraps the same arrays in the
legacy ``LaserProperties`` store. They are compatibility interfaces; normal
composed simulations use ``CrossSectionData`` directly.
