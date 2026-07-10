PumpProperties and Pumping
==========================

Pumping raises the excited-state fraction ``betaCells`` during a time-dependent
simulation.  The Python interface separates the pump description from the ASE
backend:

* ``PumpRadiationProfile`` describes the incoming beam.
* ``PumpProperties`` stores spectra, profile, solver, and custom solver values.
* A solver updates beta or provides a pump rate for ``Simulation``.

Typical Setup
-------------

.. code-block:: python

   from HASEonGPU import OneDimensionalZTraversal, PumpProperties, PumpRadiationProfile

   profile = PumpRadiationProfile(
       intensity=16e3,              # W / cm^2
       wavelengths=[940e-9],
       waist=(1.5, 1.5),
       propagationDirection=(0.0, 0.0, 1.0),
       superGaussianOrder=40,
       backReflection=True,
       reflectivity=1.0,
   )

   pump = PumpProperties(
       crossSections=spectra,
       profile=profile,
       solver=OneDimensionalZTraversal(),
       pumpSteps=50,
   )

``spectra`` is the same ``CrossSectionData`` or ``SpectralDecomposition`` object
used by ``PhiASE``.  The profile supplies pump wavelengths; the spectra supply
``sigma_a`` and ``sigma_e`` at those wavelengths.

Important Inputs
----------------

``intensity``
   Pump intensity in ``W / cm^2``.  It can be passed directly or through a
   ``PumpRadiationProfile``.

``profile``
   Optional radiation-profile object.  Dictionaries, dataclasses, and objects
   with public attributes are accepted.

``spectralProperties`` / ``crossSections``
   Material spectra.  Required for physical pump solvers.

``wavelength`` / ``wavelengths``
   One wavelength or multiple pump wavelengths.  If omitted, the first
   absorption wavelength from the spectra is used.

``radiusX`` / ``radiusY`` / ``waist``
   Beam radii for the super-Gaussian transverse profile.

``superGaussianOrder``
   Beam order.  Larger values approach a top-hat profile; ``2`` is
   Gaussian-like.

``propagationDirection``
   Direction vector for ``OneDimensionalZTraversal``.  The current built-in
   traversal uses the sign of the z component.

``backReflection`` / ``reflectivity``
   Enable and scale the reflected pump pass.

``pumpSubsteps``
   Internal substeps for ``BetaIntegrationGaussianSolver`` only.

``pumpSteps``
   Number of outer ``Simulation`` steps with active pumping when
   ``runSteps(...)`` does not receive an explicit ``pumpSteps=`` override.

``solver``
   Pump solver object.  If omitted, ``Simulation`` uses
   ``BetaIntegrationGaussianSolver``.

Built-In Solvers
----------------

``OneDimensionalZTraversal``
   Preferred continuous z-direction pump solver.  It evaluates an instantaneous
   pump rate from the current beta state, spectra, profile, and optional back
   reflection.  ``Simulation`` stores the rate in ``TimeStepState.dndtPump``.
   Use ``oneDimensionalZTraversalPumpRate(...)`` directly for diagnostics.

``BetaIntegrationGaussianSolver``
   Legacy/default super-Gaussian solver.  It advances beta analytically inside
   a pumped outer step and uses ``pumpSubsteps``.  The compatibility helpers
   ``integrateLaserPump`` and ``runLaserPumpStep`` use the same historical path.

For equations and model details, see :doc:`../theoryAndModel`.

Custom Pump Solvers
-------------------

A custom solver only needs ``step(input, pump)``.  It receives a dictionary and
the ``PumpProperties`` object, and returns beta values with the same shape as
``input["betaCell"]``.

.. code-block:: python

   class MyPumpSolver:
       def step(self, input, pump):
           beta = input["betaCell"]
           scale = pump.getProperty("scale", 1.0)
           return np.clip(beta + scale * 1e-3, 0.0, 1.0)

   pump = PumpProperties(
       spectralProperties=spectra,
       intensity=16e3,
       wavelength=940e-9,
       radiusX=1.5,
       radiusY=1.5,
       solver=MyPumpSolver(),
       scale=2.0,
   )

The input dictionary includes the current beta array, the ``GainMedium``, time
step, constants object, and optional substep override.  Custom values are read
with ``pump.getProperty(...)``.

Utilities
---------

.. code-block:: python

   pump.getProperty("myCustomVar", default)
   pump.withProperty("myCustomVar", 7)
   pump.withProperties(myCustomVar=8, anotherValue=1.0)
   pump.intensityAt(points)
   pump.toDict(timeFrame=None)
   pump.modeDict()

``toDict`` and ``modeDict`` are mainly for compatibility with the legacy pump
path.
