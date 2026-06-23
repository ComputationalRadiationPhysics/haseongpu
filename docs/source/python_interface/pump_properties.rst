PumpProperties
==============

``PumpProperties`` describes the pump contribution to the beta population
:math:`\beta` during one simulation time step.

.. code-block:: python

   from HASEonGPU import PumpProperties

   pump = PumpProperties(
       spectralProperties=spectra,
       intensity=16e3,
       pumpSubsteps=100,
       wavelength=940e-9,
       radiusX=1.5,
       radiusY=1.5,
       exponent=40,
   )

Inputs and Defaults
-------------------

``PumpProperties`` stores both generic pump information and the parameters used
by the built-in Gaussian pump solver.  The direct constructor is flexible:
custom pump solvers do not need to provide Gaussian beam parameters.

``intensity``
   Required by the constructor.  Pump intensity :math:`I` in ``W / cm^2``.

``spectralProperties`` or ``crossSections``
   Required in practice.  Pass a ``SpectralDecomposition`` object, or provide
   monochromatic ``crossSectionAbsorption`` and ``crossSectionEmission`` values.

``wavelength``
   Pump wavelength :math:`\lambda`.  It is required when monochromatic cross
   sections are supplied.  With ``spectralProperties`` or ``crossSections`` it
   can be omitted; in that case the first absorption wavelength is used.  The
   ``PumpProperties.superGaussian(...)`` helper requires ``wavelength``
   explicitly in its signature.

``radiusX``
   Required only for the built-in Gaussian pump solver and for Gaussian helper
   methods such as ``intensityAt(...)`` and ``toDict(...)``.  It is the beam
   radius parameter used in :math:`I(x, y)`.  Custom pump solvers can omit it.

``radiusY``
   Optional.  When omitted, it defaults to ``radiusX`` for the built-in
   Gaussian pump solver.

``exponent``
   Optional.  Super-Gaussian exponent in :math:`I(x, y)`.  Defaults to
   ``40.0``.

``pumpSubsteps``
   Optional.  Number of time samples used by the default pump integration
   inside one pumped simulation step.  Defaults to ``100`` and must be at
   least 2 when supplied.

   This is not the number of outer simulation steps that receive pump energy.
   Use ``Simulation.runSteps(steps, pumpSteps=...)`` or set ``pumpSteps`` on
   ``PumpProperties`` for that.  For example, ``runSteps(150, pumpSteps=50)``
   pumps the first 50 outer steps and then continues ASE and fluorescence
   without pump contribution for the remaining 100 steps.

``pumpSteps``
   Optional custom property.  Number of outer simulation steps with pump
   contribution when ``Simulation.runSteps(steps)`` is called without an
   explicit ``pumpSteps`` override.  Omit it, or set it to ``None``, to pump for
   every outer step.

``solver``
   Optional.  Custom pump solver.  If omitted, ``Simulation`` uses the built-in
   ``BetaIntegrationGaussianSolver``.

Construction Helpers
--------------------

Use the direct constructor when you want to pass custom properties freely:

.. code-block:: python

   pump = PumpProperties(
       spectralProperties=spectra,
       intensity=16e3,
       wavelength=940e-9,
       radiusX=1.5,
       myCustomVar=6,
   )

Use ``superGaussian`` for an explicit beam-shaped setup:

.. code-block:: python

   pump = PumpProperties.superGaussian(
       spectralProperties=spectra,
       intensity=16e3,
       wavelength=940e-9,
       radiusX=1.5,
       radiusY=1.5,
       exponent=40,
       backReflection=True,
       reflectivity=1.0,
   )

Default Pump Routine
--------------------

If ``pump.solver`` is not set, ``Simulation`` uses
``BetaIntegrationGaussianSolver``.  The default routine:

1. Reads the current ``betaCells`` array, the point-wise :math:`\beta_i`, with shape
   ``(numberOfPoints, numberOfLevels)``.
2. Evaluates the super-Gaussian intensity profile at every transverse topology
   point.
3. Propagates pump intensity :math:`I` along the crystal levels.
4. Optionally applies a backward reflected pump contribution when
   ``backReflection`` is true.
5. Integrates the local beta update :math:`d\beta/dt` over ``pumpSubsteps``.
6. Returns the beta array :math:`\beta` after pumping.

The intensity profile :math:`I(x, y)` is:

.. code-block:: python

   I(x, y) = intensity * exp(-r ** exponent)
   r = sqrt(x**2 / radiusY**2 + y**2 / radiusX**2)

You can evaluate it directly:

.. code-block:: python

   intensities = pump.intensityAt(medium.topology.points)

Custom Pump Solvers
-------------------

A pump solver is any object with a ``step(input, pump)`` method.  It receives a
dictionary and the ``PumpProperties`` object.  It must return updated beta
values :math:`\beta` with the same shape as ``input["betaCell"]``.

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

The ``input`` dictionary currently contains:

* ``"betaCell"``: current beta array :math:`\beta`.
* ``"_medium"``: the ``GainMedium``.
* ``"_timeStep"``: simulation time step.
* ``"_constants"``: physical constants object.
* ``"_substeps"``: optional substep override.

Properties and Utilities
------------------------

Custom values are stored in ``customProperties`` and can be accessed in three
ways:

.. code-block:: python

   pump.getProperty("myCustomVar")
   pump.withProperty("myCustomVar", 7)
   pump.withProperties(myCustomVar=8, anotherValue=1.0)

Common attributes include:

* ``crossSections`` / ``spectralProperties``
* ``radiusX`` and ``radiusY``
* ``exponent``
* ``pumpDuration`` or ``duration``
* ``temporaryFluorescence``
* ``backReflection``
* ``reflectivity``
* ``extraction``
* ``solver``

``toDict(timeFrame=None)`` produces the low-level pump dictionary used by the
default solver.  ``modeDict()`` returns the pump mode flags for back reflection,
reflectivity, and extraction.
