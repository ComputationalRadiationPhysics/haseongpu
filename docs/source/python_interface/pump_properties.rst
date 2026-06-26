PumpProperties and Pumping
==========================

Pumping is the part of a time-dependent HASEonGPU simulation that raises the
excited-state fraction :math:`\beta` in the gain medium.  The ASE backend uses
:math:`\beta` to estimate amplified spontaneous-emission flux, while the Python
``Simulation`` loop evolves :math:`\beta` by combining pump excitation, ASE
depletion, and fluorescence decay.

In the Python interface the pump is intentionally split into small objects:

* ``PumpRadiationProfile`` describes the incoming radiation field: intensity,
  wavelengths, waist, propagation direction, spectral weights, and reflection.
* ``PumpProperties`` stores the profile, material spectra, solver choice, and
  custom solver parameters.
* A pump solver converts the current ``betaCells`` into either an updated beta
  state or, through ``Simulation``, a pump rate :math:`d\beta/dt`.

A typical continuous-pump setup is:

.. code-block:: python

   from HASEonGPU import (
       OneDimensionalZTraversal,
       PumpProperties,
       PumpRadiationProfile,
   )

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
used by ``PhiASE``.  The profile supplies wavelengths; the crystal spectra
supply :math:`\sigma_a(\lambda)` and :math:`\sigma_e(\lambda)`.

Inputs and Defaults
-------------------

``PumpProperties`` accepts core arguments plus arbitrary keyword properties.
Unknown keywords are stored in ``customProperties`` so that custom solvers can
read them later.

``intensity``
   Pump intensity :math:`I_0` in ``W / cm^2``.  It may be supplied directly or
   through a ``PumpRadiationProfile``.

``profile``
   Optional radiation-profile object.  ``PumpRadiationProfile`` is the built-in
   profile container, but dictionaries, dataclasses, and simple objects with
   public attributes are also accepted.

``spectralProperties`` or ``crossSections``
   Absorption and emission spectra.  These are required in practice.  If they
   are omitted, monochromatic ``crossSectionAbsorption`` and
   ``crossSectionEmission`` values plus ``wavelength`` must be provided.

``wavelength`` / ``wavelengths``
   ``wavelength`` selects one pump wavelength.  ``wavelengths`` is a custom
   property consumed by the continuous z-traversal solver for multi-wavelength
   pumping.  If no explicit wavelength is supplied, the first absorption
   wavelength in the spectra is used by ``PumpProperties``.

``radiusX`` / ``radiusY`` or ``waist``
   Beam radii for the super-Gaussian transverse profile.  ``waist=(rx, ry)`` is
   the profile-oriented spelling; ``radiusX`` and ``radiusY`` are the lower-level
   names.  If only one radius is supplied it is used in both transverse
   directions.

``superGaussianOrder``
   Super-Gaussian order :math:`q`.  The default is ``40.0``.  Larger values
   approximate a top-hat beam more closely; ``q=2`` gives a Gaussian-like
   profile.

``propagationDirection``
   Three-component direction vector used by ``OneDimensionalZTraversal``.  Only
   the sign of the normalized z component is currently used: positive z pumps
   from the first level toward the last level, and negative z pumps in the
   opposite direction.

``spectralWeights``
   Optional weights for each wavelength sample in a multi-wavelength continuous
   pump.  When omitted, every wavelength has unit weight.

``backReflection`` and ``reflectivity``
   Enable and scale the reflected pump pass.  For the continuous z traversal,
   the reflected intensity is propagated back through the same one-dimensional
   attenuation model.

``pumpSubsteps``
   Internal time samples for the legacy ``BetaIntegrationGaussianSolver`` only.
   It defaults to ``100`` and must be at least ``2``.  The continuous
   ``OneDimensionalZTraversal`` solver evaluates an instantaneous rate once per
   outer derivative evaluation and does not use substeps.

``pumpSteps``
   Optional custom property read by ``Simulation.runSteps``.  It limits pumping
   to the first ``pumpSteps`` outer simulation steps when ``runSteps`` is called
   without an explicit ``pumpSteps=`` override.  Omit it or set it to ``None`` to
   pump on every outer step.

``solver``
   Pump solver object.  If omitted, ``Simulation`` uses
   ``BetaIntegrationGaussianSolver``.

Continuous One-Dimensional Pump
-------------------------------

``OneDimensionalZTraversal`` is the preferred built-in solver for a continuous
pump field that travels along the extruded z direction of the gain medium.  It
keeps the current beta state fixed while it computes the pump contribution for
one derivative evaluation.

First the incoming intensity is distributed over the transverse mesh with a
super-Gaussian profile centered at ``center``:

.. math::

   I_0(x,y) = I_{\mathrm{peak}}
   \exp\left[-\left(\sqrt{\frac{(x-x_c)^2}{r_y^2}
   + \frac{(y-y_c)^2}{r_x^2}}\right)^q\right].

For each wavelength sample :math:`\lambda_k`, the solver propagates intensity
between adjacent z levels using the average excited-state fraction between
those levels:

.. math::

   I_{k,z+\Delta z} = I_{k,z}
   \exp\left[-\left(\sigma_{a,k}
   - \bar\beta_z(\sigma_{a,k}+\sigma_{e,k})\right)
   N_{\mathrm{tot}}\Delta z\right].

With back reflection enabled, a reflected pass is initialized at the far end of
the crystal with ``reflectivity`` times the transmitted intensity and is
propagated back through the same attenuation factors.  The total local pump
intensity is the sum of the forward and reflected passes.

The photon-flux conversion for each wavelength is

.. math::

   \Phi_k = I_k \frac{\lambda_k}{hc},

so the pump contribution to the population equation is

.. math::

   \left.\frac{d\beta}{dt}\right|_{\mathrm{pump}}
   =
   \sum_k
   \left[\sigma_{a,k} - \beta(\sigma_{a,k}+\sigma_{e,k})\right]
   \Phi_k.

The helper ``oneDimensionalZTraversalPumpRate(...)`` exposes this rate directly
for inspection or tests.  ``Simulation`` normally calls the solver through
``PumpProperties`` and stores the resulting rate in ``TimeStepState.dndtPump``.

Legacy Analytical Gaussian Pump
-------------------------------

When no solver is set, ``Simulation`` creates ``BetaIntegrationGaussianSolver``.
This is the historical pump integrator behind ``integrateLaserPump`` and
``runLaserPumpStep``.  It also uses a super-Gaussian transverse profile and a
one-dimensional z propagation, but it advances beta with an analytical update
inside one pumped outer step and samples that update over ``pumpSubsteps``.

A compact legacy setup is:

.. code-block:: python

   pump = PumpProperties.superGaussian(
       spectralProperties=spectra,
       intensity=16e3,
       wavelength=940e-9,
       radiusX=1.5,
       radiusY=1.5,
       superGaussianOrder=40,
       backReflection=True,
       reflectivity=1.0,
   )

For a fixed local intensity during one substep, the analytical update is

.. math::

   \beta(t + \Delta t)
   =
   \frac{A}{C}\left(1 - e^{-C\Delta t}\right)
   + \beta(t)e^{-C\Delta t},

with

.. math::

   A = \sigma_a I \frac{\lambda}{hc},
   \qquad
   C = (\sigma_a + \sigma_e)I\frac{\lambda}{hc} + \frac{1}{\tau}.

This solver consumes ``pumpSubsteps``, ``pumpDuration`` or the simulation time
step, ``temporaryFluorescence`` if supplied, and the low-level mode flags from
``modeDict()``.

Custom Pump Solvers
-------------------

A pump solver is any object with a ``step(input, pump)`` method.  It receives a
dictionary and the ``PumpProperties`` object.  It must return updated beta
values with the same shape as ``input["betaCell"]``.  ``Simulation`` converts
that update into a rate by dividing by the active pump duration.

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

* ``"betaCell"``: current point-level beta array :math:`\beta`.
* ``"_medium"``: the ``GainMedium``.
* ``"_timeStep"``: simulation time step.
* ``"_constants"``: physical constants object.
* ``"_substeps"``: optional substep override, used by the legacy solver.

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
* ``profile``
* ``radiusX`` and ``radiusY``
* ``superGaussianOrder``
* ``wavelength`` / ``wavelengths``
* ``spectralWeights``
* ``propagationDirection``
* ``pumpDuration`` or ``duration``
* ``temporaryFluorescence``
* ``backReflection``
* ``reflectivity``
* ``extraction``
* ``solver``

``intensityAt(points)`` evaluates the super-Gaussian transverse profile for
``(x, y)`` points.  ``toDict(timeFrame=None)`` produces the low-level pump
dictionary used by the legacy solver.  ``modeDict()`` returns the legacy mode
flags for back reflection, reflectivity, and extraction.
