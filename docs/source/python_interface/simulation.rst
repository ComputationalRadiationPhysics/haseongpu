Simulation
==========

``Simulation`` is the high-level time loop.  It combines a ``GainMedium``,
``PumpProperties``, ``PhiASE`` configuration, and a time integration solver.

.. code-block:: python

   from HASEonGPU import Simulation, RungeKutta4

   simulation = Simulation(
       gainMedium=medium,
       pump=pump,
       phiASE=phi_ase,
       timeIntegrationSolver=RungeKutta4(),
       timeStep=1e-5,
       endTime=1e-3,
   )

Running
-------

Run a fixed number of steps:

.. code-block:: python

   simulation.runSteps(3)

Run until a target time:

.. code-block:: python

   simulation.runUntil(endTime=1e-3)

If ``endTime`` was supplied to the constructor:

.. code-block:: python

   simulation.runUntil()

Run exactly one step and inspect the returned state:

.. code-block:: python

   state = simulation.step()
   print(state.step, state.time, state.betaCells.shape, state.phiAse.shape)

Simulation Step Order
---------------------

Each call to ``step()`` performs:

1. One-time ``onInit`` callbacks.
2. ``beforeStep`` callbacks.
3. Time integration of the beta derivative :math:`d\beta/dt`.
4. Pump contribution through the configured pump solver.
5. ASE contribution through ``phiASE.run(...)``.
6. Fluorescence decay using ``crystalTFluo`` (:math:`\tau`).
7. Clipping of updated beta values to ``[0, 1]``.
8. ``betaVolume`` update from ``betaCells``.
9. Result storage and ``onStep`` callbacks.

Callbacks
---------

``onInit(callback)``
   Registers a callback that receives the ``Simulation`` object before the
   first step.

``beforeStep(callback)``
   Registers a callback that receives the ``Simulation`` object before every
   step.

``onStep(callback)``
   Registers a callback that receives the completed ``TimeStepState``.

Example:

.. code-block:: python

   def print_state(state):
       print(state.step, state.time, state.betaCells.mean())

   simulation.onStep(print_state)

Callbacks return values are ignored, so they are best used for logging,
inspection, exporting, or controlled mutation of the simulation before a step.

Results
-------

``getResults()`` returns a list of ``TimeStepState`` snapshots:

.. code-block:: python

   results = simulation.getResults()
   last = results[-1]
   last.betaCells.shape
   last.betaVolume.shape
   last.phiAse.shape

``TimeStepState`` fields are:

* ``step``: completed step index.
* ``time``: simulation time after the step.
* ``betaCells``: point and level beta values :math:`\beta_i`.
* ``betaVolume``: prism beta values :math:`\beta_j` used by ASE.
* ``phiAse``: ASE flux values :math:`\Phi_i`, or ``None`` if unavailable.
* ``dndtAse``: ASE derivative contribution to :math:`d\beta/dt`.
* ``dndtPump``: pump derivative contribution to :math:`d\beta/dt`.
* ``aseResult``: raw lower-level ASE result object.

Properties
----------

.. code-block:: python

   simulation.time
   simulation.stepIndex

``time`` and ``stepIndex`` expose the current simulation clock and completed
step count.

Spectral Properties Resolution
------------------------------

``Simulation`` requires spectral properties such as cross sections and emission
spectra. These properties may be defined at different levels, for example
globally for the whole simulation or locally for individual objects.

When multiple definitions are available, the simulation resolves them using a
fixed priority order. More specific definitions override more general ones.
This determines which spectral data is used for each object during the
simulation.

1. ``Simulation.crossSections``
2. ``phiASE.spectralProperties``
3. ``phiASE.crossSections``
4. ``pump.spectralProperties``
5. ``pump.crossSections``

The resolved spectra are also written back to ``phiASE`` if needed.

Time Integration
----------------

``timeIntegrationSolver`` must provide:

.. code-block:: python

   step(rhs, betaCells, time, timeStep)

Built-in solvers are documented in :doc:`utilities`.

Beta Volume Mapping
-------------------

``Simulation`` updates ``betaVolume`` (:math:`\beta_j`) from ``betaCells``
(:math:`\beta_i`) after each step using ``LegacyGridDataBetaVolumeMapper``.
The mapper interpolates beta values from topology points and z-levels to prism
centers.  It requires ``scipy``.

``updateTerminalLevel``
-----------------------

By default ``updateTerminalLevel=True`` and every z-level is updated.  If it is
set to ``False``, the terminal level is preserved while the other levels are
advanced.
