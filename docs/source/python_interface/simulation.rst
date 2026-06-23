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

Limit pumping to the first outer simulation steps while ASE and fluorescence
continue for the full run:

.. code-block:: python

   simulation.runSteps(150, pumpSteps=50)

``pumpSteps`` counts outer calls to ``Simulation.step()``.  You can pass it
to ``runSteps`` or store it on ``PumpProperties`` and then call
``simulation.runSteps(150)``.  When neither location provides ``pumpSteps``, the
pump is active for every step passed to ``runSteps``.  This is different from
``PumpProperties.pumpSubsteps``, which only controls the internal time
resolution of the pump integration inside one pumped simulation step.

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
9. Latest-state update and ``onStep`` callbacks.

Callbacks
---------

Lifecycle hooks always insert the simulation-provided object as the first
argument. Extra arguments passed during registration are appended after that
first argument. Callback return values are ignored, and each registration method
returns ``self`` so calls can be chained.

``onInit(callback, *args, **kwargs)``
   Runs once before the first step. The callback signature is
   ``callback(simulation, *args, **kwargs)`` and receives the live
   ``Simulation`` object first. Use it to initialize or normalize mutable
   simulation inputs.

``beforeStep(callback, *args, **kwargs)``
   Runs before every step, after ``onInit`` has run. The callback signature is
   ``callback(simulation, *args, **kwargs)`` and receives the live
   ``Simulation`` object first. Use it for controlled pre-step changes such as
   time-dependent pump settings.

``onStep(callback, *args, **kwargs)``
   Runs after every completed step. The callback signature is
   ``callback(state, *args, **kwargs)`` and receives the completed
   ``TimeStepState`` first. Use it for logging, inspection, exporting, or
   storing result snapshots.

Examples:

.. code-block:: python

   def initialize(simulation, beta0):
       simulation.gainMedium.get("betaCells").value = beta0

   def adjust_pump(simulation, scale):
       simulation.pump.withProperty("scale", scale)

   def write_state(state, output_dir):
       print(state.step, state.time, state.betaCells.mean())

   simulation.onInit(initialize, beta0)
   simulation.beforeStep(adjust_pump, 0.5)
   simulation.onStep(write_state, output_dir)

This means ``simulation.onStep(write_state, output_dir)`` calls
``write_state(state, output_dir)`` for every completed step.

Results
-------

``getLastState()`` returns the most recent ``TimeStepState`` snapshot:

.. code-block:: python

   last = simulation.getLastState()
   last.betaCells.shape
   last.betaVolume.shape
   last.phiAse.shape

``Simulation`` keeps only the latest state in memory. Use ``onStep`` to write,
inspect, or store state for every completed step. ``getResults()`` is retained
as an alias for ``getLastState()``.

``TimeStepState`` fields are:

* ``step``: completed step index.
* ``time``: simulation time after the step.
* ``betaCells``: point and level beta values :math:`\beta_i`.
* ``betaVolume``: prism beta values :math:`\beta_j` used by ASE.
* ``phiAse``: ASE flux values :math:`\Phi_i`, or ``None`` if unavailable.
* ``dndtAse``: ASE derivative contribution to :math:`d\beta/dt`.
* ``dndtPump``: pump derivative contribution to :math:`d\beta/dt`.
* ``aseResult``: raw lower-level ASE result object.
* ``topology``: static mesh topology used by geometry-aware callbacks such as VTK export.

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
