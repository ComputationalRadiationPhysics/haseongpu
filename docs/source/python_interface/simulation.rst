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
       enableASE=True,
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
``PumpProperties.pumpSubsteps``, which is a compatibility field for pump
routines with inner time integration and does not count outer simulation
steps.

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

Each call to ``step()`` or ``runSteps(...)`` serializes the initial setup to
openPMD, launches the compiled ``calcPhiASE --run-simulation`` path, and then
receives one streamed snapshot per completed step.  Pump traversal, ASE
evaluation, derivative composition, time integration, clipping, and prism beta
mapping all run in C++/Alpaka.  Python only runs ``onInit`` before launch and
``onStep`` callbacks as snapshots arrive.

The supported compiled pump routine is ``one-dimensional-z-traversal``.
Set ``enableASE=False`` to run the same pump, fluorescence, and integration
path while omitting ASE depletion; ASE result fields remain present and are
filled with zeros.

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
   Not supported by compiled simulation runs. Registering this hook causes
   ``runSteps``/``step`` to raise, because Python cannot mutate state between
   C++-owned steps.

``onStep(callback, *args, **kwargs)``
   Runs after every completed step. The callback signature is
   ``callback(state, *args, **kwargs)`` and receives the completed
   ``TimeStepState`` first. Use it for logging, inspection, exporting, or
   storing result snapshots.

Examples:

.. code-block:: python

   def initialize(simulation, beta0):
       simulation.gainMedium.get("betaCells").value = beta0

   def write_state(state, output_dir):
       print(state.step, state.time, state.betaCells.mean())

   simulation.onInit(initialize, beta0)
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

``timeIntegrationSolver`` is a compiled integrator name or a descriptor object
with a ``.name`` attribute. Built-in descriptors are documented in
:doc:`utilities`. Supported names are:

* ``explicit-euler``
* ``heun``
* ``midpoint``
* ``runge-kutta-4``
* ``frozen-phi-ase-runge-kutta-4`` (reuses one ASE evaluation across RK4 stages)
* ``frozen-phi-ase-runge-kutta-4``
* ``implicit-euler``
* ``exponential-euler``

``ImplicitEuler(iterations=8, tolerance=1e-10)`` also serializes
``implicit_iterations`` and ``implicit_tolerance`` run-control attributes.
Set ``simulation.enableAse = False`` to advance pump and fluorescence without ASE. Custom Python time integrators are not supported by compiled runs.

Beta Volume Mapping
-------------------

``Simulation`` initializes missing ``betaVolume`` values from ``betaCells`` with
``ConnectivityAverageBetaVolumeMapper``.  Each prism value is the mean of the
three lower-level and three upper-level vertex beta values, matching the
C++/Alpaka mapping kernel. ``LegacyGridDataBetaVolumeMapper`` remains as a
compatibility alias for this mapper.
