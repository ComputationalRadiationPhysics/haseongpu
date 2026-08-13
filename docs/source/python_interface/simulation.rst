Simulation
==========

``Simulation`` is the assembly and execution boundary for a time-dependent run.
It references one ``GainMedium``, one ``PhiASE`` configuration, one compiled
time-integrator descriptor, an outer-loop horizon, and the ASE cross sections.
Physical pumps are registered separately before execution and own their own
sampling and activity controls.

Assembly
--------

.. code-block:: python

   simulation = Simulation(
       gain_medium=medium,
       phi_ase=phi_ase,
       time_integrator=FrozenPhiAseRungeKutta4(),
       time_step_size=2e-5,
       simulation_steps=150,
       cross_sections=ase_spectra,
       pre_pump=True,
   )
   simulation.add_pump(
       pump,
       injection_method=SurfacePumpInjector("pump_input"),
       relays=(PlanarPumpRelay.retroreflect("pump_output"),),
   )

``gain_medium`` supplies geometry, material fields, and initial excitation.
``cross_sections`` supplies the complete spectrum used by ASE. The registered
``Pump`` carries its own physical spectrum, matching pump cross sections,
``ray_count``, ``rng_seed``, and ``pump_steps``. ``PhiASE`` carries
``ase_steps``. ``Simulation.pre_pump`` coordinates the first pump and ASE step.
See :doc:`pump_properties` for that separation.

Advance controls
----------------

``step`` advances the constructor's configured ``simulation_steps`` by default
and returns the simulation for chaining. An explicit call count overrides only
the outer-loop horizon:

.. code-block:: python

   simulation.step()
   simulation.step(150)

When ``simulation_steps`` and ``max_time`` are omitted, ``step()`` derives the
horizon as the maximum of every registered ``Pump.pump_steps`` and
``PhiASE.ase_steps``. If all are omitted or zero, supply an explicit simulation
horizon. Activity windows use cumulative steps, so splitting a run across
multiple calls does not reactivate a completed pump or ASE schedule.

``run_until(max_time)`` converts a physical end time into steps using
``time_step_size``. Constructor ``max_time`` supplies its default target.
Configure either ``simulation_steps`` or ``max_time``, not both.

Execution ownership
-------------------

The full time loop, pump evaluation, ASE evaluation, derivative composition,
time integration, and clipping run on cell-centered fields in C++/Alpaka.
Backend selection, device-mesh creation, spectrum upload, forward-ASE queues,
accumulators, reflection reservoirs, and integrator buffers happen once when a
compiled run starts. The evolving ``beta_volume`` and derived ``phi_ase`` stay
on the selected device between steps.

Execution modes
---------------

``execution_mode="autonomous"`` is the normal and fastest contract. Python
executes ``on_init`` once, writes one initialization iteration, and has no
control contact with the C++ loop until it receives a requested snapshot. An
openPMD streaming backend delivers snapshots while the run is active; a file
backend returns them after the run. Snapshot consumers can create backpressure
on a live stream, but they do not drive stepping.

``execution_mode="synchronized-debug"`` is the explicit diagnostic contract.
It requires a streaming backend, emits every completed step, and waits for a
matching Python control iteration before continuing. Register ``before_step``
callbacks to inspect the just-completed state and update fields named by
``control_fields``. This mode intentionally synchronizes Python and the backend
and should not be used for performance measurements.

Snapshot schedule
-----------------

``output_steps`` is either omitted, meaning every completed step, or is a
strictly increasing sequence of one-based completed-step indices. For example,
one snapshot after pumping and one at the end can be selected by application
code:

.. code-block:: python

   simulation = Simulation(
       # ... physical configuration ...
       output_steps=(40, 150),
   )

Final-only output is only a schedule convenience, not a separate execution
path:

.. code-block:: python

   from HASEonGPU import autonomous_final

   simulation = Simulation(
       # ...
       output_steps=autonomous_final(150),
   )

The requested indices must not exceed the number of steps passed to ``step``.
``synchronized-debug`` does not accept ``output_steps`` because every step is a
synchronization boundary.

Selectable fields
-----------------

``output_fields`` fixes the snapshot payload at initialization. Supported
cell-centered fields are:

* ``beta_volume``
* ``phi_ase``
* ``standard_error``
* ``relative_standard_error``
* ``total_rays``
* ``dndt_ase``
* ``dndt_pump``

Unselected arrays are ``None`` in ``TimeStepState``. The synchronized-debug
``control_fields`` list currently supports ``beta_volume``. Point-state names
are deliberately unsupported: the current model and transport contract are
cell/volume based throughout.

.. _simulation-callback-lifecycle:

State and callbacks
-------------------

Callbacks connect Python code to three distinct points in a compiled run:

* ``on_init`` prepares the initial simulation state before step 1.
* ``on_step`` consumes a selected completed-step snapshot.
* ``before_step`` changes selected state between two synchronized-debug steps.

Callbacks receive the additional positional and keyword arguments supplied at
registration. Their return values are ignored. Multiple callbacks for the same
hook run in registration order. An exception stops the current call to
``step`` and is raised to the caller.

Prepare the initial state
^^^^^^^^^^^^^^^^^^^^^^^^^

``on_init`` receives the live ``Simulation`` as its first argument. It runs
once, immediately before the first compiled run, regardless of execution mode,
openPMD backend, or snapshot schedule. Use it to initialize or normalize values
that step 1 must consume:

.. code-block:: python

   def initialize(simulation, initial_beta):
       simulation.gain_medium.get("betaVolume").value[...] = initial_beta

   simulation.on_init(initialize, 0.15)

Register ``on_init`` before the first call to ``step`` or ``run_until``. It does
not run again for later calls on the same ``Simulation`` object.

Consume completed snapshots
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``on_step`` receives a ``TimeStepState`` as its first argument:

.. code-block:: python

   def report(state, prefix):
       print(prefix, state.step, state.time, state.beta_volume.mean())

   simulation.on_step(report, "completed")
   simulation.step(10)
   final_state = simulation.get_last_state()

The state exposes the cumulative completed ``step``, physical ``time`` in
seconds, static ``topology``, and the cell arrays selected by ``output_fields``.
Its arrays are copies. Changing them does not update either ``GainMedium`` or
the compiled backend. ``Simulation`` stores only the latest emitted snapshot;
use ``on_step`` to retain a history or write output.

In autonomous mode, ``output_steps`` selects the snapshots that invoke
``on_step``. The indices are relative to each call to ``step``. For example,
``output_steps=(2,)`` followed by two calls to ``step(2)`` invokes ``on_step``
with cumulative state indices 2 and 4. Every requested output index must fit
within each compiled call.

The openPMD backend determines when selected snapshots reach Python:

* ``adios-sst`` delivers them while the compiled run is active. A slow
  ``on_step`` callback can delay further snapshot delivery.
* ``adios`` and ``hdf5`` store them in files. Python invokes ``on_step`` while
  reading those files after the compiled run has finished.

Neither behavior makes ``on_step`` a backend control hook.

Change state between steps
^^^^^^^^^^^^^^^^^^^^^^^^^^

``before_step`` receives the live ``Simulation`` and is available only with
``execution_mode="synchronized-debug"`` and a streaming openPMD backend. The
supported streaming backend is ``adios-sst``.

The name refers to the next step. It first runs between steps 1 and 2, not
before step 1. Use ``on_init`` to modify the state consumed by step 1. For each
nonfinal completed step, the order is:

#. The backend completes step *N* and emits its snapshot.
#. Python updates the live ``GainMedium`` when ``beta_volume`` is present in
   the snapshot, then calls ``on_step``.
#. Python calls ``before_step``.
#. Python sends the fields named by ``control_fields`` to the backend.
#. The backend begins step *N* + 1.

There is no ``before_step`` call after the final step because no next step
exists. Consequently, a compiled ``step(1)`` call never invokes this hook.

Currently, ``beta_volume`` is the only supported control field:

.. code-block:: python

   def prepare_next_step(simulation, beta):
       simulation.gain_medium.get("betaVolume").value[...] = beta

   simulation = Simulation(
       # ... physical configuration ...
       execution_mode="synchronized-debug",
       output_fields=("beta_volume", "phi_ase"),
       control_fields=("beta_volume",),
   )
   simulation.on_init(prepare_next_step, 0.15)
   simulation.on_step(report, "completed")
   simulation.before_step(prepare_next_step, 0.20)
   simulation.step(3)

This example uses ``on_init`` to set beta for step 1. After steps 1 and 2,
``before_step`` sets beta for steps 2 and 3. Include ``beta_volume`` in
``output_fields`` when the inter-step callback must inspect the beta produced
by the completed step.

Synchronized-debug emits every completed step and therefore rejects
``output_steps``. It waits for Python after every nonfinal snapshot, so use it
for inspection and control rather than performance measurements. See
:doc:`../openpmdTransport` for the transport lifecycle and backend selection.

Runtime behavior
----------------

``PhiASE.ase_steps=None`` and ``ase_steps=0`` skip ASE transport and depletion
while retaining active pumps and fluorescence terms. ``Pump.pump_steps=None``
and ``pump_steps=0`` disable that individual source. ``Simulation.pre_pump=True``
omits ASE from the first outer step. Use it when that step should seed excitation
with the active pumps before the first ASE evaluation.
``report_timings=True`` reports frontend transport/decoding, snapshot creation,
and callback time; build-time backend instrumentation is documented in
:doc:`../compilation`.

The compiled population equation and integrator behavior are described in
:ref:`pump-and-time-stepping`.
