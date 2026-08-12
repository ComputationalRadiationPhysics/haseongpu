Simulation
==========

``Simulation`` is the assembly and execution boundary for a time-dependent run.
It references one ``GainMedium``, one ``PhiASE`` configuration, one compiled
time-integrator descriptor, shared pump-solver controls, and the ASE cross
sections. Physical pumps are registered separately before execution.

Assembly
--------

.. code-block:: python

   simulation = Simulation(
       gain_medium=medium,
       phi_ase=phi_ase,
       time_integrator=FrozenPhiAseRungeKutta4(),
       time_step_size=2e-5,
       pump_solver=MonteCarloPumpSolver(
           ray_count=50_000,
           seed=5489,
           max_steps=100,
       ),
       cross_sections=ase_spectra,
       enable_ase=True,
       pre_pump=True,
   )
   simulation.add_pump(
       pump,
       injection_method=SurfacePumpInjector("pump_input"),
       relays=(PlanarPumpRelay.retroreflect("pump_output"),),
   )

``gain_medium`` supplies geometry, material fields, and initial excitation.
``cross_sections`` supplies the complete spectrum used by ASE. The registered
``Pump`` carries its own physical spectrum and matching pump cross sections.
See :doc:`pump_properties` for that separation.

Advance controls
----------------

``step`` advances an explicit number of outer time steps and returns the
simulation for chaining:

.. code-block:: python

   simulation.step(150)
   simulation.step(150, pump_steps=50)

The optional call-level ``pump_steps`` overrides
``MonteCarloPumpSolver.max_steps`` for that run. It selects how many of the
requested outer steps include pump transport; it does not change pump ray count
or Tet4 traversal length.

``run_until(max_time)`` converts a physical end time into steps using
``time_step_size``. Constructor ``max_time`` supplies its default target.
``max_steps`` is retained as an intended-run property, but ``step`` still takes
the count to execute explicitly.

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

State and callbacks
-------------------

Register output consumers before the first step:

.. code-block:: python

   def report(state, prefix):
       print(prefix, state.step, state.time, state.beta_volume.mean())

   simulation.on_step(report, "completed")
   simulation.step(10)
   final_state = simulation.get_last_state()

An ``on_step`` callback receives a copied ``TimeStepState`` followed by the
arguments supplied during registration. It exposes the completed one-based
``step``, physical ``time``, static ``topology``, and cell arrays
``beta_volume``, ``phi_ase``, ``dndt_pump``, and ``dndt_ase``. ASE uncertainty,
history counts, and reflection status remain available through ``ase_result``.
Mutating a callback snapshot does not change backend state.

``Simulation`` stores only the latest requested snapshot. Use callbacks to
retain a history or write output. ``on_init`` receives the live simulation once
before compiled execution and can finalize its inputs. Autonomous execution
does not run Python between C++-owned time steps. A ``before_step`` callback
therefore requires ``execution_mode="synchronized-debug"`` and runs after every
nonfinal snapshot, before the backend is allowed to begin the next step.
Synchronized-debug emits every step and does not accept ``output_steps``.

Runtime behavior
----------------

``enable_ase=False`` skips ASE transport and depletion while retaining pump and
fluorescence terms. ``pre_pump=True`` omits ASE from the first pumped step so
the pump can seed excitation before the first ASE evaluation.
``report_timings=True`` reports frontend transport/decoding, snapshot creation,
and callback time; build-time backend instrumentation is documented in
:doc:`../compilation`.

The compiled population equation and integrator behavior are described in
:ref:`pump-and-time-stepping`.
