Simulation
==========

``Simulation`` is the central assembly object. Its public constructor and
registration methods follow PICMI-style snake_case naming:

.. code-block:: python

   simulation = Simulation(
       gain_medium=medium,
       phi_ase=phi_ase,
       time_integrator=RungeKutta4(),
       time_step_size=1e-5,
       pump_solver=MonteCarloPumpSolver(ray_count=50000, seed=5489),
       cross_sections=spectra,
       max_steps=150,
       enable_ase=True,
       pre_pump=True,
   )
   simulation.add_pump(pump, injection_method=injector, relays=relays)
   simulation.on_step(write_state, output_directory)
   simulation.step()

As in PICMI, ``simulation.step(nsteps=1)`` advances the requested number of
steps and defaults to one. ``max_steps`` and ``max_time`` describe the intended
run limits; pass the desired count to ``step`` or use ``run_until`` for a time
limit. ``pump_steps`` can override the pump solver's ``max_steps`` for a
particular call.

``simulation.get_last_state()`` returns the latest ``TimeStepState``. The
simulation does not retain the full history, so register ``on_step`` to store or
write every snapshot.

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
matching Python control iteration before continuing. Register ``beforeStep``
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

   def pump_boundary_and_final(number_of_steps, pump_steps):
       return tuple(sorted({min(number_of_steps, pump_steps), number_of_steps}))

   simulation = Simulation(
       # ... physical configuration ...
       output_steps=pump_boundary_and_final(150, 40),
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
