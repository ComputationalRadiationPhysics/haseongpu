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
time integration, clipping, and beta mapping run in C++/Alpaka. Python executes
``on_init`` before launch and ``on_step`` callbacks as snapshots arrive.
Per-step Python mutation is not supported inside a compiled run.
