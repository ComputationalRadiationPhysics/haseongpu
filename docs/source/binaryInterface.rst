Binary Interface
================

``calcPhiASE`` is the standalone C++ executable.  Most users call it through
``PhiASE.run(...)``; use the binary directly when you already have a HASEonGPU
openPMD input series.

Build
-----

Manual compilation is required for direct command-line use.  See
:doc:`CMake Build Options <compilation>`.

.. code-block:: text

   ./build/calcPhiASE

Usage
-----

The binary reads one openPMD input series and writes one openPMD result series:

.. code-block:: bash

   ./build/calcPhiASE \
       --input-path=<openPMD-input-series> \
       --output-path=<openPMD-output-series>

The input series contains mesh topology, material fields, the dynamic
``core_beta_volume`` field, backend settings, parallel mode, sample range, and
optional RNG seed.  This is the same transport layout written by the Python
interface; see :doc:`openpmdTransport`.

For forward reflections, the request uses the schema-defined iteration
attributes ``use_reflections``, ``reflection_max_iterations``,
``reflection_tolerance``, and ``surface_reservoir_size``. ``forward_ray_length``
is retired and rejected. The local environment variables
``HASE_SRM_MAX_ITERATIONS`` and ``HASE_SRM_DIVERGENCE_STREAK`` override the
maximum reflected-pass count and the consecutive-growth divergence threshold;
both must be positive integers. The latter defaults to ``3``.

Examples
--------

Single process:

.. code-block:: bash

   ./build/calcPhiASE --input-path=./input.bp --output-path=./output.bp

MPI launch:

.. code-block:: bash

   mpiexec -npernode 4 ./build/calcPhiASE \
       --input-path=./input.bp \
       --output-path=./output.bp

Arguments
---------

``--input-path``
   Path to the HASEonGPU openPMD input series.

``--output-path``
   Path for the result series.  Results are written as
   ``core_result_phi_ase``, ``core_result_standard_error``,
   ``core_result_relative_standard_error``, ``core_result_total_rays``, and
   ``core_result_dndt_ase`` mesh records. Result iterations also report SRM
   termination through ``srm_status``, ``srm_passes``,
   ``srm_remaining_fraction``, ``srm_max_iterations``, and
   ``srm_divergence_streak`` attributes.

Compiled Simulation Mode
------------------------

``Simulation.runSteps(...)`` uses ``calcPhiASE --cpp-control`` rather than
running a Python pump loop. Python writes the initial geometry, material,
spectra, beta state, and run-control attributes; C++/Alpaka advances the
requested steps and writes a snapshot after each one. The snapshots contain the
updated beta fields, ASE results, and pump and ASE derivatives. The first
snapshot also carries the static context needed to read the series on its own.

.. code-block:: bash

   ./build/calcPhiASE \
       --input-path=./simulation-input.bp \
       --output-path=./simulation-output.bp \
       --cpp-control

The run control selects a time step, step count, pump-step limit, and one of the
compiled integrators (explicit Euler, Heun, midpoint, RK4, implicit Euler, or
exponential Euler). The supported pump routine is
``one-dimensional-z-traversal``. It is configured through ``PumpProperties``;
custom Python pump routines are not part of this execution path.

In this mode Python sends the initial mesh/material/spectra/beta state and the
binary writes one output iteration per completed time step. The output snapshot
iterations include ``core_point_beta``, ``core_beta_volume``,
``core_result_phi_ase``, ``core_result_standard_error``,
``core_result_relative_standard_error``, ``core_result_total_rays``,
``core_result_dndt_ase``, and ``core_result_dndt_pump``. Iteration 0 also
carries the static mesh/material/spectral records so the snapshot series can be
read independently.

Run-control attributes include ``time_step``, ``number_of_steps``,
``pump_steps``, ``time_integrator``, optional implicit-Euler controls
``implicit_iterations`` and ``implicit_tolerance``, and general-pump attributes: ``pump_schema_version``, ``pump_ray_count``,
``pump_rng_seed``, plus flattened source, spectrum, angular-distribution,
spatial-profile, and planar-relay arrays. Schema version 1 replaces the legacy
one-dimensional pump attributes.
