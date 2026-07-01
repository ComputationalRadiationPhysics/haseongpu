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

The input series contains mesh topology, material fields, dynamic beta fields,
backend settings, parallel mode, sample range, and optional RNG seed.  This is
the same transport layout written by the Python interface; see
:doc:`openpmdTransport`.

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
   ``core_result_phi_ase``, ``core_result_mse``, ``core_result_total_rays``,
   and ``core_result_dndt_ase`` records.

Compiled Simulation Mode
------------------------

``Simulation.runSteps(...)`` uses ``calcPhiASE --run-simulation`` rather than
running a Python pump loop. Python writes the initial geometry, material,
spectra, beta state, and run-control attributes; C++/Alpaka advances the
requested steps and writes a snapshot after each one. The snapshots contain the
updated beta fields, ASE results, and pump and ASE derivatives. The first
snapshot also carries the static context needed to read the series on its own.

The run control selects a time step, step count, pump-step limit, and one of the
compiled integrators (explicit Euler, Heun, midpoint, RK4, implicit Euler, or
exponential Euler). The supported pump routine is
``one-dimensional-z-traversal``. It is configured through ``PumpProperties``;
custom Python pump routines are not part of this execution path.
