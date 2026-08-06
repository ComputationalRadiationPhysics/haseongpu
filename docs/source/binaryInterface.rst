Binary Interface
================

``calcPhiASE`` is the standalone C++ executable behind the Python frontend.
Most users call it through ``PhiASE`` or ``Simulation``; use it directly when a
workflow already produces the HASE openPMD input series.

Build
-----

Manual CMake compilation produces ``./build/calcPhiASE``. See
:doc:`compilation` for build and provider options. A thin Python installation
also records the matching executable and normally launches it automatically.

One ASE evaluation
------------------

The default mode reads one input series and writes one result series:

.. code-block:: bash

   ./build/calcPhiASE \
       --input-path=./input.bp \
       --output-path=./output.bp

The input iteration contains Tet4 topology, gain-medium and spectral fields,
dynamic cell excitation, PhiASE controls, and compute settings. The output
contains cell-centered flux, uncertainty, history counts, depletion rate, and
reflection termination metadata. The authoritative record and attribute names
are specified in :ref:`openpmd-record-layout`.

Compiled simulation mode
------------------------

``--cpp-control`` interprets the input iteration as a complete time-stepped run
request:

.. code-block:: bash

   ./build/calcPhiASE \
       --input-path=./simulation-input.bp \
       --output-path=./simulation-output.bp \
       --cpp-control

The request adds time-step, step-count, pump-duration, integrator, ASE toggle,
pre-pump, execution mode, output schedule and field selection, control fields,
and serialized pump-source controls. C++/Alpaka owns the loop and writes one
snapshot iteration for each selected completed step. Without an output
schedule, every completed step is selected. The first emitted snapshot repeats
the static context so the output series can be read independently; later
snapshots contain the selected evolving cell fields and results.

Python ``Simulation.step`` uses this mode automatically. Physical object
composition is documented in :doc:`pythonInterface`; serialized run controls
and snapshot records are specified in :ref:`compiled-simulation-run-control`.

Arguments
---------

``--input-path=<series>``
   Required HASE openPMD input series.

``--output-path=<series>``
   Required destination for result iterations.

``--cpp-control``
   Optional compiled simulation mode. Without it, the executable performs one
   PhiASE evaluation.

No other command-line options are accepted. Physics, sampling, compute, and
transport settings belong to the openPMD request rather than to a second binary
CLI configuration surface.

MPI
---

The same executable runs under MPI and consumes the same transport layout.
Build requirements, frontend launching, rank/device distribution, and scheduler
examples are centralized in :doc:`mpi`.
