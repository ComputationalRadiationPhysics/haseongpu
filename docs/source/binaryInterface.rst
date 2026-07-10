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
