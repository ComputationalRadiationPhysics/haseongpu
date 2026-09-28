MPI Execution
=============

HASEonGPU can distribute one ASE calculation over multiple MPI ranks. The
domain scheduler assigns indivisible ``(domainId, batchId)`` work items to
rank-owned devices. It balances estimated ray/mesh work while considering
resident bytes and node locality. MPI changes execution topology, not the
physical model or openPMD record layout.

Build Requirement
-----------------

Build with MPI support before selecting ``parallelMode="mpi"`` or YAML
``parallel_mode: mpi``.  The CMake option is ``DISABLE_MPI``:

``AUTO``
   Use MPI if CMake can find it; otherwise continue without MPI.

``OFF``
   Require MPI and fail configuration if it is missing.

``ON``
   Build without MPI support.

Execution Model
---------------

* MPI ranks participate as workers with stable rank, node, and device identity.
* Each rank owns one device selected from the node-local visible device list.
* Domain/batch work is never split after assignment.
* Raw batch accumulators are reduced before one global normalization.
* GPU IDs printed by HASEonGPU are local to each node.

Direct boundary candidates stay on the device when the worker group has one owner.
With multiple rank-owned devices, candidates cross ranks through host-staged MPI
transport because device peer access is not assumed. Candidate payloads gather
only at the combing rank; selected histories are sent only to their scheduled
destination owners for the next pass.

For example, with two nodes and four visible GPUs per node, ranks on both nodes
may report GPUs ``0-3``; those are different physical devices on different
nodes.

Runtime Settings
----------------

These values are set through ``PhiASE`` or YAML. The scheduler supplies the
node allocation; for ``parallelMode="mpi"``, the Python frontend starts
``calcPhiASE`` with ``mpiexec -npernode <nPerNode>`` inside that allocation.
Alpaka compute and openPMD storage backends remain independent selections; see
:doc:`backendSelection` and :doc:`openpmdTransport`.

``parallelMode`` / ``parallel_mode``
   ``single`` runs without MPI. ``mpi`` launches the executable under MPI and
   splits samples across ranks.

``numDevices``
   Maximum number of local devices made visible to the run. In MPI mode each
   rank selects one of those devices from its node-local rank index.

``nPerNode`` / ``n_per_node``
   Number of MPI ranks per allocated node passed to ``mpiexec -npernode``.

The frontend places temporary file-based openPMD transport data below
``./IO/phiase_mpi``. For multi-node runs, launch HASE from a working directory
that is visible on every allocated node.

Direct executable invocation uses the same ``calcPhiASE`` arguments documented
in :doc:`binaryInterface`; the Python frontend normally constructs those paths
and launches the binary automatically.

Common Layouts
--------------

One rank per GPU is usually the most straightforward layout:

.. code-block:: text

   parallelMode = mpi
   numDevices = $devicesPerNode
   nPerNode = $devicesPerNode

One rank per node lets one process drive multiple GPUs, but requires enough CPU
cores for the GPU-driving host threads:

.. code-block:: text

   parallelMode = mpi
   numDevices = $devicesPerNode
   nPerNode = 1

Slurm examples:

.. code-block:: bash

   # one rank per GPU
   srun -N $numNodes --tasks-per-node=$devicesPerNode --gres=gpu:$devicesPerNode --pty bash

   # one rank per node
   srun -N $numNodes --tasks-per-node=1 --cpus-per-task=$cpusPerTask \
        --gres=gpu:$devicesPerNode --pty bash

If a scheduler binds a multi-GPU rank to too few CPU cores, the run can become
host-side limited.  Use Slurm ``--cpus-per-task`` or Open MPI
``--map-by ...:PE=<n>`` to match the number of devices driven by each rank.
``--report-bindings`` is useful for checking Open MPI CPU binding.

Limitations
-----------

A supported MPI run requires every worker rank to expose the selected Alpaka
backend and at least one compatible local device. Every rank must remain
responsive, complete its assigned batches, participate in every collective,
and exit normally. The Python launcher can validate the openPMD provider, but
its process-local Alpaka query cannot establish device visibility on remote
ranks. Verify backend and device visibility throughout the allocation before
starting a production run.

HASEonGPU does not currently provide a collective worker-health handshake,
heartbeats, collective timeouts, failed-rank recovery, or work redistribution.
If a rank exits prematurely, failure propagation depends on the MPI
implementation and launcher. If a rank becomes unresponsive, the remaining
ranks may wait indefinitely in a collective while the launcher process remains
alive. The frontend streaming watchdog checks the launched process; it does not
prove that every worker rank is making progress.

MPI reductions require exactly one producer for each domain/batch work item and
combine raw contributions before normalization. This detects inconsistent
global ray accounting when every rank reaches the collective, but it does not detect
plausible corrupted values returned by a participating worker. HASEonGPU has no
redundant computation or result-integrity mechanism for that case.

Treat an MPI result as complete only when ``calcPhiASE``/``mpiexec`` exits
successfully and all requested result snapshots have been received. Do not use
partial output from a failed or interrupted MPI run as a completed simulation
result.

Output
------

MPI-enabled runs print the active node/rank/device topology, for example:

.. code-block:: text

   [INFO] Active nodes             : 2
   [INFO] Active ranks             : 2
   [INFO] Active ranks per node    : 1 avg (min=1, max=1)
   [INFO] Active GPUs              : 8
   [INFO] GPUs per active rank     : 4 avg
   [INFO] GPUs per active node     : 4 avg (min=4, max=4)
