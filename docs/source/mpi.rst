MPI Execution
=============

HASEonGPU can distribute one ASE calculation over multiple MPI ranks.  MPI
splits the sample index range across ranks; each rank then uses the devices
assigned on its node.

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

* MPI ranks divide the global sample range.
* Ranks on the same node divide the local device list.
* GPU IDs printed by HASEonGPU are local to each node.

For example, with two nodes and four visible GPUs per node, ranks on both nodes
may report GPUs ``0-3``; those are different physical devices on different
nodes.

Runtime Settings
----------------

These values are set through ``PhiASE`` or YAML. The scheduler supplies the
node allocation; for ``parallelMode="mpi"``, the Python frontend starts
``calcPhiASE`` with ``mpiexec -npernode <nPerNode>`` inside that allocation.

``parallelMode`` / ``parallel_mode``
   ``single`` runs without MPI. ``mpi`` launches the executable under MPI and
   splits samples across ranks.

``numDevices``
   Maximum number of local devices one node should use.  In MPI mode,
   HASEonGPU divides those devices across ranks on the same node.

``nPerNode`` / ``n_per_node``
   Number of MPI ranks per allocated node passed to ``mpiexec -npernode``.

The frontend places temporary file-based openPMD transport data below
``./IO/phiase_mpi``. For multi-node runs, launch HASE from a working directory
that is visible on every allocated node.

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
