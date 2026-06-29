MPI Execution
=============

HASEonGPU can distribute one ASE calculation over multiple MPI ranks.  MPI is
used to split the sample index range across ranks; each active rank can then
use one or more devices from the node on which that rank is running.

Dependencies and Build Configuration
------------------------------------

MPI execution requires an MPI implementation (``OpenMPI >= 4.0``) and an HASEonGPU build with MPI
support enabled.  The CMake option ``DISABLE_MPI`` controls this at build time:

``DISABLE_MPI=AUTO``
   Tries to find MPI.  If MPI cannot be found, the build continues without MPI
   support.

``DISABLE_MPI=OFF``
   Requires MPI.  Configuration fails if MPI is missing.

``DISABLE_MPI=ON``
   Builds without MPI support.

An MPI-enabled build links against the MPI C++ target found by CMake.  A
runtime started from openPMD metadata with ``parallelMode = mpi`` exits with an
error if the binary was built without MPI support.

Execution Model
---------------

MPI execution has two levels of distribution:

* MPI distributes the global sample index range across active ranks.
* Each rank divides its local sample range across the GPUs assigned to that
  rank.

For example, an input with 4210 samples and 2 active ranks creates two
rank-level ranges of roughly 2105 samples each.  If each rank owns four GPUs,
each rank then splits its 2105 samples across four worker threads and four
devices.

The GPU IDs shown in HASEonGPU output are local to each node.  A report such as
``rank 0 node ga008 ... assigned GPUs 0-3`` and ``rank 1 node ga009 ...
assigned GPUs 0-3`` means GPU IDs 0-3 on two different nodes, not the same four
physical GPUs.

Important Settings
------------------
These settings can be written through ``PhiASE`` configuration and are serialized into the openPMD compute metadata. The MPI process topology is controlled by how ``calcPhiASE`` is launched, for example with ``mpiexec`` or a scheduler.

``parallelMode``
   Selects the execution path.

   ``"single"`` runs without MPI communication.  The process uses up to
   ``numDevices`` devices on the local node.

   ``"mpi"`` uses MPI communication.  The executable must be launched by an MPI
   launcher such as ``mpiexec`` or by a scheduler that starts MPI tasks.  The
   sample range is split across active ranks.

``numDevices``
   Sets the maximum number of devices visible to one HASEonGPU process on its
   node.  In MPI mode, HASEonGPU first limits the local device list to
   ``numDevices`` and then distributes that local list across the ranks on the
   same node.

   With four visible GPUs per node:

   * 4 ranks per node and ``numDevices=4`` gives each rank one GPU.
   * 1 rank per node and ``numDevices=4`` gives that rank four GPUs.
   * 2 ranks per node and ``numDevices=4`` gives each rank two GPUs.

``nPerNode``
   Controls the number of MPI ranks launched per node in Python helper paths
   that explicitly call ``mpiexec`` for the ``calcPhiASE`` binary.

   ``nPerNode`` is a launcher setting, not a GPU count inside the C++ compute
   kernel.  After MPI has started the processes, the C++ executable detects the
   actual ranks per node from MPI and divides ``numDevices`` across those local
   ranks.

Common Usage Modes
------------------

In general, any variation of MPI ranks per node is supported. The two most
common usage modes are described below.

One MPI rank per device
^^^^^^^^^^^^^^^^^^^^^^^

This mode starts one MPI rank for each device. It is usually the faster option,
because each rank is responsible for driving a single device. For example, on
two nodes with four devices per node, this results in eight ranks in total and
one device per rank.

.. code-block:: text

   parallelMode = mpi
   numDevices = $devicesPerNode
   nPerNode = $devicesPerNode

Slurm allocation example:

.. code-block:: bash

   srun -N $numNodes \
        --tasks-per-node=$devicesPerNode \
        --gres=gpu:$devicesPerNode \
        --pty bash


One MPI rank per node with multiple devices per rank
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This mode starts one MPI rank per node and lets each rank use multiple devices.
Depending on the workload and host-side scheduling, this can be slower due to
thread contention inside a single rank. For example, on two nodes with four
devices per node, this results in two ranks in total and four devices per rank.

.. code-block:: text

   parallelMode = mpi
   numDevices = $devicesPerNode
   nPerNode = 1

Slurm allocation example:

.. code-block:: bash

   srun -N $numNodes \
        --tasks-per-node=1 \
        --cpus-per-task=$cpusPerTask \
        --gres=gpu:$devicesPerNode \
        --pty bash

For this mode, ``--cpus-per-task`` should be chosen significantly larger than
``$devicesPerNode``. Each MPI rank starts multiple host threads to drive multiple
devices. If the rank is bound to too few CPU cores, the run can become host-side
limited even though all devices are allocated.


Interpreting Output
-------------------

MPI-enabled HASEonGPU prints topology information in the ``[INFO]`` output:

.. code-block:: text

   [INFO] Active nodes             : 2
   [INFO] Active ranks             : 2
   [INFO] Active ranks per node    : 1 avg (min=1, max=1)
   [INFO] Active GPUs              : 8
   [INFO] GPUs per active rank     : 4 avg
   [INFO] GPUs per active node     : 4 avg (min=4, max=4)

.. code-block:: text

Runtime Considerations
----------------------

The two common layouts use the same total GPU count but do not always have the
same runtime:

* Many ranks with one GPU each create one host process per GPU.
* One rank with multiple GPUs creates multiple host threads inside one process.

If a scheduler binds one MPI rank to one CPU core, the multi-GPU-per-rank mode
can become CPU limited because several GPU-driving threads share the same core.
Slurm ``--cpus-per-task`` or Open MPI ``--map-by ...:PE=<n>`` should match the
number of GPUs driven by each rank.  Open MPI ``--report-bindings`` is useful
for checking the effective CPU binding.
