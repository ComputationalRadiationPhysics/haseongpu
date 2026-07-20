openPMD Transport
=================

HASEonGPU uses openPMD as the transport boundary between the Python frontend and
the C++ ``calcPhiASE`` backend.  Users normally work with Python objects such
as ``GainMedium`` and ``PhiASE``; the transport converts those objects into the
records and attributes consumed by the backend.

Storage Backends
----------------

The openPMD storage backend is independent from the Alpaka compute backend.
``auto`` is the default: it chooses the first backend supported by both the
compiled and Python openPMD providers in this order: ``adios-sst``, ``adios``,
then ``hdf5``. Explicit runtime values are:

``adios-sst``
   ADIOS2 SST streaming series.  This is the default when supported.

``adios``
   ADIOS2 file-backed ``.bp`` series.

``hdf5``
   HDF5 ``.h5`` series.  Requires HDF5 support in the selected openPMD-api
   provider.

Select it in Python or YAML:

.. code-block:: python

   phi_ase = PhiASE(..., openpmdBackend="auto")

.. code-block:: yaml

   compute:
     backend: Host_Cpu_CpuSerial       # Alpaka compute backend
     openpmd_backend: auto              # choose a compatible backend

Set ``PhiASE.openpmdBackend`` or the YAML ``openpmd_backend`` value to override
automatic selection for a particular run.

Streaming Sessions
------------------

``PhiASE.run(...)`` defaults to one open/write/read/close cycle per call.  For
repeated ``adios-sst`` calls, keep a stream open:

.. code-block:: python

   session = phi_ase.openStream()
   try:
       for _ in range(steps):
           phi_ase.run(gainMedium=medium, crossSections=spectra, openpmdSession=session)
           result = phi_ase.getResults()
   finally:
       phi_ase.closeStream()

Use ``openpmdSession="persistent"`` to let ``PhiASE`` own a reusable stream, or
``openpmdSession="interval"`` to force one-shot behavior.  ``Simulation`` uses a
persistent stream automatically for ``adios-sst`` unless another mode or session
is supplied.

Provider Compatibility
----------------------

The Python ``openpmd_api`` module and the C++ ``openPMD::openPMD`` provider
must be compatible and must both support the selected runtime backend.  The
guided setup checks this for common installs:

.. code-block:: bash

   python3 utils/configure_hase.py

For manual checks against an existing provider:

.. code-block:: bash

   python3 utils/check_openpmd_compatibility.py \
     --backend adios-sst \
     --cmake-prefix-path /path/to/openpmd/prefix

Then point installation or CMake configuration at the same provider, for
example with ``CMAKE_PREFIX_PATH`` or ``openPMD_DIR``.  If the matching Python
package is not on the normal Python path, set ``HASE_OPENPMD_PYTHON_PACKAGE_DIR``
at build time or ``HASE_OPENPMD_PYTHONPATH`` before importing HASEonGPU.

The HASEonGPU wheel does not vendor openPMD runtime libraries or generated
``openpmd_api`` bindings.  The runtime environment must provide compatible
openPMD libraries and Python bindings.

Record Layout
-------------

All array data is written as openPMD ``Mesh`` records below each iteration's
``meshes`` group.  Scalar simulation settings such as ``number_of_points``,
``backend``, ``openpmd_backend``, and ``parallel_mode`` are iteration
attributes, not mesh records.

Topology records use a VTK-compatible unstructured wedge layout inside openPMD:

* ``core_points`` stores point coordinates as ``x``, ``y``, and ``z``
  components.
* ``core_cells_connectivity`` stores VTK wedge point ids.
* ``core_cells_offsets`` stores offsets into the connectivity array.
* ``core_cells_types`` stores the VTK cell type id; wedge cells use type ``13``.

Main input field records are:

* ``core_beta_volume`` and ``core_point_beta`` for dynamic excited-state data
* ``core_cladding_cell_type``, ``core_refractive_index``, and
  ``core_reflectivity`` for static material/surface data
* ``core_lambda_absorption``, ``core_lambda_emission``,
  ``core_sigma_absorption``, and ``core_sigma_emission`` for spectra

The backend writes result records named ``core_result_phi_ase``,
``core_result_mse``, ``core_result_total_rays``, and
``core_result_dndt_ase``.

Custom fields declared with ``GainMedium.defineField(...)`` or
``PrimitiveFieldSpec`` are serialized as additional openPMD mesh records with
their unit metadata.  They are available to downstream readers; the current ASE
backend ignores them unless a future backend explicitly consumes them.

Iteration Updates
-----------------

The first Python-written iteration contains the full static context: topology,
material records, spectra, compute attributes, and dynamic beta fields.  Later
iterations normally contain only ``core_beta_volume`` and ``core_point_beta``
and reuse the cached static context from iteration 0.

Changing topology, spectra, material constants, or compute settings requires a
new input series whose first iteration carries a complete static update.  This
keeps repeated ASE evaluations and streaming runs small while preserving a
stable backend contract.

MPI Launching
-------------

The standalone binary reads the same transport layout under MPI:

.. code-block:: bash

   mpiexec -npernode 4 ./build/calcPhiASE \
       --input-path=input.sst \
       --output-path=output.sst

The high-level Python frontend launches the binary automatically when MPI mode
is selected:

.. code-block:: python

   import HASEonGPU

   phi_ase = HASEonGPU.PhiASE(
       parallelMode="mpi",
       nPerNode=4,
       openpmdBackend="adios-sst",
   )
   phi_ase.run(gainMedium=medium, crossSections=spectra)

The scheduler controls the node allocation, while ``nPerNode`` controls the
number of ranks launched on each allocated node. File-based transport data is
created below ``./IO/phiase_mpi`` so the launch directory must be shared for a
multi-node run.

Artifact Retention
------------------

Temporary transport artifacts are normally removed when a session exits.  These
environment variables help with debugging:

``HASE_OPENPMD_KEEP_ARTIFACTS=1``
   Keep artifacts below ``./hase-openpmd-artifacts``.

``HASE_OPENPMD_ARTIFACT_DIR=/path``
   Write artifacts to an explicit directory.

``HASE_OPENPMD_ARTIFACT_PREFIX=name``
   Prefix generated artifact names.

``HASE_OPENPMD_ARTIFACT_RUN_ID=id``
   Use a stable run id instead of a timestamped id.

``HASE_OPENPMD_WATCHDOG_INTERVAL=30``
   Watchdog interval while the result receiver waits.  Use ``0`` or ``none`` to
   disable the watchdog.

``HASE_OPENPMD_THREAD_JOIN_TIMEOUT=10``
   Time allowed for streaming helper threads to stop during session close.

``HASE_CALCPHIASE=/path/to/calcPhiASE``
   Force the Python transport to use a specific binary.
