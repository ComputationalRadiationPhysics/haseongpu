openPMD Transport
=================

HASEonGPU uses openPMD as the transport boundary between the Python frontend and
the C++ ``calcPhiASE`` backend.  Users normally work with Python objects such
as ``GainMedium`` and ``PhiASE``; the transport converts those objects into the
records and attributes consumed by the backend. This page owns storage-backend,
provider-compatibility, and record-layout details; Alpaka compute selection is
documented separately in :doc:`backendSelection`.

Storage Backends
----------------

The openPMD storage backend is independent from the Alpaka compute backend.
``auto`` is the default: it chooses the first backend supported by both the
compiled and Python openPMD providers in this order: ``adios``, ``adios-sst``,
then ``hdf5``. Explicit runtime values are:

``adios-sst``
   ADIOS2 SST streaming series. Select it explicitly when a live stream is
   preferable to file-backed exchange.

``adios``
   ADIOS2 file-backed ``.bp`` series. This is the automatic default when
   supported because it is currently more robust and usually faster than SST
   for HASEonGPU's frontend/backend exchange.

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
``openpmdSession="interval"`` to force one-shot behavior.  ``Simulation`` owns a separate transport session for each compiled run;
caller-managed simulation sessions are not supported.

``Simulation.runSteps(...)`` and ``Simulation.runUntil(...)`` launch the
compiled ``calcPhiASE --cpp-control`` path. Python writes one initial input
iteration with run-control attributes, then reads the snapshot series produced
by the C++ time loop. For streaming backends, Python starts a dedicated
snapshot receiver thread before sending the input iteration. The autonomous
backend owns stepping; Python only consumes the completed-step indices and
fields selected by ``output_steps`` and ``output_fields`` in the initial
iteration. A bounded handoff keeps memory use finite, so a slow callback can
apply ordinary stream backpressure without turning Python into the step
controller. Caller-managed simulation openPMD sessions are not supported; the
compiled run owns its transport lifetime.

With ``execution_mode="synchronized-debug"``, the input series remains open.
After output step *N*, Python writes dynamic input iteration *N* and the backend
waits for it before starting step *N+1*. Only records listed in
``control_fields`` are written in these later iterations; currently that list
may contain ``beta_volume``. Static topology, spectra, backend selection, and
all other initialization records are transferred only in iteration zero.

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
openPMD libraries and Python bindings. Provider build options are listed in
:ref:`openpmd-provider-options`.

.. _openpmd-record-layout:

openPMD Record Layout
---------------------

The Python frontend has its own object model: ``MeshTopology``,
``GainMedium``, ``CrossSectionData``, and ``PhiASE``. Those names are not an
openPMD schema. The transport boundary is the openPMD series written for the
C++ backend.

All array data at that boundary is written as openPMD ``Mesh`` records below
each ``Iteration``'s ``meshes`` group. Scalar arrays are named openPMD records
with the scalar ``SCALAR`` record component. Component records, currently
``core_points``, use named components such as ``x``, ``y``, and ``z``. The
record names are HASE-owned, which is allowed by openPMD, but the records carry
the normal openPMD mesh and component metadata: ``geometry``,
``geometryParameters``, ``dataOrder``, ``axisLabels``, ``gridSpacing``,
``gridGlobalOffset``, ``gridUnitSI``, ``unitDimension``, component ``unitSI``,
and component ``position``.

Scalar simulation and backend settings are not openPMD field records. Values
such as ``number_of_points``, ``thickness``, ``rng_seed``, ``backend``, and
``parallel_mode`` are stored as attributes on the openPMD iteration.
These values configure the HASE backend and do
not represent sampled mesh data. They therefore are not part of ``/meshes``
and do not carry record metadata such as ``axisLabels`` or component
``position``.

Forward-reflection request attributes follow the same HASE openPMD extension
schema: ``use_reflections``, ``reflection_max_iterations``,
``reflection_tolerance``, and ``surface_reservoir_size``. The parser rejects
the retired ``forward_ray_length`` attribute; forward rays now traverse to a
physical boundary. Runtime environment overrides ``HASE_SRM_MAX_ITERATIONS``
and ``HASE_SRM_DIVERGENCE_STREAK`` are deliberately not serialized because
they are local execution policy rather than portable request data.

The topology convention inside ``/meshes`` follows VTK's unstructured-grid
model and tetrahedral cell. openPMD provides the mesh-record model, but it
does not standardize VTK-style Tet4 connectivity itself. HASEonGPU therefore
stores a VTK-compatible unstructured-cell layout in openPMD records:

* ``core_points`` stores VTK ``POINTS`` as ``x``, ``y``, and ``z`` components.
* ``core_cells_connectivity`` stores the VTK cell connectivity point ids.
* ``core_cells_offsets`` stores offsets into the connectivity array.
* ``core_cells_types`` stores the VTK cell type id; Tet4 cells use type ``10``.

Main input field records are:

* ``core_beta_volume`` for the authoritative dynamic time-integrator state
* ``core_cladding_cell_type``, ``core_refractive_index``, and
  ``core_reflectivity`` for static material/surface data
* ``core_lambda_absorption``, ``core_lambda_emission``,
  ``core_sigma_absorption``, and ``core_sigma_emission`` for spectra

The C++ backend writes result records under ``core_result_``:
``phi_ase``, ``standard_error``, ``relative_standard_error``, ``total_rays``,
and ``dndt_ase``. ``standard_error`` has the same flux unit as ``phi_ase``;
``relative_standard_error`` is dimensionless. Result records use record-C
layout.

Custom fields declared with ``GainMedium.defineField(...)`` or
``PrimitiveFieldSpec`` are serialized as additional openPMD mesh records with
their unit metadata.  They are available to downstream readers; the current ASE
backend ignores them unless a future backend explicitly consumes them.

Result iterations also carry registered HASE extension attributes for SRM
termination: ``srm_status``, ``srm_passes``, ``srm_remaining_fraction``,
``srm_max_iterations``, and ``srm_divergence_streak``. They are scalar
iteration metadata, not mesh records. Python readers expose them as
``Result.srmStatus``, ``srmPasses``, ``srmRemainingFraction``,
``srmMaxIterations``, and ``srmDivergenceStreak``.

.. _compiled-simulation-run-control:

Compiled Simulation Run Control
--------------------------------

For compiled ``Simulation`` runs, iteration attributes also include run-control
metadata:

* ``time_step`` and ``number_of_steps``
* ``pump_steps``
* ``enable_ase`` and ``pre_pump``
* ``execution_mode`` (``autonomous`` or ``synchronized-debug``)
* optional ``output_steps`` containing one-based completed-step indices; when
  omitted, every completed step is emitted
* ``output_fields_string`` and ``control_fields_string`` as scalar JSON arrays
  of field names; this scalar encoding is identical for ADIOS BP, ADIOS SST,
  and HDF5
* ``time_integrator`` (``explicit-euler``, ``heun``, ``midpoint``,
  ``runge-kutta-4``, ``frozen-phi-ase-runge-kutta-4``,
  ``implicit-euler``, or ``exponential-euler``)
* ``implicit_iterations`` and ``implicit_tolerance`` for implicit Euler
* ``pump_schema_version`` (currently ``1``), ``pump_ray_count``, and ``pump_rng_seed``
* flattened source, spectrum, angular, profile, and planar-relay arrays

Readers continue to accept the former ``output_fields`` and ``control_fields``
string-vector attributes for existing series. New writers use only the scalar
JSON attributes so that run control has one backend-independent wire format.

The C++ backend writes one output iteration for each selected completed step.
Snapshot iterations contain the records selected by ``output_fields_string``;
by default these are the cell-centered ``core_beta_volume`` record plus
``core_result_phi_ase``, ``core_result_standard_error``,
``core_result_relative_standard_error``, ``core_result_total_rays``,
``core_result_dndt_ase``, and ``core_result_dndt_pump``. Every dynamic and
result field has the ``cell`` axis; point beta, point PhiASE, and point
derivatives are not part of the compiled simulation contract. The first emitted
snapshot also includes the static canonical mesh/material/spectral records.


Iteration Updates
-----------------

The first Python-written iteration contains the full static context: topology,
material records, spectra, compute attributes, and ``core_beta_volume``.
Synchronized-debug control iterations contain only ``core_beta_volume`` and
reuse the cached static context from iteration 0.

Changing topology, spectra, material constants, or compute settings requires a
new input series whose first iteration carries a complete static update.  This
keeps repeated ASE evaluations and streaming runs small while preserving a
stable backend contract.

MPI uses the same records and attributes; it changes execution topology, not
the transport schema. Rank/device layout, automatic frontend launching, shared
working-directory requirements, and scheduler examples are documented in
:doc:`mpi`.

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
