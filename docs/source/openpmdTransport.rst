openPMD Transport
=================

HASEonGPU uses an openPMD series as the transport boundary between the Python
frontend and the C++ ``calcPhiASE`` backend. The high-level Python objects
remain the user-facing model; the transport serializes the data needed by the
backend into openPMD records and attributes.

The transport is used by ``PhiASE.run(...)`` and by the standalone
``calcPhiASE`` binary.  Advanced users can also import
``pyInclude.openpmd.transport`` to write input series, run the backend, and read
result series explicitly.

Storage Backends
----------------

The openPMD storage backend controls the on-disk or streaming format. It is
separate from the alpaka compute backend selected by the ``PhiASE.backend``
setting.

Accepted transport backend names are:

``adios``
   ``.bp`` series with an explicit ADIOS2 openPMD configuration.

``hdf5``
   ``.h5`` series. Requires an openPMD-api build with HDF5 support.

``adios-sst``
   ``.sst`` ADIOS2 SST streaming series. Use this backend only when a producer
   and consumer are intended to run concurrently.

Select the backend on ``PhiASE`` in Python or YAML:

.. code-block:: python

   phi_ase = PhiASE(..., openpmdBackend="adios")

.. code-block:: yaml

   compute:
     backend: Host_Cpu_CpuSerial       # Alpaka compute backend
     openpmd_backend: adios-sst       # openPMD storage/streaming backend

For a process-wide default, use ``HASE_OPENPMD_BACKEND``:

.. code-block:: bash

   HASE_OPENPMD_BACKEND=adios python3 my_simulation.py

For direct helper calls, pass ``transport=``:

.. code-block:: python

   import pyInclude.openpmd.transport as transport

   result = transport.runPhiASE(
       phi_ase,
       medium,
       spectra,
       transport="adios",
   )

For repeated calls over a streaming backend, keep one session open and pass it
back into each run:

.. code-block:: python

   phi_ase = PhiASE(..., openpmdBackend="adios-sst")

   openpmdSession = phi_ase.openStream()
   try:
       for _ in range(steps):
           phi_ase.run(
               gainMedium=medium,
               crossSections=spectra,
               openpmdSession=openpmdSession,
           )
           result = phi_ase.getResults()
   finally:
       phi_ase.closeStream()

``PhiASE.run(...)`` defaults to interval ownership: one open/write/read/close
cycle per call. Pass ``openpmdSession="persistent"`` to let the ``PhiASE``
object open and reuse its own stream, or pass an existing session object when
the caller owns the stream lifetime.

``Simulation.runSteps(...)`` and ``Simulation.runUntil(...)`` keep one
persistent stream automatically for ``adios-sst``. Pass
``openpmdSession="interval"`` to force one-shot open/write/read/close behavior,
or pass an existing session object to share caller-managed ownership.

The CMake build also has ``HASE_OPENPMD_BACKEND``. That option selects which
openPMD dependencies and default file extension are built into the C++ test and
binary configuration. Runtime Python selection still requires the matching
openPMD-api Python module and backend support to be available.

The Python transport requires the frontend ``openpmd_api`` module and the
compiled ``calcPhiASE`` reader to use the same openPMD-api provider. With the
default bundled build, HASEonGPU builds and installs that Python module next to
the package automatically. With ``HASE_USE_SYSTEM_OPENPMD=ON``, the system
openPMD C++ package must be paired with the Python ``openpmd_api`` package from
the same system installation. CMake records the matching package directory in
``HASE_OPENPMD_PYTHON_PACKAGE_DIR`` and the runtime imports from that directory
before accepting a pre-imported module.

For a pip editable install, either make the matching module importable in the
same Python environment used by the build, or pass the package directory
explicitly:

.. code-block:: bash

   CMAKE_ARGS="-DHASE_USE_SYSTEM_OPENPMD=ON \
     -DHASE_OPENPMD_PYTHON_PACKAGE_DIR=/path/to/site-packages" \
     python -m pip install -e .

For a source-tree CMake build with Python enabled, use the same cache variable:

.. code-block:: bash

   cmake -S . -B build \
     -DHASE_ENABLE_PYTHON=ON \
     -DHASE_USE_SYSTEM_OPENPMD=ON \
     -DHASE_OPENPMD_PYTHON_PACKAGE_DIR=/path/to/site-packages

``HASE_OPENPMD_PYTHON_PACKAGE_DIR`` must point to the directory containing the
``openpmd_api`` package, not to an unrelated PyPI wheel. As an explicit runtime
escape hatch, ``HASE_OPENPMD_PYTHONPATH`` may be set to the matching package
directory before importing HASEonGPU. Do not use it to mix unrelated frontend
and backend openPMD builds.

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

The topology convention inside ``/meshes`` follows VTK's unstructured-grid
model and ``vtkWedge`` cell. openPMD provides the mesh-record model, but it
does not standardize VTK-style wedge connectivity itself. HASEonGPU therefore
stores a VTK-compatible unstructured-cell layout in openPMD records:

* ``core_points`` stores VTK ``POINTS`` as ``x``, ``y``, and ``z`` components.
* ``core_cells_connectivity`` stores the VTK cell connectivity point ids.
* ``core_cells_offsets`` stores offsets into the connectivity array.
* ``core_cells_types`` stores the VTK cell type id for each cell.

For the current extruded-triangular-prism topology, every cell is a
``VTK_WEDGE`` with cell type ``13`` and six point ids. Connectivity uses the VTK
``vtkWedge`` point ordering/orientation.

Dynamic fields are also openPMD ``Mesh`` records. They are flattened in backend
order and stored with ``axisLabels`` ``["flatIndex"]``; the original primitive
axes and shape are retained in HASE metadata as ``haseAxes`` and
``hasePrimitiveShape``:

* ``core_beta_volume`` for prism-centered excited-state fraction.
* ``core_point_beta`` for point-level excited-state fraction.

Static material and spectral records include ``core_cladding_cell_type``,
``core_refractive_index``, ``core_reflectivity``,
``core_lambda_absorption``, ``core_lambda_emission``,
``core_sigma_absorption``, and ``core_sigma_emission``.

The C++ backend writes result records under ``core_result_``:
``phi_ase``, ``mse``, ``total_rays``, and ``dndt_ase``. Result records use
record-C layout with openPMD ``axisLabels`` ``["point", "level"]``.

User-defined fields from ``GainMedium.defineField(...)`` or
``PrimitiveFieldSpec`` inheritance are serialized through the same scalar
mesh-record path. They do not change the base openPMD contract: extra record
names are valid openPMD records as long as the normal record metadata is
present. The current ASE backend ignores user-defined records unless a future
backend explicitly opts in.

Iteration Updates
-----------------

An openPMD input series may contain multiple iterations.

The first Python-written iteration includes ``haseStaticUpdate = true`` and
contains topology, static material records, spectra, and dynamic fields.
Later iterations default to ``haseStaticUpdate = false`` and contain only
``core_beta_volume`` and ``core_point_beta``. The C++ parser caches the
iteration-0 topology, material, spectra, and compute settings, then applies
each dynamic-only iteration to that cached context.

After iteration 0, the parser's contract is intentionally dynamic-only:
``core_beta_volume`` and ``core_point_beta`` are the only records that may
change. Additional non-dynamic mesh records, ``haseStaticUpdate = true``, or
changed HASE-known topology/material/spectral/compute attributes are rejected
with an openPMD validation error. Repeating unchanged scalar attributes is
accepted for compatibility, but does not update the cached static context.

This split is important for streaming or repeated ASE evaluations because the
large static topology does not need to be rewritten for every dynamic state.
Changing topology, spectra, material constants, or compute settings requires a
new input series whose first iteration carries a complete static update.

Unsupported Options
-------------------

The openPMD transport does not preserve explicit device lists and does not
write VTK output from the backend. ``PhiASE.devices`` and ``PhiASE.writeVtk``
are therefore rejected by the Python transport when they request unsupported
behavior. The C++ parser also rejects ``write_vtk = true`` and explicit
``devices`` metadata.

Use the Python-side VTK helpers documented in :doc:`python_interface/utilities`
for visualization output.

MPI Launching
-------------

The standalone ``calcPhiASE`` binary can be launched under MPI and reads the
same openPMD transport layout:

.. code-block:: bash

   mpiexec -npernode 4 ./build/calcPhiASE \
       --input-path=input.bp \
       --output-path=output.bp

For Python-controlled runs that need an MPI launcher, pass a command prefix to
the transport helper:

.. code-block:: python

   result = transport.runPhiASE(
       phi_ase,
       medium,
       spectra,
       command_prefix=["mpiexec", "-npernode", "4"],
       workspace_dir="IO/phiase_mpi",
   )

The ``parallel_mode`` metadata is still written into the input series and
consumed by the backend compute configuration. The process topology itself is
controlled by how ``calcPhiASE`` is launched.

Artifact Retention
------------------

By default, temporary input and output series are removed when the Python
transport session exits. Use these environment variables when debugging the
transport files:

``HASE_OPENPMD_KEEP_ARTIFACTS=1``
   Keep artifacts below ``./hase-openpmd-artifacts``.

``HASE_OPENPMD_ARTIFACT_DIR=/path/to/dir``
   Write artifacts to an explicit directory.

``HASE_OPENPMD_ARTIFACT_PREFIX=name``
   Prefix generated artifact names.

``HASE_OPENPMD_ARTIFACT_RUN_ID=name``
   Use a stable run id instead of a timestamped id.

``HASE_OPENPMD_WATCHDOG_INTERVAL=seconds``
   Set the Python-side watchdog interval for streaming backends. The watchdog
   checks that the launched backend process is still alive while the result
   receiver waits. The default is ``30`` seconds. Use ``0`` or ``none`` to
   disable the watchdog.

``HASE_OPENPMD_THREAD_JOIN_TIMEOUT=seconds``
   Set how long the Python transport waits for streaming helper threads to stop
   during session close before reporting a cleanup error. The default is ``10``
   seconds.

``HASE_CALCPHIASE=/path/to/calcPhiASE``
   Force the Python transport to use a specific openPMD ``calcPhiASE`` binary.
