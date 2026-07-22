:doc:`<- Back to Getting Started <gettingStarted>`

:doc:`<- Back to Binary Interface <binaryInterface>`

CMake Build Options
===================

This page summarizes manual CMake configuration for HASEonGPU.  Most Python
users should start with :doc:`Getting Started <gettingStarted>` and let
``hase-configure`` print the matching ``CMAKE_ARGS`` command.

Manual Build
------------

.. code-block:: bash

   git clone https://github.com/computationalradiationphysics/haseongpu.git
   cd haseongpu
   cmake -S . -B build
   cmake --build build

The standalone binary is then available as:

.. code-block:: text

   ./build/calcPhiASE

For Python source installs, ``python3 -m pip install .`` configures or reuses
this same ``build/`` directory and installs only the thin frontend. Pass custom
CMake cache options through ``CMAKE_ARGS``:

.. code-block:: bash

   CMAKE_ARGS="-DHASE_OPENPMD_PROVIDER=system -DCMAKE_PREFIX_PATH=/path/to/openpmd" \
     python3 -m pip install -v .

Common Configurations
---------------------

CPU-only first build:

.. code-block:: bash

   cmake -S . -B build -DHASE_NATIVE_OPTIMIZATIONS=OFF
   cmake --build build

Require MPI support:

.. code-block:: bash

   cmake -S . -B build -DDISABLE_MPI=OFF
   cmake --build build

Use a system openPMD-api provider:

.. code-block:: bash

   cmake -S . -B build \
     -DHASE_OPENPMD_PROVIDER=system \
     -DCMAKE_PREFIX_PATH=/path/to/openpmd/prefix

Use HASE-managed bundled openPMD dependencies:

.. code-block:: bash

   cmake -S . -B build -DHASE_OPENPMD_PROVIDER=bundled

Manually select Alpaka backends, for example CUDA only:

.. code-block:: bash

   cmake -S . -B build \
     -DHASE_SELECT_BACKEND_ALPAKA=ON \
     -Dalpaka_DEP_CUDA=ON \
     -Dalpaka_DEP_HIP=OFF \
     -Dalpaka_DEP_TBB=OFF \
     -Dalpaka_EXEC_CpuSerial=OFF

Core HASE Options
-----------------

``HASE_BUILD_RELEASE``
   Default ``ON``.  Forces the project release configuration, including
   ``CMAKE_BUILD_TYPE=Release`` and release optimization flags.  Set ``OFF``
   when you need a custom debug or profiling build type.

``HASE_NATIVE_OPTIMIZATIONS``
   Default ``OFF``. Adds host-specific
   CPU tuning such as ``-march=native``.  Enable for local performance builds;
   disable for redistributable wheels or unknown target CPUs.

``HASE_ENABLE_PYTHON``
   Default ``ON``. Generates the Python frontend metadata describing the
   matching native runtime. Set ``OFF`` for command-line-only C++ builds.

``HASE_BUILD_RUNTIME``
   Default ``ON`` for ordinary CMake builds. The thin-wheel build sets this to
   ``OFF`` and builds the runtime separately in the durable directory selected
   by ``HASE_RUNTIME_DIR``.

``HASE_RUNTIME_DIR``
   Default ``<source>/build`` for a thin frontend install. It is the single
   root used for the executable, probe libraries, and generated openPMD Python
   provider metadata. Prefer ``hase-configure --runtime-dir <path>`` when
   selecting a nonstandard location.

``HASE_TESTING``
   Default ``OFF``.  Enables test targets.

``HASE_BENCHMARK``
   Default ``OFF``.  Enables scoped PhiASE benchmark CSV output.

``HASE_FORWARD_LOGGING``
   Default ``OFF``.  Forwards captured ``calcPhiASE`` stdout/stderr through the
   Python openPMD launcher.

MPI Option
----------

``DISABLE_MPI``
   Default ``AUTO``.

   * ``AUTO``: use MPI if CMake can find it, otherwise build without MPI
   * ``OFF``: require MPI; configuration fails if MPI is unavailable
   * ``ON``: disable MPI support

The runtime ``parallel_mode``/``parallelMode`` setting must agree with how the
binary was built and launched.  See :doc:`mpi` for execution guidance.

Alpaka and Accelerator Options
------------------------------

``HASE_SELECT_BACKEND_ALPAKA``
   Default ``OFF``.  ``OFF`` lets HASEonGPU auto-detect supported Alpaka
   dependencies such as CUDA, HIP, and TBB.  ``ON`` delegates backend selection
   to Alpaka CMake options such as ``alpaka_DEP_CUDA`` and
   ``alpaka_EXEC_CpuSerial``.

``HASE_USE_SYSTEM_ALPAKA``
   Default ``OFF``.  Uses an installed Alpaka package found through
   ``alpaka_DIR`` or ``CMAKE_PREFIX_PATH`` instead of fetching the pinned
   version.

``HASE_CUDA_ARCHITECTURES``
   Default ``native`` with fallback architectures when no local NVIDIA GPU is
   visible.  Set explicit values such as ``80`` or ``90`` for reproducible CUDA
   builds on machines different from the target system.

``HASE_CUDA_FLUSHTOZERO``
   Default ``OFF``.  Enables CUDA flush-to-zero behavior.

CUDA debugging options
   ``HASE_CUDA_SHOW_REGISTER``, ``HASE_CUDA_KEEP_FILES``, and
   ``HASE_CUDA_SHOW_CODELINES`` default to ``OFF``.  They are debugging aids
   for CUDA compilation and diagnostics.

openPMD Provider Options
------------------------

CMake chooses the openPMD-api provider used to build HASEonGPU.  The storage
backend used by a simulation is selected at runtime with
``PhiASE.openpmdBackend`` or YAML ``openpmd_backend``.

``HASE_OPENPMD_PROVIDER``
   Default ``auto``.

   * ``auto``: use a system openPMD-api C++ package if found, otherwise build a
     bundled provider
   * ``system``: require an installed ``openPMD::openPMD`` package; the runtime
     Python environment must separately provide the matching ``openpmd_api``
     module (it is intentionally not an unconditional Python dependency)
   * ``bundled``: build the pinned openPMD-api provider and selected
     dependencies and matching Python bindings into a build-local prefix; the
     HASE frontend selects those bindings automatically at process startup

``HASE_OPENPMD_USE_ADIOS2``
   Default ``ON`` for the bundled provider.  Enables ADIOS2-backed runtime
   backends.

``HASE_OPENPMD_USE_SST``
   Default ``ON`` when ADIOS2 is enabled.  Enables the ``adios-sst`` runtime
   backend.

``HASE_OPENPMD_USE_HDF5``
   Default ``OFF`` for the bundled provider.  Enable when the ``hdf5`` runtime
   backend is required.

``HASE_OPENPMD_FETCH_ADIOS2`` / ``HASE_OPENPMD_FETCH_HDF5``
   Default ``ON``.  With the bundled provider, fetch and build the pinned
   dependency.  Set ``OFF`` to use a system package found through
   ``CMAKE_PREFIX_PATH``, ``ADIOS2_DIR``, or ``HDF5_DIR``.

``HASE_OPENPMD_SUPERBUILD``
   Default ``ON`` unless inherited from ``openPMD_SUPERBUILD``.  Allows
   openPMD-api to handle its helper dependencies.  This is not the ADIOS2/HDF5
   selection switch.

``HASE_OPENPMD_BUILD_PYTHON_BINDINGS``
   Default ``OFF`` in CMake; the configurator may enable it for bundled local
   installs.  Builds matching ``openpmd_api`` Python bindings in the build tree;
   they are not vendored into the HASEonGPU wheel.

``HASE_OPENPMD_PYTHON_PACKAGE_DIR``
   Default empty.  Path to a matching ``openpmd_api`` package directory when it
   is not found on the normal Python path.

``HASE_OPENPMD_RUNTIME_RPATH``
   Default empty.  Extra semicolon-separated runtime library directories to
   encode into installed HASEonGPU targets.

``HASE_OPENPMD_BUNDLED_PREFIX`` / ``HASE_OPENPMD_BUNDLED_BUILD_DIR``
   Default below the CMake build directory.  Override only when you need a
   custom location for bundled-provider artifacts.

``HASE_OPENPMD_BUNDLED_REBUILD``
   Default ``OFF``.  Forces rebuilding the HASE-managed bundled provider.

Deprecated aliases ``HASE_BUILD_OPENPMD_FROM_SOURCE`` and
``HASE_USE_SYSTEM_OPENPMD`` are still recognized for compatibility, but new
commands should use ``HASE_OPENPMD_PROVIDER``.

Runtime Backend Reminder
------------------------

Build options decide which providers and compute backends are available.
Runtime selection is still done through user configuration:

.. code-block:: python

   phi_ase = PhiASE(..., backend="Host_Cpu_CpuSerial", openpmdBackend="auto")

.. code-block:: yaml

   compute:
     backend: Host_Cpu_CpuSerial
     openpmd_backend: auto

Use ``python3 utils/configure_hase.py`` or ``hase-configure`` to generate a
matching YAML snippet and install command for common setups.
