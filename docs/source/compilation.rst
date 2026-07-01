:doc:`<- Back to Getting Started <gettingStarted>`

:doc:`<- Back to Binary Interface <binaryInterface>`

Compilation
===========

This page describes how to build HASEonGPU manually from source.

For most users, manual compilation is only required when using the standalone
binary directly or when adjusting build options. For general setup and
dependency information, see :doc:`Getting Started <gettingStarted>`.

Basic Build
-----------

Clone the repository and build HASEonGPU with CMake:

.. code-block:: bash

   git clone https://github.com/computationalradiationphysics/haseongpu.git
   cd haseongpu
   cmake -S . -B build
   cmake --build build

After compilation, the ``calcPhiASE`` binary is available under:

.. code-block:: text

   ./build/calcPhiASE

When Python bindings are enabled, they are built alongside the C++ backend.

Typical Build Variants
----------------------

Minimal default build:

.. code-block:: bash

   cmake -S . -B build
   cmake --build build

Build with MPI support:

.. code-block:: bash

   cmake -S . -B build -DDISABLE_MPI=OFF
   cmake --build build

CMake Options
-------------

The following CMake variables control important build options.

``DISABLE_MPI``
^^^^^^^^^^^^^^^

* Default: ``AUTO``
* Description:
  Controls whether MPI support is required, disabled, or auto-detected.

* Values:

  * ``AUTO``: try to find MPI and continue without MPI if it is unavailable
  * ``OFF``: require MPI support; configuration fails if MPI is missing
  * ``ON``: disable MPI support

``HASE_BUILD_PhiAse``
^^^^^^^^^^^^^^^^^^^^^

* Default: ``ON``
* Description:
  Builds the standalone ``calcPhiASE`` command-line executable.  Disable this
  only when a build needs the libraries or Python package but not the binary.

``HASE_USE_SYSTEM_ALPAKA``
^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: ``OFF``
* Description:
  Uses an existing alpaka package from ``alpaka_DIR`` or ``CMAKE_PREFIX_PATH``
  instead of fetching the pinned alpaka version during CMake configuration.

``HASE_CUDA_ARCHITECTURES``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: ``native`` with fallback to ``80`` when no local NVIDIA GPU is visible
* Description:
  Selects the CUDA target architecture used for compilation.

* Typical values:

  * ``native``: detect the local GPU architecture automatically
  * explicit CUDA architectures such as ``80``, ``86``, ``89``, ``90``

Using ``native`` is convenient for local builds when CMake can query a local GPU.
For reproducible builds on different systems, specifying the CUDA architecture is recommended.

``HASE_ENABLE_PYTHON``
^^^^^^^^^^^^^^^^^^^^^^

* Default: ``ON``
* Description:
  Controls whether the Python extension and package helpers are built. Turn
  this off only for command-line-only C++ builds. For normal Python
  installation and usage, see :doc:`Python Interface Guide <pythonInterface>`.

* Values:

  * ``OFF``: build only the C++ project and binary interface
  * ``ON``: build the Python interface

``HASE_BUILD_RELEASE``
^^^^^^^^^^^^^^^^^^^^^^

* Default: ``ON``
* Description:
  Controls whether HASEonGPU applies its release build configuration.  When
  enabled, CMake forces ``CMAKE_BUILD_TYPE=Release`` and enables the release
  optimization options used by the project, including CUDA/HIP fast-math
  related flags where applicable.

  Important: ``HASE_BUILD_RELEASE=ON`` overwrites user-provided
  ``CMAKE_BUILD_TYPE`` values and related optimization settings during
  configuration.  Set ``HASE_BUILD_RELEASE=OFF`` if you need a custom build
  type, debug flags, or manually controlled compiler optimization options.

* Values:

  * ``OFF``: keep user-provided build type and optimization settings
  * ``ON``: force the project release configuration

``HASE_NATIVE_OPTIMIZATIONS``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: ``ON``
* Description:
  Enables host-specific CPU tuning with ``-march=native`` and
  ``-mtune=native``. Keep this enabled for local source builds when peak
  performance on the build machine is desired. Disable it for redistributable
  binaries or wheels that need to run on unknown CPUs.

* Values:

  * ``OFF``: do not add host-specific native CPU tuning flags
  * ``ON``: build for the local host CPU

For Python source installs, pass the option through ``CMAKE_ARGS``:

.. code-block:: bash

   CMAKE_ARGS="-DHASE_NATIVE_OPTIMIZATIONS=ON" python3 -m pip install .

For redistributable wheels or binaries, configure with:

.. code-block:: bash

   CMAKE_ARGS="-DHASE_NATIVE_OPTIMIZATIONS=OFF" python3 -m pip install .

``HASE_SELECT_BACKEND_ALPAKA``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: ``OFF``

* Description:
  Controls whether HASEonGPU selects available alpaka backends automatically or
  whether backend selection is delegated to alpaka's CMake options.

  For general information about backend names and runtime backend selection,
  see :doc:`Backend Selection <backendSelection>`.

* Values:

  * ``OFF``:
    HASEonGPU automatically searches for supported backend dependencies such as
    CUDA, HIP, and TBB and enables the corresponding alpaka backends when
    possible.

    If both HIP and CUDA are installed on the same system, automatic detection
    may cause configuration issues in alpaka. In that case, manual backend
    selection should be used to explicitly disable one of the conflicting
    backends.

  * ``ON``:
    Enables manual backend selection using alpaka's existing CMake options.

    The relevant alpaka CMake options are documented in the
    `alpaka CMake argument documentation <https://alpaka3.readthedocs.io/en/latest/advanced/cmake.html#arguments>`__.

    Example: configure HASEonGPU only for an NVIDIA GPU backend:

    .. code-block:: bash

       cmake -S . -B build \
         -DHASE_SELECT_BACKEND_ALPAKA=ON \
         -Dalpaka_DEP_CUDA=ON \
         -Dalpaka_DEP_HIP=OFF \
         -Dalpaka_DEP_TBB=OFF \
         -Dalpaka_EXEC_CpuSerial=OFF


openPMD Provider and Runtime Backend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CMake build selects the openPMD provider, not the storage backend used by a
simulation. Runtime backend selection belongs to Python or YAML through
``PhiASE.openpmdBackend`` or ``openpmd_backend``. The default runtime backend is
``adios-sst``; accepted values are ``adios-sst``, ``adios``, and ``hdf5``.

Provider modes are separate from runtime backend selection:

* External provider, ``HASE_BUILD_OPENPMD_FROM_SOURCE=OFF``: HASEonGPU uses an
  already installed C++ ``openPMD::openPMD`` package and a compatible Python
  ``openpmd_api`` module. This mode does not build ADIOS2. For source-built
  external openPMD-api providers, install or build ADIOS2 first when the
  provider must support ``adios`` or ``adios-sst``, then expose it through
  openPMD-api's ``CMAKE_PREFIX_PATH`` or ``ADIOS2_DIR``. For ADIOS2 source
  installs used this way, pass ``-DADIOS2_INSTALL_GENERATE_CONFIG=OFF`` unless
  you need the legacy ``adios2-config`` shell helper; HASEonGPU and openPMD-api
  use ADIOS2's CMake package config. openPMD-api 0.17.x does not FetchContent
  ADIOS2 through ``openPMD_SUPERBUILD``; it only finds an existing ADIOS2
  package.
* Bundled source provider, ``HASE_BUILD_OPENPMD_FROM_SOURCE=ON``: HASEonGPU
  uses FetchContent to build ADIOS2 before configuring the pinned openPMD-api
  provider, then enables ADIOS2, ADIOS2 SST, and HDF5 support for the CMake
  build.

Use ``utils/check_openpmd_compatibility.py`` before installing when in doubt.
The preflight prints the selected backend, the backends confirmed by the active
Python and CMake providers, and the built-in HASE source-build backend set
(``adios-sst`` by default, plus ``adios`` and ``hdf5``).

The CI matrix keeps most jobs on the bundled source-build provider and includes
a dedicated Ubuntu job that runs an ordered external-provider smoke matrix in a
single container. For each pinned openPMD-api version, it builds and
installs the external provider prerequisites, builds an openPMD-api C++ prefix
against those prerequisites, installs the matching Python package, installs
HASEonGPU against that prefix, and then checks ``adios``, ``adios-sst``, and
``hdf5`` smoke scenarios. The current pinned external versions are ``0.17.0``
and ``0.17.1``.
This keeps external-versus-bundled provider selection covered without
multiplying every compiler and accelerator combination.

``HASE_BUILD_OPENPMD_FROM_SOURCE``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: ``OFF``
* Description:
  Controls the openPMD install contract.

  * ``OFF``: use an externally installed openPMD-api C++ package found through
    ``find_package(openPMD CONFIG REQUIRED)``. The Python ``openpmd_api``
    package installed in the runtime environment must come from the same
    openPMD provider family and must support the runtime backend used by your
    workflow.
  * ``ON``: use the FetchContent source-build provider. HASEonGPU builds ADIOS2
    before configuring the pinned openPMD-api provider and enables ADIOS2,
    ADIOS2 SST, and HDF5 support for the CMake build.

The HASE wheel does not vendor openPMD runtime libraries or generated
``openpmd_api`` bindings in either mode. The target runtime environment must
provide compatible openPMD shared libraries and a compatible Python
``openpmd_api`` package. HASEonGPU records the detected openPMD provider
library directory in the installed targets' RPATH; set
``HASE_OPENPMD_RUNTIME_RPATH`` when the provider needs additional runtime
library directories.

``HASE_OPENPMD_SUPERBUILD``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: ``ON`` unless inherited from ``openPMD_SUPERBUILD``
* Description:
  Applies only when ``HASE_BUILD_OPENPMD_FROM_SOURCE=ON``. Allows the pinned
  openPMD-api build to fetch or build its helper dependencies. This is not an
  ADIOS2 install switch; HASE's source-build path uses FetchContent to build
  ADIOS2 before openPMD-api is configured. Disable this only when the required
  openPMD dependencies are provided externally.

``HASE_USE_SYSTEM_OPENPMD``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: deprecated
* Description:
  Deprecated compatibility alias. Use ``HASE_BUILD_OPENPMD_FROM_SOURCE``
  instead. ``HASE_USE_SYSTEM_OPENPMD=ON`` maps to
  ``HASE_BUILD_OPENPMD_FROM_SOURCE=OFF``; ``HASE_USE_SYSTEM_OPENPMD=OFF`` maps
  to ``HASE_BUILD_OPENPMD_FROM_SOURCE=ON``.

``HASE_OPENPMD_BUILD_PYTHON_BINDINGS``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: ``OFF``
* Description:
  Builds openPMD-api Python bindings as part of the HASE CMake build tree.
  These bindings are not installed into the HASE wheel. The runtime Python
  module must be supplied by the selected external openPMD installation.

``HASE_OPENPMD_PYTHON_PACKAGE_DIR``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: empty
* Description:
  Directory containing the matching ``openpmd_api`` Python package. Normally
  leave this empty so the installed Python environment supplies
  ``openpmd_api``. Set it only when the runtime must prefer a specific
  site-packages directory from the same external openPMD installation as
  ``openPMD::openPMD``.

``HASE_OPENPMD_RUNTIME_RPATH``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: empty
* Description:
  Additional semicolon-separated runtime library directories to encode into
  the installed HASE targets. Use this for non-vendored openPMD providers when
  ``openPMD::openPMD`` does not expose every required runtime directory through
  its CMake target metadata.

``HASE_TESTING``
^^^^^^^^^^^^^^^^

* Default: ``OFF``
* Description:
  Enables the test suite during configuration and build.

* Values:

  * ``OFF``: tests are not built
  * ``ON``: test targets are enabled

This option is primarily useful for development and validation work.

Notes
-----

Manual compilation is usually not required for every workflow.

For example, when using the Python interface, the backend may be built as part
of the Python installation process. However, a manual build can be useful for
adjusting CMake options, debugging, or working directly with the standalone
binary.
