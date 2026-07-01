:doc:`<- Back to overview <index>`

Getting Started
===============

This guide covers setting up HASEonGPU and choosing an interface for a
simulation workflow. For a compact project overview, see
:doc:`HASEonGPU Documentation <index>`. For the ASE model and estimator, see
:doc:`Theory and Model <theoryAndModel>`.

HASEonGPU supports multiple usage paths:

* :doc:`Python Interface Guide <pythonInterface>` for workflow-oriented Python usage
* :doc:`Binary Interface <binaryInterface>` for command-line execution
* :doc:`openPMD Transport <openpmdTransport>` for the Python-to-C++ transport contract

Repository Setup
----------------

Clone the repository from GitHub:

.. code-block:: bash

   git clone https://github.com/computationalradiationphysics/haseongpu.git
   cd haseongpu

Dependencies
------------

Required software:

* ``cmake``
* ``ninja``
* ``Python >= 3.10`` with ``pip``
* openPMD-api C++ package with the backends used by HASEonGPU
* A C++20-capable compiler, tested with:

  * ``gcc >= 12``
  * ``clang >= 17``

Optional backend/runtime dependencies:

* ``cuda >= 12.5`` tested up to ``cuda 13.2``
* ``hip/rocm >= 6.2.4`` tested up to ``rocm 7.2``
* ``hwloc`` unless disabled

Optional software and tools:

* ``OpenMPI >= 4.0``
* ``ParaView`` for visualization of ``.vtk`` output
* ``matplotlib`` for the optional plotting helper in :doc:`scripts`

Hardware requirements:

* CPU execution works through host alpaka backends.
* GPU execution requires a backend-supported device, for example NVIDIA/CUDA or
  AMD/HIP hardware, and a build configured for that backend.

Additional Notes
----------------

For Windows-specific limitations, see :doc:`windows`.

Compilation Notes
-----------------

A manual compilation step is available, but is not required for every workflow.

For example, when using the Python interface, the C++ backend is built under
the hood during installation. For details on manual compilation, see
:doc:`compilation`.

Recommended Python Source Install
---------------------------------

For performance-sensitive Python use, install from source so the C++ backend is
compiled on the target machine. The default workflow uses a PyPI
``openpmd-api`` Python wheel and a local external openPMD-api C++ prefix that
CMake can find.

1. Create and activate a Python environment:

.. code-block:: bash

   python3 -m venv .venv
   source .venv/bin/activate
   python3 -m pip install -U pip

2. Provide openPMD-api for both Python and CMake.

For Conda environments, use the same environment for Python and CMake:

.. code-block:: bash

   conda install -c conda-forge openpmd-api
   python3 utils/check_openpmd_compatibility.py \
     --backend adios-sst \
     --cmake-prefix-path "$CONDA_PREFIX"

For Spack or module-provided openPMD installations, load the provider first and
use its prefix. The Python interpreter and every ``openpmd_api`` entry on
``PYTHONPATH`` must use the same Python ABI; for example, do not run Python
3.14 with a module path that injects ``lib/python3.11/site-packages``.

.. code-block:: bash

   spack load openpmd-api
   python3 -c "import sys; print(sys.version); print(sys.path)"
   python3 utils/check_openpmd_compatibility.py \
     --backend adios-sst \
     --cmake-prefix-path "$OPENPMD_API_ROOT"

For a manually installed source build of openPMD-api, install the matching
Python package in the active environment and point CMake to the C++ install
prefix:

.. code-block:: bash

   python3 -m pip install openpmd-api
   python3 utils/check_openpmd_compatibility.py \
     --backend adios-sst \
     --cmake-prefix-path /path/to/openpmd/prefix

3. Install HASEonGPU:

.. code-block:: bash

   CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/path/to/openpmd/prefix -DHASE_NATIVE_OPTIMIZATIONS=ON" \
     python3 -m pip install .

Use ``$CONDA_PREFIX`` or ``$OPENPMD_API_ROOT`` in place of
``/path/to/openpmd/prefix`` for the Conda or Spack/module paths above.

This builds and installs a release-mode backend with host-specific CPU tuning.
Use ``HASE_NATIVE_OPTIMIZATIONS=OFF`` only when building redistributable wheels
or binaries for unknown CPUs. See :doc:`compilation` for the full CMake option
reference.

4. Select the runtime openPMD backend from Python or YAML when needed:

.. code-block:: python

   phi_ase = PhiASE(..., openpmdBackend="adios-sst")

.. code-block:: yaml

   compute:
     openpmd_backend: adios-sst

Bundled Build Provider
^^^^^^^^^^^^^^^^^^^^^^

If no compatible external C++ provider is available during the HASE build, use
the FetchContent source-build provider:

.. code-block:: bash

   CMAKE_ARGS="-DHASE_BUILD_OPENPMD_FROM_SOURCE=ON" python3 -m pip install .

This path fetches and builds the pinned openPMD-api provider with ADIOS2,
ADIOS2 SST, and HDF5 support for the HASE CMake build. The HASE wheel does not
vendor the resulting openPMD runtime libraries or generated ``openpmd_api``
Python bindings. The target runtime environment must still provide compatible
openPMD shared libraries and a compatible Python ``openpmd_api`` package.

For redistributable wheels, prefer a normal external provider and make the
runtime dependency explicit in the deployment environment.

Optional Provider Notes
^^^^^^^^^^^^^^^^^^^^^^^

If you already have a matching external openPMD-api installation, install the
Python package in the active environment and point CMake to the C++ provider:

.. code-block:: bash

   python3 -m pip install openpmd-api
   python3 utils/check_openpmd_compatibility.py \
     --backend adios-sst \
     --cmake-prefix-path /path/to/openpmd/prefix
   CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/path/to/openpmd/prefix" python3 -m pip install .

Troubleshooting openPMD Discovery
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If CMake cannot find ``openPMD::openPMD``, pass either the installation prefix
or the exact config directory:

.. code-block:: bash

   CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/path/to/openpmd/prefix" python3 -m pip install .
   CMAKE_ARGS="-DopenPMD_DIR=/path/to/openPMDConfig.cmake-directory" python3 -m pip install .

If ``utils/check_openpmd_compatibility.py`` reports missing backend support,
use one of these paths:

* provide a Python ``openpmd_api`` package and CMake prefix that support the
  runtime backend you want, then select that backend in Python or YAML
* use ``HASE_BUILD_OPENPMD_FROM_SOURCE=ON`` to build the bundled provider with
  all built-in HASE openPMD backends

Choose the same MPI setting for the Python and C++ providers. With bundled
openPMD, HASEonGPU derives openPMD MPI support from ``DISABLE_MPI``.

Choose an Interface
-------------------

Python Interface Guide
^^^^^^^^^^^^^^^^^^^^^^

The Python interface guide is the recommended starting point for most new users.
It explains the high-level workflow and shows how to assemble geometry, material
data, pump settings, ASE configuration, and time stepping in Python.

Continue with :doc:`Python Interface Guide <pythonInterface>`.  When you need
exact signatures or generated member lists, use the :doc:`Python API Reference <pythonAPI>`.

Binary Interface
^^^^^^^^^^^^^^^^

The binary interface is useful for running HASEonGPU directly from the command
line or as an entry point for constructing a custom API call.

Continue with :doc:`Binary Interface <binaryInterface>`. For the openPMD record layout and transport backends used by both Python and the binary, see :doc:`openPMD Transport <openpmdTransport>`.

Typical Workflow
----------------

A typical HASEonGPU setup consists of the following steps:

#. Clone the repository
#. Install the required dependencies
#. Decide which interface to use
#. Follow installation steps provided on the page of the chosen interface
#. Verify the setup with one of the provided interface-specific examples
#. Inspect generated VTK output directly or with the helper scripts
#. Use HASEonGPU in the target workflow

Next Step
---------

Proceed to one of the interface pages listed above, depending on the intended
usage path.
