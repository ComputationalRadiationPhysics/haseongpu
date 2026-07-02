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

Guided Setup
------------

For a first installation, use the guided configurator instead of guessing the
openPMD provider, transport backend, MPI mode, and Alpaka backend names. From a
source checkout run:

.. code-block:: bash

   python3 utils/configure_hase.py

After installation the same helper is available as:

.. code-block:: bash

   hase-configure

The first guide question only chooses between the bundled source-build provider
and an external openPMD-api provider. For an external provider, pass the
openPMD CMake prefix or ``openPMD_DIR`` and let the helper validate the active
Python ``openpmd_api`` module against that C++ provider. If the external
openPMD-api was built with ADIOS2 support and CMake cannot find ADIOS2
transitively, the guide can add an ADIOS2 prefix or ``ADIOS2_DIR`` hint.
HASEonGPU does not install a matching external-provider Python ``openpmd_api``
module for you; install/load one from the same provider family.

For the bundled provider, the guide asks separately how to handle ADIOS2:
fetch/build pinned ADIOS2, build an HDF5-only provider without ADIOS2, or use a
system ADIOS2 installation. By default it also builds matching openPMD Python
bindings and records their build-tree ``site-packages`` path in the installed
HASE configuration. Those bindings are not copied into the HASE wheel, but
HASE can use them from outside the source checkout as long as the build-tree
provider path remains available. Fetching the bundled dependencies is the
easiest starting point, but it can take noticeably longer than using already
installed providers.

After the openPMD setup, the guide points out the optional MPI path
(``parallel_mode: mpi`` and ``n_per_node``) and explains that Python and C++
openPMD providers must use compatible MPI settings. It then lists installed
Alpaka compute backends when the backend-name helper is available. Use a CPU
host backend for first validation; CUDA and HIP backends require matching
hardware, compiler/toolkit setup, and Alpaka build options.

The command writes a small PhiASE YAML run-control file under
``config/hase-phiase.yaml`` by default and prints a ready-to-run install snippet
using ``HASE_CONFIGURE_CMAKE_ARGS`` and ``CMAKE_ARGS``. At the
end of an interactive run it asks whether to install immediately; the default
answer is yes. Answering no keeps the generated files and repeats the exact
install command to run later. If pip reports an externally managed Python
environment, prefer a virtual environment. If you intentionally install into
that environment anyway, launch the guide with ``--break-system-packages`` so
the printed command and optional install add pip's ``--break-system-packages``
flag; the guide does not prompt for this option interactively.
The native-optimization prompt defaults to on for local performance. It enables host-specific CPU tuning,
so turn it off for redistributable wheels or unknown CPUs. The final guidance
says where the configuration file is present and that it can be modified,
lists the supported ``openpmd_backend`` values for the selected provider,
summarizes MPI keys, and explains that if no GPU alpaka backend is listed,
alpaka probably did not find a usable GPU toolchain/backend or matching device.
It points CUDA users to ``CUDACXX`` or ``nvcc``/CUDA toolkit paths and HIP
users to ``HIPCXX``,
``hipcc``, or ``ROCM_PATH``. Pass the YAML to ``PhiASE.fromYaml(...)`` and keep
physics objects such as geometry, spectra, and pump settings in Python. The
laser-pump cladding examples use the same ``config/hase-phiase.yaml`` path;
running ``hase-configure`` overwrites that file with the selected setup.

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

Manual external source builds need one extra dependency step. Build or install
ADIOS2 before configuring openPMD-api when the provider must support ``adios``
or ``adios-sst``, and include the ADIOS2 prefix in openPMD-api's
``CMAKE_PREFIX_PATH``. For ADIOS2 source installs used as an openPMD provider
dependency, pass ``-DADIOS2_INSTALL_GENERATE_CONFIG=OFF`` unless you need the
legacy ``adios2-config`` shell helper; HASEonGPU and openPMD-api use ADIOS2's
CMake package config. openPMD-api 0.17.x does not fetch ADIOS2 through
``openPMD_SUPERBUILD``; that option only covers openPMD helper dependencies.
After installing openPMD-api, install the matching Python package in the active
environment and point CMake to the C++ install prefix.

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

This path is the only HASE install mode that uses FetchContent to build ADIOS2
before configuring openPMD-api. It builds the pinned openPMD-api provider with
ADIOS2, ADIOS2 SST, and HDF5 support for the HASE CMake build.

The HASE wheel does not vendor the resulting ADIOS2/openPMD runtime libraries
or generated ``openpmd_api`` Python bindings. The target runtime environment
must still provide compatible openPMD shared libraries and a compatible Python
``openpmd_api`` package; HASE records provider library directories through
RPATH rather than copying provider libraries into the wheel.

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
  runtime backend you want, then select that backend in Python or YAML; for a
  manual external openPMD-api source build, install ADIOS2 first when selecting
  ``adios`` or ``adios-sst``
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
