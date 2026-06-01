:doc:`<- Back to Getting Started <gettingStarted>`

:doc:`<- Back to Matlab Interface <MATLABInterface>`

:doc:`<- Back to Binaray Interface <binaryInterface>`

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
   mkdir build
   cd build
   cmake ..
   cmake --build .

After compilation, the ``calcPhiASE`` binary is available under:

.. code-block:: text

   ./build/calcPhiASE

When Python bindings are enabled, they are built alongside the C++ backend.

Typical Build Variants
----------------------

Minimal default build:

.. code-block:: bash

   cmake ..
   cmake --build .

Build with MPI support:

.. code-block:: bash

   cmake .. -DDISABLE_MPI=OFF
   cmake --build .

CMake Options
-------------

The following CMake variables control important build options.

``DISABLE_MPI``
^^^^^^^^^^^^^^^

* Default: ``AUTO``
* Description:
  Enabling allows compilation without requiring MPI or BoostMPI as a dependency.

* Values:

  * ``AUTO``: CMAKE tries to detect whether MPI exists and sets this knob in correspondence
  * ``OFF``: MPI support remains - dependencies are required
  * ``ON``: MPI support is disabled

``HASE_CUDA_ARCHITECTURES``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Default: ``native``
* Description:
  Selects the CUDA target architecture used for compilation.

* Typical values:

  * ``native``: detect the local GPU architecture automatically
  * explicit CUDA architectures such as ``75``, ``80``, ``86``

Using ``native`` is convenient for local builds. For reproducible and performance on
different systems, specifying the CUDA architecture is recommended.

``HASE_ENABLE_PYTHON``
^^^^^^^^^^^^^^^^^^^^^^

* Default: ``ON``
* Description:

  If python as a dependency is missing on your system this knob can be turned off in order to use HASEonGPU from the command-line only.
  For normal Python installation and usage, please refer to :doc:`Python Interface Guide <pythonInterface>`.


* Values:

  * ``OFF``: build only the C++ project and binary interface
  * ``ON``: build the Python interface

``HASE_BUILD_RELEASE``
^^^^^^^^^^^^^^^^^^^^^^

* Default: ``ON``
* Description:
  Controls whether HASEonGPU applies its release build configuration.  When
  enabled, CMake forces ``CMAKE_BUILD_TYPE=Release`` and enables the release
  optimization options used by the project, including native host compiler
  tuning flags and CUDA/HIP fast-math related flags where applicable.

  Important: ``HASE_BUILD_RELEASE=ON`` overwrites user-provided
  ``CMAKE_BUILD_TYPE`` values and related optimization settings during
  configuration.  Set ``HASE_BUILD_RELEASE=OFF`` if you need a custom build
  type, debug flags, or manually controlled compiler optimization options.

* Values:

  * ``OFF``: keep user-provided build type and optimization settings
  * ``ON``: force the project release configuration

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
