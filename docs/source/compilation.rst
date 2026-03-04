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
   make

After compilation, the ``calcPhiASE`` binary is available under:

.. code-block:: text

   ./build/calcPhiASE

When Python bindings are enabled, they are built alongside the C++ backend.

Typical Build Variants
----------------------

Minimal default build:

.. code-block:: bash

   cmake ..
   make

Build with MPI support:

.. code-block:: bash

   cmake .. -DDISABLE_MPI=OFF
   make

Build with GrayBat support:

.. code-block:: bash

   cmake .. -DUSE_GRAYBAT=ON
   make

CMake Options
-------------

The following CMake variables control important build options.

``DISABLE_MPI``
^^^^^^^^^^^^^^^

* Default: ``OFF``
* Description:
  Enabling allows compilation without requiring MPI or BoostMPI as a dependency.

* Values:

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

* Default: ``OFF``
* Description:

  This option usually does not need to be enabled manually unless you want to
  customize the Python interface build. For normal Python installation and
  usage, please refer to :doc:`Python Interface <pythonInterface>`.


* Values:

  * ``OFF``: build only the C++ project and binary interface
  * ``ON``: build the Python interface

``HASE_RELEASE``
^^^^^^^^^^^^^^^^

* Default: ``ON``
* Description:
  Controls whether HASEonGPU is built in release mode.

* Values:

  * ``OFF``: non-release build
  * ``ON``: release-oriented build

  This option is enabled by default for performance-oriented builds. Disabling
  it allows internal runtime assertions and is therefore mainly useful for
  debugging and development.

``HASE_TESTING``
^^^^^^^^^^^^^^^^

* Default: ``OFF``
* Description:
  Enables the test suite during configuration and build.

* Values:

  * ``OFF``: tests are not built
  * ``ON``: test targets are enabled

This option is primarily useful for development and validation work.

``USE_GRAYBAT``
^^^^^^^^^^^^^^^

* Default: ``OFF``
* Description:
  Enables compilation with the experimental GrayBat support.
  When deactivated graybat is not required as a dependency.

* Values:

  * ``OFF``: GrayBat support is disabled
  * ``ON``: build with GrayBat support

For more information, see :doc:`graybat`.

Notes
-----

Manual compilation is usually not required for every workflow.

For example, when using the Python interface, the backend may be built as part
of the Python installation process. However, a manual build can be useful for
adjusting CMake options, debugging, or working directly with the standalone
binary.