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
* :doc:`MATLAB Interface <MATLABInterface>` for integration into existing
  MATLAB or Octave workflows

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
* ``Python >= 3.11`` with ``pip``
* A C++20-capable compiler, tested with:

  * ``gcc >= 12``
  * ``clang >= 17``

Optional backend/runtime dependencies:

* ``cuda >= 12.5`` tested up to ``cuda 13.2``
* ``hip/rocm >= 6.2.4`` tested up to ``rocm 7.2``
* ``hwloc`` unless disabled

Optional software and tools:

* ``OpenMPI >= 4.0``
* MATLAB or Octave
* ``ParaView`` for visualization of ``.vtk`` output

Hardware requirements:

* CPU execution works through host alpaka backends.
* GPU execution requires a backend-supported device, for example NVIDIA/CUDA or
  AMD/HIP hardware, and a build configured for that backend.

Additional Notes
----------------

For Windows-specific installation notes [deprecated], see :doc:`windows`.

Compilation Notes
-----------------

A manual compilation step is available, but is not required for every workflow.

For example, when using the Python interface, the C++ backend is built under
the hood during installation. For details on manual compilation, see
:doc:`compilation`.

Recommended Python Source Install
---------------------------------

For performance-sensitive Python use, install from source so the C++ backend is
compiled on the target machine:

.. code-block:: bash

   CMAKE_ARGS="-DHASE_NATIVE_OPTIMIZATIONS=ON" python3 -m pip install .

This builds and installs a release-mode backend with host-specific CPU tuning.
Use ``HASE_NATIVE_OPTIMIZATIONS=OFF`` only when
building redistributable wheels or binaries for unknown CPUs. See
:doc:`compilation` for the full CMake option reference.

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

Continue with :doc:`Binary Interface <binaryInterface>`.

MATLAB Interface
^^^^^^^^^^^^^^^^

The MATLAB-compatible interface is mainly intended for existing MATLAB or Octave
workflows. For new workflows, the Python interface is usually the better
choice.

Continue with :doc:`MATLAB Interface <MATLABInterface>`.

Typical Workflow
----------------

A typical HASEonGPU setup consists of the following steps:

#. Clone the repository
#. Install the required dependencies
#. Decide which interface to use
#. Follow installation steps provided on the page of the chosen interface
#. Verify the setup with one of the provided interface-specific examples
#. Use HASEonGPU in the target workflow

Next Step
---------

Proceed to one of the interface pages listed above, depending on the intended
usage path.
