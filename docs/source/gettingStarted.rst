:doc:`<- Back to overview <index>`

Getting Started
===============

This page provides a short guide to setting up HASEonGPU and choosing the
interface that best fits your workflow.

HASEonGPU supports multiple usage paths:

* :doc:`Python Interface <pythonInterface>` for direct usage from Python
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
-------------------

Required software:

* ``cmake >= 3.0.1``
* ``gcc >= 11``
* ``cuda >= 11.0``

Optional software and tools:

* ``OpenMPI >= 4.0``
* ``Python >= 3.10``
* MATLAB or Octave
* ``ParaView`` for visualization of ``.vtk`` output

Hardware requirements:

* NVIDIA GPU with CUDA support
* For GPU execution, a CUDA-capable device is required

Additional Notes
----------------

For experimental GrayBat support, see :doc:`graybat`.

For Windows-specific installation notes [deprecated], see :doc:`windows`.

Compilation Notes
-----------------

A manual compilation step is available, but is not required for every workflow.

For example, when using the Python interface, the C++ backend is built under
the hood during installation. For details on manual compilation, see
:doc:`compilation`.

Choose an Interface
-------------------

Python Interface
^^^^^^^^^^^^^^^^

The Python interface is the recommended starting point for most new users.
It provides a library call that can be easily integrated into custom workflows.
Additionally, the provided examples demonstrate this usage.


Continue with :doc:`Python Interface <pythonInterface>`.

Binary Interface
^^^^^^^^^^^^^^^^

The binary interface is useful if you want to run HASEonGPU directly from the
command line, or use it as an entry point for constructing a custom API call.

Continue with :doc:`Binary Interface <binaryInterface>`.

MATLAB Interface
^^^^^^^^^^^^^^^^

The MATLAB-compatible interface is mainly intended for existing MATLAB or Octave
workflows. If you are starting a new workflow, the Python interface is usually
the better choice.

Continue with :doc:`MATLAB Interface <MATLABInterface>`.

Typical Workflow
----------------

A typical HASEonGPU setup consists of the following steps:

#. Clone the repository
#. Install the required dependencies
#. Decide which interface to use
#. Follow installation steps provided on the page of the chosen interface
#. In order to verify your current setup run one of the provided interface-specific examples
#. Use HASEonGPU in your own workflow

Next Step
---------

Proceed to one of the interface pages listed above, depending on how you want
to use HASEonGPU.