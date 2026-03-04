:doc:`<- Back to Getting Started <gettingStarted>`

GrayBat
=======

This page documents the (highly) experimental GrayBat support in HASEonGPU.

GrayBat is an optional communication layer for distributed execution and is
primarily relevant when using ``parallelMode="graybat"``.

Overview
--------

GrayBat is intended as an alternative communication framework for distributed
parallel execution. In HASEonGPU, it is used similarly to MPI-based execution,
but provides a more flexible communication topology.

GrayBat support in HASEonGPU is currently experimental.

Build Support
-------------

To enable GrayBat support during compilation, configure HASEonGPU with:

.. code-block:: bash

   cmake .. -DUSE_GRAYBAT=ON

Usage
-----

When GrayBat support is enabled, it can be selected through the interface-level
``parallelMode`` parameter or the corresponding binary command-line option.

For example:

* Python interface: ``parallelMode="graybat"``
* MATLAB interface: ``parallelMode="graybat"``
* Binary interface: ``--parallel-mode=graybat``

Notes
-----

GrayBat is not required for normal HASEonGPU usage.

If you do not specifically need GrayBat, the standard threaded or MPI execution
modes are usually the simpler choice.