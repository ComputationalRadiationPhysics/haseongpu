:doc:`<- Back to Getting Started <gettingStarted>`

Windows Notes
=============

This page contains Windows-specific notes for HASEonGPU.

Windows support is experimental.

Overview
--------

HASEonGPU is primarily developed and documented for Linux environments.
Windows builds may require additional manual setup and are not the main target
platform of the project.

Build Notes
-----------

A Windows build typically requires:

* Visual Studio
* CMake
* CUDA
* Boost
* optionally an MPI implementation

Depending on the local setup, it may be necessary to configure Boost-related
paths manually in the CMake configuration.

Notes
-----

If possible, Linux is the recommended platform for building and using
HASEonGPU.

For standard setup instructions, see :doc:`Getting Started <gettingStarted>`.