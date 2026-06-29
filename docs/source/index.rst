HASEonGPU Documentation
=======================

.. image:: _static/Logo_Lightmode.png
   :alt: HASEonGPU logo
   :width: 400px
   :align: center

HASEonGPU (**H**\ igh performance **A**\ mplified **S**\ pontaneous **E**\ mission on **GPU**) is an
open-source HPC software for calculating amplified spontaneous emission (ASE)
flux in laser gain media.

HASEonGPU 2.x releases ship the Python frontend and compiled C++ backend together. The
Python package is the primary user-facing interface, and the C++ backend,
headers, command-line executable, and CMake package export use the same release
version for downstream integration.

It is intended to support the design and analysis of high-power laser systems,
where ASE is an important limiting effect for stored energy, gain distribution,
and overall amplifier performance.

Start with :doc:`Getting Started <gettingStarted>` for setup and interface
selection.  The :doc:`Theory and Model <theoryAndModel>` page explains
the ASE model, estimator, pump coupling, and implementation restrictions.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   gettingStarted
   theoryAndModel
   backendSelection
   compilation
   scripts
   windows
   binaryInterface
   openpmdTransport
   mpi
   pythonInterface
   pythonAPI
   pythonInterfaceLegacy
