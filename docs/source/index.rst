HASEonGPU Documentation
=======================

.. image:: _static/Logo_Lightmode.png
   :alt: HASEonGPU logo
   :width: 400px
   :align: center

HASEonGPU (**H**\ igh performance **A**\ mplified **S**\ pontaneous **E**\ mission on **GPU**) is an
open-source HPC software for calculating amplified spontaneous emission (ASE)
flux in laser gain media.

HASEonGPU 2.1 ships the Python frontend and compiled C++ backend together. The
frontend uses openPMD to communicate with the standalone backend; legacy
in-process Python bindings are not supported. The C++ backend, headers,
command-line executable, and CMake package export use the same release version
for downstream integration.

It is intended to support the design and analysis of high-power laser systems,
where ASE is an important limiting effect for stored energy, gain distribution,
and overall amplifier performance.

Start with :doc:`Getting Started <gettingStarted>` for installation and
interface selection.  Use :doc:`CMake Build Options <compilation>` only when
you need manual build configuration.  The :doc:`Theory and Model
<theoryAndModel>` page explains the ASE model and pump coupling.

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
