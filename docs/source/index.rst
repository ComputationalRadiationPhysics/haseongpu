HASEonGPU Documentation
=======================

.. image:: _static/Logo_Lightmode.png
   :alt: HASEonGPU logo
   :width: 400px
   :align: center

HASEonGPU (**H**\ igh performance **A**\ mplified **S**\ pontaneous **E**\ mission on **GPU**) is an
open-source HPC software for calculating amplified spontaneous emission (ASE)
flux in laser gain media.

The Python frontend uses openPMD to communicate with a standalone C++ backend;
legacy in-process Python bindings are not supported. A Python source install
keeps that backend in a resident CMake runtime directory and installs a thin
frontend that refers to it.

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
