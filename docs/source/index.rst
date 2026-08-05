HASEonGPU Documentation
=======================

.. image:: _static/Logo_Lightmode.png
   :alt: HASEonGPU logo
   :width: 400px
   :align: center

HASEonGPU (**H**\ igh performance **A**\ mplified **S**\ pontaneous **E**\ mission on **GPU**) is an
open-source HPC software for calculating amplified spontaneous emission (ASE)
flux in laser gain media.

It is intended to support the design and analysis of high-power laser systems,
where ASE is an important limiting effect for stored energy, gain distribution,
and overall amplifier performance.

Start with :doc:`Getting Started <gettingStarted>` to install HASEonGPU, then
follow the :doc:`laserPumpCladding tutorial <laserPumpCladding>` for a complete
Tet4 pump-and-ASE simulation. The remaining pages separate reusable modeling
concepts, execution and transport details, and reference material.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Start here

   gettingStarted
   laserPumpCladding

.. toctree::
   :maxdepth: 2
   :caption: Modeling

   pythonInterface
   theoryAndModel

.. toctree::
   :maxdepth: 2
   :caption: Execution and data

   backendSelection
   openpmdTransport
   mpi
   binaryInterface
   scripts

.. toctree::
   :maxdepth: 2
   :caption: Build and platforms

   compilation
   additional_deps
   windows

.. toctree::
   :maxdepth: 2
   :caption: Reference

   pythonAPI
