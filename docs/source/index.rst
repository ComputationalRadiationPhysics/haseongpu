.. HASEonGPU documentation master file, created by
   sphinx-quickstart on Wed Mar  4 13:00:19 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HASEonGPU Documentation
=======================

HASEonGPU (High performance Amplified Spontaneous Emission on GPU) is a
high-performance simulation framework for calculating **amplified spontaneous
emission (ASE)** in laser gain media.

The software was developed at the Helmholtz-Zentrum Dresden-Rossendorf (HZDR)
and is designed to accelerate Monte-Carlo based simulations of ASE using modern
GPU hardware and distributed computing techniques.

Overview
--------

Amplified spontaneous emission (ASE) occurs when spontaneous emission in a
laser gain medium is amplified through stimulated emission. In high-power laser
systems, ASE can significantly reduce the achievable energy storage in the gain
medium and therefore must be carefully analyzed during laser design.

HASEonGPU presents a scalable Monte-Carlo simulation to study these effects in
large gain media by combining, which is suitable for running in HPC-Environments.

:contentReference[oaicite:0]{index=0}

Features
--------

* Monte-Carlo simulation of amplified spontaneous emission
* GPU acceleration for large scale simulations
* Distributed multi-GPU execution using MPI
* Adaptive load balancing
* Interfaces for MATLAB and C-based applications
* Visualization support through VTK output (e.g. with ParaView)

Scientific Background
---------------------

ASE is generated when a population-inverted gain medium amplifies spontaneous
emission photons. In high-gain or large-volume laser systems this process can
deplete stored energy and degrade amplifier performance. :contentReference[oaicite:1]{index=1}

Numerical modeling of ASE is therefore essential for the design and optimization
of high-power laser systems, including diode-pumped solid-state lasers and
large slab amplifiers.

HASEonGPU implements an adaptive Monte-Carlo algorithm for tracking photon
trajectories and computing ASE flux in two-dimensional gain media.
The geometry is extruded along the z direction and represented as a stack of layers.


Reference
---------

If you use HASEonGPU for scientific work, please cite:

Eckert, C. H. J., Zenker, E., Bussmann, M., and Albach, D.
"HASEonGPU — An adaptive, load-balanced MPI/GPU code for calculating
the amplified spontaneous emission in high power laser media."
Computer Physics Communications (2016).

Contents
--------

.. toctree::
   :maxdepth: 2

   installation
   usage
   theory
   examples

.. toctree::
   :maxdepth: 2
   :caption: Contents:

