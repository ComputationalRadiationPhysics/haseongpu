HASEonGPU
=========

<p style="font-size: 1.5em; margin-top: -0.5em;">
 <strong>H</strong>igh performance <strong>A</strong>mplified <strong>S</strong>pontaneous <strong>E</strong>mission on <strong>GPU</strong>
</p>

<p align="left">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/source/_static/Logo_Darkmode.png">
    <source media="(prefers-color-scheme: light)" srcset="docs/source/_static/Logo_Lightmode.png">
    <img alt="HASEonGPU logo" src="docs/source/_static/Logo_Lightmode.png" width="400">
  </picture>
</p>


Description
-----------

HASEonGPU is an open-source HPC software for calculating amplified spontaneous emission (ASE) flux in laser gain media. It is intended to support the design and analysis of high-power laser systems by providing guidance on the expected ASE flux for a given gain-medium geometry, pumping configuration, and material setup. The code implements an adaptive Monte Carlo ray-tracing approach and was developed to accelerate ASE simulations using GPU hardware and distributed execution.

Documentation
----------------

The full project documentation is available at [https://haseongpu.readthedocs.io/en/latest/](https://haseongpu.readthedocs.io/en/latest/).

Referencing
-----------

HASEonGPU is a *scientific project*. If you present and/or publish scientific results that used HASEonGPU, you should cite the associated publication:

C.H.J. Eckert, E. Zenker, M. Bussmann, and D. Albach,  
*HASEonGPU—An adaptive, load-balanced MPI/GPU-code for calculating the amplified spontaneous emission in high power laser media*,  
Computer Physics Communications, **207** (2016), 362–374.  
DOI: `10.1016/j.cpc.2016.05.019`

Software License
----------------

*HASEonGPU* is licensed under the **GPLv3+**.

Active Team
-----------

### Scientific Supervision

- Dr. Michael Bussmann
- Dr. Daniel Albach

### Maintainers and core developers

- Erik Zenker
- Carlchristian Eckert
- Tim Hanel

### Participants, Former Members and Thanks

- Marius Melzer
- Frank Liebold
- Daniel Höhne