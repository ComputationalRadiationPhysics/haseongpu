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

Release Scope
-------------

HASEonGPU 2.x releases ship the Python frontend and the compiled C++ backend as one supported release unit. The Python package is the primary user-facing interface, while the C++ backend, command-line executable, headers, and CMake package export are versioned with the same release number for downstream builds and reproducible integration.

Documentation
----------------

The full project documentation is available at [https://haseongpu.readthedocs.io/en/latest/](https://haseongpu.readthedocs.io/en/latest/).

Installation Notes
------------------

HASEonGPU uses openPMD as the Python-to-C++ transport. By default, a source
install expects an external openPMD-api C++ package discoverable by CMake
(`openPMD::openPMD`) and a matching Python `openpmd_api` package in the active
Python environment:

```bash
python3 -m pip install openpmd-api
CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/path/to/openpmd/prefix" python3 -m pip install .
```

The C++ and Python openPMD providers must be built with the same required
backends, for example ADIOS2 for `adios`/`adios-sst` or HDF5 for `hdf5`. See the
[openPMD-api installation guide](https://openpmd-api.readthedocs.io/en/latest/install/install.html)
for Spack, Conda, pip, and source-install options. For ADIOS2-capable installs,
use provider options such as `spack install openpmd-api +python +adios2` or
`openPMD_CMAKE_openPMD_USE_ADIOS2=ON python3 -m pip install --no-binary openpmd-api openpmd-api`.

If no matching external openPMD installation is available, HASEonGPU can use the
legacy self-contained source-build path:

```bash
CMAKE_ARGS="-DHASE_BUILD_OPENPMD_FROM_SOURCE=ON" python3 -m pip install .
```

Citation
--------

If you use **HASEonGPU** in research, please cite it using the metadata in [CITATION.cff](CITATION.cff).

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
- René Widera

### Maintainers and core developers

- Erik Zenker
- Carlchristian Eckert
- Tim Hanel

### Participants, Former Members and Thanks

- Marius Melzer
- Frank Liebold
- Daniel Höhne
