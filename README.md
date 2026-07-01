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

HASEonGPU uses openPMD as the Python-to-C++ transport. Provide an openPMD-api
installation that supplies both the Python `openpmd_api` module and the CMake
`openPMD::openPMD` package:

```bash
conda install -c conda-forge openpmd-api
python3 utils/check_openpmd_compatibility.py --backend adios-sst --cmake-prefix-path "$CONDA_PREFIX"
CMAKE_ARGS="-DCMAKE_PREFIX_PATH=$CONDA_PREFIX" python3 -m pip install .
```

The default runtime openPMD backend is `adios-sst`. Select `adios` or `hdf5`
from Python/YAML only when the Python and C++ openPMD providers support that
backend. Spack, system modules, or manual external openPMD installs can also be
used by passing `CMAKE_PREFIX_PATH` or `openPMD_DIR`.

If no matching external openPMD installation is available, HASEonGPU can use the
FetchContent source-build path:

```bash
CMAKE_ARGS="-DHASE_BUILD_OPENPMD_FROM_SOURCE=ON" python3 -m pip install .
```

That path builds openPMD-api for the HASE CMake build, but the HASE wheel does
not vendor openPMD runtime libraries or generated `openpmd_api` bindings.

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
