# Changelog

## HASEonGPU 2.1.1

Patch release that moves pump time stepping and frozen-RK4 integration into the
compiled C++/Alpaka backend, exposes compiled ASE controls, consumes all pending
openPMD simulation output streams, and expands runtime-backend regression
coverage. Python-side fallback stepping and its retired compatibility surfaces
are removed. A plain ``python3 -m pip install -v .`` continues to build and
update the durable native runtime without preliminary CMake configuration,
while installing only the thin Python frontend. Direct CMake builds retain an
explicit option for installing the public C++ artifacts.

## HASEonGPU 2.1.0

Breaking release following 2.0.2 (released 2026-07-03) that makes openPMD the
only Python-to-backend transport.
The legacy pybind11 package, direct ``calcPhiASE(...)`` Python API, and retired
text-file compatibility paths were removed. The Python frontend now packages
the standalone backend and discovery libraries privately, selects a compatible
openPMD backend automatically (SST, ADIOS, then HDF5), and documents the
provider/runtime contract, migration path, and corrected transport diagnostics.
In MPI mode, the frontend again uses ``nPerNode`` to launch ``calcPhiASE``
through ``mpiexec`` automatically while exchanging data through openPMD.
The bundled openPMD-api provider is updated to 0.17.1; CI retains coverage for
both 0.17.0 and 0.17.1 system-provider installations.

## HASEonGPU 2.0.2

Patch release for the 2.0 series with a Julia 1D SSG regression check for
disabled-ASE laser-pump cladding output, synchronized release metadata, and
CodeRabbit review configuration.

## HASEonGPU 2.0.1

Patch release for the 2.0 series with pump-workflow fixes, Python API
refinements, restored VTK pump-rate output, and updated package and publication
metadata.

## HASEonGPU 2.0 ⸻ alpaka3

This release ports HASEonGPU to **alpaka 3**, introduces a more generic Python
interface, and adds regression testing and CI coverage for the updated workflow.

By adopting alpaka 3 as its backend abstraction, HASEonGPU becomes the
first scientific projects to use the new alpaka 3 library for high-performance,
platform-independent code.

### Highlights

- Ported the backend implementation to alpaka 3 with support for CPU serial,
  OpenMP, CUDA, and HIP backends.
- Improved CMake backend discovery, backend selection, and related
  documentation.
- Added a generic Python interface with reusable simulation setup helpers.
- Added geometry and mesh utilities for STL, gmsh, and VTK workflows.
- Updated Python examples and interface documentation.
- Added C++ and Python regression tests, including analytical comparison tests.
- Added GitHub Actions CI for selected compiler and CPU backend configurations.
