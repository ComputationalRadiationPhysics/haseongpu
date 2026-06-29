# Changelog

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
