:doc:`<- Back to Getting Started <gettingStarted>`

Additional Dependencies
=======================

Python package installation provides only the Python-level dependencies listed
in ``pyproject.toml``. Source builds and tests may still need development tools
and runtime libraries from the operating system, especially in minimal
containers.

On Ubuntu-like systems, a CPU-only source install or test environment commonly
starts with:

.. code-block:: bash

   sudo apt-get install -y --no-install-recommends \
     build-essential cmake ninja-build pkg-config python3-dev \
     libgl1 libglu1 libomp-dev \
     libxcursor1 libxft2 libxinerama1 libxrender1

``libgl1`` provides ``libGL.so.1`` and ``libglu1`` provides ``libGLU.so.1``.
The ``libx*`` packages provide X11 runtime libraries loaded by the Python
``gmsh`` wheel, including ``libXrender.so.1``. These libraries are used by
Python packages such as ``gmsh`` and ``vtk`` but are not installed by pip.
Depending on the selected build options, add:

* ``libhwloc-dev`` for Alpaka host-memory support when it is enabled
* ``openmpi-bin libopenmpi-dev`` for MPI builds
* ``libtbb-dev`` when enabling TBB-backed host execution
* ``libhdf5-dev`` when using a system HDF5/openPMD provider
* vendor CUDA or HIP/ROCm toolkit and runtime packages for accelerator builds

Distribution package names vary; use the equivalent packages for your operating
system or cluster module environment.
