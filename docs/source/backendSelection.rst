Backend Selection
=================

HASEonGPU uses alpaka to support different hardware backends.

Backend selection happens in two steps:

* build-time backend selection
* runtime backend selection

The build-time selection defines which alpaka backends are compiled into
HASEonGPU. The runtime selection chooses one of these compiled and available
backends for a specific simulation run.

HASEonGPU also builds a small backend-name library,
``HaseAlpakaBackendNames``, which exposes the backend names detected by the
compiled alpaka configuration.  The Python interface uses this library to ask
the installed build which backend names are valid on the current machine.


Build-time backend selection
----------------------------

The backends available at runtime depend on which alpaka backends were enabled
when configuring HASEonGPU with CMake.

By default, HASEonGPU tries to detect supported backend dependencies such as
CUDA, HIP, and TBB automatically.
In the case of manual compilation: discrete/manual backend selection can be enabled with
the CMake option ``HASE_SELECT_BACKEND_ALPAKA``.

When manual selection is enabled, the existing alpaka CMake options can be used
directly to select APIs, device kinds, and executors.

The relevant alpaka CMake options are documented in the
`alpaka CMake argument documentation <https://alpaka3.readthedocs.io/en/latest/advanced/cmake.html#arguments>`__.

For example, a manual CUDA-focused backend configuration can use:

.. code-block:: bash

   cmake -S . -B build \
     -DHASE_SELECT_BACKEND_ALPAKA=ON \
     -Dalpaka_DEP_CUDA=ON \
     -Dalpaka_DEP_HIP=OFF \
     -Dalpaka_DEP_TBB=OFF \
     -Dalpaka_EXEC_CpuSerial=OFF

For more information about manually compiling HASEonGPU, see
:doc:`Compilation <compilation>`.

The exact set of available alpaka CMake options may change between alpaka
versions. Therefore, the alpaka documentation should be treated as the primary
reference for the supported backend configuration flags.
It is possible to compile HASEonGPU with all backends provided by alpaka,
including SYCL CPU/GPU backends. However, HASEonGPU is currently tested mainly
with CUDA, HIP and host backends. Other backends should therefore be treated
as experimental.

Backend-name helper library
---------------------------

During the CMake build, HASEonGPU also builds the shared library target
``HaseAlpakaBackendNames``.  This library links against the same HASEonGPU core
and alpaka configuration as the simulation code.  It enumerates alpaka's
enabled APIs and executors, filters them to backends for which a device can be
created, and exposes the resulting names through a small C ABI.

For Python builds, CMake copies this library into the build package next to the
Python extension:

.. code-block:: text

   build/python/HASEonGPU_Bindings/libHaseAlpakaBackendNames.so

On macOS the filename ends in ``.dylib``.  On Windows it is
``HaseAlpakaBackendNames.dll``.

This library is useful because backend availability is not only a static CMake
setting.  A backend must be compiled in and usable on the machine where the
code runs.  For example, a CUDA backend name should only be reported when that
backend was enabled in the build and alpaka can create a CUDA device.

If the Python package cannot find this helper library, build the target
explicitly:

.. code-block:: bash

   cmake --build build --target HaseAlpakaBackendNames

When the Python interface is built, the regular Python extension target also
copies the helper library into the Python runtime package.


Runtime backend selection
-------------------------

At runtime, the ``--backend=`` command-line option, or the equivalent
``backend`` option in the MATLAB and Python interfaces, selects which of the
compiled backends should be used for the simulation.

If the command-line binary is given a backend string that does not match any
compiled backend with an available device, it prints the list of currently
available runtime backends before exiting.  Use one of those names as the
``--backend=`` value.

Backend names are constructed from alpaka's ``api``, ``deviceKind``, and
``executor``:

.. code-block:: text

   api_deviceKind_executor

Only backends that were enabled at build time can be selected at runtime.

Query available backends from Python
------------------------------------

The Python front end exposes the CMake-built backend-name library through
``AlpakaBackends``:

.. code-block:: python

   from HASEonGPU import AlpakaBackends

   print(AlpakaBackends.all())

``AlpakaBackends.known()`` is an alias for ``AlpakaBackends.all()``:

.. code-block:: python

   available = AlpakaBackends.known()
   backend = available[0]

If a backend name is a valid Python identifier, it is also available as a class
attribute:

.. code-block:: python

   backend = AlpakaBackends.Host_Cpu_CpuSerial

The returned strings are the values to pass to the new Python interface:

.. code-block:: python

   from HASEonGPU import PhiASE

   phi_ase = PhiASE(
       spectralProperties=spectra,
       backend=AlpakaBackends.Host_Cpu_CpuSerial,
       parallelMode="single",
   )

or to the command-line binary:

.. code-block:: bash

   ./build/calcPhiASE --backend=Host_Cpu_CpuSerial ...

If importing ``AlpakaBackends`` raises an error about the backend-name library,
the Python package can see the Python files but not the compiled helper
library.  Build ``HaseAlpakaBackendNames`` and make sure the generated library
is either installed with the Python package or present under
``build/python/HASEonGPU_Bindings``.

The available APIs, device kinds, and executors are documented in the
`alpaka cheatsheet <https://alpaka3.readthedocs.io/en/latest/basic/cheatsheet.html#accelerator-platform-and-device>`__.

For a more detailed explanation of alpaka devices and backend selection, see
the `alpaka device documentation <https://alpaka3.readthedocs.io/en/latest/tutorial/device.html>`__.
