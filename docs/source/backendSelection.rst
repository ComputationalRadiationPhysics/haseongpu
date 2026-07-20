Backend Selection
=================

HASEonGPU has two independent backend choices:

* **Alpaka compute backend**: where kernels run, for example a CPU, CUDA, or HIP
  backend.
* **openPMD storage backend**: how Python and ``calcPhiASE`` exchange data.
  The default ``auto`` selects ``adios-sst``, then ``adios``, then ``hdf5``
  from the backends supported by both runtime providers.

Do not use an openPMD backend name as the Alpaka ``backend`` value.

Build-Time Compute Backends
---------------------------

Only Alpaka backends enabled at build time can be selected at runtime.  By
default, HASEonGPU tries to detect supported dependencies such as CUDA, HIP,
and TBB automatically.

For manual selection, configure with ``HASE_SELECT_BACKEND_ALPAKA=ON`` and pass
Alpaka's CMake options directly:

.. code-block:: bash

   cmake -S . -B build \
     -DHASE_SELECT_BACKEND_ALPAKA=ON \
     -Dalpaka_DEP_CUDA=ON \
     -Dalpaka_DEP_HIP=OFF \
     -Dalpaka_DEP_TBB=OFF \
     -Dalpaka_EXEC_CpuSerial=OFF

The exact Alpaka option names are defined by Alpaka; see the
`Alpaka CMake argument documentation <https://alpaka3.readthedocs.io/en/latest/advanced/cmake.html#arguments>`__.
HASEonGPU is mainly validated with host, CUDA, and HIP backends.  Other Alpaka
backends should be treated as experimental unless validated in your environment.

Runtime Compute Backend
-----------------------

The Python ``PhiASE.backend`` value, or the equivalent openPMD metadata used by
``calcPhiASE``, selects one compiled-and-available compute backend.  Backend
names have the form:

.. code-block:: text

   api_deviceKind_executor

Use ``hase-configure`` for an interactive list of backends available in the
installed build.  It also keeps this compute choice separate from the openPMD
storage backend.

From Python, query the same list with ``AlpakaBackends``:

.. code-block:: python

   from HASEonGPU import AlpakaBackends, PhiASE

   available = AlpakaBackends.all()
   phi_ase = PhiASE(backend=available[0])

``AlpakaBackends.known()`` is an alias for ``all()``.  Names that are valid
Python identifiers are also exposed as attributes, for example
``AlpakaBackends.Host_Cpu_CpuSerial``.

Backend-Name Helper Library
---------------------------

The Python query uses the small CMake-built helper library
``HaseAlpakaBackendNames``.  If importing ``AlpakaBackends`` reports that this
library is missing, rebuild/install the Python package or build the target
explicitly:

.. code-block:: bash

   cmake --build build --target HaseAlpakaBackendNames

For CMake option details, see :doc:`CMake Build Options <compilation>`.
