:doc:`<- Back to overview <index>`

Getting Started
===============

This page is a compact installation guide for a source checkout of HASEonGPU.
For modeling concepts, see :doc:`Theory and Model <theoryAndModel>`.  For the
main user workflow, continue with :doc:`Python Interface Guide <pythonInterface>`
after installation.

1. Clone the Repository
-----------------------

.. code-block:: bash

   git clone https://github.com/computationalradiationphysics/haseongpu.git
   cd haseongpu

2. Install Prerequisites
------------------------

Required tools are:

* ``Python >= 3.10`` with ``pip``
* ``cmake`` and ``ninja``
* a C++20 compiler, tested with ``gcc >= 12`` and ``clang >= 17``
* an openPMD-api provider for the storage backend you want to use, or the
  bundled provider selected by ``hase-configure``

Optional dependencies depend on the run mode:

* CUDA or HIP/ROCm for GPU builds
* OpenMPI for MPI runs
* ParaView for VTK visualization
* ``matplotlib`` for helper plotting scripts

Windows support is experimental; see :doc:`windows`.

3. Create a Python Environment
------------------------------

Use a virtual environment unless your site already provides a managed Python
module or Conda environment:

.. code-block:: bash

   python3 -m venv .venv
   source .venv/bin/activate
   python3 -m pip install -U pip

If you use Conda, Spack, or environment modules for openPMD-api, activate/load
that environment before configuring HASEonGPU so Python and CMake see the same
provider.

4. Run the Guided Configurator
------------------------------

From the source checkout, run:

.. code-block:: bash

   python3 utils/configure_hase.py

After HASEonGPU is installed, the same helper is available as:

.. code-block:: bash

   hase-configure

The configurator asks only for choices that affect installation or runtime
selection:

* openPMD provider: auto, bundled, or system
* ADIOS2/HDF5 handling for the selected provider
* runtime openPMD backend: automatic selection, or ``adios-sst``, ``adios``, or ``hdf5``
* Alpaka compute backend
* single-process or MPI mode
* native CPU optimizations for the local machine
* whether to run the printed install command immediately

The script writes a small PhiASE YAML run-control file to
``config/hase-phiase.yaml`` by default, prints the exact install command, and
finishes with guidance for the selected openPMD backend, MPI setting, and
available compute backends.  The generated YAML contains compute settings only;
physics inputs such as geometry, spectra, pump settings, and material state are
still constructed in Python.

Useful non-interactive options include ``--autoinstall``, ``--reinstall``,
``--use-ccache``, ``--provider``, ``--openpmd-backend``, ``--runtime-dir``, and
``--output``.  Run
``python3 utils/configure_hase.py --help`` for the complete list.

5. Install
----------

For the default configuration, installation remains one command:

.. code-block:: bash

   python3 -m pip install .

This configures or reuses the durable native runtime under ``build/``, builds
it, and installs a thin Python frontend that remembers that runtime.  Rebuilding
``build/`` later changes the binary and provider metadata used by the installed
frontend without reinstalling the wheel.

The configurator prints a command of this form when options are selected:

.. code-block:: bash

   CMAKE_ARGS="<selected CMake options>" python3 -m pip install -v .

Run the printed command if you did not let the configurator install
immediately. ``CMAKE_ARGS`` is how you pass build options such as the openPMD
provider, MPI mode, Alpaka choices, and native CPU optimization setting.

To keep the runtime somewhere other than ``build/``, let the configurator
select and remember it:

.. code-block:: bash

   python3 utils/configure_hase.py --runtime-dir /opt/hase/runtime-2.1 --autoinstall

Simulation scripts still use ``import HASEonGPU``; they do not need to set the
runtime path themselves.

If pip reports an externally managed Python environment, prefer a virtual
environment.  Use ``--break-system-packages`` with the configurator only when
you intentionally install into such an environment.

6. Verify and Continue
----------------------

Check that the package imports:

.. code-block:: bash

   python3 -c "import HASEonGPU; print(HASEonGPU.__version__)"

For the recommended user workflow, continue with
:doc:`Python Interface Guide <pythonInterface>`.  Use
:doc:`Binary Interface <binaryInterface>` only when running ``calcPhiASE``
directly, and :doc:`CMake Build Options <compilation>` when you need manual
CMake configuration.
