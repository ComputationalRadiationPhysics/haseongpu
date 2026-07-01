Python Interface Guide
======================

The Python interface is the recommended way to build HASEonGPU simulations.  It
keeps the physical setup in Python objects and uses the compiled C++ backend for
the ASE calculation.

The usual workflow is:

#. build a :doc:`topology <python_interface/topology>`
#. attach material/state arrays with :doc:`GainMedium <python_interface/gain_medium>`
#. provide :doc:`spectra <python_interface/spectral_decomposition>`
#. define :doc:`pump properties <python_interface/pump_properties>`
#. configure :doc:`PhiASE <python_interface/phi_ase>`
#. run a :doc:`Simulation <python_interface/simulation>`

For generated signatures and member lists, use the
:doc:`Python API Reference <pythonAPI>`.

Installation
------------

Install from the repository root after following
:doc:`Getting Started <gettingStarted>`:

.. code-block:: bash

   python3 utils/configure_hase.py
   # run the printed command, typically:
   CMAKE_ARGS="<selected CMake options>" python3 -m pip install -v .

``pip install -v .`` configures or reuses the native runtime in ``build/`` and
installs a thin Python frontend that remembers it. The configurator chooses
compatible openPMD provider settings, prints the install command, and writes the optional
``config/hase-phiase.yaml`` compute-settings file.

If you change CMake options or native sources, reconfigure or rebuild the
durable runtime. The installed frontend uses the updated binary and provider
metadata directly, so it does not need to be reinstalled:

.. code-block:: bash

   cmake -S . -B build <updated options>
   cmake --build build

Concept Pages
-------------

.. toctree::
   :maxdepth: 2
   :caption: Python interface guide

   python_interface/topology
   python_interface/gain_medium
   python_interface/spectral_decomposition
   python_interface/pump_properties
   python_interface/phi_ase
   python_interface/simulation
   python_interface/utilities

Minimal Example Tutorial
------------------------

The snippets below come from ``example/minimalExampleNewInterface.py`` and stay
in sync with that runnable example.

Geometry
^^^^^^^^

Describe the transverse mesh and z-levels.  HASEonGPU currently uses a 2D
triangular mesh extruded into prism layers.

.. literalinclude:: ../../example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: topology
   :end-before: # docs:end: topology

``Grid`` is the shortest path for rectangular media.  ``MeshTopology`` can also
be created from point clouds, planar STL files, legacy VTK wedge files, and
gmsh triangle meshes; see :doc:`python_interface/topology`.

Material and State
^^^^^^^^^^^^^^^^^^

Attach physical arrays and scalar material properties to the topology:

.. literalinclude:: ../../example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: gain-medium
   :end-before: # docs:end: gain-medium

Use ``medium.get("...").expectedShape`` when allocating mesh-dependent arrays.
The most important fields are ``betaCells`` for point-level excited-state
fraction, ``betaVolume`` for prism-centered beta, cladding labels,
reflectivities, active-ion density ``nTot``, and fluorescence lifetime
``crystalTFluo``.  See :doc:`python_interface/gain_medium` for field shapes and
openPMD metadata.

Spectra
^^^^^^^

Provide absorption and emission cross sections for pump and ASE calculations:

.. literalinclude:: ../../example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: spectral-decomposition
   :end-before: # docs:end: spectral-decomposition

Each wavelength array must match the length of its cross-section array.  The
``resolution`` value is passed to the ASE backend as spectral interpolation
resolution.

Pump
^^^^

``PumpProperties`` describes the physical pump sent to the compiled simulation:
intensity, wavelength, cross sections, transverse radii, super-Gaussian
exponent, and optional reflected pass. The backend implements the pump; a
custom Python pump solver is not supported.

.. literalinclude:: ../../example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: pump-properties
   :end-before: # docs:end: pump-properties

See :doc:`python_interface/pump_properties` for parameters and the model.

PhiASE
^^^^^^

``PhiASE`` configures the ASE backend: Monte Carlo ray limits, adaptive
sampling, reflections, compute backend, openPMD backend, and parallel mode.

.. literalinclude:: ../../example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: phi-ase
   :end-before: # docs:end: phi-ase

Use a backend name reported by the installed build:

.. code-block:: python

   from HASEonGPU import AlpakaBackends

   backend = AlpakaBackends.all()[0]

``backend`` is the Alpaka compute backend.  ``openpmdBackend`` is the openPMD
storage/streaming backend such as ``adios-sst``.  See
:doc:`Backend Selection <backendSelection>` and
:doc:`openPMD Transport <openpmdTransport>`.

Simulation
^^^^^^^^^^

``Simulation`` sends the setup to the C++/Alpaka backend, which owns the pump,
ASE, and time-integration loop. Python receives one ``TimeStepState`` snapshot
per completed step and invokes ``onStep`` callbacks for each snapshot.

.. literalinclude:: ../../example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: simulation
   :end-before: # docs:end: simulation

Use ``runSteps(150, pumpSteps=50)`` to pump only the first 50 outer steps, or
``runUntil(endTime=1e-3)`` for a time target. ``pumpSubsteps`` instead controls
the internal resolution of one pumped step. ``onInit`` is available before the
backend launch; per-step Python mutation through ``beforeStep`` is unsupported.

YAML Compute Settings
^^^^^^^^^^^^^^^^^^^^^

``PhiASE.fromYaml(...)`` can load run-control settings while keeping geometry,
material arrays, spectra, and pump setup in Python:

.. code-block:: yaml

   phiASE:
     min_rays_per_sample: 100000
     max_rays_per_sample: 1000000
     mse_threshold: 0.05
     repetitions: 2
     adaptive_steps: 4
     use_reflections: true

   compute:
     backend: Host_Cpu_CpuSerial
     openpmd_backend: auto
     parallel_mode: single
     numDevices: 1

.. code-block:: python

   phi_ase = PhiASE.fromYaml("config/hase-phiase.yaml", spectralProperties=spectra)

Constructor keyword arguments override YAML values.  ``hase-configure`` writes a
small YAML file with only the compute settings.

Results
^^^^^^^

.. literalinclude:: ../../example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: results
   :end-before: # docs:end: results

``simulation.getLastState()`` returns the latest ``TimeStepState`` with step,
time, ``betaCells``, ``betaVolume``, ``phiAse``, pump derivative, ASE
derivative, and the raw ASE result object.  Use callbacks to store or export
every step.
