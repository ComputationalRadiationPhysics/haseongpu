Python Interface Guide
======================

This guide is the recommended starting point for new HASEonGPU simulations in
Python.  It explains the workflow and concepts: mesh construction, material
properties, pump propagation, ASE execution, and time stepping with the high-level
Python objects.

For generated signatures, class members, and direct object lookup, use the
:doc:`Python API Reference <pythonAPI>`.  The older low-level
``calcPhiASE(...)`` interface is still supported and is documented separately in
:doc:`pythonInterfaceLegacy`.

The current interface allows users to describe the physical problem using
separate high-level abstraction objects. These objects define the
``MeshTopology`` and the physical properties of the ``GainMedium``.

Absorption and emission cross sections are described through a
``SpectralDecomposition``.

The pump setup (``PumpProperties``) and the settings for the C++ ``PhiASE``
computation can be configured independently.

Together, these configuration objects can be used to run a time-stepped
physical ``Simulation`` of the laser crystal. During the simulation, the
evolution of gain propagation and population inversion can be monitored.

Installation
------------

Install the haseongpu Python package from the repository root:

.. code-block:: bash

   python3 -m pip install -e .

This installs the Python package in editable mode and also installs the Python
dependencies declared in ``pyproject.toml``.

Source-Tree Import
------------------

If the wheel installation above does not work, you can use the source tree directly
as a temporary workaround.

Place or run your script inside the repository root. From there, Python should be
able to import ``HASEonGPU`` directly without installing the package.

For this workaround the dependencies have to be installed manually:

.. code-block:: bash

   python3 -m pip install -r requirements.txt


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
   pythonInterfaceLegacy

Minimal Example Tutorial
------------------------

The Python interface guide is built around the objects that appear in a physical ASE
simulation.  First describe the crystal geometry.  Then attach material and
state data to that geometry.  Then describe the spectra and pump.  Finally,
configure the ASE solver and let ``Simulation`` advance the system in time.

The code snippets in this section are taken from the minimal new-interface
example so that the tutorial code and the runnable example stay in sync.  The
snippets are included by named code markers instead of fixed line numbers, so
the documentation follows the example when code is moved inside the file.

Describe the Crystal Geometry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to simulate laser pump and ASE behavior, HASEonGPU first needs a
geometry for the laser crystal.  The new interface calls this geometry a
:doc:`python_interface/topology`: it contains the transverse mesh, the triangle
connectivity, and the z-levels that are used to extrude the 2D mesh into prism
cells.

.. literalinclude:: ../../example/python_example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: topology
   :end-before: # docs:end: topology

Here ``Grid`` describes a rectangular crystal.  ``xExtent`` and ``yExtent`` are
the transverse size of the mesh.  ``zExtent`` is the crystal length in the
propagation direction.  ``tileSizeX`` and ``tileSizeY`` control the transverse
mesh spacing, while ``tileSizeZ`` controls the spacing between z-levels.

``MeshTopology.fromGrid(...)`` triangulates that rectangular grid and keeps the
z-level information needed later by the gain medium. Since HASEonGPU currently does not support 3 dimensional input data,
but rather infers three dimensionality by extruding a 2 dimensional mesh given the number of layers and the physical-distance
called "thickness" between the layers. The same topology object
can also be created from 2 dimensional point clouds, planar STL files, or gmsh triangle
meshes; see :doc:`python_interface/topology`.

Attach Material and State Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A topology only describes the dimensions and geometrical properties of the crystal. The :doc:`python_interface/gain_medium`
adds what the crystal is made of and what state it is currently in.  This is
where the population inversion, cladding, surface data, doping density, and
fluorescence lifetime are assigned.

.. literalinclude:: ../../example/python_example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: gain-medium
   :end-before: # docs:end: gain-medium

The shape of most arrays depends on the mesh, so the example asks the medium
for the expected shape before allocating data.  For example,
``medium.get("betaCells").expectedShape`` returns
``(numberOfPoints, numberOfLevels)``.  ``reflectivities`` returns
``(2, numberOfTriangles)`` because there is one value for the bottom surface
and one for the top surface of every triangle.
In general the data-layout is still in-line with the pythonInterfaceLegacy in :doc:`pythonInterfaceLegacy`.

The properties in this minimal setup are:

``betaCells``
   Excited-state fraction :math:`\beta_i` at mesh points and z-levels.  The
   time integration updates this array.

``claddingCellTypes``
   Triangle-wise cladding labels.  In an extruded 2D mesh this is how side
   regions can be assigned to cladding groups.

``refractiveIndices``
   Refractive indices for the material transitions at the lower and upper
   crystal surfaces.  The layout is
   ``[bottomInside, bottomOutside, topInside, topOutside]`` and is used when
   reflections are enabled.

``reflectivities``
   Surface reflectivity per triangle.  Row 0 describes the bottom surface and
   row 1 describes the top surface.

``nTot``
   Active-ion concentration :math:`N_{\mathrm{tot}}` in the gain medium.  Pump
   absorption, emission, and ASE depletion are scaled by this density.

``crystalTFluo``
   Fluorescence lifetime :math:`\tau`.  The time loop uses it for spontaneous
   decay of the excited population.

``claddingNumber`` and ``claddingAbsorption``
   Select which triangle label is treated as cladding and how strongly that
   cladding absorbs.

See :doc:`python_interface/gain_medium` and the low-level argument reference in
:doc:`pythonInterfaceLegacy` for the exact layouts passed to HASEonGPU.

Provide Absorption and Emission Spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The pump and ASE calculation need material cross sections.  The pump uses the
absorption cross section :math:`\sigma_a` at the pump wavelength
:math:`\lambda` to update the excited population.  The ASE calculation uses
the absorption and emission spectra, :math:`\sigma_a(\lambda)` and
:math:`\sigma_e(\lambda)`, to compute wavelength-dependent amplification and
loss. This enables mult-chromatic ASE calculation.

.. literalinclude:: ../../example/python_example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: spectral-decomposition
   :end-before: # docs:end: spectral-decomposition

``SpectralDecomposition`` stores these tables in one object.  Each wavelength
array must have the same length as its matching cross-section array:
``wavelengthsAbsorption`` with ``crossSectionAbsorption`` and
``wavelengthsEmission`` with ``crossSectionEmission``.  The ``resolution`` value
is passed to the ASE backend as the spectral interpolation resolution.

Describe the Pump
^^^^^^^^^^^^^^^^^

The pump defines how energy is deposited into ``betaCells``
(:math:`\beta_i`) during one time step.  In the minimal setup the pump is a
super-Gaussian beam with wavelength :math:`\lambda`, intensity :math:`I`,
radii, and exponent.  It uses the same spectral data object created above.

.. literalinclude:: ../../example/python_example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: pump-properties
   :end-before: # docs:end: pump-properties

``PumpProperties`` is intentionally flexible.  Standard pump parameters such as
``intensity``, ``wavelength``, ``radiusX``, ``radiusY``, and ``exponent`` are
read by the default pump routine.  Additional keyword arguments are stored as
custom properties and can be read by user-defined pump solvers.

The default solver is ``BetaIntegrationGaussianSolver``.  It evaluates the
super-Gaussian intensity profile :math:`I(x, y)` at each transverse topology
point, propagates the pump through the z-levels (one-dimensional pump), optionally adds back
reflection, and returns the beta distribution :math:`\beta` after pumping.

The example passes a custom solver to show the extension point:

.. literalinclude:: ../../example/python_example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: custom-pump-solver
   :end-before: # docs:end: custom-pump-solver

A custom pump solver only needs a ``step(input, pump)`` method.  It receives
the current beta array :math:`\beta` as ``input["betaCell"]`` and returns the
updated beta array with the same shape.  ``pump.getProperty(...)`` is how
custom solver parameters are read.  More realistic pump construction is covered in
:doc:`python_interface/pump_properties`.

Configure PhiASE
^^^^^^^^^^^^^^^^

At this point the geometry, material state, spectra, and pump are known.  The
remaining question is how the ASE calculation should be executed.
To answer that question, the interface provides a ``PhiASE`` object, which is a configuration object for the HASEonGPU ASE
c++ backend, which sets sampling limits, adaptive convergence settings, reflection handling,
backend name, and parallel execution mode.

.. literalinclude:: ../../example/python_example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: phi-ase
   :end-before: # docs:end: phi-ase

``minRaysPerSample`` and ``maxRaysPerSample`` bound the Monte Carlo ray count
per sample, the :math:`N` used in the estimator in the scientific background.
``mseThreshold``, ``repetitions``, and ``adaptiveSteps`` control how adaptive
sampling can increase the accuracy of the Monte Carlo integration.
``useReflections`` enables the surface reflection model that
uses ``reflectivities`` from the gain medium.

``backend``, ``parallelMode``, and ``numDevices`` describe how HASEonGPU should
run.  When ``phi_ase.run(...)`` is called, ``PhiASE`` builds the low-level host
mesh and compute parameter objects and forwards them to the compiled
HASEonGPU implementation.  ``PhiASE`` can also be used to execute a one-shot ASE call without a time loop.
This is shown in :doc:`python_interface/phi_ase`.  Generated signatures and member lists are available in the :doc:`Python API Reference <pythonAPI>`.

The backend value is a runtime selection string, not a fixed Python enum.
The minimal example uses ``"Host_Cpu_CpuSerial"`` because that backend is
available in a plain CPU build, but GPU and threaded CPU builds can expose
additional names.  The Python interface exposes the CMake-built backend-name
library through ``AlpakaBackends`` so scripts can query the names supported by
the installed build:

.. code-block:: python

   from HASEonGPU import AlpakaBackends

   print(AlpakaBackends.all())
   backend = AlpakaBackends.Host_Cpu_CpuSerial

Pass one of these strings to ``PhiASE(..., backend=...)``.  The same backend
names are used by the low-level ``calcPhiASE(...)`` interface and by the
``--backend=`` option of the command-line binary.  See
:doc:`Backend Selection <backendSelection>` for build-time backend selection,
runtime backend naming, and troubleshooting the backend-name helper library.

Assemble the Time Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Simulation`` connects the objects above into a time-dependent calculation.
It owns the current time, runs a pump routine, calls ``PhiASE`` for ASE,
combines pump gain, ASE depletion, and fluorescence decay into a derivative
called ``dndtASE`` (:math:`d\beta/dt`), and advances ``betaCells`` with the
selected time integration method.

.. literalinclude:: ../../example/python_example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: simulation
   :end-before: # docs:end: simulation

``RungeKutta4`` is the time integration solver in this example.  Other built-in
solvers and the custom solver protocol are listed in
:doc:`python_interface/utilities`.

The three callback registrations are optional but useful:

``onInit``
   Runs once before the first step and receives the ``Simulation`` object.

``onStep(printState)``
   Runs after every completed step and receives a ``TimeStepState`` snapshot.

``onStep(vtkWedge(...))``
   Writes ``phiAse`` (:math:`\Phi_i`) to VTK files for visualization.  The
   filename template can use values from the state, such as ``{step:03d}``.

One simulation step performs the physical update in this order: pump
contribution, ASE calculation, fluorescence decay, time integration, beta
volume update, result storage, and step callbacks.

``simulation.runSteps(3)`` runs exactly three steps.  For a time-based run use
``simulation.runUntil(endTime=1e-3)`` or set ``endTime`` in the constructor and
call ``simulation.runUntil()``.

Passing PhiASE Settings from YAML
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Python interface can load the ASE compute settings from a YAML file.  This
is useful when the physical setup is assembled in Python, but sampling,
backend, MPI, and output settings should stay in a small run configuration.

.. code-block:: yaml

   phiASE:
     min_rays_per_sample: 100000
     max_rays_per_sample: 1000000
     mse_threshold: 0.05
     repetitions: 2
     adaptive_steps: 4
     use_reflections: true
     monochromatic: false

   compute:
     backend: Host_Cpu_CpuSerial
     parallel_mode: single
     numDevices: 1
     write_vtk: false

.. code-block:: python

   phi_ase = PhiASE.fromYaml(
       "phi_ase.yaml",
       spectralProperties=spectra,
   )

   simulation = Simulation(
       gainMedium=medium,
       pump=pump,
       phiASE=phi_ase,
       timeIntegrationSolver=RungeKutta4(),
       timeStep=1e-5,
   )

``PhiASE.fromYaml(...)`` accepts the same keyword overrides as the constructor.
Overrides are applied after reading the file, so they are the right place for
objects that cannot be represented by the YAML settings, such as
``spectralProperties``.  Keys can be written in the Python attribute style
(``minRaysPerSample``) or in common snake-case style
(``min_rays_per_sample``).  The same settings may be placed at the YAML top
level or grouped under ``phiASE``, ``phi_ase``, ``experiment``, or ``compute``.

The YAML file configures ``PhiASE`` only.  Geometry, gain-medium arrays,
spectral tables, pump properties, and the time-integration solver are still
normal Python objects.  This keeps mesh parsing and array validation in the
typed Python interface while letting run-control parameters live in a
configuration file.  See :doc:`python_interface/phi_ase` for the complete list
of accepted keys and command-line helper support.

Inspect Results
^^^^^^^^^^^^^^^

.. literalinclude:: ../../example/python_example/minimalExampleNewInterface.py
   :language: python
   :start-after: # docs:start: results
   :end-before: # docs:end: results

``simulation.getResults()`` returns the stored ``TimeStepState`` objects.  Each
state contains the completed step index, physical time, ``betaCells``,
``betaVolume``, ``phiAse`` (:math:`\Phi_i`), pump derivative, ASE derivative,
and the raw ASE result object.  These states are the main Python-side output
of the time loop.
