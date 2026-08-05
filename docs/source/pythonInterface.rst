Python Interface Guide
======================

The Python interface describes a simulation as a graph of physical and
numerical objects. Python assembles that graph; the compiled C++/Alpaka runtime
performs pump transport, ASE transport, and time integration. Current forward
simulations use an explicit Tet4 ``VolumeTopology`` and one evolving excitation
value per cell.

The :doc:`laserPumpCladding tutorial <laserPumpCladding>` shows these objects in
one runnable workflow. This guide explains how the pieces relate and links to
the page that owns each concept. Generated signatures and complete member lists
are in the :doc:`Python API Reference <pythonAPI>`.

Object model
------------

The setup proceeds from geometry to physical state and then to solvers:

#. ``VolumeTopology`` owns points, Tet4 connectivity, derived geometry, and
   named cell and surface domains.
#. ``GainMedium`` combines that topology with excitation, material constants,
   cladding fields, and boundary optics.
#. ``CrossSectionData`` supplies wavelength-dependent absorption and emission
   cross sections.
#. ``Pump`` describes a physical source. ``SurfacePumpInjector`` and optional
   ``PlanarPumpRelay`` objects place that source in the geometry, while
   ``MonteCarloPumpSolver`` controls its numerical sampling.
#. ``PhiASE`` controls the ASE Monte Carlo estimator and its compute and
   transport settings.
#. ``Simulation`` combines those objects with a time integrator and returns one
   ``TimeStepState`` snapshot after each completed step.

A compact assembly has this shape:

.. code-block:: python

   import numpy as np
   from HASEonGPU import (
       CrossSectionData, FrozenPhiAseRungeKutta4, GainMedium,
       MonteCarloPumpSolver, PhiASE, Pump, PumpSpectrum, Simulation,
       SurfacePumpInjector, VolumeTopology,
   )

   pump_cross_sections = ...
   write_state = ...

   topology = VolumeTopology.fromFile("crystal.msh")
   medium = GainMedium(topology).withPhysicalProperties(
       betaVolume=np.zeros(topology.numberOfCells),
       nTot=2.776e20,
       crystalTFluo=9.41e-4,
   )

   spectra = CrossSectionData.fromDirectory("spectra", resolution=1000)
   phi_ase = PhiASE.fromYaml("phiase.yaml", spectralProperties=spectra)
   pump = Pump(
       total_power=16_000.0,
       spectrum=PumpSpectrum.monochromatic(940e-9),
       cross_sections=pump_cross_sections,
   )

   simulation = Simulation(
       gain_medium=medium,
       phi_ase=phi_ase,
       time_integrator=FrozenPhiAseRungeKutta4(),
       time_step_size=2e-5,
       pump_solver=MonteCarloPumpSolver(ray_count=50_000, seed=5489),
       cross_sections=spectra,
   ).add_pump(
       pump,
       injection_method=SurfacePumpInjector("pump_input"),
   )

   simulation.on_step(write_state)
   simulation.step(10)

The names passed to ``SurfacePumpInjector`` and ``SurfaceOptics`` resolve
against surface domains on the topology. The ASE spectrum passed to ``PhiASE``
describes spontaneous-emission transport; a pump may use a separate
single-wavelength ``CrossSectionData`` object. ``Simulation`` owns neither a
second geometry nor a second excitation field: all solvers operate on the
cell-centered state in ``GainMedium``.

Configuration boundaries
------------------------

Geometry, material state, spectra, and physical pumps are Python objects.
``PhiASE.fromYaml`` loads numerical ASE, compute-backend, openPMD, and MPI run
controls; keyword arguments override values from the file. The Alpaka compute
backend and openPMD storage backend are independent choices. See
:doc:`Backend Selection <backendSelection>` and :doc:`openPMD Transport
<openpmdTransport>` rather than duplicating those choices in simulation code.

Concept pages
-------------

.. toctree::
   :maxdepth: 2
   :caption: Python modeling concepts

   python_interface/topology
   python_interface/gain_medium
   python_interface/spectral_decomposition
   python_interface/pump_properties
   python_interface/phi_ase
   python_interface/simulation
   python_interface/utilities
