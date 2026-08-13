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
   ``PlanarPumpRelay`` objects place that source in the geometry. Each pump owns
   its ray count, RNG seed, and active outer steps.
#. ``PhiASE`` controls the ASE Monte Carlo estimator and its compute and
   transport settings.
#. ``Simulation`` combines those objects with a time integrator and returns one
   ``TimeStepState`` snapshot for each selected output step. With no explicit
   ``output_steps`` schedule, every completed step is selected.

A compact assembly has this shape:

.. code-block:: python

   import numpy as np
   from HASEonGPU import (
       CrossSectionData, FrozenPhiAseRungeKutta4, GainMedium,
       PhiASE, Pump, PumpSpectrum, Simulation,
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
   phi_ase = PhiASE.fromYaml(
       "phiase.yaml", spectralProperties=spectra, ase_steps=150
   )
   pump = Pump(
       total_power=16_000.0,
       spectrum=PumpSpectrum.monochromatic(940e-9),
       cross_sections=pump_cross_sections,
       ray_count=50_000,
       pump_steps=50,
       rng_seed=5489,
   )

   simulation = Simulation(
       gain_medium=medium,
       phi_ase=phi_ase,
       time_integrator=FrozenPhiAseRungeKutta4(),
       time_step_size=2e-5,
       simulation_steps=150,
       cross_sections=spectra,
   ).add_pump(
       pump,
       injection_method=SurfacePumpInjector("pump_input"),
   )

   simulation.on_step(write_state)
   simulation.step()

The names passed to ``SurfacePumpInjector`` and ``SurfaceOptics`` resolve
against surface domains on the topology. The ASE spectrum passed to ``PhiASE``
describes spontaneous-emission transport; a pump may use a separate
single-wavelength ``CrossSectionData`` object. ``Simulation`` owns neither a
second geometry nor a second excitation field: all solvers operate on the
cell-centered state in ``GainMedium``.

Configuration boundaries
------------------------

Geometry, material state, spectra, and physical pumps are Python objects.
``Simulation.from_yaml`` constructs a complete schema-v2 simulation object
graph. ``PhiASE.fromYaml`` can read the ``simulation.phi_ase`` subsection for
applications that assemble the remaining objects in Python. The Alpaka compute
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
   forwardAseRse
   python_interface/simulation
   python_interface/configuration
   python_interface/utilities
