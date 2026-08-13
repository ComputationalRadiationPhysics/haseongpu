laserPumpCladding Tutorial
==========================

This tutorial follows ``example/laserPumpCladdingApi.py`` from material data to
a time-dependent pump-and-ASE result. The equivalent
``example/laserPumpCladdingYaml.py`` constructs the same public objects from
``config/laserPumpCladding.yaml``. The configuration includes comments for
units, supported choices, scheduling, and reflection behavior.

The example models a cylindrical Tet4 gain medium. A collimated 940 nm pump
enters through the bottom surface, traverses the volume, and is returned once
from the top. The compiled C++/Alpaka simulation evolves one excited-state
fraction per Tet4 cell from pump excitation, ASE depletion, and fluorescence
decay.

Run a short simulation
----------------------

Install HASEonGPU as described in :doc:`gettingStarted`, then run from the
repository root:

.. code-block:: bash

   python3 example/laserPumpCladdingApi.py \
       --backend Host_Cpu_CpuSerial \
       --simulation-steps 3 \
       --pump-steps 3 \
       --ase-steps 3 \
       --pump-ray-count 1000 \
       --rng-seed 1234 \
       --pump-rng-seed 5489 \
       --vtk-output-dir output

``--backend`` selects an Alpaka compute backend compiled into the installed
runtime. ``--openpmd-backend`` is a separate storage/streaming choice. See
:doc:`backendSelection` for the distinction and :doc:`openpmdTransport` for the
storage backends.

To use only declarative configuration, edit a copy of
``config/laserPumpCladding.yaml`` and run:

.. code-block:: bash

   python3 example/laserPumpCladdingYaml.py \
       --config config/laserPumpCladding.yaml \
       --vtk-output-dir output

The YAML entry point accepts only the configuration path and output locations.
Backend, time-step, pump-ray, seed, and physics values are set in the YAML, not
as launcher options. The API entry point exposes command-line controls for the
model assembled directly in Python.

The default output schedule writes every completed step, producing
``output/laserPumpCladding_001.vtk`` through
``output/laserPumpCladding_003.vtk``. The example also prints the final PhiASE
and excitation array shapes. The small ray count makes this a workflow check,
not a statistically converged production calculation.

Load the material spectra
-------------------------

The four files in ``example/input`` contain measured absorption and emission
wavelengths and their matching cross sections. The helper
``laserPumpCladdingSpectralProperties`` loads them into one
``CrossSectionData`` object:

.. code-block:: python

   spectralProperties = CrossSectionData(
       wavelengthsAbsorption=raw_wavelengths_absorption,
       crossSectionAbsorption=raw_cross_section_absorption,
       wavelengthsEmission=raw_wavelengths_emission,
       crossSectionEmission=raw_cross_section_emission,
       resolution=spectralResolution,
   )

``PhiASE`` uses this complete spectrum when sampling spontaneous-emission
wavelengths. The ``resolution`` value controls the backend's numerical spectral
table; it does not add information to the measured input. See
:doc:`Cross-section data <python_interface/spectral_decomposition>` for the
field units, constructors, and interpolation rules.

The physical pump is monochromatic, so it needs one absorption and one emission
cross section at 940 nm. ``runExample`` obtains them from the same material
tables:

.. code-block:: python

   pumpWavelength = 940e-9
   pumpCrossSections = CrossSectionData.monochromatic(
       wavelength=pumpWavelength,
       crossSectionAbsorption=np.interp(
           pumpWavelength * 1e9,
           raw_wavelengths_absorption,
           raw_cross_section_absorption,
       ),
       crossSectionEmission=np.interp(
           pumpWavelength * 1e9,
           raw_wavelengths_emission,
           raw_cross_section_emission,
       ),
   )

The wavelength tables use nanometers, hence the conversion from meters before
``np.interp``. In the shipped tables, 940 nm is itself a sample and the call
returns that tabulated value. For a pump wavelength between samples, the same
code would linearly interpolate between the two neighboring rows.

Build the Tet4 gain medium
--------------------------

``VolumeTopology`` owns geometry: points, tetrahedral connectivity, derived
cell volumes and face geometry, and domain labels. The example loads
``example/data/ptTet4.vtk`` and identifies the cylinder's exterior regions:

.. code-block:: python

   topology = VolumeTopology.fromVtk(materialPath)
   topology = topology.withCellDomains(
       domain=1,
       name="gain_medium",
       where="all",
   ).withSurfaceDomains([
       {"domain": 1, "name": "ase_bottom", "faceIndices": np.argwhere(bottom)},
       {"domain": 2, "name": "ase_top", "faceIndices": np.argwhere(top)},
       {"domain": 3, "name": "cladding", "faceIndices": np.argwhere(side)},
   ])

``gain_medium`` labels all volume cells. ``ase_bottom`` is later used by the
pump injector, ``ase_top`` by the return relay, and all three surface names may
carry ASE boundary optics. New gmsh inputs should normally provide these names
as physical groups; the geometry-specific face selection exists because this
converted VTK fixture has no named surface groups. The general domain API is
described in :doc:`Volume topology <python_interface/topology>`.

``GainMedium`` attaches physical state and material properties to that topology:

.. code-block:: python

   medium = GainMedium(topology).withPhysicalProperties(
       betaVolume=np.zeros(topology.numberOfCells),
       claddingCellTypes=np.zeros(topology.numberOfCells, dtype=np.uint32),
       nTot=2 * 1.388e20,
       crystalTFluo=9.41e-4,
       claddingNumber=1,
       claddingAbsorption=5.5,
   ).with_surface_optics({
       "ase_bottom": SurfaceOptics(
           reflectivity=0.0, n_inside=1.83, n_outside=1.0
       ),
       "ase_top": SurfaceOptics(
           reflectivity=0.0, n_inside=1.83, n_outside=1.0
       ),
       "cladding": SurfaceOptics(
           reflectivity=0.0, n_inside=1.0, n_outside=1.0
       ),
   })

``betaVolume`` is the evolving excited-state fraction :math:`\beta_j` and has
one value per Tet4 cell. ``nTot`` is the active-ion concentration and
``crystalTFluo`` the fluorescence lifetime. ``claddingCellTypes`` is initialized
to zero here, so the default fixture's cells use the active-medium gain model;
the named ``cladding`` region is its cylindrical side boundary. A cell whose
type equals ``claddingNumber`` instead uses ``claddingAbsorption`` during
transport. See :doc:`GainMedium <python_interface/gain_medium>` for the full
field contract and :doc:`theoryAndModel` for the gain equation.

Configure ASE
-------------

The API example constructs adaptive sampling, reflection, scheduling, compute,
openPMD, and MPI controls directly:

.. code-block:: python

   phiAse = PhiASE(
       spectralProperties=spectralProperties,
       minRays=10_000,
       maxRays=1_000_000,
       useReflections=True,
       ase_steps=150,
       backend="Host_Cpu_CpuOmpBlocks",
       openpmdBackend="auto",
   )

``ase_steps=None`` and ``ase_steps=0`` both disable ASE in a time-stepped
simulation. A positive value enables ASE through that outer step. The shipped
API and YAML examples both enable ASE surface reflections for 150 steps.
The estimator and its uncertainty controls are described in :doc:`PhiASE
<python_interface/phi_ase>`.

Describe the physical pump
--------------------------

The pump profile weights launch positions over ``ase_bottom``. Its integral
converts the example's legacy peak-power-density convention into the
aperture-integrated power expected by ``Pump.total_power``:

.. code-block:: python

   pumpProfile = SuperGaussianPumpProfile(
       radius_u=1.5,
       radius_v=1.5,
       exponent=40,
   )
   profileArea = integrate_pump_profile(
       medium.topology,
       "ase_bottom",
       pumpProfile,
   )
   pump = Pump(
       total_power=16e3 * profileArea,
       spectrum=PumpSpectrum.monochromatic(pumpWavelength),
       cross_sections=pumpCrossSections,
       ray_count=50_000,
       pump_steps=50,
       rng_seed=5489,
       angular_distribution=PumpAngularDistribution.collimated(),
       profile=pumpProfile,
   )

Each ``Pump`` owns its sampling budget, reproducibility seed, and active outer
steps. This permits independently sampled pumps with different durations.
``pump_steps=None`` and ``pump_steps=0`` disable that pump. The injector still
owns the geometric launch location. See :doc:`Pump configuration
<python_interface/pump_properties>` for other spectra, angular distributions,
profiles, and relay transforms.

Assemble the simulation
-----------------------

``Simulation`` combines the medium, ASE estimator, spectra, and time integrator,
and owns only the outer-loop horizon:

.. code-block:: python

   simulation = Simulation(
       gain_medium=medium,
       phi_ase=phiAse,
       time_integrator=FrozenPhiAseRungeKutta4(),
       time_step_size=2e-5,
       simulation_steps=150,
       cross_sections=spectralProperties,
       pre_pump=True,
       output_steps=(
           None
           if outputSteps is None
           else tuple(int(step) for step in outputSteps)
       ),
   ).add_pump(
       pump,
       injection_method=SurfacePumpInjector("ase_bottom"),
       relays=(PlanarPumpRelay.retroreflect("ase_top"),),
   )

``SurfacePumpInjector`` launches the pump on the named bottom faces.
``PlanarPumpRelay.retroreflect("ase_top")`` represents one explicit return pass
through the top aperture. Relay ``transmission`` is the retained power fraction
through that return mapping; it is not transmission through the gain-medium
boundary. A value of one reinjects every exiting pump ray inward without relay
loss. The relay is not an unlimited resonator, coating model, or Fresnel
calculation. ASE reflection is a separate ``PhiASE`` process using the
``SurfaceOptics`` fields.

``Simulation.simulation_steps`` may be longer or shorter than any pump or ASE
activity window. If it and ``max_time`` are omitted, ``simulation.step()`` uses
the maximum ``Pump.pump_steps`` and ``PhiASE.ase_steps``. An explicitly longer
run continues fluorescence evolution and ``on_step``/``before_step`` callbacks
after pump and ASE stop. With
``output_steps`` omitted, the backend emits every completed step
and the registered callbacks write one VTK snapshot per step. Use
``--output-steps`` to select explicit one-based completed-step indices when
only part of the trajectory is needed. For example, ``--output-steps 150``
writes only the final snapshot of a 150-step run.

``FrozenPhiAseRungeKutta4`` evaluates PhiASE once per RK4 step and reuses it for
the remaining RK stages, while the pump term is still evaluated at each stage.
The resulting population equation and the behavior of ``pre_pump`` are
described in :doc:`Pump and Time Stepping
<theoryAndModel>`.

Read and write results
----------------------

The example registers callbacks before starting the compiled run:

.. code-block:: python

   simulation.on_step(printState)
   simulation.on_step(
       writeVtkFields,
       vtkOutputDir,
       absorption,
       spectralProperties,
       medium.get("nTot").value,
   )
   simulation.step()
   state = simulation.get_last_state()

``on_step`` runs once for each snapshot selected by ``output_steps``. It does
not control or synchronize backend stepping.

Each callback receives a ``TimeStepState`` for an emitted step, with cell-centered
``beta_volume``, ``phi_ase``, ``dndt_pump``, and ``dndt_ase`` arrays. The VTK
callback additionally writes ``cladAbs`` and ``localGain``. ``Simulation`` keeps
only the last snapshot; callbacks are the place to retain a history.

Printed array means are quick diagnostics. For nonuniform Tet4 cells, use the
cell volumes for a physical integral or volume-weighted mean:

.. code-block:: python

   volumes = state.topology.cellVolumes
   phi_integral = np.sum(state.phi_ase * volumes)
   phi_volume_mean = np.average(state.phi_ase, weights=volumes)

Use :doc:`scripts` for the small-signal-gain plotting helper and
:doc:`Utilities <python_interface/utilities>` for VTK and openPMD/ParaView
output.

Evaluate a prepared state without pumping
-----------------------------------------

``--phiase-only`` skips pump construction and time integration. It loads a Tet4
VTK file containing geometry, ``betaVolume``, and supported material fields,
then performs one forward PhiASE evaluation:

.. code-block:: bash

   python3 example/laserPumpCladdingApi.py \
       --phiase-only \
       --tet4-input prepared-state.vtk \
       --backend Host_Cpu_CpuSerial \
       --rng-seed 1234

Use ``VolumeTopology.fromVtk`` for geometry-only input and
``GainMedium.fromVtk`` for a prepared geometry-plus-state file.

Command-line options
--------------------

Execution and configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``-h``, ``--help``
   Print the generated option synopsis and exit.

``--backend BACKEND``
   Alpaka compute backend. Default: ``Host_Cpu_CpuOmpBlocks``.

``--openpmd-backend BACKEND``
   openPMD storage backend. Default: ``auto``.

``--timings``
   Print frontend time spent in compiled transport and decoding, snapshot
   materialization, and callbacks. This is distinct from build-time benchmark
   instrumentation documented in :doc:`compilation`.

Time and physics controls
^^^^^^^^^^^^^^^^^^^^^^^^^

``--simulation-steps N``
   Total outer simulation steps. Default: ``150``.

``--output-steps STEP [STEP ...]``
   One-based completed-step indices to emit. If omitted, emit every completed
   step.

``--pump-steps N``
   Initial outer steps that include this pump. Default: ``50``. Zero disables
   the pump.

``--ase-steps N``
   Initial outer steps that include ASE. Default: ``150``. Zero disables ASE.

``--disable-pre-pump``
   Enable ASE during the first pumped step. By default, the first step seeds
   excitation with pump and fluorescence before the first ASE evaluation.

``--spectral-resolution N``
   Numerical spectral-table resolution passed to the ASE backend. Default:
   ``1000``.

Monte Carlo and reproducibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--rng-seed N``
   ASE-history seed. By default no explicit seed is supplied.

``--pump-ray-count N``
   Equal-power Monte Carlo histories per registered pump source. Default:
   ``50000``.

``--pump-rng-seed N``
   Pump-history seed. Default: ``5489``. Set both pump and ASE seeds when
   comparing complete runs.

Input, output, and alternate mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--vtk-output-dir PATH``
   Directory for one ``laserPumpCladding_<step>.vtk`` file per emitted step.
   Defaults to the example directory.

``--openpmd-output-dir PATH``
   Optional directory for callback snapshots written as an openPMD series and
   ParaView handle. No openPMD visualization output is written by default.

``--phiase-only``
   Run one forward PhiASE evaluation instead of the time-stepped pump workflow.

``--tet4-input PATH``
   Prepared Tet4 VTK state used by ``--phiase-only``. The parser rejects
   ``--phiase-only`` when this option is missing.
