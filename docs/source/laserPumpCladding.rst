laserPumpCladding Tutorial
==========================

This tutorial follows ``example/laserPumpCladding.py`` from material data to a
time-dependent pump-and-ASE result. It introduces one concept at a time and
links to the corresponding concept or reference page for details.

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

   python3 example/laserPumpCladding.py \
       --backend Host_Cpu_CpuSerial \
       --timeSteps 3 \
       --pumpSteps 3 \
       --output-steps 1 2 3 \
       --pump-ray-count 1000 \
       --rng-seed 1234 \
       --pump-rng-seed 5489 \
       --vtk-output-dir output

``--backend`` selects an Alpaka compute backend compiled into the installed
runtime. ``--openpmd-backend`` is a separate storage/streaming choice. See
:doc:`backendSelection` for the distinction and :doc:`openpmdTransport` for the
storage backends.

The explicit output schedule makes the short run write
``output/laserPumpCladding_001.vtk`` through
``output/laserPumpCladding_003.vtk``. Without ``--output-steps``, the example
emits only the pump-boundary and final snapshots. It also prints the final
PhiASE and excitation array shapes. The small ray count makes this a workflow
check, not a statistically converged production calculation.

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
   ).withSurfaceOptics({
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

The example reads adaptive sampling, reflection, compute, openPMD, and MPI
controls from ``config/hase-phiase.yaml``:

.. code-block:: python

   phiAse = PhiASE.fromYaml(
       phiAseConfigPath,
       spectralProperties=spectralProperties,
       **AseOverride,
   )

Command-line backend and seed values override the corresponding YAML settings.
Geometry and excitation are not duplicated in the YAML; they remain in
``medium`` and are supplied when the simulation runs. The estimator and its
uncertainty controls are described in :doc:`PhiASE
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
       angular_distribution=PumpAngularDistribution.collimated(),
       profile=pumpProfile,
   )

``Pump`` contains physical source properties. It does not decide where rays
enter or how many rays approximate the source. Those numerical and geometrical
choices are attached when the simulation is assembled. See :doc:`Pump
configuration <python_interface/pump_properties>` for other spectra, angular
distributions, profiles, and relay transforms.

Assemble the simulation
-----------------------

``Simulation`` combines the medium, ASE estimator, spectra, time integrator, and
pump sampling controls:

.. code-block:: python

   simulation = Simulation(
       gain_medium=medium,
       phi_ase=phiAse,
       time_integrator=FrozenPhiAseRungeKutta4(),
       time_step_size=2e-5,
       pump_solver=MonteCarloPumpSolver(
           ray_count=pumpRayCount,
           seed=pumpRngSeed,
           max_steps=pumpSteps,
       ),
       cross_sections=spectralProperties,
       enable_ase=enableASE,
       pre_pump=prePump,
       output_steps=(
           laserPumpCladdingOutputSteps(timeSlices, pumpSteps)
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
through the top aperture. It is not an unlimited resonator, coating model, or
Fresnel calculation. ASE reflection is a separate ``PhiASE`` process using the
``SurfaceOptics`` fields.

``MonteCarloPumpSolver.ray_count`` controls equal-power pump histories per
source evaluation. Its ``max_steps`` limits the outer simulation steps that
include the pump. ``simulation.step(timeSlices)`` controls the total number of
outer steps, so the run can continue with ASE and fluorescence after pumping
stops. By default, ``laserPumpCladdingOutputSteps`` emits the pump boundary and
the final step; ``--output-steps`` replaces that schedule with explicit
one-based completed-step indices.

``FrozenPhiAseRungeKutta4`` evaluates PhiASE once per RK4 step and reuses it for
the remaining RK stages, while the pump term is still evaluated at each stage.
The resulting population equation and the behavior of ``pre_pump`` and
``enable_ase`` are described in :doc:`Pump and Time Stepping
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
   simulation.step(timeSlices)
   state = simulation.get_last_state()

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

   python3 example/laserPumpCladding.py \
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
   Alpaka compute backend. The default ``UseConfig`` keeps the YAML value.

``--openpmd-backend BACKEND``
   openPMD storage backend. The default ``UseConfig`` keeps the YAML value; use
   ``auto`` to select a compatible installed backend automatically.

``--phi-ase-config PATH``
   PhiASE run-control YAML. Defaults to ``config/hase-phiase.yaml``.

``--timings``
   Print frontend time spent in compiled transport and decoding, snapshot
   materialization, and callbacks. This is distinct from build-time benchmark
   instrumentation documented in :doc:`compilation`.

Time and physics controls
^^^^^^^^^^^^^^^^^^^^^^^^^

``--timeSteps N``
   Total outer simulation steps. Default: ``150``.

``--output-steps STEP [STEP ...]``
   One-based completed-step indices to emit. By default, the example emits the
   pump boundary, when reached, and the final step.

``--pumpSteps N``
   Initial outer steps that include pump transport. Default: ``100``. Set it to
   ``timeSteps`` to pump for the complete run.

``--disable-ase``
   Evolve pump excitation and fluorescence without ASE transport or depletion.

``--disable-pre-pump``
   Enable ASE during the first pumped step. By default, the first step seeds
   excitation with pump and fluorescence before the first ASE evaluation.

``--disable-reflections``
   Override the PhiASE configuration and disable ASE surface reflections.

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
