Simulation YAML
===============

Schema version 2 constructs the same public objects used by the Python
frontend. It does not expose openPMD records or C++ backend structures.
``Simulation.from_yaml`` constructs the object graph without executing it:

.. code-block:: python

   simulation = Simulation.from_yaml("simulation.yaml")
   simulation.step()

Paths are resolved relative to the YAML file. Unknown keys are errors.
``config/laserPumpCladding.yaml`` is a complete commented example with units,
supported choices, activity schedules, and boundary/relay semantics.

Document structure
------------------

.. code-block:: yaml

   schema_version: 2

   cross_sections:
     ase:
       from_directory:
         path: spectra/ase
         resolution: 1000
     pump:
       monochromatic:
         wavelength: 9.40e-7
         cross_section_absorption: 1.2e-21
         cross_section_emission: 0.0

   simulation:
     gain_medium: {}
     phi_ase:
       ase_steps: 150
     pumps:
       - total_power: 16000.0
         ray_count: 50000
         pump_steps: 50
         rng_seed: 5489
         # spectrum, cross_sections, profile, injection, and optional relays
         # are described below.
     time_integrator: {}
     time_step_size: 2.0e-5
     simulation_steps: 150
     pre_pump: true

Cross-section data
------------------

Each named entry constructs one ``CrossSectionData``. Select exactly one form:

``from_directory``
   ``path`` and optional ``resolution``. The directory contains
   ``lambda_a.txt``, ``sigma_a.txt``, ``lambda_e.txt``, and ``sigma_e.txt``.

``monochromatic``
   ``wavelength``, ``cross_section_absorption``, and
   ``cross_section_emission``.

``inline``
   ``wavelengths_absorption``, ``cross_section_absorption``,
   ``wavelengths_emission``, ``cross_section_emission``, and optional
   ``resolution``.

Gain medium
-----------

Select ``from_vtk`` to construct a complete ``GainMedium`` from a Tet4 VTK
file, or provide ``topology``. A topology selects either ``from_file`` or
``from_tetrahedra``. File formats are ``vtk``, ``gmsh``, and ``stl``;
``boundary_default`` and ``mesh_size`` apply to importers that support them.

Optional ``cell_domains`` and ``surface_domains`` contain the same assignments
accepted by ``VolumeTopology.withCellDomains`` and ``withSurfaceDomains``.
Snake-case selectors are ``cell_indices``, ``face_indices``, ``gmsh_name``,
``gmsh_tag``, ``where``, and ``allow_internal``.

``properties`` accepts only these physical ``GainMedium`` properties:

* ``beta_volume``
* ``cladding_cell_types``
* ``n_tot``
* ``fluorescence_lifetime``
* ``cladding_number``
* ``cladding_absorption``

Boundary optics use structured values exclusively:

.. code-block:: yaml

   surface_optics:
     pump_input:
       reflectivity: 0.0
       n_inside: 1.83
       n_outside: 1.0

Raw reflectivity and refractive-index arrays are private transport data and are
rejected by schema version 2. ``custom_fields`` maps to
``GainMedium.defineField`` and accepts ``name``, ``entity``, ``values``,
``dtype``, ``unit``, ``unit_si``, ``unit_dimension``, ``dynamic``, and
``backend_required``.

PhiASE
------

``phi_ase.cross_sections`` names one cross-section entry. Remaining supported
options are:

* ``propagation_mode``
* ``min_rays``, ``max_rays``, and ``forward_ray_count``
* ``relative_standard_error_threshold``, ``repetitions``, and
  ``adaptive_steps``
* ``rng_seed`` and ``monochromatic``
* ``ase_steps``
* ``use_reflections``, ``reflection_max_iterations``,
  ``reflection_tolerance``, and ``surface_reservoir_size``
* ``min_sample_range`` and ``max_sample_range``
* ``backend`` and ``openpmd_backend``
* ``parallel_mode``, ``num_devices``, and ``n_per_node``

``propagation_mode`` currently accepts only ``forward``. Compute backend names
come from ``AlpakaBackends.all()``. The openPMD choice is ``auto``, ``adios``,
``adios-sst``, or ``hdf5``, subject to provider support.

Pump configuration
------------------

Each ``pumps`` entry accepts ``name``, ``total_power``, ``ray_count``,
``pump_steps``, ``rng_seed``, a named ``cross_sections`` reference, ``spectrum``,
``angular_distribution``, ``profile``, ``injection``, and ``relays``.
``ray_count`` is required. ``pump_steps`` omitted or zero disables that pump;
``rng_seed`` defaults to ``5489``.

``spectrum``
   Either ``monochromatic`` or matching ``wavelengths`` and ``weights``.

``angular_distribution``
   Exactly one of ``collimated``, ``uniform_cone``, or ``discrete``. A uniform
   cone accepts ``half_angle``, ``polar_samples``, and ``azimuthal_samples``.
   Discrete data supplies ``polar_angles``, ``azimuthal_angles``, and
   ``weights``.

``profile``
   ``kind`` is ``uniform`` or ``super_gaussian``. A super-Gaussian accepts
   ``radius_u``, ``radius_v``, ``exponent``, ``center``, ``axis_u``, and
   ``axis_v``.

``injection``
   ``surface_domains`` configures the sole public ``SurfacePumpInjector``.

``relays``
   Each ``PlanarPumpRelay`` accepts ``exit_domains``, ``entry_domains``,
   ``flip_u``, ``flip_v``, ``rotation``, ``offset``, ``tilt``,
   ``magnification``, and ``transmission``. ``transmission`` is the retained
   power fraction through the return mapping, not boundary transmission.

Simulation controls
-------------------

``time_integrator.method`` is one of ``explicit_euler``, ``heun``,
``midpoint``, ``runge_kutta4``, ``frozen_phi_ase_runge_kutta4``,
``implicit_euler``, or ``exponential_euler``. ``implicit_euler`` additionally
accepts ``iterations`` and ``tolerance``.

Configure at most one of ``simulation_steps`` and ``max_time``. If both are
omitted, the run length is the maximum ``phi_ase.ase_steps`` and
``pumps[*].pump_steps``. ``ase_steps`` and each ``pump_steps`` are independent;
omitting either value or setting it to zero disables that component. An
explicitly longer simulation keeps integrating fluorescence and invoking
callbacks after all transport sources stop.

``pre_pump`` is a simulation sequencing control. When true, ASE is suppressed
on the first outer step. Use it when the first pump step should seed excitation
before ASE starts. It is not a ``phi_ase`` or per-pump property.

The remaining controls are ``report_timings``, ``execution_mode``,
``output_steps``, ``output_fields``, and ``control_fields``. Execution modes are
``autonomous`` and ``synchronized-debug``. Supported output fields are
``beta_volume``, ``phi_ase``, ``standard_error``,
``relative_standard_error``, ``total_rays``, ``dndt_ase``, and ``dndt_pump``.
The only control field is currently ``beta_volume``.

Callbacks, logging, and progress rendering are not YAML data. Register Python
callbacks in Python. ``HASE_DEBUG_LOGGING`` and ``HASE_FORWARD_LOGGING`` remain
runtime/build behavior.

``example/laserPumpCladdingYaml.py`` accepts only ``--config`` and output paths.
Edit model, backend, sampling, and schedule values in YAML; the launcher does
not duplicate them as parser arguments.
