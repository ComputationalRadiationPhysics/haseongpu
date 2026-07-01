PhiASE
======

``PhiASE`` configures and runs the ASE calculation from domain-level Python
objects: a ``GainMedium``, spectral data, and solver/backend settings. The backend transport data is built internally and is not part of the
frontend modeling interface.

.. code-block:: python

   from HASEonGPU import PhiASE

   phi_ase = PhiASE(
       spectralProperties=spectra,
       minRaysPerSample=1000,
       maxRaysPerSample=1000,
       repetitions=1,
       adaptiveSteps=1,
       mseThreshold=0.005,
       useReflections=True,
       backend="Host_Cpu_CpuSerial",
       parallelMode="single",
       numDevices=1,
       rngSeed=1234,
   )

Run One PhiASE Step
-------------------

ASE can be run once without a ``Simulation`` time loop:

.. code-block:: python

   phi_ase.run(gainMedium=medium, crossSections=spectra)
   result = phi_ase.getResults()

   phi = np.asarray(result.phiAse).reshape(
       medium.get("betaCells").expectedShape,
       order="F",
   )

``run(...)`` canonicalizes the domain objects for the openPMD transport,
launches the compiled ``calcPhiASE`` backend, stores the raw result, and
returns ``self``. The returned ``result.phiAse`` values correspond to the ASE
flux :math:`\Phi_i` described in the scientific background.

Sampling and Physics Settings
-----------------------------

``minRaysPerSample``
   Minimum number of Monte Carlo rays :math:`N` used for each sample point.

``maxRaysPerSample``
   Maximum number of Monte Carlo rays used for adaptive sampling.

``mseThreshold``
   Target mean squared error threshold for the ASE estimate.

``repetitions``
   Maximum number of repeated ``PhiASE`` runs using the same number of rays if
   the MSE target is not reached.  Since the importance-sampling distribution
   assigns rays to prisms stochastically, repeating the calculation can improve
   the estimate without increasing the number of rays per run.

``adaptiveSteps``
   If the MSE threshold is not reached for one backend run, HASEonGPU increases
   the rays per sample up to ``maxRaysPerSample``.  This parameter controls how
   many ray-count steps are available between ``minRaysPerSample`` and
   ``maxRaysPerSample``.

``useReflections``
   Enables or disables reflections at the top and bottom surfaces.

``monochromatic``
   Forces the ASE computation to use only the first absorption and emission
   cross-section samples instead of spectral interpolation.

``rngSeed``
   Optional unsigned RNG seed for reproducible Monte Carlo ray sampling.  Set
   this explicitly for reproducible runs.  If omitted, the Python wrapper
   initializes a process-local NumPy seed stream from ``np.random.SeedSequence()``
   and draws one unsigned 32-bit backend seed for each ASE invocation.

Backend and Parallel Settings
-----------------------------

``backend``
   Alpaka backend name.  The minimal example uses ``"Host_Cpu_CpuSerial"``
   because it is available in a plain CPU build.  Query the installed build
   with ``AlpakaBackends.all()`` and pass one of the returned strings here.
   See :doc:`../backendSelection`.

``parallelMode``
   Backend compute mode written to the openPMD ``parallel_mode`` metadata.
   ``"single"`` runs without MPI communication inside one launched process.
   ``"mpi"`` is meaningful when ``calcPhiASE`` is launched under MPI, for
   example through the transport ``command_prefix`` helper or directly with
   ``mpiexec``. Setting this value alone does not create MPI ranks.

``numDevices``
   Maximum number of devices made available on each node for the compute run.
   In MPI execution, HASEonGPU distributes those devices across the MPI ranks
   that are active on the same node.

``nPerNode``
   Launcher setting for Python paths that explicitly call ``mpiexec``. It is
   not serialized as a HASE openPMD transport attribute. See :doc:`../mpi` and
   :doc:`../openpmdTransport` for the interaction between process launching,
   ``parallelMode``, and ``numDevices``.

``minSampleRange`` and ``maxSampleRange``
   Optional inclusive sample-index range.  When omitted, all beta samples
   :math:`\beta_i` are processed.

``writeVtk``
   Not supported by the openPMD transport. Keep this false and use the
   Python-side VTK helpers in :doc:`utilities` for visualization output.

openPMD Transport Options
-------------------------

The openPMD storage backend is selected separately from ``PhiASE.backend``.
The default is ``adios-sst``. Set ``PhiASE.openpmdBackend`` in Python, use
``openpmd_backend`` in YAML, or pass ``--openpmd-backend`` through the
command-line helper to choose a different runtime backend. Lower-level helper
calls also accept ``pyInclude.openpmd.transport.runPhiASE(..., transport=...)``.
Accepted values are ``adios-sst``, ``adios``, and ``hdf5``.

For repeated or streaming use, the transport can keep a session open and write
only dynamic fields after the first iteration. See :doc:`../openpmdTransport`
for the openPMD record layout, storage backend options, artifact-retention environment
variables, and MPI command-prefix examples.

Configuration Helpers
---------------------

``PhiASE`` can read settings from a dictionary or YAML file.  This is intended
for run-control values: sampling, convergence, reflection flags, Alpaka compute
backend selection, openPMD storage backend selection, MPI launcher settings,
optional sample ranges, and the disabled ``writeVtk`` compatibility key.
Objects such as ``GainMedium``, ``SpectralDecomposition``, and pump solvers are
still passed from Python.

.. code-block:: python

   phi_ase = PhiASE({"minRaysPerSample": 1000, "backend": "Host_Cpu_CpuSerial"})
   phi_ase = PhiASE.fromYaml(
       "phi_ase.yaml",
       spectralProperties=spectra,
       gainMedium=medium,
   )

A YAML file can keep experiment and compute settings together:

.. code-block:: yaml

   experiment:
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
     n_per_node: 1
     write_vtk: false  # must remain false for the openPMD transport
     min_sample_range: 0
     max_sample_range: 999
     rng_seed: 1234

YAML keys may be placed at the top level or under ``phiASE``, ``phi_ase``,
``experiment``, or ``compute``.  If the same setting appears more than once,
``PhiASE`` applies sections in this order: ``phiASE``, ``phi_ase``,
``experiment``, ``compute``, then the top-level mapping.  Explicit keyword
overrides passed to ``fromYaml(...)`` are applied after the file is read.

Accepted setting names are the ``PhiASE`` attribute names plus these aliases:
``minRays`` -> ``minRaysPerSample``, ``maxRays`` ->
``maxRaysPerSample``, ``min_rays_per_sample``, ``max_rays_per_sample``,
``mse_threshold``, ``adaptive_steps``, ``use_reflections``,
``parallel_mode``, ``max_gpus`` -> ``numDevices``, ``n_per_node``,
``write_vtk``, ``min_sample_range``, ``max_sample_range``, and
``rng_seed``. ``write_vtk`` is accepted as a configuration key, but the
openPMD transport rejects true values.

Loading YAML requires ``PyYAML``. The package installation installs this
dependency from ``pyproject.toml``; source-tree usage must provide it in the
Python environment.

For command-line tools:

.. code-block:: python

   parser = PhiASE.addArguments(parser)
   args = parser.parse_args()
   phi_ase = PhiASE.fromArgs(args, spectralProperties=spectra)

The command-line helper accepts ``--phi-ase-config`` first and then applies
explicit command-line options such as ``--backend`` or
``--min-rays-per-sample`` as overrides.  It also accepts ``--rng-seed`` for
reproducible Monte Carlo sampling.

Inspection After a Run
----------------------

After ``run(...)``, inspect the simulation result and the original domain
objects rather than backend adapter containers:

.. code-block:: python

   result = phi_ase.getResults()
   points = medium.getPoints()
   prisms = medium.getPrisms()

``getResults()`` raises ``RuntimeError`` if the object has not been run yet.
