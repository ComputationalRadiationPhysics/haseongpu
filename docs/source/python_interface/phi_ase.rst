PhiASE
======

``PhiASE`` configures and runs the ASE calculation from domain-level Python
objects: a ``GainMedium``, spectral data, and solver/backend settings. The backend transport data is built internally and is not part of the
frontend modeling interface.

.. code-block:: python

   from HASEonGPU import PhiASE

   phi_ase = PhiASE(
       spectralProperties=spectra,
       forwardRayCount=1000,
       repetitions=1,
       relativeStandardErrorThreshold=0.1,
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

``minRays``
   Initial number of globally launched Monte Carlo rays :math:`N`.

``maxRays``
   Hard cap on the cumulative forward-ray histories.  If ``forwardRayCount`` is
   unset and this exceeds ``minRays``, the backend adds geometrically sized
   batches until every cell meets the RSE target or this cap is reached.

``forwardRayCount``
   Explicit fixed number of globally launched forward rays. When nonzero it
   overrides the adaptive ``minRays``/``maxRays`` range.

``relativeStandardErrorThreshold``
   Target dimensionless relative standard error; ``0.1`` requests 10%.

``repetitions``
   Retained transport compatibility field. The forward backend currently uses
   the adaptive ray range instead; it does not schedule extra fixed-count
   repetitions.

``adaptiveSteps``
   Maximum number of geometric ray-count increases between ``minRays`` and
   ``maxRays``. It is ignored for an explicit ``forwardRayCount``.

How direct rays are sampled
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each direct ray samples a spectral bin, an emitting Tet4 volume and point, and
an isotropic direction. The volume probability is proportional to
``betaVolume * cellVolume``; this matches the source strength and uses a unit
source weight. A ray then deposits its gain-weighted track-length contribution
in each Tet4 volume it crosses. The next Tet4 boundary is selected from
precomputed barycentric face planes instead of testing every triangular face.

Within every adaptive batch, direct rays are distributed as evenly as possible
over the discrete ASE spectral bins. They also use one randomly shifted
systematic point per interval of the beta-volume CDF. Both strategies preserve
the original sampling distributions while reducing spectral and source-cell
selection noise.

The stratifiers use global batch indices, so their coverage is preserved when a
batch is split across devices or MPI ranks. They apply to direct volume sources
only. Reflected rays retain the surface-reservoir sampler because their source
distribution is determined by incident boundary weight.

Directions remain isotropic. HASE does not infer a preferred axis from an
arbitrary Tet4 mesh; directional importance sampling requires a separate,
explicitly validated proposal model.

``useReflections``
   Enables or disables the surface-reservoir method (SRM) for reflected ASE
   sources. Direct and relaunched rays always propagate from their starting
   point to a physical mesh boundary; there is no forward ray-length cutoff.

``reflectionMaxIterations``
   Hard cap on reflected SRM passes after the direct pass. It is serialized as
   the openPMD ``reflection_max_iterations`` request attribute and defaults to
   ``40``. ``HASE_SRM_MAX_ITERATIONS`` overrides it at runtime and must be a
   positive integer.

``reflectionTolerance``
   Fraction used for the SRM convergence and stability checks. A pass is
   converged when its remaining reflected source weight divided by the direct
   pass's reflected source weight is below this value.

``surfaceReservoirSize``
   Number of weighted SRM source records retained per physical boundary face.
   It is serialized as ``surface_reservoir_size``.

``HASE_SRM_DIVERGENCE_STREAK``
   Runtime-only positive-integer override for the number of consecutive growing
   reflected passes required to report ``diverged``. It defaults to ``3``.
   This setting is intentionally not an openPMD request attribute: it is a
   runtime safety policy, like ``HASE_SRM_MAX_ITERATIONS``.

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
   ``"mpi"`` makes the Python frontend launch ``calcPhiASE`` through
   ``mpiexec``. The number of ranks per allocated node comes from ``nPerNode``.

``numDevices``
   Maximum number of devices made available on each node for the compute run.
   In MPI execution, HASEonGPU distributes those devices across the MPI ranks
   that are active on the same node.

``nPerNode``
   Number of MPI ranks launched per allocated node when ``parallelMode`` is
   ``"mpi"``. It is not serialized as a HASE openPMD transport attribute. See
   :doc:`../mpi` and :doc:`../openpmdTransport` for the interaction between
   process launching, ``parallelMode``, and ``numDevices``.

``minSampleRange`` and ``maxSampleRange``
   Optional inclusive sample-index range.  When omitted, all beta samples
   :math:`\beta_i` are processed.

openPMD Transport Options
-------------------------

The openPMD storage backend is selected separately from ``PhiASE.backend``.
The default is ``auto``: it selects ``adios-sst``, then ``adios``, then
``hdf5`` from backends supported by both the compiled and Python openPMD
providers. Set ``PhiASE.openpmdBackend`` in Python, use
``openpmd_backend`` in YAML, or pass ``--openpmd-backend`` through the
command-line helper to choose a different runtime backend.
Accepted values are ``auto``, ``adios-sst``, ``adios``, and ``hdf5``.

For repeated or streaming use, the transport can keep a session open and write
only dynamic fields after the first iteration. See :doc:`../openpmdTransport`
for the openPMD record layout, storage backend options, artifact-retention environment
variables, and MPI command-prefix examples.

Configuration Helpers
---------------------

``PhiASE`` can read settings from a dictionary or YAML file.  This is intended
for run-control values: sampling, convergence, reflection flags, Alpaka compute
backend selection, openPMD storage backend selection, MPI launcher settings,
and optional sample ranges.
Objects such as ``GainMedium``, ``SpectralDecomposition``, and pump solvers are
still passed from Python.

.. code-block:: python

   phi_ase = PhiASE({"forwardRayCount": 1000, "backend": "Host_Cpu_CpuSerial"})
   phi_ase = PhiASE.fromYaml(
       "phi_ase.yaml",
       spectralProperties=spectra,
       gainMedium=medium,
   )

A YAML file can keep experiment and compute settings together:

.. code-block:: yaml

   experiment:
     minRays: 100000
     maxRays: 1000000
     relativeStandardErrorThreshold: 0.05
     repetitions: 2
     adaptive_steps: 4
     use_reflections: true
     reflection_max_iterations: 40
     reflection_tolerance: 1.0e-4
     surface_reservoir_size: 32
     monochromatic: false

   compute:
     backend: Host_Cpu_CpuSerial
     parallel_mode: single
     numDevices: 1
     n_per_node: 1
     min_sample_range: 0
     max_sample_range: 999
     rng_seed: 1234

YAML keys may be placed at the top level or under ``phiASE``, ``phi_ase``,
``experiment``, or ``compute``.  If the same setting appears more than once,
``PhiASE`` applies sections in this order: ``phiASE``, ``phi_ase``,
``experiment``, ``compute``, then the top-level mapping.  Explicit keyword
overrides passed to ``fromYaml(...)`` are applied after the file is read.

Accepted setting names are the ``PhiASE`` attribute names plus the legacy aliases
``minRaysPerSample`` -> ``minRays``, ``maxRaysPerSample`` -> ``maxRays``,
``min_rays_per_sample`` -> ``minRays``, ``max_rays_per_sample`` -> ``maxRays``,
``relative_standard_error_threshold``, ``adaptive_steps``, ``use_reflections``,
``parallel_mode``, ``max_gpus`` -> ``numDevices``, ``n_per_node``,
``min_sample_range``, ``max_sample_range``, and ``rng_seed``.

Loading YAML requires ``PyYAML``. Package installation installs that
dependency; the openPMD Python bindings come from the runtime provider selected
at build time. Source-tree usage must provide both in the Python environment.

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
For reflected runs, the returned result also exposes ``srmStatus``,
``srmPasses``, ``srmRemainingFraction``, ``srmMaxIterations``, and
``srmDivergenceStreak``. Terminal statuses are ``converged`` (the reflected
source decayed), ``stable`` (non-growing equilibrium), ``diverged`` (the
configured consecutive-growth streak), and ``max_iterations`` (hard cap
reached). ``disabled`` is reported when reflections were not requested.
