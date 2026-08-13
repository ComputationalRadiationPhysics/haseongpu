PhiASE
======

``PhiASE`` configures the forward, source-driven ASE estimator. It consumes a
Tet4 ``GainMedium`` and cross-section data and returns one flux estimate per
cell. It owns numerical sampling, reflection, compute, transport, and parallel
controls; it does not own geometry or the evolving excitation state.
``propagationMode="forward"`` is the only supported mode.

.. code-block:: python

   phi_ase = PhiASE(
       spectralProperties=spectra,
       minRays=100_000,
       maxRays=1_000_000,
       adaptiveSteps=4,
       relativeStandardErrorThreshold=0.05,
       useReflections=True,
       backend="Host_Cpu_CpuSerial",
       openpmdBackend="auto",
       rngSeed=1234,
       ase_steps=150,
   )

Call it directly for one state or pass it to ``Simulation``:

.. code-block:: python

   phi_ase.run(gainMedium=medium, crossSections=spectra)
   result = phi_ase.getResults()
   phi = np.asarray(result.phiAse)

``getResults`` raises ``RuntimeError`` before a successful run.
The result includes ``phiAse``, ``standardError``,
``relativeStandardError``, ``totalRays``, and ``dndtAse`` plus surface-reservoir
termination information when reflections are enabled. A time-stepped
``Simulation`` exposes the same raw object as ``TimeStepState.ase_result``.

Sampling controls
-----------------

``ase_steps``
   Initial outer simulation steps that include ASE. ``None`` and zero disable
   ASE in ``Simulation``. Direct one-state ``PhiASE.run`` calls are unaffected.

``minRays`` and ``maxRays``
   Initial and maximum global history counts. Adaptive execution adds
   geometrically growing batches until every cell reaches the requested RSE or
   the maximum is reached.

``adaptiveSteps``
   Maximum geometric count increases between the two ray limits.

``forwardRayCount``
   Fixed global history count. Setting it disables adaptive count selection.

``relativeStandardErrorThreshold``
   Target one-sigma uncertainty relative to each cell's estimated mean. ``0.05``
   requests 5%. It measures sampling uncertainty, not discretization or model
   error.

``rngSeed``
   Unsigned seed for reproducible ASE histories. If omitted, each invocation
   draws a process-local seed.

``monochromatic``
   Use the first absorption and emission samples instead of integrating the
   spectrum.

Each direct history samples a spectral bin, a source cell with probability
proportional to ``betaVolume * cellVolume``, a uniform point in that Tet4 cell,
and an isotropic direction. It then deposits a gain-weighted track-length score
in every traversed cell. Spectral bins and source cells are stratified over the
global batch, including when devices or MPI ranks split that batch. See
:ref:`forward-ase-model` for normalization and uncertainty equations.

Reflections
-----------

``useReflections`` enables specular ASE reflection on domain-assigned
``SurfaceOptics``. Direct and reflected rays travel to a physical mesh boundary;
there is no configurable forward ray-length cutoff.

``reflectionMaxIterations``
   Maximum surface-reservoir passes after the direct pass. The positive integer
   environment override is ``HASE_SRM_MAX_ITERATIONS``.

``reflectionTolerance``
   Stop when remaining reflected source weight, relative to the direct pass,
   falls below this fraction.

``surfaceReservoirSize``
   Maximum weighted source records retained per physical boundary face.

The runtime reports ``srmStatus``, ``srmPasses``, ``srmRemainingFraction``,
``srmMaxIterations``, and ``srmDivergenceStreak``. Terminal status can be
``converged``, ``stable``, ``diverged``, or ``max_iterations``; ``disabled``
means reflections were not requested. ``HASE_SRM_DIVERGENCE_STREAK`` controls
how many consecutive growing passes report divergence.

The current boundary model uses configured constant reflectivity and total
internal reflection. It does not calculate Fresnel coefficients or launch a
transmitted/refracted ray. See :ref:`ase-surface-reflections` for the physical
model and termination criteria.

Compute and transport
---------------------

``backend`` selects an Alpaka compute backend reported by
``AlpakaBackends.all()``. ``openpmdBackend`` independently selects the transport
format/engine. See :doc:`../backendSelection` for selection syntax and
:doc:`../openpmdTransport` for available storage backends and provider
compatibility.

``parallelMode="single"`` runs one process. ``parallelMode="mpi"`` asks the
frontend to launch through MPI; ``nPerNode`` selects ranks per node and
``numDevices`` limits devices available on each node. See :doc:`../mpi`.

``minSampleRange`` and ``maxSampleRange`` optionally restrict the inclusive
flattened cell range. Normal full-volume runs leave them unset.

YAML and CLI helpers
--------------------

``fromYaml`` accepts schema-v2 PhiASE settings under ``simulation.phi_ase``:

.. code-block:: yaml

   schema_version: 2
   simulation:
     phi_ase:
       min_rays: 100000
       max_rays: 1000000
       relative_standard_error_threshold: 0.05
       adaptive_steps: 4
       ase_steps: 150
       use_reflections: true
       reflection_max_iterations: 40
       backend: Host_Cpu_CpuSerial
       openpmd_backend: auto
       parallel_mode: single
       rng_seed: 1234

.. code-block:: python

   phi_ase = PhiASE.fromYaml(
       "config/hase-phiase.yaml",
       spectralProperties=spectra,
       maxRays=2_000_000,
   )

Keyword arguments override file values. ``addArguments`` and ``fromArgs`` add
the same controls to an ``argparse`` command. Geometry, state, spectra, and pump
objects remain Python inputs and are not YAML run-control data.
