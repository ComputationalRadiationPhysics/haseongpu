PhiASE
======

``PhiASE`` configures the forward, source-driven ASE estimator. ``Simulation``
supplies its gain domains and the authoritative spectra from each component's
``Material``. ``PhiASE`` owns numerical sampling, reflection, compute,
transport, and parallel controls; it does not own geometry, material, or the
evolving excitation state.
``propagationMode="forward"`` is the only supported mode.

.. code-block:: python

   phi_ase = PhiASE(
       minRays=100_000,
       maxRays=1_000_000,
       adaptiveSteps=4,
       relativeStandardErrorThreshold=0.05,
       enableDiagnostics=True,
       useReflections=True,
       reflectionMode="direct",
       backend="Host_Cpu_CpuSerial",
       openpmdBackend="auto",
       rngSeed=1234,
       ase_steps=150,
   )

Normal applications pass it to ``Simulation``. The direct ``run`` entry point
accepts the same physical ``GainMedium`` graph. Pass ``opticalComponents`` when
the trace also traverses passive components. Spectra are obtained from the
``Material`` referenced by each ``OpticalComponent``; ``PhiASE`` does not
create a second transport copy:

.. code-block:: python

   phi_ase.run(
       gainMedium=medium,
       opticalComponents=components,
       initialExcitation=0.25,
   )
   result = phi_ase.getResults()
   phi = np.asarray(result.phiAse)

``getResults`` raises ``RuntimeError`` before a successful run.
The result includes ``phiAse``, ``standardError``,
``relativeStandardError``, ``totalRays``, and ``dndtAse`` plus boundary-pass
termination information when reflection or inter-component transmission is
active. A time-stepped
``Simulation`` exposes the same raw object as ``TimeStepState.aseResult``.

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
   Maximum geometric count increases between the two ray limits. Zero performs
   one evaluation at ``minRays`` (unless ``forwardRayCount`` is set).

``forwardRayCount``
   Fixed global history count. Setting it disables adaptive count selection.

``relativeStandardErrorThreshold``
   Target one-sigma uncertainty relative to each cell's estimated mean. ``0.05``
   requests 5%. It measures sampling uncertainty, not discretization or model
   error.

``enableDiagnostics``
   Enable optional per-cell ray visits in ``totalRays``. Diagnostics use 128
   threads per block; the performance specialization uses 512 and omits visit
   counters. Essential failed-ray accounting remains enabled in both modes.
   With diagnostics disabled, ``totalRays`` retains its cell-shaped layout but
   contains zeros; dropped-ray counts still report failures.
   Visits are not the global sampling budget: zero-weight primary samples are
   not transported and make no visits, even with diagnostics enabled. In an
   entirely unpumped solve, adaptive budget accounting can still reach
   ``maxRays`` while ``totalRays`` remains zero and RSE is undefined.
   ``trackRayVisits`` remains accepted as a deprecated constructor and
   transport input alias, but new output contains only ``enableDiagnostics``.

``rngSeed``
   Unsigned seed for reproducible ASE histories. If omitted, each invocation
   draws a process-local seed.

``monochromatic``
   Use the first absorption and emission samples instead of integrating the
   spectrum.

Each direct history samples a spectral bin, a source cell with probability
proportional to its spontaneous-source strength, a uniform point in that Tet4
cell, and an isotropic direction. The global count is divided among optical
components from their current total source strengths. Setting
``OpticalComponent.aseRays`` reserves that component's exact final primary-ray
count; unspecified positive-source components share the remainder. It then
deposits a gain-weighted track-length score in every traversed cell. Spectral
bins and source cells are stratified within each domain and statistical batch,
with an independently keyed permutation separating the two dimensions.
Every independent ray population represents all positive-source components,
using population-specific importance weights. ``numIndependentRayPopulations``
is a runtime integer (default 8), independent of worker count. YAML accepts
``num_independent_ray_populations`` and the CLI accepts
``--num-independent-ray-populations``. Each emitting component must have enough
rays to occur in every population; insufficient quotas are rejected. Adaptive
increases can be coalesced while preserving exact final quotas.
One population cannot estimate RSE and reports the maximum error sentinel. See
:ref:`forward-ase-model` for normalization and uncertainty equations.

Domain-local source CDFs preserve very weak sources beside strong components.
Zero-source domains receive no automatically allocated primary rays, but remain
eligible for transmitted rays and their resulting field and RSE estimates.

Forward traversal has no fixed cell-crossing limit. A valid ray continues until
it reaches a physical boundary or a cell policy terminates it, so increasing
mesh resolution cannot make a ray fail merely because it requires more
crossings. Invalid geometric transitions and non-finite contributions are
always counted as dropped rays. ``Simulation`` rejects dropped histories and
non-finite or unrepresentable ASE output before updating excitation.

Domain boundaries
-----------------

``useReflections`` enables specular ASE reflection on domain-assigned
``SurfaceOptics``. Direct and reflected rays travel to a physical mesh boundary;
there is no configurable forward ray-length cutoff.

``reflectionMode``
   Selects the boundary-history representation. ``"direct"`` (the default)
   stores up to two compact exact-intersection children per boundary hit and
   performs systematic particle combing. ``"srm"`` instead offers both children
   to a bounded weighted sample per mesh face. Both policies route transmitted
   histories between adjacent optical components and keep their large buffers,
   scans, and selections on the accelerator.

``surfaceReservoirSize``
   Number of statistically retained ray records per boundary face when
   ``reflectionMode="srm"``, per logical batch. Reflected and transmitted weight
   is accumulated independently of this bounded record count.

   SRM divides each source domain's ray population into logical batches of at
   most 65,536 primary rays before assigning workers. Source and wavelength
   strata belong to the complete domain population, not to individual logical
   batches. Each batch keeps its own reservoirs, random stream, and boundary
   stopping state throughout transport. Reservoirs are reused between completed
   batches on a worker; they are not merged merely because batches share a
   worker. Transmitted rays can score domains with no primary emission.

   Raw batch contributions are summed under their ``rayPopulationId`` before
   population normalization and RSE estimation. Worker count changes ownership,
   not these statistical boundaries. Floating-point reduction order can still
   cause small numerical differences. Changing the logical batch cap can change
   sampling variance; it must not be tuned automatically from worker count.
   More workers than logical batches can leave workers idle.

``srmPositionMode``
   Selects where retained SRM records are relaunched. ``"exact"`` retains each
   sampled boundary intersection, while ``"centroid"`` allocates no device
   position buffers and relaunches every record at the centroid of its owning
   face. Directions, weights, wavelengths, and selection keys remain bounded by
   faces times ``surfaceReservoirSize``. This setting affects only
   ``reflectionMode="srm"``; direct mode always uses exact intersections.

``boundaryMaxPasses``
   Optional explicit limit for direct or SRM boundary passes. ``None`` uses
   ``reflectionMaxIterations``.

``reflectionMaxIterations``
   Default limit for direct or SRM boundary passes after the direct volume-source
   pass. For ``reflectionMode="srm"``, the positive integer environment override
   is ``HASE_SRM_MAX_ITERATIONS``.

``reflectionTolerance``
   Stop when remaining reflected source weight, relative to the direct pass,
   falls below this fraction.

The runtime reports ``boundaryStatus``, ``boundaryPasses``,
``boundaryRemainingFraction``, ``boundaryMaxPasses``,
``boundaryDivergenceStreak``, ``boundaryGamma``,
``boundaryGammaStandardError``, ``boundaryTailFactor``, and
``boundaryTailClosure``. Terminal status can be ``converged``, ``stable``,
``diverged``, or ``maxPasses``; ``disabled`` means neither reflections nor
inter-component routing required boundary passes.
``HASE_SRM_DIVERGENCE_STREAK`` controls how many consecutive growing SRM passes
report divergence.

For a truncated series, HASE fits the recent reflected population as
:math:`W_p \simeq W_0\Gamma^p`. The scalar fit and candidate tail factor are
diagnostics only: geometric decay of total weight does not establish a
geometric spatial or spectral field. No analytical tail is added and the fit
does not promote an unresolved result to ``converged``.
``Simulation`` stops before updating excitation when the status is ``stable``,
``diverged`` or ``maxPasses``; standalone ``PhiASE.run`` returns those partial
tallies and diagnostics for analysis. Increasing the pass limit can provide
more evidence, but does not by itself establish a finite steady-state field.

With ``useReflections`` enabled, each eligible interface hit creates both
histories: the reflected child has weight ``R W`` and the transmitted child has
weight ``(1-R) W``. Total internal reflection creates only a reflected child
with weight ``W``. Disabling reflections discards the reflected contribution.
Direct particle combing groups spatially sorted candidates into strata, selects
one candidate proportionally to weight in each stratum, and assigns that
candidate the stratum's total weight. This preserves each candidate's expected
contribution without duplicating candidates. If the surviving population is
smaller than the number of live destination domains, global weighted combing
samples destinations instead of requiring an impossible slot per domain. The model
uses configured constant reflectivity and Snell refraction; it does not
calculate Fresnel or polarization-dependent coefficients. See
:ref:`ase-surface-reflections`.

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

``fromYaml`` accepts schema-v3 PhiASE settings under ``simulation.phi_ase``:

.. code-block:: yaml

   schema_version: 3
   simulation:
     phi_ase:
       min_rays: 100000
       max_rays: 1000000
       relative_standard_error_threshold: 0.05
       enable_diagnostics: true
       adaptive_steps: 4
       ase_steps: 150
       use_reflections: true
       reflection_mode: srm
       surface_reservoir_size: 256
       srm_position_mode: centroid
       reflection_max_iterations: 40
       boundary_max_passes: 256
       backend: Host_Cpu_CpuSerial
       openpmd_backend: auto
       parallel_mode: single
       rng_seed: 1234

.. code-block:: python

   phi_ase = PhiASE.fromYaml(
       "config/hase-phiase.yaml",
       maxRays=2_000_000,
   )

Keyword arguments override file values. Material spectra are attached later by
``Simulation``. ``addArguments`` and ``fromArgs`` add the same controls to an
``argparse`` command. Boolean pairs allow either a constructor default or a
loaded YAML value to be overridden explicitly:

.. code-block:: console

   --use-reflections | --no-reflections
   --enable-diagnostics | --disable-diagnostics
   --monochromatic | --polychromatic
   --write-vtk | --no-write-vtk

Explicit device IDs use ``--devices ID [ID ...]``. Inclusive sample bounds use
``--min-sample-range`` and ``--max-sample-range``; ``--ase-steps`` controls the
outer ASE-active step count.
The deprecated ``track_ray_visits`` YAML key and
``--track-ray-visits``/``--no-track-ray-visits`` CLI spellings remain accepted
as input aliases.
