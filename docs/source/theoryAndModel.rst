Theory and Model
================

HASEonGPU builds on earlier ASE modeling work by D. Albach et al. [2],
where ray-tracing techniques and Monte Carlo integration were used to calculate
ASE in laser gain media in a single-threaded CPU-centered context. Based on
this scientific foundation, HASEonGPU [1] extends the approach with portable
Alpaka execution, adaptive sampling, and distributed multi-node execution.

Nature of the Problem
---------------------

Accurate ASE simulations require a spatially resolved estimate of radiation
transport through an inhomogeneous gain medium. A useful solver must sample the
emission spectrum and source volume, integrate gain and absorption along many
paths, and optionally account for repeated surface reflections. The cost grows
with the mesh size and requested statistical accuracy, which makes the
transport algorithm and its mapping to parallel hardware central to practical
simulations.

Frontend quantities and model symbols
-------------------------------------

The Python objects name the same quantities used in the equations below. This
mapping is the bridge between a physical model and a runnable input:

.. list-table::
   :header-rows: 1
   :widths: 18 32 50

   * - Symbol
     - Python field or object
     - Meaning
   * - :math:`\beta_j`
     - ``Simulation.initialExcitation`` / ``TimeStepState.betaVolume``
     - Cell-centered excited-state fraction, evaluated only in components
       whose material has ``active=True`` and which belong to ``GainMedium``.
   * - :math:`N_{\mathrm{tot}}`
     - ``component.material.activeIonDensity``
     - Active-ion concentration used in ASE and pump gain.
   * - :math:`\tau`
     - ``component.material.fluorescenceLifetime``
     - Fluorescence lifetime.
   * - :math:`\sigma_a(\lambda)`, :math:`\sigma_e(\lambda)`
     - ``component.material.crossSections``
     - Wavelength-dependent absorption and emission cross sections.
   * - :math:`\alpha_{\mathrm{bulk}}`
     - ``passiveComponent.material.bulkAttenuation``
     - Constant passive intensity attenuation coefficient.
   * - :math:`\Phi_j`
     - ``TimeStepState.phiAse``
     - ASE flux estimate after backend physical scaling.
   * - :math:`d\beta/dt|_{\mathrm{ASE}}`
     - ``TimeStepState.dndtAse``
     - ASE depletion contribution.
   * - :math:`d\beta/dt|_{\mathrm{pump}}`
     - ``TimeStepState.dndtPump``
     - Pump-induced population contribution.
   * - :math:`P`
     - ``Pump.total_power``
     - Pump power integrated over the injection aperture.
   * - boundary reflectivity and indices
     - ``component.surfaceOptics`` / ``SurfaceOptics``
     - Domain-assigned ASE boundary properties.

Object construction and units are documented in the :doc:`Python Interface
Guide <pythonInterface>`. This page owns the estimator equations, physical
normalization, and model limits.

The frontend obtains :math:`N_{\mathrm{tot}}`, :math:`\tau`, and both cross
sections from the resolved gain ``Material``. Lowering converts these
unit-bearing values to the backend fields used in the equations. Pump and ASE
therefore cannot select inconsistent material spectra in a high-level
``Simulation``.

.. _forward-ase-model:

Forward ASE Model
-----------------

HASEonGPU 2.2.0 uses a source-driven forward Monte Carlo estimator on explicit
Tet4 volume meshes. This replaces the former target-driven algorithm, which
traced a separate population of source-to-target paths for every observation
point. A forward history starts at an emitting volume and contributes to every
cell it traverses. Geometric work performed for one history is therefore reused
across all cells on that path.

Tet4 State and Emission Source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`V_j` be the volume of Tet4 cell :math:`j` and :math:`\beta_j` its
cell-centered excited-state fraction. The compiled model has no
vertex-centered excitation state.

The integrated source strength used for spatial sampling is

.. math::

   B = \sum_j \beta_j V_j.

A direct history selects cell :math:`j` with probability
:math:`\beta_j V_j / B`, samples a point uniformly inside that tetrahedron,
selects a spectral bin, and samples an isotropic direction. Sampling the source
with its physical :math:`\beta V` density absorbs that factor into the source
probability and leaves a unit importance weight. If :math:`B` is zero, the ASE
estimate is zero.

The spectral-bin and source-cell selections are stratified within each domain
and statistical batch. An independently keyed permutation of wavelength strata
prevents the two dimensions from sharing a monotone history ordering.
The stratification uses global history indices, so splitting a batch
over Alpaka devices or MPI ranks preserves the intended coverage. Directions
remain isotropic; HASEonGPU does not infer a preferred direction from an
arbitrary Tet4 mesh.

Gain Along a Ray
^^^^^^^^^^^^^^^^

Within a cell, the local gain coefficient for wavelength-dependent absorption
and emission cross sections is

.. math::

   g_j = N_{\mathrm{tot}}
         \left[\beta_j(\sigma_e + \sigma_a) - \sigma_a\right].

Cells belonging to a passive optical component instead use

.. math::

   g_j = -\alpha_{\mathrm{bulk}}.

The tracing kernels read this coefficient directly from the receiving cell's
material-owned ``bulkAttenuation``. The material coefficient is constant and
wavelength-independent in the current model; it
is distinct from the active-ion absorption cross section
:math:`\sigma_a(\lambda)`. An omitted coefficient contributes
:math:`\alpha_{\mathrm{bulk}}=0`; surface reflection remains a separate
``SurfaceOptics`` interaction. For a segment of length :math:`\ell` in cell
:math:`j`, the ray weight changes by

.. math::

   G_{\mathrm{out}} = G_{\mathrm{in}}\exp(g_j\ell).

The contribution deposited in that cell uses the exact gain-weighted
track-length integral

.. math::

   T(g_j, \ell)
   = \int_0^\ell \exp(g_j s)\,ds
   = \begin{cases}
       \dfrac{\exp(g_j\ell)-1}{g_j}, & g_j \ne 0,\\
       \ell, & g_j = 0.
     \end{cases}

The implementation evaluates the :math:`\ell` limit directly when
:math:`|g_j\ell|` is small. The score contributed by a history to a visited cell
is its gain on entering the cell multiplied by :math:`T(g_j,\ell)`. A history
that does not visit the cell contributes zero.

Direct Estimator and Physical Scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`X_{rj}` be the score deposited by direct history :math:`r` in cell
:math:`j`, including a zero score when the history does not visit that cell.
For :math:`M` active independent statistical batches with :math:`n_b` histories
each, the unscaled volume estimator is

.. math::

   Y_{bj} = \frac{1}{n_b}\sum_{r\in b} X_{rj}, \qquad
   \Phi_j^0 = \frac{B}{V_j}\frac{1}{M}\sum_{b=1}^{M}Y_{bj}.

Every batch represents the complete source. A history in domain :math:`d`
receives importance weight :math:`n_b B_d/(n_{db}B)`, where :math:`n_{db}` is
that domain's history count in the batch. Unequal batch sizes use an equal-weight
mean of batch estimates, matching the uncertainty calculation below.

Uniform sampling over the sphere already represents the angular
:math:`1/(4\pi)` average; no additional inverse-square target factor is applied.
This is a track-length estimator of the radiation field in a receiving volume,
not the old point-target estimator.

The high-level backend scales the result by active-ion density and fluorescence
lifetime,

.. math::

   \Phi_j = \frac{N_{\mathrm{tot}}}{\tau}\Phi_j^0.

``phiAse`` therefore already contains the :math:`N_{\mathrm{tot}}/\tau`
factor. It must not be applied a second time in the population derivative.

Statistical Uncertainty and Adaptation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The backend retains raw score sums and history counts per independent batch.
For each cell it calculates the sample variance of normalized batch estimates
using centered online moments and reports

.. math::

   \bar Y = \frac{1}{M}\sum_b Y_b, \qquad
   s_Y^2 = \frac{1}{M-1}\sum_b (Y_b-\bar Y)^2, \qquad
   \mathrm{RSE} = \frac{\sqrt{s_Y^2/M}}{|\bar Y|}.

At least two active batches are required; otherwise RSE is the maximum error
sentinel. A zero mean has undefined RSE (NaN). Statistical batches, not
correlated rays within a stratified or combed population, are the replicates.

Unlike the former mean-squared-error threshold, RSE is dimensionless: an RSE of
``0.1`` represents an estimated one-standard-error uncertainty of 10% relative
to the mean. ``standardError`` carries the same physical unit as ``phiAse``;
``relativeStandardError`` is dimensionless.

Adaptive execution starts with ``minRays`` and adds geometrically increasing
global batches until all cells satisfy the configured RSE threshold or
``maxRays`` is reached. ``forwardRayCount`` selects a fixed global history count
instead. ``numIndependentRayPopulations`` selects the runtime population count
(default 8), independently of workers. Every positive-source domain must occur
in every population; insufficient quotas are rejected. Adaptive increments
are coalesced when necessary to preserve complete sources and exact final
quotas. With ``adaptiveSteps=0``, evaluation stops after ``minRays``.
Dropped or non-finite histories prevent the affected cell from being
reported as converged. RSE measures Monte Carlo sampling uncertainty; it does
not include mesh discretization or model error.

Tet4 Traversal
^^^^^^^^^^^^^^

For each cell face, HASEonGPU precomputes the affine barycentric coordinate

.. math::

   \lambda_f(x) = a_f x + b_f y + c_f z + d_f.

Along a ray :math:`x(t)=x_0+t v`, the next decreasing coordinate to reach zero
identifies the exit face. Cell adjacency then provides the next tetrahedron and
the local face that must be excluded from the following intersection. This
avoids rebuilding triangle-intersection data in the hot traversal loop.

Bounded recovery handles rays that meet several faces at a shared edge or
vertex. Invalid connectivity, non-finite contributions, or a traversal that
cannot be recovered are counted as dropped histories rather than silently
contributing an invalid value.

.. _ase-surface-reflections:

Domain-boundary transport
-------------------------

Every optical component is an ASE scheduling domain. The global primary-ray
count is divided among domains from their integrated spontaneous-source
strength, except for exact counts reserved by ``OpticalComponent.aseRays``.
Source CDFs restart at zero in each domain. Their totals are local endpoints,
not differences between large global prefixes. A zero-source domain receives no
automatic primary allocation but remains a receiving transport domain.
Workers receive ``(rayPopulationId, domainId, batchId)`` work items. SRM logical
batch boundaries are independent of the number of workers; direct execution
chunks may vary with worker count because its resampling remains population-wide.

The current execution context still keeps the complete prepared trace on every
worker device. Source sampling, CDFs, scans, and reservoir operations execute
on that device. SRM retains boundary records on the owner of each logical batch
and downloads only aggregate weights for pass control. Direct transport gathers
compact boundary records in canonical order through host memory for
population-wide combing, then distributes the selected records to workers.
Both paths download combined raw scores for population aggregation.

At a domain boundary, configured constant reflectivity :math:`R` splits incoming
weight :math:`W` deterministically. The specularly reflected child receives
:math:`R W`, and the Snell-refracted transmitted child receives
:math:`(1-R)W`. Total internal reflection produces only a reflected child with
weight :math:`W`. Disabling reflections discards the reflected contribution but
does not change the transmitted child's weight.

The direct policy stores both boundary children in persistent
structure-of-arrays device buffers. A per-domain systematic particle comb
restores the requested population before the next pass. The comb preserves the
sum of candidate weights, transfers the weight represented by discarded
candidates to retained histories, and duplicates histories when too few
candidates reach a domain. Thus population control changes the number of
histories, not the represented boundary weight. Positive boundary routes receive
relaunch slots when the surviving population can represent them all; otherwise
a global comb samples the routes without deterministically dropping a domain.

The surface-reservoir method (SRM) retains a bounded weighted sample of the
direction, spectral bin, and weight from both children arriving at their target
faces. Each source domain's independent ray population is divided into logical
batches of at most 65,536 primary rays, independently of worker count. Source
and wavelength sampling is prepared for the complete domain population before
tracing; splitting the prepared rays does not restart its strata. Each logical
batch executes its complete SRM transport on one owner, including transmission
into other domains, and retains separate reservoirs and stopping state. The
next pass samples the batch's occupied target faces in proportion to aggregate
face weight, including faces in domains with no primary source.

The new pass contributes through the same track-length estimator and fills the
batch's reservoir for the following pass. Reservoir randomness is keyed by
ray population, source domain, logical batch, and pass, never by worker identity.
Raw scores from all logical batches in a population are summed before
normalization. Only the complete population estimates enter the RSE calculation;
logical batches are not additional independent RSE samples. Changing worker
ownership preserves this model, apart from floating-point reduction order.
Changing the logical batch cap can change its sampling variance and effective
reservoir capacity. This cap is not an automatic worker-load tuning parameter.

Boundary propagation converges when the remaining source weight is zero or its
fraction of the initial boundary source weight is below
``reflectionTolerance``. Otherwise, tracing continues until the configured
divergence streak or ``boundaryMaxPasses`` is reached. Small changes between
successive pass weights do not establish convergence. The result exposes
boundary status, pass count, remaining fraction, and active safety limits for
either boundary policy.

A pass-limited exit is a truncated Neumann series, not by itself a
finite steady field. HASE fits the residual reflected weights to a per-pass
multiplier :math:`\Gamma`, but retains the candidate tail factor only as a
diagnostic. Even perfect scalar decay does not establish stationarity of the
spatial or spectral distribution, so the final-pass field is not extrapolated.
A confidently supercritical multiplier reports ``diverged``.
``maxPasses`` remains an unresolved partial tally and is rejected by material
integration, as are dropped histories and non-finite ASE outputs.

The ASE boundary model does not calculate angle- or polarization-dependent
Fresnel coefficients. Transmission is available only across prepared,
conforming component interfaces; the transmitted fraction at an exterior
surface leaves the simulated assembly.

ASE Population Derivative
-------------------------

The ASE depletion term for the excited-state fraction is

.. math::

   \left.\frac{d\beta_j}{dt}\right|_{\mathrm{ASE}}
   = \left[\beta_j(\sigma_e+\sigma_a)-\sigma_a\right]\Phi_j.

The compiled simulation evaluates and integrates this term directly on Tet4
cells. It does not maintain a point-centered beta or PhiASE representation.

.. _general-monte-carlo-pump:

General Monte Carlo Pump
------------------------

The compiled pump no longer assumes a super-Gaussian profile propagated only
through ordered z levels. It launches equal-power Monte Carlo rays from tagged
exterior Tet4 faces. Each source defines an aperture-integrated total power, a
normalized spatial profile, a discrete wavelength spectrum, and an angular
distribution.

Pump entry positions use randomized systematic stratification over spatial
subregions of the injector aperture. The subregion CDF is weighted by the
integrated pump profile, rather than only by boundary-face area, and sampling
within a selected subregion retains the configured continuous profile. Global
ray indices preserve the same coverage when a ray batch is partitioned.

Within a cell, pump power follows

.. math::

   P_{\mathrm{out}} = P_{\mathrm{in}}\exp(g_p\ell),
   \qquad
   g_p = N_{\mathrm{tot}}
         \left[\beta(\sigma_a+\sigma_e)-\sigma_a\right].

The corresponding net photon exchange is distributed barycentrically to the
Tet4 vertices, normalized with lumped vertex volumes, and averaged back to the
cells. This temporary projection smooths the pump rate; it does not introduce
an evolving point-centered beta field. The time integrator receives one
:math:`d\beta_j/dt` value per cell.

``SurfacePumpInjector`` selects the tagged launch faces. Finite
``PlanarPumpRelay`` stages can map rays from coplanar exit domains to entry
domains with flips, rotation, offset, tilt, magnification, scalar transmission,
and aperture vignetting. These relays are explicit affine return paths, not a
Fresnel or unlimited cavity model.

.. _pump-and-time-stepping:

Pump and Time Stepping
----------------------

The compiled C++/Alpaka simulation advances one authoritative beta value per
Tet4 cell:

.. math::

   \frac{d\beta}{dt}
   = \left.\frac{d\beta}{dt}\right|_{\mathrm{pump}}
     - \left.\frac{d\beta}{dt}\right|_{\mathrm{ASE}}
     - \frac{\beta}{\tau}.

Standard RK4 reevaluates ASE and pump transport at every stage.
``FrozenPhiAseRungeKutta4`` reuses its first ASE calculation for the remaining
stages, while still evaluating the pump contribution.
``FrozenSourcesRungeKutta4`` evaluates pump and ASE transport once at the
pre-step beta field, then holds both resulting rate fields fixed during all
four RK4 stages. Only the fluorescence term :math:`-\beta/\tau` changes with
the intermediate stage field. Consequently, a snapshot from this solver pairs
the post-step beta field with the pump and ASE rate fields frozen at the
pre-step beta field. Setting
``PhiASE(ase_steps=0, ...)`` advances active pump excitation and fluorescence
without an ASE calculation. Setting one pump's ``pump_steps`` to zero disables
only that source.

Model Limits
------------

The current transport model does not include polarization, Fresnel
coefficients, detailed coating stacks, non-conforming internal optical interfaces,
unlimited pump-cavity recirculation, or custom Python transport callbacks
inside compiled time steps. Volume transport supports Tet4 cells. Runtime and
available device memory limit practical ray and boundary-buffer counts.

References
----------

[1] C.H.J. Eckert, E. Zenker, M. Bussmann, and D. Albach,
    *HASEonGPU-An adaptive, load-balanced MPI/GPU-code for calculating the
    amplified spontaneous emission in high power laser media*,
    Computer Physics Communications, 207, 2016, 362-374.
    DOI: `10.1016/j.cpc.2016.05.019`

[2] D. Albach, J.-C. Chanteloup, G. l. Touze,
    *Influence of ASE on the gain distribution in large size, high gain
    Yb3+:YAG slabs*,
    Optics Express, 17(5), 2009, 3792-3801.
