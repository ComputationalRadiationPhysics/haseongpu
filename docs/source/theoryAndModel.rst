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
     - ``medium.get("betaVolume")`` / ``TimeStepState.beta_volume``
     - Cell-centered excited-state fraction.
   * - :math:`N_{\mathrm{tot}}`
     - ``medium.get("nTot")``
     - Active-ion concentration used in ASE and pump gain.
   * - :math:`\tau`
     - ``medium.get("crystalTFluo")``
     - Fluorescence lifetime.
   * - :math:`\sigma_a(\lambda)`, :math:`\sigma_e(\lambda)`
     - ``CrossSectionData``
     - Wavelength-dependent absorption and emission cross sections.
   * - :math:`\Phi_j`
     - ``TimeStepState.phi_ase``
     - ASE flux estimate after backend physical scaling.
   * - :math:`d\beta/dt|_{\mathrm{ASE}}`
     - ``TimeStepState.dndt_ase``
     - ASE depletion contribution.
   * - :math:`d\beta/dt|_{\mathrm{pump}}`
     - ``TimeStepState.dndt_pump``
     - Pump-induced population contribution.
   * - :math:`P`
     - ``Pump.total_power``
     - Pump power integrated over the injection aperture.
   * - boundary reflectivity and indices
     - ``SurfaceOptics``
     - Domain-assigned ASE boundary properties.

Object construction and units are documented in the :doc:`Python Interface
Guide <pythonInterface>`. This page owns the estimator equations, physical
normalization, and model limits.

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

The spectral-bin and source-cell selections are stratified within each global
batch. The stratification uses global history indices, so splitting a batch
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

Cladding cells instead use the configured cladding absorption. For a segment of
length :math:`\ell` in cell :math:`j`, the ray weight changes by

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
For :math:`N` globally launched histories, the unscaled volume estimator is

.. math::

   \Phi_j^0 = \frac{B}{N V_j}\sum_{r=1}^{N} X_{rj}.

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

For each cell, the backend accumulates the score sum :math:`S_1=\sum_r X_r`
and squared-score sum :math:`S_2=\sum_r X_r^2`. It reports the relative
standard error

.. math::

   \mathrm{RSE}
   = \sqrt{\frac{N S_2/S_1^2 - 1}{N}}.

Unlike the former mean-squared-error threshold, RSE is dimensionless: an RSE of
``0.1`` represents an estimated one-standard-error uncertainty of 10% relative
to the mean. ``standardError`` carries the same physical unit as ``phiAse``;
``relativeStandardError`` is dimensionless.

Adaptive execution starts with ``minRays`` and adds geometrically increasing
global batches until all cells satisfy the configured RSE threshold or
``maxRays`` is reached. ``forwardRayCount`` selects a fixed global history count
instead. Dropped or non-finite histories prevent the affected cell from being
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

Surface Reflections
-------------------

When reflections are enabled, direct rays still propagate to a physical mesh
boundary. Domain-assigned ``SurfaceOptics`` determines the reflected fraction.
The current model applies specular reflection with either total internal
reflection or a configured constant reflectivity.

The surface-reservoir method (SRM) retains a bounded weighted sample of the
direction, spectral bin, and weight arriving at each reflective boundary face.
The face weights define the source distribution for a new reflected pass. The
new pass contributes through the same track-length estimator and fills the
reservoir for the following pass.

Reflection terminates when the remaining source weight is below
``reflectionTolerance``, reaches a non-growing stable state, exceeds the
configured divergence streak, or reaches ``reflectionMaxIterations``. The
result exposes the SRM status, pass count, remaining fraction, and the active
safety limits.

The ASE boundary model does not yet launch a transmitted or refracted ray and
does not calculate angle- or polarization-dependent Fresnel coefficients. A
non-reflected fraction leaves the simulated domain. Refractive indices are
currently used to detect total internal reflection only.

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
stages, while still evaluating the pump contribution. Setting
``Simulation(enable_ase=False, ...)`` advances pump excitation and fluorescence
without an ASE calculation.

Model Limits
------------

The current transport model does not include polarization, Fresnel
transmission, detailed coating stacks, general internal optical interfaces,
unlimited pump-cavity recirculation, or custom Python transport callbacks
inside compiled time steps. Volume transport supports Tet4 cells. Runtime and
available device memory limit practical ray and surface-reservoir counts.

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
