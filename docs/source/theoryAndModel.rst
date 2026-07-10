Theory and Model
================

HASEonGPU builds on earlier ASE modeling work by D. Albach et. al [2],
where ray-tracing techniques and Monte Carlo integration were used to calculate
ASE in laser gain media in a single-threaded CPU-centered context. Based on
this scientific foundation, HASEonGPU [1] extends the approach with multi-GPU
acceleration, adaptive sampling, and distributed multi-node execution.

Nature of the Problem
---------------------

Accurate ASE simulations require a large number of sampling points as well as a
high number of rays. In addition, reflections on the upper and lower surfaces of
the gain medium can substantially increase the computational workload. On
traditional CPU-based systems, this makes detailed ASE simulations
time-consuming and limits the number of practical simulation runs.

Solution Method
---------------

HASEonGPU addresses this with non-uniform spatial sampling, Monte Carlo
integration, importance sampling, adaptive sampling, and GPU parallelization
[1].  MPI execution distributes the sample range across ranks; each rank can
use one or more local devices.

ASE Estimator
-------------

The central calculation is the ASE contribution at each sample point of the
gain medium.  Let

.. math::

   p_i

be sample point ``i``.  The mesh is represented as an extruded triangular mesh.
For a prism index

.. math::

   j = t + z N_{\mathrm{tri}},

:math:`t` is the 2D triangle index, :math:`z` is the layer index, and
:math:`N_{\mathrm{tri}}` is the number of base triangles.  The prism volume is

.. math::

   V_j = A_t \Delta z,

where :math:`A_t` is the triangle area and :math:`\Delta z` is the layer thickness.

For a source point ``q`` inside prism ``j``, HASEonGPU propagates a ray from
``q`` to ``p_i``.  Along that ray the path is split into prism segments
``s``.  The serial implementation evaluates the gain as the product

.. math::

   G(q \rightarrow p_i)
   =
   \frac{1}{\lVert p_i - q \rVert^2}
   \prod_s
   \exp\left(
   N_{\mathrm{tot}}
   \left[
   \beta_s(\sigma_e + \sigma_a) - \sigma_a
   \right]
   \Delta l_s
   \right).

Equivalently, in continuous notation this is

.. math::

   G(q \rightarrow p_i)
   =
   \frac{1}{\lVert p_i - q \rVert^2}
   \exp\left(
   \int_q^{p_i}
   N_{\mathrm{tot}}
   \left[
   \beta(x)(\sigma_e + \sigma_a) - \sigma_a
   \right]
   \,dl
   \right).

The explicit inverse-square factor is applied in the implementation as
``1 / ray.length^2`` after the segment product has been accumulated.

Importance sampling is used only to decide where rays should start.  For each
sample point, HASEonGPU first evaluates a raw importance from prism centers:

.. math::

   I_j^{\mathrm{raw}}
   =
   \beta_j G(c_j \rightarrow p_i),
   \qquad
   S = \sum_j I_j^{\mathrm{raw}}.

The requested ray count for a prism is assigned approximately proportional to
this raw importance:

.. math::

   N_j =
   \left\lfloor
   \frac{I_j^{\mathrm{raw}}}{S} N
   \right\rfloor,

where :math:`N` is ``ctx.NumRays``.  Any leftover rays caused by integer rounding
are distributed randomly over the prisms.  If ``N_j`` rays are drawn from prism
``j``, then ``N_j / N`` approximates the sampling probability.  The final
Monte Carlo weight is not the raw importance.  It is the prism volume corrected
by the allocated ray count:

.. math::

   W_j = \frac{V_j}{P_j}
   \approx
   \frac{N V_j}{N_j}.

In the serial code this is the assignment

.. math::

   W_j =
   \begin{cases}
   \dfrac{N A_t \Delta z}{N_j}, & N_j > 0,\\
   0, & N_j = 0,
   \end{cases}

with :math:`A_t` = ``surface_new[t]`` and :math:`\Delta z` = ``z_mesh``.

The unscaled ASE flux estimator at sample point ``i`` is therefore

.. math::

   \Phi_i^{0}
   =
   \frac{1}{4\pi N_{\mathrm{real}}}
   \sum_j
   \sum_{k=1}^{N_j}
   \beta_j
   G(q_{j,k} \rightarrow p_i)
   W_j,

where :math:`N_{\mathrm{real}}` is the actual number of rays used after integer ray allocation.
For the serial implementation, :math:`N_{\mathrm{real}}` is the accumulated sum of allocated
:math:`N_j` values for that sample.

The value returned by the high-level HASEonGPU backend as ``phiAse`` is already
scaled by active-ion density and fluorescence lifetime:

.. math::

   \Phi_i
   =
   \frac{N_{\mathrm{tot}}}{\tau}
   \Phi_i^{0}.

This convention matters for time stepping: ``phiAse`` already contains the
:math:`N_{\mathrm{tot}} / \tau` factor.

ASE Population Derivative
-------------------------

The ASE depletion term for the excited-state fraction is computed from the
returned ``phiAse`` as

.. math::

   \left.\frac{d\beta_i}{dt}\right|_{\mathrm{ASE}}
   =
   \left[
   \beta_i(\sigma_e + \sigma_a) - \sigma_a
   \right]
   \Phi_i.

Because :math:`\Phi_i` already includes :math:`N_{\mathrm{tot}} / \tau`, this derivative must not be
multiplied by :math:`N_{\mathrm{tot}}` or divided by :math:`\tau` again.  In the
Python ``Simulation`` time loop, the cross sections used for this scalar
conversion are the maximum emission cross section and the absorption cross
section at the emission-peak wavelength.

Pump and Time Stepping
----------------------

The ASE calculation can be used directly, but the Python interface also
provides a time-stepped pump model.  In that model ``betaCells`` is advanced by
combining pump excitation, ASE depletion, and fluorescence decay:

.. math::

   \frac{d\beta}{dt}
   =
   \left.\frac{d\beta}{dt}\right|_{\mathrm{pump}}
   -
   \left.\frac{d\beta}{dt}\right|_{\mathrm{ASE}}
   -
   \frac{\beta}{\tau}.

HASEonGPU provides two built-in pump solvers.  Both use the same local pump
physics but package the time update differently.  The transverse input field is
represented as a super-Gaussian profile,

.. math::

   I_0(x,y) = I_{\mathrm{peak}}
   \exp\left[-\left(\sqrt{\frac{(x-x_c)^2}{r_y^2}
   + \frac{(y-y_c)^2}{r_x^2}}\right)^q\right].

The continuous ``OneDimensionalZTraversal`` solver propagates each wavelength
sample along the z levels.  Between two adjacent levels it uses the average
population :math:`\bar\beta_z` and the active-ion density
:math:`N_{\mathrm{tot}}`:

.. math::

   I_{k,z+\Delta z} = I_{k,z}
   \exp\left[-\left(\sigma_{a,k}
   - \bar\beta_z(\sigma_{a,k}+\sigma_{e,k})\right)
   N_{\mathrm{tot}}\Delta z\right].

The local photon flux is

.. math::

   \Phi_k = I_k\frac{\lambda_k}{hc},

and the frozen-state pump rate is

.. math::

   \left.\frac{d\beta}{dt}\right|_{\mathrm{pump}}
   =
   \sum_k\left[\sigma_{a,k}
   - \beta(\sigma_{a,k}+\sigma_{e,k})\right]\Phi_k.

If back reflection is enabled, the transmitted intensity at the far crystal
surface is multiplied by the configured reflectivity and propagated back through
the same layer factors before the forward and reflected intensities are summed.

The legacy/default ``BetaIntegrationGaussianSolver`` instead advances beta with
a local analytical update over internal pump substeps:

.. math::

   \beta(t + \Delta t)
   =
   \frac{A}{C}\left(1 - e^{-C\Delta t}\right)
   +
   \beta(t)e^{-C\Delta t},

with

.. math::

   A = \sigma_a I \frac{\lambda}{hc},
   \qquad
   C = (\sigma_a + \sigma_e) I \frac{\lambda}{hc} + \frac{1}{\tau}.

The pump model is deliberately separate from the ASE estimator. Python workflows
can replace it with a custom pump solver while keeping the same HASEonGPU ASE
calculation.

Restrictions
------------

The number of rays used for the Monte Carlo integration of a single sampling
point is limited by the available GPU memory.

User-Relevant Features
----------------------

HASEonGPU can run on a single workstation or on GPU clusters with MPI.  The
ASE setup can include polychromatic spectra, cladding, surface coatings,
refractive-index changes, upper/lower-surface reflections, and adaptive
sampling up to the configured ray limit.

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
