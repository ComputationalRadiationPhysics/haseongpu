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
multiplied by :math:`N_{\mathrm{tot}}` or divided by :math:`\tau` again.  In the compiled ``Simulation`` time loop, the cross sections used for this scalar
conversion are the maximum emission cross section and the absorption cross
section at the emission-peak wavelength.

Pump and Time Stepping
----------------------

The high-level Python interface constructs the physical state, then sends it to
the C++/Alpaka backend, which advances ``betaCells``.  It combines pump
excitation, ASE depletion, and fluorescence decay:

.. math::

   \frac{d\beta}{dt}
   =
   \left.\frac{d\beta}{dt}\right|_{\mathrm{pump}}
   -
   \left.\frac{d\beta}{dt}\right|_{\mathrm{ASE}}
   - \frac{\beta}{\tau}.

The compiled pump is ``one-dimensional-z-traversal``. At each transverse mesh
point it uses the super-Gaussian inlet profile

.. math::

   I_0(x,y) = I_{\mathrm{peak}}
   \exp\left[-\left(\sqrt{\frac{x^2}{r_y^2} + \frac{y^2}{r_x^2}}\right)^q\right].

For the configured pump wavelength, it propagates intensity through adjacent z
levels using the mean population :math:`\bar\beta_z`:

.. math::

   I_{z+\Delta z} = I_z
   \exp\left[-\left(\sigma_a
   - \bar\beta_z(\sigma_a+\sigma_e)\right)
   N_{\mathrm{tot}}\Delta z\right].

The local pump rate is

.. math::

   \left.\frac{d\beta}{dt}\right|_{\mathrm{pump}}
   = \left[\sigma_a - \beta(\sigma_a+\sigma_e)\right]
   I\frac{\lambda}{hc}.

With back reflection, the transmitted far-end intensity is scaled by the chosen
reflectivity and propagated in the reverse direction before both intensities
are added. ``extraction`` suppresses the inlet pump. The configured
``pumpSubsteps`` value is retained in the transport format for compatibility;
the current compiled traversal evaluates this rate directly and does not use it
to subdivide a time step.

The compiled integrator then advances the combined derivative. Standard RK4
re-evaluates ASE at each stage; ``FrozenPhiAseRungeKutta4`` reuses its first ASE
calculation for the remaining stages when that approximation is appropriate.
``Simulation(enableASE=False, ...)`` disables ASE and advances only pump and
fluorescence. Custom Python pump or time-integration functions cannot run
inside the compiled loop.

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
