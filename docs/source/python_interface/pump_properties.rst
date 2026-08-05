Pump configuration
==================

The pump API separates four concerns that are easy to conflate:

* ``Pump`` describes a physical source;
* ``SurfacePumpInjector`` selects the boundary through which it enters;
* ``PlanarPumpRelay`` describes an optional finite return path;
* ``MonteCarloPumpSolver`` controls numerical ray sampling and its seed.

This composition lets several physical pumps share one numerical solver while
using different spectra, profiles, directions, and launch surfaces.

Physical source
---------------

``Pump`` combines aperture-integrated power with wavelength, material, spatial,
and angular distributions:

.. code-block:: python

   profile = SuperGaussianPumpProfile(
       radius_u=1.5,
       radius_v=1.5,
       exponent=40,
   )
   pump = Pump(
       total_power=16_000.0,
       spectrum=PumpSpectrum.monochromatic(940e-9),
       cross_sections=pump_cross_sections,
       angular_distribution=PumpAngularDistribution.collimated(),
       profile=profile,
   )

``PumpSpectrum`` is a normalized discrete wavelength distribution. The pump's
``cross_sections`` must describe absorption and emission at those wavelengths;
they are distinct from the full spectrum commonly supplied to ``PhiASE``. See
:doc:`spectral_decomposition`.

``PumpAngularDistribution`` stores normalized polar/azimuthal direction samples
in the injector's inward-local frame. ``collimated`` supplies one normal-incidence
direction; ``uniform_cone`` creates a discrete cone.

``UniformPumpProfile`` weights every position equally.
``SuperGaussianPumpProfile`` evaluates an elliptical world-space profile from
its center, two orthogonal axes, radii, and exponent. ``GaussianPump`` is a
convenience ``Pump`` subclass that constructs such a profile from ``waist``.

Injection surfaces
------------------

An injector refers to one or more named or numeric exterior surface domains on
the medium's topology:

.. code-block:: python

   injector = SurfacePumpInjector("pump_input")
   simulation.add_pump(pump, injection_method=injector)

The domain assignment belongs to ``VolumeTopology``; the injector gives that
geometry a pump-specific role. Repeated ``add_pump`` calls register multiple
sources before the simulation is initialized.

Numerical sampling
------------------

The solver controls are shared by all registered pump sources:

.. code-block:: python

   pump_solver = MonteCarloPumpSolver(
       ray_count=50_000,
       seed=5489,
       max_steps=100,
   )

``ray_count`` is the number of equal-power histories launched for each source
evaluation. ``seed`` makes pump sampling reproducible. ``max_steps`` limits the
initial outer simulation steps that include pump transport; it is not a limit
on how many Tet4 cells one ray may traverse. ``Simulation.step`` independently
selects the total number of outer time steps.

Power normalization
-------------------

``Pump.total_power`` is the power integrated over the selected aperture. The
profile changes where rays launch but is normalized so it does not change that
total.

When adapting a legacy peak-power-density value, integrate the dimensionless
profile over the actual tagged surface before constructing the pump:

.. code-block:: python

   profile_area = integrate_pump_profile(
       topology,
       "pump_input",
       profile,
   )
   pump = Pump(
       total_power=peak_power_density * profile_area,
       spectrum=pump_spectrum,
       cross_sections=pump_cross_sections,
       profile=profile,
   )

The surface area uses the topology's geometry unit, so the power-density input
must use the corresponding inverse-area unit.

Finite planar relays
--------------------

``PlanarPumpRelay`` maps rays that leave coplanar exit domains back to entry
domains. A simple return from the same aperture is:

.. code-block:: python

   relay = PlanarPumpRelay.retroreflect(
       "pump_output",
       transmission=0.95,
   )
   simulation.add_pump(
       pump,
       injection_method=injector,
       relays=(relay,),
   )

General relays support flips, in-plane rotation, offset, tilt, magnification,
scalar transmission, and aperture vignetting. The ordered tuple represents a
finite number of explicit passes. It does not create an unlimited cavity.

Transport model
---------------

The compiled solver samples launch position, wavelength, and direction, then
traverses neighboring Tet4 cells to a physical boundary. Material interaction
uses the local excited-state fraction and the pump wavelength's cross sections.
The deposited pump rate is temporarily projected through Tet4 vertices for
smoothing and returned to one value per cell; the evolving excitation remains
cell-centered.

The gain equation, deposited population rate, temporary projection, and model
limits are specified once in :ref:`general-monte-carlo-pump`. Relay transmission
is a configured scalar: the current pump model does not calculate polarization,
coating response, refraction, or Fresnel transmission.
