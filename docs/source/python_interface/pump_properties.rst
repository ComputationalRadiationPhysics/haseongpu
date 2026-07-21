PICMI-style pump configuration
==============================

The public pump API separates the physical pump, numerical injection method,
and Monte Carlo solver controls. This follows the composition used by PICMI's
laser and injection objects without claiming that HASEonGPU implements the PIC
model itself.

.. code-block:: python

   from HASEonGPU import (
       MonteCarloPumpSolver,
       PlanarPumpRelay,
       Pump,
       PumpAngularDistribution,
       PumpSpectrum,
       SuperGaussianPumpProfile,
       SurfacePumpInjector,
   )

   profile = SuperGaussianPumpProfile(
       radius_u=1.5,
       radius_v=1.5,
       exponent=40,
       center=(0, 0, 0),
       axis_u=(1, 0, 0),
       axis_v=(0, 1, 0),
   )
   pump = Pump(
       total_power=100000.0,
       spectrum=PumpSpectrum.monochromatic(940e-9),
       cross_sections=pump_cross_sections,
       angular_distribution=PumpAngularDistribution.collimated(),
       profile=profile,
   )
   injector = SurfacePumpInjector(surface_domains=("pump_input",))
   pump_solver = MonteCarloPumpSolver(
       ray_count=100000,
       seed=5489,
       max_steps=50,
   )

   simulation = Simulation(
       gain_medium=medium,
       phi_ase=phi_ase,
       time_integrator=FrozenPhiAseRungeKutta4(),
       time_step_size=2e-5,
       pump_solver=pump_solver,
   )
   simulation.add_pump(
       pump,
       injection_method=injector,
       relays=(PlanarPumpRelay.retroreflect("pump_output"),),
   )
   simulation.step(150)

Physical and numerical objects
------------------------------

``Pump`` contains total power, spectrum, material cross sections, spatial
profile, and angular distribution. ``GaussianPump`` is a convenience class
that creates a super-Gaussian physical pump from a scalar or two-component
``waist``.

``SurfacePumpInjector`` selects tagged exterior triangle domains. It describes
where the physical pump is introduced, rather than duplicating those domains
inside the pump itself. Multiple pumps can be registered through repeated
``Simulation.add_pump`` calls.

``MonteCarloPumpSolver`` owns the numerical ``ray_count`` and reproducibility
``seed``. Its optional ``max_steps`` limits the pump contribution while the
outer simulation may continue.

Power normalization
-------------------

``total_power`` is integrated over the selected injector aperture. The profile
changes the launch distribution but not total power.
``integrate_pump_profile(topology, domains, profile)`` converts a legacy peak
intensity into total power on a tagged triangular aperture.

Planar relays
-------------

``PlanarPumpRelay`` maps rays leaving tagged, coplanar ``exit_domains`` to
``entry_domains``. It supports ``flip_u``, ``flip_v``, in-plane rotation,
offset, tilt, magnification, transmission, and aperture vignetting. Ordered
relays passed to ``add_pump`` represent finite return passes.

Pump polarization, coating interactions, residual cavity recirculation, and
arbitrary Python transport callbacks are outside the general pump core.
