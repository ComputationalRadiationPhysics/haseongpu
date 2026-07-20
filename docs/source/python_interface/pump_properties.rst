PumpProperties
==============

The compiled pump is a boundary-launched Monte Carlo transport model. Pump
power is normalized explicitly: every source supplies total power, a normalized
spatial profile, a normalized discrete spectrum, and a normalized angular
distribution.

A monochromatic collimated source is configured as follows:

.. code-block:: python

   from HASEonGPU import (
       PumpAngularDistribution, PumpProperties, PumpSource, PumpSpectrum,
       SuperGaussianPumpProfile,
   )

   profile = SuperGaussianPumpProfile(
       radiusU=1.5, radiusV=1.5, exponent=40,
       center=(0, 0, 0), axisU=(1, 0, 0), axisV=(0, 1, 0),
   )
   source = PumpSource(
       surfaceDomains=("pump_input",),
       totalPower=100000.0,
       spectrum=PumpSpectrum.monochromatic(940e-9),
       crossSections=pump_cross_sections,
       angularDistribution=PumpAngularDistribution.collimated(),
       profile=profile,
   )
   pump = PumpProperties(
       sources=(source,), rayCount=100000, rngSeed=5489, pumpSteps=50,
   )

``surfaceDomains`` are tagged exterior triangle domains. Multiple independent
sources may be supplied. ``UniformPumpProfile`` distributes power uniformly;
``SuperGaussianPumpProfile`` defines a world-space elliptical profile.
``PumpAngularDistribution.uniformCone`` and an explicit ``PumpSpectrum`` enable
finite-divergence and multi-wavelength pumps. Cross sections are interpolated
at every discrete pump wavelength before transport.

Power normalization
-------------------

``totalPower`` is the integral over the selected source aperture. The profile
only changes where that power is launched; it does not change total power.
``integratePumpProfile(topology, domains, profile)`` can convert a legacy peak
intensity to total power on a tagged triangular aperture.

Planar relays and passes
------------------------

``PlanarPumpRelay`` maps rays leaving tagged, coplanar exit faces to tagged,
coplanar entry faces. It supports flips, in-plane rotation, offset, tilt,
magnification, and transmission. Relays are evaluated in source order and may
represent return passes. ``PlanarPumpRelay.retroreflect("pump_output")`` is the
unit-magnification same-aperture return used to approximate the former
back-reflected z traversal.

Rays that miss the mapped entry aperture are vignetted. Coating physics,
polarization, residual cavity recirculation, and arbitrary Python callbacks are
not part of this core model.

Sampling controls
-----------------

``rayCount`` is the number of equal-power launch rays per source and is
independent of PhiASE ray counts. ``rngSeed`` makes host-side source sampling
reproducible. ``pumpSteps`` limits pump contribution to the first outer time
steps; it may also be overridden by ``Simulation.runSteps(..., pumpSteps=...)``.
