PumpProperties
==============

``PumpProperties`` defines the pump contribution to a compiled
``Simulation``. The Python object is a compact physical description; the
C++/Alpaka backend performs the one-dimensional z traversal and beta update.

.. code-block:: python

   from HASEonGPU import PumpProperties

   pump = PumpProperties(
       spectralProperties=spectra,
       intensity=16e3,       # W / cm^2
       wavelength=940e-9,    # m
       radiusX=1.5,
       radiusY=1.5,
       exponent=40,
       pumpSubsteps=100,
       backReflection=True,
       reflectivity=1.0,
   )

Required Inputs
---------------

``intensity`` is the incident intensity in ``W / cm^2``. Supply
``spectralProperties`` (or ``crossSections``) for the absorption and emission
cross sections. With spectra, the wavelength defaults to the first absorption
wavelength; otherwise provide a wavelength and monochromatic absorption and
emission cross sections. ``radiusX`` is required and ``radiusY`` defaults to
``radiusX``.

The transverse inlet profile is a super-Gaussian,

.. math::

   I(x,y) = I_0 \exp[-r^p], \qquad
   r = \sqrt{x^2/r_y^2 + y^2/r_x^2},

where ``exponent`` is :math:`p` (default 40). ``intensityAt(points)`` evaluates
this profile in Python for diagnostics.

Compiled Pump Model
-------------------

The ``one-dimensional-z-traversal`` routine evaluates the inlet profile at
each transverse point and propagates it through the z levels using the current
excited-state fraction. It exposes the resulting local pump rate to the compiled
time integrator. With ``backReflection=True``, a second pass starts at the far
end with intensity scaled by ``reflectivity``. ``extraction=True`` suppresses
the inlet pump.

``pumpSteps`` limits the number of *outer* simulation steps that receive pump
energy; set it on the object or pass it to ``Simulation.runSteps``. The retained
``pumpSubsteps`` setting is currently ignored by the compiled traversal and does
not count outer simulation steps.

Limits and Utilities
--------------------

Compiled simulations accept the built-in routine only. A ``solver`` custom
property is retained to detect legacy configurations, but causes a clear error
before launch rather than running Python code between backend steps.

Use ``PumpProperties.superGaussian(...)`` for an explicit shaped-beam
constructor. ``getProperty``, ``withProperty``, and ``withProperties`` manage
additional metadata; ``toDict`` and ``modeDict`` expose the serialized
low-level values for advanced diagnostics.
