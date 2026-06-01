Utilities
=========

This page lists supporting public objects that are commonly used with the new
Python interface.

Time Integration Solvers
------------------------

Import built-in solvers from ``HASEonGPU``:

.. code-block:: python

   from HASEonGPU import (
       ExplicitEuler,
       Heun,
       Midpoint,
       RungeKutta4,
       ImplicitEuler,
       ExponentialEuler,
   )

Available solvers:

* ``ExplicitEuler()``
* ``Heun()``
* ``Midpoint()``
* ``RungeKutta4()``
* ``ImplicitEuler(iterations=8, tolerance=1e-10)``
* ``ExponentialEuler()``

All solvers implement:

.. code-block:: python

   step(rhs, betaCells, time, timeStep)

The return value is a ``TimeIntegrationResult`` with:

* ``betaCells``: updated beta array :math:`\beta_i`.
* ``evaluation``: the ``TimeDerivative`` used by the solver.

Custom Time Integration
-----------------------

Custom solvers only need the same ``step`` method:

.. code-block:: python

   class MyEuler:
       def step(self, rhs, betaCells, time, timeStep):
           evaluation = rhs(betaCells, time)
           return TimeIntegrationResult(
               betaCells=betaCells + timeStep * evaluation.derivative,
               evaluation=evaluation,
           )

``TimeDerivative`` contains:

* ``betaCells``
* ``dndtPump``
* ``dndtAse``
* ``derivative``
* ``tau``
* ``phiAse``
* ``aseResult``

Here ``dndtPump``, ``dndtAse``, and ``derivative`` are contributions to
:math:`d\beta/dt`; ``tau`` is the fluorescence lifetime :math:`\tau`, and
``phiAse`` is the ASE flux :math:`\Phi_i`.

VTK Export
----------

``vtkWedge`` writes point or cell data on the wedge mesh to a legacy ASCII
VTK file.

Callback use:

.. code-block:: python

   simulation.onStep(vtkWedge("minimal_phi_ase_{step:03d}.vtk", medium))

Direct use:

.. code-block:: python

   state = simulation.step()
   vtkWedge("phi.vtk", state, medium)
   vtkWedge("beta.vtk", state.betaCells, medium, scalar_name="betaCells")

By default it reads ``state.phiAse`` (:math:`\Phi_i`).  Use ``field`` to export
another field:

.. code-block:: python

   simulation.onStep(vtkWedge("beta_{step:03d}.vtk", medium, field="betaCells"))

Use ``every`` to reduce output frequency:

.. code-block:: python

   simulation.onStep(vtkWedge("phi_{step:03d}.vtk", medium, every=10))

The data shape must match either:

* point data: ``(numberOfPoints, numberOfLevels)``
* cell data: ``(numberOfTriangles, numberOfLevels - 1)``

Physical Constants
------------------

``Constants`` stores the predefined physical constants used by the default pump
routine.  ``Simulation`` creates this object automatically, so most users do
not need to pass constants explicitly.

.. code-block:: python

   from HASEonGPU import Constants

   constants = Constants()
   constants.speedOfLight
   constants.planckConstant

These correspond to :math:`c` and :math:`h` in the pump equations.

The short names used by the low-level pump routine remain available as aliases:

.. code-block:: python

   constants.c
   constants.h
   constants.toDict()
   constants.describeConstant("speedOfLight")

Backend Names
-------------

``AlpakaBackends`` can list backend names discovered from the installed
HASEonGPU backend-name library.  Use these strings wherever the Python
interface accepts a ``backend`` option, for example in ``PhiASE`` or the
low-level ``calcPhiASE(...)`` call:

.. code-block:: python

   from HASEonGPU import AlpakaBackends, PhiASE

   available = AlpakaBackends.all()
   backend = available[0]

   phi_ase = PhiASE(backend=backend)

``AlpakaBackends.known()`` is an alias for ``AlpakaBackends.all()``.  Backend
names that are valid Python identifiers are also exposed as class attributes,
for example ``AlpakaBackends.Host_Cpu_CpuSerial``.

For details on how the helper library is built and how backend names are
formed, see :doc:`../backendSelection`.
