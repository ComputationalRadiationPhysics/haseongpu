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
VTK file. In a ``Simulation.onStep`` callback, pass the ``TimeStepState`` to
``vtkWedge``; the state carries the static topology and the dynamic arrays.

Callback use:

.. code-block:: python

   def write_vtk(state, output_dir, cladding_absorption):
       vtkWedge(
           output_dir / "fields_{step:03d}.vtk",
           state,
           fields={
               "betaCells": state.betaCells,
               "phiASE": state.phiAse,
               "dndtAse": state.dndtAse,
               "cladAbs": state.phiAse * cladding_absorption,
           },
       )

   simulation.onStep(write_vtk, output_dir, 5.5)

Direct use after one step:

.. code-block:: python

   state = simulation.step()
   vtkWedge("phi.vtk", state)
   vtkWedge("fields.vtk", state, field=["phiAse", "dndtAse"])
   vtkWedge("named.vtk", state, field={"phi": "phiAse", "dn": "dndtAse"})

For standalone array exports outside a simulation state, pass ``geometry`` as a
``GainMedium`` or ``MeshTopology``:

.. code-block:: python

   vtkWedge("fields.vtk", geometry=medium, fields={"phi": phi, "dn": dndt})

The older callback-factory form is still accepted and can use ``every`` to
reduce output frequency:

.. code-block:: python

   simulation.onStep(vtkWedge("phi_{step:03d}.vtk", medium, every=10))

For new code, prefer an explicit callback when output frequency or derived
fields are needed:

.. code-block:: python

   def write_every_tenth(state, output_dir):
       if state.step % 10 == 0:
           vtkWedge(output_dir / "phi_{step:03d}.vtk", state)

   simulation.onStep(write_every_tenth, output_dir)

The data shape must match either:

* point data: ``(numberOfPoints, numberOfLevels)``
* cell data: ``(numberOfTriangles, numberOfLevels - 1)``

Gain Field Export
-----------------

``calcGainFromState`` calculates small-signal laser gain from a
``TimeStepState`` and returns a point-shaped array that can be written directly
with ``vtkWedge``:

.. code-block:: python

   vtkWedge(
       output_path,
       state,
       fields={
           "gain": calcGainFromState(state, spectra, nTot),
       },
   )


Physical Constants
------------------

``Constants`` stores predefined physical constants.
``Simulation`` creates this object automatically, so most users do not need to pass constants explicitly.

.. code-block:: python

   from HASEonGPU import Constants

   constants = Constants()
   constants.speedOfLight
   constants.planckConstant

These correspond to :math:`c` and :math:`h` in the pump equations.

The short names used by the pump implementation remain available as aliases:

.. code-block:: python

   constants.c
   constants.h
   constants.toDict()
   constants.describeConstant("speedOfLight")

Backend Names
-------------

``AlpakaBackends`` can list backend names discovered from the installed
HASEonGPU backend-name library.  Use these strings wherever the Python
interface accepts a ``backend`` option, for example in ``PhiASE``:

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
