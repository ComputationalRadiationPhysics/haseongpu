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
       FrozenPhiAseRungeKutta4,
       ImplicitEuler,
       ExponentialEuler,
       FrozenPhiAseRungeKutta4,
   )

Available solvers:

* ``ExplicitEuler()``
* ``Heun()``
* ``Midpoint()``
* ``RungeKutta4()``
* ``FrozenPhiAseRungeKutta4()``: reuses one ASE evaluation across RK4 stages.
* ``FrozenPhiAseRungeKutta4()``
* ``ImplicitEuler(iterations=8, tolerance=1e-10)``
* ``ExponentialEuler()``

These objects are lightweight descriptors. Python serializes their ``name`` to
the openPMD run-control record and the compiled C++/Alpaka backend performs the
actual time integration. You can also pass one of the names directly as a
string.

Custom Python time integrators are not supported by compiled simulation runs.

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
