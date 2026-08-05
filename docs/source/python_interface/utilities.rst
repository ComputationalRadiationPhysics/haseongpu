Utilities
=========

Time integration
----------------

The public solver descriptors are:

* ``ExplicitEuler()``
* ``Heun()``
* ``Midpoint()``
* ``RungeKutta4()``
* ``FrozenPhiAseRungeKutta4()``
* ``ImplicitEuler(iterations=8, tolerance=1e-10)``
* ``ExponentialEuler()``

Python serializes the descriptor name and controls; the C++/Alpaka backend
performs the integration on cell-centered fields. Standard RK4 evaluates ASE
at each stage. ``FrozenPhiAseRungeKutta4`` reuses the first ASE result for the
remaining stages while pump transport is still evaluated. Custom Python time
integrators cannot run inside the compiled loop. See
:ref:`pump-and-time-stepping` for how these evaluations enter the population
equation.

Tet4 VTK
--------

Write a static explicit-volume medium, including supported physical fields,
with either form:

.. code-block:: python

   medium.toVtk("prepared-state.vtk")
   writeGainMediumVtk("prepared-state.vtk", medium)

``GainMedium.fromVtk`` restores that prepared state. For dynamic
``TimeStepState`` output, use a cell-data writer in an ``on_step`` callback; the
:doc:`laserPumpCladding tutorial <../laserPumpCladding>` contains the current
Tet4 VTK callback used by the complete example.

``vtkWedge`` is the compatibility writer for ``MeshTopology``'s extruded
triangle layout. It accepts point arrays shaped
``(numberOfPoints, numberOfLevels)`` and cell arrays shaped
``(numberOfTriangles, numberOfLevels - 1)``. It is not the writer for an
explicit ``VolumeTopology``.

Local gain
----------

``calcGainFromState`` returns one small-signal gain value per volume cell:

.. code-block:: python

   gain = calcGainFromState(state, spectra, nTot)
   assert gain.shape == (state.topology.numberOfCells,)

By default it evaluates at the emission-spectrum maximum. Pass ``wavelength``
or explicit ``sigmaAbsorption`` and ``sigmaEmission`` for another wavelength.
The function requires ``nTot`` because a ``TimeStepState`` deliberately
contains evolving state and topology, not all medium constants.

openPMD/ParaView output
-----------------------

``writeParaviewState`` appends callback snapshots to an openPMD series and
writes a small ``.pmd`` handle for ParaView. Pass an output directory and,
optionally, the cladding absorption used for the derived field:

.. code-block:: python

   simulation.on_step(writeParaviewState, "openpmd-output", 5.5)

The output contains beta, PhiASE, pump/ASE derivatives, and topology records.
Its storage engine comes from the installed ``openpmd_api`` provider. See
:doc:`../openpmdTransport` for provider compatibility, backend selection, and
the record layout.

Backend names
-------------

``AlpakaBackends.all()`` returns compute backend names from the installed
runtime:

.. code-block:: python

   available = AlpakaBackends.all()
   phi_ase = PhiASE(backend=available[0])

Valid identifier-like names are also class attributes, such as
``AlpakaBackends.Host_Cpu_CpuSerial``. These names choose compute execution;
they do not choose the openPMD storage backend. Backend discovery, build-time
availability, and storage-backend selection are owned by
:doc:`../backendSelection`.
