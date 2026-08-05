Python API Reference
====================

This page is generated from the public Python objects exposed by ``HASEonGPU``.
It is a reference for signatures and members, not a workflow description. Start
with the :doc:`laserPumpCladding tutorial <laserPumpCladding>` for a complete
example and use the :doc:`Python Interface Guide <pythonInterface>` to navigate
the reusable concepts.

Public API
----------

.. currentmodule:: HASEonGPU

Geometry and gain media
^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   VolumeTopology
   GainMedium
   Gmsh
   DomainMap
   SurfaceDomainMap
   SurfaceOptics
   writeGainMediumVtk

Legacy planar geometry
^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   Grid
   MeshTopology
   GainMediumGeometry
   vtkWedge

Spectra, pump, and ASE
^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   CrossSectionData
   LaserProperties
   SpectralDecomposition
   PumpSpectrum
   PumpAngularDistribution
   UniformPumpProfile
   SuperGaussianPumpProfile
   Pump
   GaussianPump
   SurfacePumpInjector
   PlanarPumpRelay
   MonteCarloPumpSolver
   integrate_pump_profile
   PhiASE

Simulation and time integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   Simulation
   TimeStepState
   TimeSteppedSimulation
   TransportResult
   ExplicitEuler
   ExponentialEuler
   FrozenPhiAseRungeKutta4
   Heun
   ImplicitEuler
   Midpoint
   RungeKutta4
   TimeIntegrationSolver

Utilities
^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   AlpakaBackends
   OpenPmdBackends
   backendFlat
   calcGainFromState
   writeParaviewState

Low-level transport schemas
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   BaseGroup
   BaseSchema
   GroupFieldSpec
   PointSchema
   TriangleSchema
   PrismSchema
   PrimitiveFieldSpec

The ``unitDimension`` namespace is exported from ``HASEonGPU`` and contains
predefined openPMD unit-dimension tuples for HASE transport variables and common
dimensions.
