Python API Reference
====================

This page is generated from the Python objects exposed by ``HASEonGPU``.  It is
the reference page for signatures, class members, functions, and direct object
lookup.  For workflow-oriented explanations and complete simulation examples,
start with the :doc:`Python Interface Guide <pythonInterface>`.

Public API
----------

.. currentmodule:: HASEonGPU

Geometry and gain media
^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   Grid
   MeshTopology
   GainMedium
   GainMediumGeometry
   Gmsh

Spectra, pump, and ASE configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   CrossSectionData
   LaserProperties
   SpectralDecomposition
   PumpRadiationProfile
   PumpProperties
   PhiASE

Simulation and time integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   Simulation
   TimeStepState
   TimeSteppedSimulation
   TimeDerivative
   TimeIntegrationResult
   TimeIntegrationSolver
   ExplicitEuler
   ExponentialEuler
   FrozenPhiAseRungeKutta4
   Heun
   ImplicitEuler
   Midpoint
   RungeKutta4

Pump solvers, field helpers, and utilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   BetaInt3PumpSolver
   BetaIntegrationSolver
   BetaIntegrationGaussianSolver
   OneDimensionalZTraversal
   oneDimensionalZTraversalPumpRate
   Constants
   LegacyGridDataBetaVolumeMapper
   calcGainFromState
   AlpakaBackends
   BaseSchema
   PointSchema
   TriangleSchema
   PrismSchema
   PrimitiveFieldSpec
   backendFlat
   vtkWedge
   writeGainMediumVtk
   writeParaviewState

The ``unitDimension`` namespace is exported from ``HASEonGPU`` and contains
predefined openPMD unit-dimension tuples for HASE transport variables and common
dimensions.

Compatibility helpers
^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   beta_int3Main
   integrateLaserPump
   runLaserPumpStep
