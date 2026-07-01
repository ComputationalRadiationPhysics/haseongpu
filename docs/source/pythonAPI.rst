Python API Reference
====================

This page is generated from the public Python objects exposed by ``HASEonGPU``.
It is a reference for signatures and members. For the user workflow and a
complete example, start with the :doc:`Python Interface Guide <pythonInterface>`.

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

Spectra, pump, and ASE
^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   CrossSectionData
   LaserProperties
   SpectralDecomposition
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
   ExplicitEuler
   ExponentialEuler
   Heun
   ImplicitEuler
   Midpoint
   RungeKutta4
   TimeIntegrationSolver

Utilities and transport schemas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   AlpakaBackends
   BaseSchema
   PointSchema
   TriangleSchema
   PrismSchema
   PrimitiveFieldSpec
   backendFlat
   calcGainFromState
   vtkWedge
   writeGainMediumVtk
   writeParaviewState

The ``unitDimension`` namespace is exported from ``HASEonGPU`` and contains
predefined openPMD unit-dimension tuples for HASE transport variables and common
dimensions.
