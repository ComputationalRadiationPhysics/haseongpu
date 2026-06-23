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
   Heun
   ImplicitEuler
   Midpoint
   RungeKutta4

Pump solvers and utilities
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   BetaInt3PumpSolver
   BetaIntegrationSolver
   BetaIntegrationGaussianSolver
   Constants
   LegacyGridDataBetaVolumeMapper
   calcGainFromState
   AlpakaBackends
   vtkWedge
   writeGainMediumVtk

Lower-level compatibility helpers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:

   calcPhiASE
   beta_int3Main
   integrateLaserPump
   runLaserPumpStep

Module Summaries
----------------

These modules back the public API.  They are documented for users who need to
inspect implementation-level helpers or extension points.

.. currentmodule:: pyInclude

.. autosummary::
   :toctree: generated
   :recursive:

   geometry.core
   geometry.msh
   geometry.vtk
   laser
   pumping
   simulation
   gainMap
   timeIntegration
   vtkWedge
   calcPhiASE
