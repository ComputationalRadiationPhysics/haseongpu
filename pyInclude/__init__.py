# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Public Python convenience exports for HASEonGPU's openPMD frontend."""

__version__ = "2.1.0"

from ._runtime import activate_openpmd_python_provider as _activate_openpmd_python_provider

_activate_openpmd_python_provider()
del _activate_openpmd_python_provider

from .alpakaUtils import AlpakaBackends
from .openpmd import (
    BaseGroup,
    BaseSchema,
    GroupFieldSpec,
    OpenPmdBackends,
    PointSchema,
    PrimitiveFieldSpec,
    PrismSchema,
    TriangleSchema,
    backendFlat,
    unitDimension,
)
from .geometry import GainMedium, GainMediumGeometry, Gmsh, Grid, MeshTopology, writeGainMediumVtk
from .laser import CrossSectionData, LaserProperties, PumpProperties, SpectralDecomposition, PumpRadiationProfile
from .pumping import (
    BetaInt3PumpSolver,
    BetaIntegrationSolver,
    BetaIntegrationGaussianSolver,
    Constants,
    OneDimensionalZTraversal,
    oneDimensionalZTraversalPumpRate,
    beta_int3Main,
    integrateLaserPump,
    runLaserPumpStep,
)
from .simulation import (
    LegacyGridDataBetaVolumeMapper,
    PhiASE,
    Simulation,
    TimeStepState,
    TimeSteppedSimulation,
)
from .gainMap import calcGainFromState
from .openpmd.paraview import writeParaviewState
from .vtkWedge import vtkWedge
from .timeIntegration import (
    ExplicitEuler,
    ExponentialEuler,
    FrozenPhiAseRungeKutta4,
    Heun,
    ImplicitEuler,
    Midpoint,
    RungeKutta4,
    TimeDerivative,
    TimeIntegrationResult,
    TimeIntegrationSolver,
)
