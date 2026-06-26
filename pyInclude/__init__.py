# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Public Python convenience exports for HASEonGPU."""

import HASEonGPU_Bindings

__version__ = "2.0"

from .alpakaUtils import AlpakaBackends
from .calcPhiASE import calcPhiASE
from .openpmd import (
    BaseGroup,
    BaseSchema,
    GroupFieldSpec,
    PointSchema,
    PrimitiveFieldSpec,
    PrismSchema,
    TriangleSchema,
    backendFlat,
    unitDimension,
)
from .geometry import GainMedium, GainMediumGeometry, Gmsh, Grid, MeshTopology, writeGainMediumVtk
from .laser import CrossSectionData, LaserProperties, PumpProperties, SpectralDecomposition
from .pumping import (
    BetaInt3PumpSolver,
    BetaIntegrationSolver,
    BetaIntegrationGaussianSolver,
    Constants,
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
from .structures import Result as TransportResult
from .vtkWedge import vtkWedge
from .timeIntegration import (
    ExplicitEuler,
    ExponentialEuler,
    Heun,
    ImplicitEuler,
    Midpoint,
    RungeKutta4,
    TimeDerivative,
    TimeIntegrationResult,
    TimeIntegrationSolver,
)

Mesh = HASEonGPU_Bindings.HostMesh
