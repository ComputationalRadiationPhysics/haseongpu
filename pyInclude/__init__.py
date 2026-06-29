# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Public Python convenience exports for HASEonGPU."""

import HASEonGPU_Bindings

__version__ = "2.0.1"

from .alpakaUtils import AlpakaBackends
from .calcPhiASE import calcPhiASE
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
Mesh=HASEonGPU_Bindings.HostMesh
