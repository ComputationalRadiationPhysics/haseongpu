# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Public Python convenience exports for HASEonGPU's openPMD frontend."""

__version__ = "2.1.1"

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
from .geometry import (
    DomainMap,
    GainMedium,
    GainMediumGeometry,
    Gmsh,
    Grid,
    MeshTopology,
    SurfaceDomainMap,
    SurfaceOptics,
    VolumeTopology,
    writeGainMediumVtk,
)
from .laser import (
    CrossSectionData,
    GaussianPump,
    LaserProperties,
    MonteCarloPumpSolver,
    PlanarPumpRelay,
    Pump,
    PumpAngularDistribution,
    PumpSpectrum,
    SpectralDecomposition,
    SuperGaussianPumpProfile,
    SurfacePumpInjector,
    UniformPumpProfile,
    integrate_pump_profile,
)
from .simulation import (
    ConnectivityAverageBetaVolumeMapper,
    LegacyGridDataBetaVolumeMapper,
    PhiASE,
    Simulation,
    TimeStepState,
    TimeSteppedSimulation,
)
from .structures import Result as TransportResult
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
    TimeIntegrationSolver,
)
