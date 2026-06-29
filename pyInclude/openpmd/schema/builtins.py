from __future__ import annotations

import numpy as np

from .core import (
    ExtensionAttributeSpec,
    PrimitiveFieldSpec,
    PrimitiveSchemaDefinition,
)
from .unit_dimension import unitDimension


HASE_TRANSPORT_ATTRIBUTES = {
    "haseTopologyConvention": "vtkWedgeExtrudedTriangularPrism",
    "simulationType": "laserCrystalASE",
    "geometryType": "extrudedTriangularPrism",
}


SIMULATION_ATTRIBUTE_SPECS = (
    ExtensionAttributeSpec("numberOfPoints", "number_of_points", "int", unit="count", unitDimension=unitDimension.numberOfPoints),
    ExtensionAttributeSpec("numberOfTriangles", "number_of_cells", "int", unit="count", unitDimension=unitDimension.numberOfTriangles),
    ExtensionAttributeSpec("numberOfLevels", "number_of_levels", "int", unit="count", unitDimension=unitDimension.numberOfLevels),
    ExtensionAttributeSpec("thickness", "thickness", "float", unit="m", unitDimension=unitDimension.thickness),
    ExtensionAttributeSpec("nTot", "n_tot", "float", unit="cm^-3", unitSI=1.0e6, unitDimension=unitDimension.nTot),
    ExtensionAttributeSpec("crystalTFluo", "crystal_t_fluo", "float", unit="s", unitDimension=unitDimension.crystalTFluo),
    ExtensionAttributeSpec("claddingNumber", "cladding_number", "int", unit="count", unitDimension=unitDimension.claddingNumber),
    ExtensionAttributeSpec(
        "claddingAbsorption", "cladding_absorption", "float", unit="cm^-1", unitSI=100.0, unitDimension=unitDimension.claddingAbsorption
    ),
    ExtensionAttributeSpec("minRaysPerSample", "min_rays_per_sample", "int", unit="count", unitDimension=unitDimension.minRaysPerSample),
    ExtensionAttributeSpec("maxRaysPerSample", "max_rays_per_sample", "int", unit="count", unitDimension=unitDimension.maxRaysPerSample),
    ExtensionAttributeSpec("mseThreshold", "mse_threshold", "float", unitDimension=unitDimension.mseThreshold),
    ExtensionAttributeSpec("repetitions", "repetitions", "int", unit="count", unitDimension=unitDimension.repetitions),
    ExtensionAttributeSpec("adaptiveSteps", "adaptive_steps", "int", unit="count", unitDimension=unitDimension.adaptiveSteps),
    ExtensionAttributeSpec("useReflections", "use_reflections", "bool", unitDimension=unitDimension.useReflections),
    ExtensionAttributeSpec("spectralResolution", "spectral_resolution", "int", unit="count", unitDimension=unitDimension.spectralResolution),
    ExtensionAttributeSpec("monochromatic", "monochromatic", "bool", unitDimension=unitDimension.monochromatic),
    ExtensionAttributeSpec(
        "maxSigmaAbsorption", "max_sigma_absorption", "float", unit="cm^2", unitSI=1.0e-4, unitDimension=unitDimension.crossSection
    ),
    ExtensionAttributeSpec(
        "maxSigmaEmission", "max_sigma_emission", "float", unit="cm^2", unitSI=1.0e-4, unitDimension=unitDimension.crossSection
    ),
    ExtensionAttributeSpec("backend", "backend", "str", unitDimension=unitDimension.backend),
    ExtensionAttributeSpec("maxGpus", "max_gpus", "int", unit="count", unitDimension=unitDimension.maxGpus),
    ExtensionAttributeSpec("parallelMode", "parallel_mode", "str", unitDimension=unitDimension.parallelMode),
    ExtensionAttributeSpec("minSampleRange", "min_sample_range", "int", unit="index", unitDimension=unitDimension.minSampleRange),
    ExtensionAttributeSpec("maxSampleRange", "max_sample_range", "int", unit="index", unitDimension=unitDimension.maxSampleRange),
    ExtensionAttributeSpec("rngSeed", "rng_seed", "int", unitDimension=unitDimension.rngSeed),
)


class PointSchema(PrimitiveSchemaDefinition):
    primitiveName = "point"
    axes = ("point",)
    shapeField = "position"

    position = PrimitiveFieldSpec(
        "position",
        "vertices",
        np.float64,
        axes=("coordinate", "point"),
        shape=lambda context: (2, context.numberOfPoints),
        unit="m",
        unitDimension=unitDimension.position,
    )
    pointBeta = PrimitiveFieldSpec(
        "pointBeta", np.float64, axes=("point", "level"), unitDimension=unitDimension.pointBeta, dynamic=True
    )
    phiAse = PrimitiveFieldSpec(
        "phiAse",
        np.float32,
        axes=("point", "level"),
        unit="cm^-2 s^-1",
        unitSI=1.0e4,
        unitDimension=unitDimension.phiAse,
        dynamic=True,
        backendRequired=False,
        schemaRole="result",
    )
    mse = PrimitiveFieldSpec(
        "mse", np.float64, axes=("point", "level"), unitDimension=unitDimension.mse, dynamic=True, backendRequired=False, schemaRole="result"
    )
    totalRays = PrimitiveFieldSpec(
        "totalRays",
        np.uint32,
        axes=("point", "level"),
        unit="count",
        unitDimension=unitDimension.totalRays,
        dynamic=True,
        backendRequired=False,
        schemaRole="result",
    )
    dndtAse = PrimitiveFieldSpec(
        "dndtAse",
        np.float64,
        axes=("point", "level"),
        unit="s^-1",
        unitDimension=unitDimension.dndtAse,
        dynamic=True,
        backendRequired=False,
        schemaRole="result",
    )


class TriangleSchema(PrimitiveSchemaDefinition):
    primitiveName = "triangle"
    axes = ("cell",)
    shapeField = "connectivity"
    connectivity = PrimitiveFieldSpec("connectivity", np.uint32, axes=("cell", "local_vertex"), unitDimension=unitDimension.connectivity)
    neighbors = PrimitiveFieldSpec("neighbors", np.int32, axes=("cell", "local_side"), unitDimension=unitDimension.neighbors)
    forbiddenEdges = PrimitiveFieldSpec("forbiddenEdges", np.int32, axes=("cell", "local_side"), unitDimension=unitDimension.forbiddenEdges)
    normalPoints = PrimitiveFieldSpec("normalPoints", np.uint32, axes=("cell", "local_side"), unitDimension=unitDimension.normalPoints)
    center = PrimitiveFieldSpec("center", "cell_center", np.float64, axes=("coordinate", "cell"), unit="m", unitDimension=unitDimension.center)
    normal = PrimitiveFieldSpec("normal", np.float64, axes=("cell", "local_side", "coordinate"), unitDimension=unitDimension.normal)
    surface = PrimitiveFieldSpec("surface", np.float32, axes=("cell",), unit="m^2", unitDimension=unitDimension.surface)
    claddingGroup = PrimitiveFieldSpec("claddingGroup", np.uint32, axes=("cell",), unitDimension=unitDimension.claddingGroup)
    refractiveIndex = PrimitiveFieldSpec("refractiveIndex", np.float32, axes=("interface",), unitDimension=unitDimension.refractiveIndex)
    reflectivity = PrimitiveFieldSpec("reflectivity", np.float32, axes=("cell", "interface"), unitDimension=unitDimension.reflectivity)


class PrismSchema(PrimitiveSchemaDefinition):
    primitiveName = "prism"
    axes = ("cell", "layer")
    shapeField = "betaVolume"

    betaVolume = PrimitiveFieldSpec("betaVolume", np.float64, axes=("cell", "layer"), unitDimension=unitDimension.betaVolume, dynamic=True)


BACKEND_FIELD_SPECS = (
    PrimitiveFieldSpec("cellCenterX", "cell_center_x", np.float64, axes=("cell",), unit="m", unitDimension=unitDimension.cellCenterX),
    PrimitiveFieldSpec("cellCenterY", "cell_center_y", np.float64, axes=("cell",), unit="m", unitDimension=unitDimension.cellCenterY),
    PrimitiveFieldSpec("cellNormalX", "cell_normal_x", np.float64, axes=("cell", "local_side"), unitDimension=unitDimension.cellNormalX),
    PrimitiveFieldSpec("cellNormalY", "cell_normal_y", np.float64, axes=("cell", "local_side"), unitDimension=unitDimension.cellNormalY),
    PrimitiveFieldSpec("claddingCellType", "cladding_cell_type", np.uint32, axes=("cell",), unitDimension=unitDimension.claddingCellType),
)


COMPONENT_FIELD_SPECS = (
    PrimitiveFieldSpec(
        "lambdaAbsorption", "lambda_absorption", np.float64, axes=("wavelength",), unit="m", unitDimension=unitDimension.lambdaAbsorption, backendRequired=False
    ),
    PrimitiveFieldSpec(
        "lambdaEmission", "lambda_emission", np.float64, axes=("wavelength",), unit="m", unitDimension=unitDimension.lambdaEmission, backendRequired=False
    ),
    PrimitiveFieldSpec(
        "sigmaAbsorption",
        "sigma_absorption",
        np.float64,
        axes=("wavelength",),
        unit="cm^2",
        unitSI=1.0e-4,
        unitDimension=unitDimension.sigmaAbsorption,
        backendRequired=False,
    ),
    PrimitiveFieldSpec(
        "sigmaEmission",
        "sigma_emission",
        np.float64,
        axes=("wavelength",),
        unit="cm^2",
        unitSI=1.0e-4,
        unitDimension=unitDimension.sigmaEmission,
        backendRequired=False,
    ),
)


PRIMITIVE_SCHEMA_CLASSES = {
    "point": PointSchema,
    "triangle": TriangleSchema,
    "prism": PrismSchema,
}


FIELD_ALIASES = {
    "points": "position",
}
