from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

from .FieldSpec import FieldSpec
from .schema import (
    BACKEND_FIELD_SPECS,
    COMPONENT_FIELD_SPECS,
    FIELD_ALIASES,
    HASE_TRANSPORT_ATTRIBUTES,
    PRIMITIVE_SCHEMA_CLASSES,
    SIMULATION_ATTRIBUTE_SPECS,
    BaseGroup,
    BaseSchema,
    ExtensionAttributeSpec,
    GroupFieldSpec,
    PointSchema,
    PrimitiveFieldSpec,
    PrimitiveSchema,
    PrimitiveSchemaDefinition,
    PrismSchema,
    TriangleSchema,
    unitDimension,
)


BACKEND_FLAT = "backendFlat"


@dataclass(frozen=True)
class BackendFlatArray:
    values: object


def backendFlat(values):
    return BackendFlatArray(values)


haseTransportAttributes = HASE_TRANSPORT_ATTRIBUTES
simulationAttributeSpecs = SIMULATION_ATTRIBUTE_SPECS
primitiveSchemaClasses = PRIMITIVE_SCHEMA_CLASSES
primitiveSchemas = {name: schema_class.primitiveSchema() for name, schema_class in primitiveSchemaClasses.items()}
componentFieldSpecs = tuple(spec.toFieldSpec(tuple(spec.axes)) for spec in COMPONENT_FIELD_SPECS)
backendFieldSpecs = tuple(spec.toFieldSpec(tuple(spec.axes)) for spec in BACKEND_FIELD_SPECS)
HASE_TRANSPORT_VERSION = "0.1"
globals().update({schema_class.__name__: schema_class for schema_class in primitiveSchemaClasses.values()})
PointSchema = primitiveSchemaClasses["point"]
TriangleSchema = primitiveSchemaClasses["triangle"]
PrismSchema = primitiveSchemaClasses["prism"]


def _fieldSpecsFromPrimitiveSchemas(schemas=primitiveSchemas):
    fields = {
        spec.name: spec
        for schema in schemas.values()
        for spec in schema.fieldSpecs()
    }
    fields.update({spec.name: spec for spec in componentFieldSpecs})
    fields.update({spec.name: spec for spec in backendFieldSpecs})
    for alias, canonical in FIELD_ALIASES.items():
        fields[alias] = fields[canonical]
    return fields


schemaFields = _fieldSpecsFromPrimitiveSchemas()


def primitiveSchema(name: str) -> PrimitiveSchema:
    return primitiveSchemas[name]


def primitiveFieldSpecs(name: str) -> tuple[FieldSpec, ...]:
    return primitiveSchema(name).fieldSpecs()


def simulationAttributeSpec(name: str) -> ExtensionAttributeSpec:
    for spec in simulationAttributeSpecs:
        if spec.name == name:
            return spec
    raise KeyError(name)


def fieldSpec(name: str) -> FieldSpec:
    return schemaFields[name]


def resultFieldSpecs():
    return tuple(spec for spec in schemaFields.values() if spec.schemaRole == "result")


def spectralContext(values):
    return SimpleNamespace(spectral=np.asarray(values).size)


def flatEntityLabel(spec: FieldSpec) -> str:
    return "flatIndex"


def _isBackendFlat(values, layoutOrder):
    return isinstance(values, BackendFlatArray) or layoutOrder == BACKEND_FLAT


def _unwrap(values):
    if isinstance(values, BackendFlatArray):
        return values.values
    return values


def backendFlatArray(values, spec: FieldSpec, context, *, layoutOrder=None):
    expectedShape = spec.expectedShape(context)
    expectedSize = int(np.prod(expectedShape, dtype=np.int64))
    arr = np.asarray(_unwrap(values), dtype=spec.dtypeObject)
    if arr.size != expectedSize:
        raise ValueError(
            f"{spec.name} expects {expectedSize} values for entity {spec.entity} "
            f"with primitive shape {expectedShape}, got shape {arr.shape}"
        )

    if _isBackendFlat(values, layoutOrder):
        if arr.ndim != 1:
            raise ValueError(
                f"{spec.name} marked backend-flat must be a 1-D array with "
                f"{expectedSize} values for entity {spec.entity}, got shape {arr.shape}"
            )
        return arr

    if arr.shape == expectedShape:
        return arr.reshape(-1, order="F")

    if arr.ndim == 1:
        raise ValueError(
            f"{spec.name} got ambiguous flat array with {arr.size} values for entity "
            f"{spec.entity}; pass backendFlat(values) or layoutOrder='backendFlat' "
            f"to declare canonical backend-flat order, or pass primitive shape "
            f"{expectedShape}"
        )

    raise ValueError(
        f"{spec.name} expects primitive shape {expectedShape} for entity "
        f"{spec.entity}, got shape {arr.shape}"
    )


def primitiveArray(values, spec: FieldSpec, context, *, layoutOrder=None):
    arr = backendFlatArray(values, spec, context, layoutOrder=layoutOrder)
    return arr.reshape(spec.expectedShape(context), order="F")


def primitiveView(values, spec: FieldSpec, context, *, layoutOrder=None):
    return primitiveArray(values, spec, context, layoutOrder=layoutOrder)
