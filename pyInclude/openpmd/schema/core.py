from __future__ import annotations

from dataclasses import dataclass

from ..FieldSpec import FieldSpec
from .unit_dimension import DIMENSIONLESS


BACKEND_FLAT = "backendFlat"


def _camel_to_snake(name: str) -> str:
    out = []
    for index, char in enumerate(name):
        if char.isupper() and index > 0 and (not name[index - 1].isupper()):
            out.append("_")
        out.append(char.lower())
    return "".join(out)


def _shapePoint(context):
    return (context.numberOfPoints,)


def _shapeCell(context):
    return (context.numberOfTriangles,)


def _shapeCellSide(context):
    return (context.numberOfTriangles, 3)


def _shapeCellLayer(context):
    return (context.numberOfTriangles, context.numberOfLevels - 1)


def _shapePointLevel(context):
    return (context.numberOfPoints, context.numberOfLevels)


def _shapeInterface(context):
    return (4,)


def _shapeSurface(context):
    return (context.numberOfTriangles, 2)


def _shapeWavelength(context):
    return (context.spectral,)



_SHAPES_BY_AXES = {
    ("point",): _shapePoint,
    ("cell",): _shapeCell,
    ("coordinate", "point"): lambda context: (2, context.numberOfPoints),
    ("point", "coordinate"): lambda context: (context.numberOfPoints, 2),
    ("coordinate", "cell"): lambda context: (2, context.numberOfTriangles),
    ("cell", "coordinate"): lambda context: (context.numberOfTriangles, 2),
    ("cell", "local_vertex"): _shapeCellSide,
    ("cell", "local_side"): _shapeCellSide,
    ("cell", "local_side", "coordinate"): lambda context: (context.numberOfTriangles, 3, 2),
    ("cell", "layer"): _shapeCellLayer,
    ("point", "level"): _shapePointLevel,
    ("interface",): _shapeInterface,
    ("cell", "interface"): _shapeSurface,
    ("wavelength",): _shapeWavelength,
}


@dataclass(frozen=True, init=False)
class PrimitiveFieldSpec:
    name: str
    recordName: str
    dtype: object
    axes: tuple[str, ...] | None = None
    shape: object | None = None
    unit: str = "1"
    unitSI: float = 1.0
    unitDimension: tuple[float, float, float, float, float, float, float] = DIMENSIONLESS
    dynamic: bool = False
    backendRequired: bool = True
    userDefined: bool = False
    schemaRole: str = "input"

    def __init__(
        self,
        name: str,
        recordNameOrDtype: object,
        dtype: object | None = None,
        axes: tuple[str, ...] | None = None,
        shape: object | None = None,
        unit: str = "1",
        unitSI: float = 1.0,
        unitDimension: tuple[float, float, float, float, float, float, float] = DIMENSIONLESS,
        dynamic: bool = False,
        backendRequired: bool = True,
        userDefined: bool = False,
        schemaRole: str = "input",
        *,
        recordName: str | None = None,
    ):
        if dtype is None:
            dtype = recordNameOrDtype
            resolved_record_name = recordName or _camel_to_snake(name)
        else:
            resolved_record_name = str(recordName or recordNameOrDtype)

        object.__setattr__(self, "name", name)
        object.__setattr__(self, "recordName", resolved_record_name)
        object.__setattr__(self, "dtype", dtype)
        object.__setattr__(self, "axes", axes)
        object.__setattr__(self, "shape", shape)
        object.__setattr__(self, "unit", unit)
        object.__setattr__(self, "unitSI", unitSI)
        object.__setattr__(self, "unitDimension", unitDimension)
        object.__setattr__(self, "dynamic", dynamic)
        object.__setattr__(self, "backendRequired", backendRequired)
        object.__setattr__(self, "userDefined", userDefined)
        object.__setattr__(self, "schemaRole", schemaRole)

    def toFieldSpec(self, primitiveAxes: tuple[str, ...]) -> FieldSpec:
        axes = primitiveAxes if self.axes is None else tuple(self.axes)
        shape = self.shape if self.shape is not None else _SHAPES_BY_AXES[axes]
        return FieldSpec(
            self.name,
            self.recordName,
            axes,
            self.dtype,
            shape,
            unit=self.unit,
            unitSI=self.unitSI,
            unitDimension=self.unitDimension,
            dynamic=self.dynamic,
            backendRequired=self.backendRequired,
            userDefined=self.userDefined,
            schemaRole=self.schemaRole,
        )


@dataclass(frozen=True)
class PrimitiveSchema:
    name: str
    axes: tuple[str, ...]
    fields: tuple[PrimitiveFieldSpec, ...]
    shapeField: str | None = None

    def fieldSpecs(self) -> tuple[FieldSpec, ...]:
        return tuple(field.toFieldSpec(self.axes) for field in self.fields)

    def extend(self, *fields: PrimitiveFieldSpec) -> "PrimitiveSchema":
        return PrimitiveSchema(self.name, self.axes, self.fields + tuple(fields), self.shapeField)

    def fieldSpec(self, name: str) -> FieldSpec:
        for spec in self.fieldSpecs():
            if spec.name == name:
                return spec
        raise KeyError(name)

    def expectedShape(self, context, parentShape: tuple[int, ...] | None = None) -> tuple[int, ...]:
        if self.shapeField is not None:
            spec = self.fieldSpec(self.shapeField)
            field_shape = spec.expectedShape(context)
            return tuple(size for size, axis in zip(field_shape, spec.axes) if axis in self.axes)
        if parentShape is not None:
            return tuple(parentShape)
        return tuple(_SHAPES_BY_AXES[self.axes](context))


@dataclass(frozen=True)
class ExtensionAttributeSpec:
    name: str
    attribute: str
    dtype: str
    unit: str = "1"
    unitSI: float = 1.0
    unitDimension: tuple[float, float, float, float, float, float, float] = DIMENSIONLESS

    def cast(self, value):
        if self.dtype == "bool":
            return bool(value)
        if self.dtype == "int":
            return int(value)
        if self.dtype == "float":
            return float(value)
        if self.dtype == "str":
            return str(value)
        raise ValueError(f"unknown openPMD transport attribute dtype '{self.dtype}'")


_GROUP_FIELD_MISSING = object()


@dataclass(frozen=True)
class GroupFieldSpec:
    name: str
    dtype: object
    default: object = _GROUP_FIELD_MISSING
    targetName: str | None = None
    unit: str = "1"
    unitSI: float = 1.0
    unitDimension: tuple[float, float, float, float, float, float, float] = DIMENSIONLESS

    @property
    def target(self) -> str:
        return self.targetName or self.name


class BaseGroup:
    groupName: str | None = None

    def __init__(self, **values):
        specs = {spec.name: spec for spec in self.fieldSpecs()}
        unknown = set(values) - set(specs)
        if unknown:
            raise TypeError(f"unknown group field '{sorted(unknown)[0]}'")
        resolved = {}
        for name, spec in specs.items():
            if name in values:
                resolved[name] = values[name]
            elif spec.default is not _GROUP_FIELD_MISSING:
                resolved[name] = spec.default
            else:
                raise TypeError(f"missing value for group field '{name}'")
        self._values = resolved
        self._groupToken = object()

    @classmethod
    def declaredFields(cls) -> tuple[GroupFieldSpec, ...]:
        fields = []
        seen = set()
        for group_cls in reversed(cls.mro()):
            for value in group_cls.__dict__.values():
                if isinstance(value, GroupFieldSpec) and value.name not in seen:
                    fields.append(value)
                    seen.add(value.name)
        return tuple(fields)

    @classmethod
    def fieldSpecs(cls) -> tuple[GroupFieldSpec, ...]:
        return cls.declaredFields()

    @property
    def name(self) -> str:
        return self.groupName or self.__class__.__name__

    def fieldItems(self):
        for spec in self.fieldSpecs():
            yield spec, self._values[spec.name]

    def add(self, primitive):
        if hasattr(primitive, "_applyGroup"):
            primitive._applyGroup(self)
            return primitive
        for spec, value in self.fieldItems():
            current = getattr(primitive, spec.target, _GROUP_FIELD_MISSING)
            if current is not _GROUP_FIELD_MISSING:
                raise ValueError(f"field '{spec.target}' is already assigned on this primitive")
            setattr(primitive, spec.target, value)
        return primitive


class BaseSchema:
    primitiveName: str | None = None
    axes: tuple[str, ...] = ()
    shapeField: str | None = None

    @classmethod
    def declaredFields(cls) -> tuple[PrimitiveFieldSpec, ...]:
        fields = []
        seen = set()
        for schema_cls in reversed(cls.mro()):
            for value in schema_cls.__dict__.values():
                if isinstance(value, PrimitiveFieldSpec) and value.name not in seen:
                    fields.append(value)
                    seen.add(value.name)
        return tuple(fields)

    @classmethod
    def primitiveSchema(cls) -> PrimitiveSchema:
        name = cls.primitiveName or cls.__name__
        return PrimitiveSchema(name, tuple(cls.axes), cls.declaredFields(), cls.shapeField)

    @classmethod
    def fieldSpecs(cls) -> tuple[FieldSpec, ...]:
        return cls.primitiveSchema().fieldSpecs()


class PrimitiveSchemaDefinition(BaseSchema):
    pass
