# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Laser spectra and pump-beam configuration used by the Python interface."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .geometry import OpenPmdScalarField
from .openpmd import fieldSpec


@dataclass(frozen=True)
class LaserPropertySpec:
    """Metadata for one low-level laser property.

    The names match the historical ``calcPhiASE`` input fields. User-facing
    code usually works through ``CrossSectionData`` or ``LaserProperties``.
    """

    name: str
    dtype: object
    shape: tuple
    description: str
    required: bool = True


class LaserProperty:
    """Handle returned by ``LaserProperties.get(...)``.

    It exposes the property's physical description, expected dtype/shape, and
    a validated ``value`` setter. The handle writes back to its parent
    ``LaserProperties`` object.
    """

    def __init__(self, laser, spec):
        self._laser = laser
        self.spec = spec

    @property
    def name(self):
        return self.spec.name

    @property
    def description(self):
        return self.spec.description

    @property
    def dtype(self):
        return np.dtype(self.spec.dtype)

    @property
    def expectedShape(self):
        return self.spec.shape

    @property
    def value(self):
        return self._laser.values.get(self.name)

    @value.setter
    def value(self, values):
        self._laser.set(self.name, values)

    def meta(self):
        """Return serializable metadata for documentation or validation UIs."""
        return {
            "name": self.name,
            "description": self.description,
            "dtype": str(self.dtype),
            "expectedShape": self.expectedShape,
            "required": self.spec.required,
            "isSet": self.name in self._laser.values,
        }


LASER_PROPERTY_SPECS = {
    "l_abs": LaserPropertySpec(
        name="l_abs",
        dtype=np.float64,
        shape=("nAbsorptionSamples",),
        description="Wavelength values for the absorption spectrum in nm.",
    ),
    "s_abs": LaserPropertySpec(
        name="s_abs",
        dtype=np.float64,
        shape=("nAbsorptionSamples",),
        description="Absorption cross-section values in cm^2, corresponding to l_abs.",
    ),
    "l_ems": LaserPropertySpec(
        name="l_ems",
        dtype=np.float64,
        shape=("nEmissionSamples",),
        description="Wavelength values for the emission spectrum in nm.",
    ),
    "s_ems": LaserPropertySpec(
        name="s_ems",
        dtype=np.float64,
        shape=("nEmissionSamples",),
        description="Emission cross-section values in cm^2, corresponding to l_ems.",
    ),
    "l_res": LaserPropertySpec(
        name="l_res",
        dtype=np.uint32,
        shape=(),
        description="Spectral interpolation resolution used by calcPhiASE.",
    )
}

LASER_ALIASES = {
    "lambdaA": "l_abs",
    "lambdaE": "l_ems",
    "sigmaA": "s_abs",
    "sigmaE": "s_ems",
    "crossSectionAbsorption": "s_abs",
    "crossSectionEmission": "s_ems",
    "spectral": "l_res",
    "spectralResolution": "l_res",
}


_CROSS_SECTION_FIELD_ATTRS = {
    "lambdaAbsorption": "wavelengthsAbsorption",
    "lambdaEmission": "wavelengthsEmission",
    "sigmaAbsorption": "crossSectionAbsorption",
    "sigmaEmission": "crossSectionEmission",
}
_CROSS_SECTION_FIELD_ALIASES = {
    "wavelengthsAbsorption": "lambdaAbsorption",
    "wavelengthsEmission": "lambdaEmission",
    "crossSectionAbsorption": "sigmaAbsorption",
    "crossSectionEmission": "sigmaEmission",
    "lambdaA": "lambdaAbsorption",
    "lambdaE": "lambdaEmission",
    "sigmaA": "sigmaAbsorption",
    "sigmaE": "sigmaEmission",
}
_MISSING_SPECTRAL_FIELD_VALUE = object()


class SpectralField:
    def __init__(self, crossSections, name):
        self._crossSections = crossSections
        self.name = name
        self.spec = fieldSpec(name)

    def value(self, newValue=_MISSING_SPECTRAL_FIELD_VALUE):
        attr = _CROSS_SECTION_FIELD_ATTRS[self.name]
        if newValue is _MISSING_SPECTRAL_FIELD_VALUE:
            return getattr(self._crossSections, attr)
        self._crossSections._setField(self.name, newValue)
        return self

    def meta(self):
        values = self.value()
        return {
            "name": self.name,
            "recordName": self.spec.recordName,
            "entity": self.spec.entity,
            "axes": self.spec.axes,
            "dtype": str(self.spec.dtypeObject),
            "unit": self.spec.unit,
            "unitSI": self.spec.unitSI,
            "expectedShape": values.shape,
            "isSet": True,
        }

    def __repr__(self):
        meta = self.meta()
        return (
            f"SpectralField(name={self.name!r}, dtype={meta['dtype']!r}, "
            f"unit={meta['unit']!r}, shape={meta['expectedShape']!r})"
        )


@dataclass
class CrossSectionData:
    r"""Absorption and emission spectra for ASE and pump calculations.

    Wavelength arrays store :math:`\lambda`; matching cross-section arrays
    store :math:`\sigma_a` and :math:`\sigma_e` in ``cm^2``. The wavelength
    unit is kept as supplied, with interpolation helpers handling the common
    ``nm`` table versus ``m`` query mismatch.
    """

    wavelengthsAbsorption: object
    """Wavelength samples for the absorption spectrum."""
    crossSectionAbsorption: object
    r"""Absorption cross sections :math:`\sigma_a` in ``cm^2``."""
    wavelengthsEmission: object
    """Wavelength samples for the emission spectrum."""
    crossSectionEmission: object
    r"""Emission cross sections :math:`\sigma_e` in ``cm^2``."""
    resolution: int = 1
    """Spectral interpolation resolution passed to ``calcPhiASE``."""

    def __post_init__(self):
        for attr in _CROSS_SECTION_FIELD_ATTRS.values():
            setattr(self, attr, np.asarray(getattr(self, attr), dtype=np.float64).reshape(-1))
        self.resolution = int(self.resolution)
        self._validate()

    def _validate(self):
        if self.wavelengthsAbsorption.size != self.crossSectionAbsorption.size:
            raise ValueError("wavelengthsAbsorption and crossSectionAbsorption must have the same length")
        if self.wavelengthsEmission.size != self.crossSectionEmission.size:
            raise ValueError("wavelengthsEmission and crossSectionEmission must have the same length")
        if self.resolution < 1:
            raise ValueError("resolution must be positive")

    def _canonicalFieldName(self, name):
        canonical = _CROSS_SECTION_FIELD_ALIASES.get(name, name)
        if canonical not in _CROSS_SECTION_FIELD_ATTRS:
            known = ", ".join(_CROSS_SECTION_FIELD_ATTRS)
            raise KeyError(f"unknown spectral field '{name}'. Known fields: {known}")
        return canonical

    def _setField(self, name, values):
        canonical = self._canonicalFieldName(name)
        attr = _CROSS_SECTION_FIELD_ATTRS[canonical]
        old = getattr(self, attr)
        setattr(self, attr, np.asarray(values, dtype=fieldSpec(canonical).dtypeObject).reshape(-1))
        try:
            self._validate()
        except Exception:
            setattr(self, attr, old)
            raise

    def getField(self, name):
        return SpectralField(self, self._canonicalFieldName(name))

    def getFields(self):
        return [SpectralField(self, name) for name in _CROSS_SECTION_FIELD_ATTRS]

    @classmethod
    def monochromatic(cls, *, wavelength, crossSectionAbsorption, crossSectionEmission):
        """Build a single-wavelength spectrum for monochromatic workflows."""
        return cls(
            wavelengthsAbsorption=[wavelength],
            crossSectionAbsorption=[crossSectionAbsorption],
            wavelengthsEmission=[wavelength],
            crossSectionEmission=[crossSectionEmission],
            resolution=1,
        )

    @classmethod
    def fromDirectory(cls, path, resolution=1000):
        """Load ``lambda_a``, ``sigma_a``, ``lambda_e``, and ``sigma_e`` text files."""
        root = Path(path)
        return cls(
            wavelengthsAbsorption=np.loadtxt(root / "lambda_a.txt"),
            crossSectionAbsorption=np.loadtxt(root / "sigma_a.txt"),
            wavelengthsEmission=np.loadtxt(root / "lambda_e.txt"),
            crossSectionEmission=np.loadtxt(root / "sigma_e.txt"),
            resolution=resolution,
        )

    def toLaserProperties(self):
        """Wrap the same spectra in the lower-level ``LaserProperties`` store."""
        return LaserProperties(crossSections=self)

    def absorptionAt(self, wavelength):
        r"""Interpolate :math:`\sigma_a` at ``wavelength``."""
        return self._interpolate(self.wavelengthsAbsorption, self.crossSectionAbsorption, wavelength)

    def emissionAt(self, wavelength):
        r"""Interpolate :math:`\sigma_e` at ``wavelength``."""
        return self._interpolate(self.wavelengthsEmission, self.crossSectionEmission, wavelength)

    @staticmethod
    def _interpolate(wavelengths, values, wavelength):
        wavelengths = np.asarray(wavelengths, dtype=np.float64).reshape(-1)
        values = np.asarray(values, dtype=np.float64).reshape(-1)
        query = float(wavelength)
        if wavelengths.size == 1:
            return float(values[0])

        # Existing material files use nm, while pump wavelengths are commonly
        # specified in m. Convert only when the magnitude makes that unambiguous.
        scale = np.nanmax(np.abs(wavelengths))
        if scale > 1e-6 and abs(query) < 1e-6:
            query *= 1e9
        elif scale < 1e-6 and abs(query) > 1e-6:
            query *= 1e-9

        order = np.argsort(wavelengths)
        return float(np.interp(query, wavelengths[order], values[order]))

    def toDict(self):
        """Return the dictionary layout expected by the low-level bindings."""
        return {
            "l_abs": self.wavelengthsAbsorption,
            "l_ems": self.wavelengthsEmission,
            "s_abs": self.crossSectionAbsorption,
            "s_ems": self.crossSectionEmission,
            "l_res": int(self.resolution),
        }

    def openPmdAttributes(self):
        laser = LaserProperties(crossSections=self)
        return {
            "spectralResolution": self.resolution,
            "maxSigmaAbsorption": laser.maxSigmaA,
            "maxSigmaEmission": laser.maxSigmaE,
        }

    def openPmdFields(self, spectralContext):
        for field in self.getFields():
            yield OpenPmdScalarField(
                field.name,
                field.value(),
                spectralContext(field.value()),
                spec=field.spec,
            )


SpectralDecomposition = CrossSectionData


def _positive_normalized(values, name):
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    if values.size == 0 or np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError(f"{name} must contain finite non-negative values")
    total = float(values.sum())
    if total <= 0.0:
        raise ValueError(f"{name} must have positive total weight")
    return values / total


def _unit_vector(value, name):
    value = np.asarray(value, dtype=np.float64).reshape(-1)
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be a finite three-vector")
    length = float(np.linalg.norm(value))
    if length <= 0.0:
        raise ValueError(f"{name} must be non-zero")
    return tuple(value / length)


@dataclass(frozen=True)
class PumpSpectrum:
    """Normalized discrete spectrum for a physical pump."""

    wavelengths: object
    weights: object

    def __post_init__(self):
        wavelengths = np.asarray(self.wavelengths, dtype=np.float64).reshape(-1)
        weights = _positive_normalized(self.weights, "pump spectrum weights")
        if wavelengths.size != weights.size:
            raise ValueError("pump spectrum wavelengths and weights must have the same length")
        if np.any(~np.isfinite(wavelengths)) or np.any(wavelengths <= 0.0):
            raise ValueError("pump spectrum wavelengths must be finite and positive")
        object.__setattr__(self, "wavelengths", wavelengths)
        object.__setattr__(self, "weights", weights)

    @classmethod
    def monochromatic(cls, wavelength):
        return cls([wavelength], [1.0])


@dataclass(frozen=True)
class PumpAngularDistribution:
    """Normalized discrete directions in the injector's inward-local frame."""

    polar_angles: object
    azimuthal_angles: object
    weights: object

    def __post_init__(self):
        polar = np.asarray(self.polar_angles, dtype=np.float64).reshape(-1)
        azimuthal = np.asarray(self.azimuthal_angles, dtype=np.float64).reshape(-1)
        weights = _positive_normalized(self.weights, "pump angular weights")
        if polar.size != azimuthal.size or polar.size != weights.size:
            raise ValueError("pump angular angles and weights must have the same length")
        if np.any(~np.isfinite(polar)) or np.any((polar < 0.0) | (polar >= 0.5 * np.pi)):
            raise ValueError("pump polar angles must be finite and in [0, pi/2)")
        if np.any(~np.isfinite(azimuthal)):
            raise ValueError("pump azimuthal angles must be finite")
        object.__setattr__(self, "polar_angles", polar)
        object.__setattr__(self, "azimuthal_angles", azimuthal)
        object.__setattr__(self, "weights", weights)

    @property
    def polarAngles(self):
        return self.polar_angles

    @property
    def azimuthalAngles(self):
        return self.azimuthal_angles

    @classmethod
    def collimated(cls):
        return cls([0.0], [0.0], [1.0])

    @classmethod
    def uniform_cone(cls, half_angle, *, polar_samples=8, azimuthal_samples=16):
        if half_angle <= 0.0 or half_angle >= 0.5 * np.pi:
            raise ValueError("half_angle must be in (0, pi/2)")
        cos_edges = np.linspace(1.0, np.cos(float(half_angle)), int(polar_samples) + 1)
        polar = np.arccos(0.5 * (cos_edges[:-1] + cos_edges[1:]))
        azimuthal = (np.arange(int(azimuthal_samples)) + 0.5) * (2.0 * np.pi / int(azimuthal_samples))
        theta, phi = np.meshgrid(polar, azimuthal, indexing="ij")
        return cls(theta.reshape(-1), phi.reshape(-1), np.ones(theta.size))


@dataclass(frozen=True)
class UniformPumpProfile:
    """Uniform power density over the selected injector aperture."""

    kind: str = field(default="uniform", init=False)


@dataclass(frozen=True)
class SuperGaussianPumpProfile:
    """Normalized world-space elliptical super-Gaussian profile."""

    radius_u: float
    radius_v: float | None = None
    exponent: float = 40.0
    center: object = (0.0, 0.0, 0.0)
    axis_u: object = (1.0, 0.0, 0.0)
    axis_v: object = (0.0, 1.0, 0.0)
    kind: str = field(default="super-gaussian", init=False)

    def __post_init__(self):
        radius_v = self.radius_u if self.radius_v is None else self.radius_v
        if self.radius_u <= 0.0 or radius_v <= 0.0 or self.exponent <= 0.0:
            raise ValueError("super-Gaussian radii and exponent must be positive")
        center = np.asarray(self.center, dtype=np.float64).reshape(-1)
        if center.shape != (3,) or np.any(~np.isfinite(center)):
            raise ValueError("pump profile center must be a finite three-vector")
        axis_u = np.asarray(_unit_vector(self.axis_u, "axis_u"))
        axis_v = np.asarray(_unit_vector(self.axis_v, "axis_v"))
        if abs(float(np.dot(axis_u, axis_v))) > 1.0e-10:
            raise ValueError("axis_u and axis_v must be orthogonal")
        object.__setattr__(self, "radius_v", float(radius_v))
        object.__setattr__(self, "center", tuple(center))
        object.__setattr__(self, "axis_u", tuple(axis_u))
        object.__setattr__(self, "axis_v", tuple(axis_v))

    @property
    def radiusU(self):
        return self.radius_u

    @property
    def radiusV(self):
        return self.radius_v

    @property
    def axisU(self):
        return self.axis_u

    @property
    def axisV(self):
        return self.axis_v

    def weightAt(self, points):
        points = np.asarray(points, dtype=np.float64)
        relative = points - np.asarray(self.center)
        u = relative @ np.asarray(self.axis_u) / self.radius_u
        v = relative @ np.asarray(self.axis_v) / self.radius_v
        return np.exp(-((u * u + v * v) ** (0.5 * self.exponent)))

    def weight_at(self, points):
        """Evaluate the dimensionless spatial profile at world coordinates."""
        return self.weightAt(points)


@dataclass(frozen=True)
class Pump:
    """Physical pump definition, independent of its numerical injection method."""

    total_power: float
    spectrum: PumpSpectrum
    cross_sections: CrossSectionData
    profile: object = field(default_factory=UniformPumpProfile)
    angular_distribution: PumpAngularDistribution = field(default_factory=PumpAngularDistribution.collimated)
    name: str | None = None

    def __post_init__(self):
        if not np.isfinite(self.total_power) or self.total_power <= 0.0:
            raise ValueError("pump total_power must be finite and positive")
        if not isinstance(self.spectrum, PumpSpectrum):
            raise TypeError("pump spectrum must be PumpSpectrum")
        if not isinstance(self.cross_sections, CrossSectionData):
            raise TypeError("pump cross_sections must be CrossSectionData")
        if not isinstance(self.angular_distribution, PumpAngularDistribution):
            raise TypeError("pump angular_distribution must be PumpAngularDistribution")
        if not isinstance(self.profile, (UniformPumpProfile, SuperGaussianPumpProfile)):
            raise TypeError("pump profile must be UniformPumpProfile or SuperGaussianPumpProfile")


class GaussianPump(Pump):
    """Convenience physical pump with an elliptical super-Gaussian profile."""

    def __init__(
        self,
        *,
        total_power,
        spectrum,
        cross_sections,
        waist,
        exponent=2.0,
        center=(0.0, 0.0, 0.0),
        axis_u=(1.0, 0.0, 0.0),
        axis_v=(0.0, 1.0, 0.0),
        angular_distribution=None,
        name=None,
    ):
        radii = (float(waist), float(waist)) if np.isscalar(waist) else tuple(waist)
        if len(radii) != 2:
            raise ValueError("waist must be a scalar or a two-vector")
        super().__init__(
            total_power=total_power,
            spectrum=spectrum,
            cross_sections=cross_sections,
            profile=SuperGaussianPumpProfile(
                radius_u=radii[0],
                radius_v=radii[1],
                exponent=exponent,
                center=center,
                axis_u=axis_u,
                axis_v=axis_v,
            ),
            angular_distribution=(
                PumpAngularDistribution.collimated()
                if angular_distribution is None
                else angular_distribution
            ),
            name=name,
        )


@dataclass(frozen=True)
class SurfacePumpInjector:
    """Numerical injection method that launches a pump from tagged surfaces."""

    surface_domains: object

    def __post_init__(self):
        domains = (
            (self.surface_domains,)
            if isinstance(self.surface_domains, (str, int))
            else tuple(self.surface_domains)
        )
        if not domains:
            raise ValueError("SurfacePumpInjector requires at least one surface domain")
        object.__setattr__(self, "surface_domains", domains)


@dataclass(frozen=True)
class PlanarPumpRelay:
    """Affine re-imaging link between tagged planar boundary surfaces."""

    exit_domains: object
    entry_domains: object
    flip_u: bool = False
    flip_v: bool = False
    rotation: float = 0.0
    offset: tuple[float, float] = (0.0, 0.0)
    tilt: tuple[float, float] = (0.0, 0.0)
    magnification: float = 1.0
    transmission: float = 1.0

    def __post_init__(self):
        exit_domains = (self.exit_domains,) if isinstance(self.exit_domains, (str, int)) else tuple(self.exit_domains)
        entry_domains = (
            (self.entry_domains,)
            if isinstance(self.entry_domains, (str, int))
            else tuple(self.entry_domains)
        )
        if not exit_domains or not entry_domains:
            raise ValueError("relay requires exit and entry surface domains")
        object.__setattr__(self, "exit_domains", exit_domains)
        object.__setattr__(self, "entry_domains", entry_domains)
        if self.magnification <= 0.0:
            raise ValueError("relay magnification must be positive")
        if self.transmission < 0.0 or self.transmission > 1.0:
            raise ValueError("relay transmission must be in [0, 1]")
        if len(tuple(self.offset)) != 2 or len(tuple(self.tilt)) != 2:
            raise ValueError("relay offset and tilt must be two-vectors")

    @property
    def exitDomains(self):
        return self.exit_domains

    @property
    def entryDomains(self):
        return self.entry_domains

    @property
    def flipU(self):
        return self.flip_u

    @property
    def flipV(self):
        return self.flip_v

    @classmethod
    def retroreflect(cls, domains, *, transmission=1.0):
        return cls(domains, domains, transmission=transmission)


@dataclass(frozen=True)
class MonteCarloPumpSolver:
    """Numerical controls shared by pumps registered on a simulation."""

    ray_count: int = 100_000
    seed: int = 5489
    max_steps: int | None = None

    def __post_init__(self):
        if self.ray_count <= 0:
            raise ValueError("MonteCarloPumpSolver.ray_count must be positive")
        if self.seed < 0 or self.seed >= 2**32:
            raise ValueError("MonteCarloPumpSolver.seed must fit uint32")
        if self.max_steps is not None and self.max_steps < 0:
            raise ValueError("MonteCarloPumpSolver.max_steps must be non-negative")


@dataclass(frozen=True)
class _PumpSource:
    """CamelCase transport adapter assembled by ``Simulation.add_pump``."""

    surfaceDomains: tuple
    totalPower: float
    spectrum: PumpSpectrum
    crossSections: CrossSectionData
    angularDistribution: PumpAngularDistribution
    profile: object
    relays: tuple[PlanarPumpRelay, ...]


@dataclass(frozen=True)
class _PumpProperties:
    """CamelCase compiled-transport pump graph."""

    sources: tuple[_PumpSource, ...]
    rayCount: int
    rngSeed: int
    pumpSteps: int | None


def integrate_pump_profile(topology, surface_domains, profile):
    """Integrate a normalized profile over tagged triangular boundary faces."""
    domains = (surface_domains,) if isinstance(surface_domains, (str, int)) else tuple(surface_domains)
    domain_map = topology.surfaceDomainMap()
    domain_ids = {domain_map.resolve(value) if isinstance(value, str) else int(value) for value in domains}
    mask = (np.asarray(topology.neighborCells) < 0) & np.isin(
        np.asarray(topology.faceBoundaries), list(domain_ids)
    )
    cell_faces = np.argwhere(mask)
    if cell_faces.size == 0:
        raise ValueError("pump profile integration selected no exterior faces")
    indices = np.asarray(topology.facePointIndices)[mask]
    triangles = np.asarray(topology.points, dtype=np.float64)[indices]
    area = 0.5 * np.linalg.norm(
        np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]), axis=1
    )
    barycentric = np.asarray([
        [1 / 3, 1 / 3, 1 / 3],
        [0.059715871789770, 0.470142064105115, 0.470142064105115],
        [0.470142064105115, 0.059715871789770, 0.470142064105115],
        [0.470142064105115, 0.470142064105115, 0.059715871789770],
        [0.797426985353087, 0.101286507323456, 0.101286507323456],
        [0.101286507323456, 0.797426985353087, 0.101286507323456],
        [0.101286507323456, 0.101286507323456, 0.797426985353087],
    ])
    weights = np.asarray([
        0.225,
        0.132394152788506, 0.132394152788506, 0.132394152788506,
        0.125939180544827, 0.125939180544827, 0.125939180544827,
    ])
    points = np.einsum("qb,tbc->tqc", barycentric, triangles)
    values = np.ones(points.shape[:2]) if isinstance(profile, UniformPumpProfile) else profile.weightAt(points)
    return float(np.sum(area * (values @ weights)))


@dataclass
class LaserProperties:
    """Mutable low-level laser-property store.

    Prefer ``CrossSectionData`` for new simulations. This class remains useful
    when code needs the historical ``l_abs``, ``s_abs``, ``l_ems``, ``s_ems``,
    and ``l_res`` handles or aliases used by ``calcPhiASE``.
    """

    crossSections: CrossSectionData | None = None
    """Optional spectral data used to populate the property store."""
    values: dict = field(default_factory=dict)
    """Canonical low-level property values."""

    def __post_init__(self):
        if self.crossSections is not None:
            self.withProperties(**self.crossSections.toDict())

    @classmethod
    def spectral(cls, **kwargs):
        """Create ``LaserProperties`` from explicit spectral arrays."""
        if "absorption" in kwargs and "crossSectionAbsorption" not in kwargs:
            kwargs["crossSectionAbsorption"] = kwargs.pop("absorption")
        if "emission" in kwargs and "crossSectionEmission" not in kwargs:
            kwargs["crossSectionEmission"] = kwargs.pop("emission")
        return cls(crossSections=CrossSectionData(**kwargs))

    @classmethod
    def monochromatic(cls, *, absorption, emission, wavelengthAbsorption=0.0, wavelengthEmission=0.0):
        """Create a single-sample absorption/emission data set."""
        return cls.spectral(
            wavelengthsAbsorption=[wavelengthAbsorption],
            crossSectionAbsorption=[absorption],
            wavelengthsEmission=[wavelengthEmission],
            crossSectionEmission=[emission],
            resolution=1,
        )

    @classmethod
    def fromDirectory(cls, path):
        """Load spectral text files and wrap them as ``LaserProperties``."""
        return cls(crossSections=CrossSectionData.fromDirectory(path))

    def withProperties(self, **properties):
        """Set multiple laser properties and return ``self`` for chaining."""
        for name, value in properties.items():
            self.set(name, value)
        return self

    def get(self, name):
        """Return a ``LaserProperty`` handle by canonical name or alias."""
        canonical = LASER_ALIASES.get(name, name)
        if canonical not in LASER_PROPERTY_SPECS:
            known = ", ".join(LASER_PROPERTY_SPECS)
            raise KeyError(f"unknown laser property '{name}'. Known properties: {known}")
        return LaserProperty(self, LASER_PROPERTY_SPECS[canonical])

    def set(self, name, value):
        """Validate and store one laser property by canonical name or alias."""
        prop = self.get(name)
        if prop.expectedShape == ():
            self.values[prop.name] = prop.dtype.type(value).item()
            return self

        arr = np.asarray(value, dtype=prop.dtype).reshape(-1)
        if arr.size == 0:
            raise ValueError(f"{prop.name} must not be empty")
        self.values[prop.name] = arr
        self._validate_pairs()
        return self

    def listProperties(self):
        """Return metadata for all known laser properties."""
        return [self.get(name).meta() for name in LASER_PROPERTY_SPECS]

    def toDict(self):
        """Return the complete low-level laser dictionary after validation."""
        self.validate(requiredOnly=True)
        return {
            "l_abs": self.values["l_abs"],
            "l_ems": self.values["l_ems"],
            "s_abs": self.values["s_abs"],
            "s_ems": self.values["s_ems"],
            "l_res": int(self.values["l_res"]),
        }

    @property
    def maxSigmaA(self):
        r"""Maximum absorption cross section :math:`\max(\sigma_a)`."""
        self.validate(requiredOnly=True)
        return float(np.max(self.values["s_abs"]))

    @property
    def maxSigmaE(self):
        r"""Maximum emission cross section :math:`\max(\sigma_e)`."""
        self.validate(requiredOnly=True)
        return float(np.max(self.values["s_ems"]))

    @property
    def emissionPeakIndex(self):
        """Index of the largest emission cross-section sample."""
        self.validate(requiredOnly=True)
        return int(np.argmax(self.values["s_ems"]))

    @property
    def absorptionAtEmissionPeak(self):
        """Absorption cross section sampled at the emission peak index."""
        self.validate(requiredOnly=True)
        idx = min(self.emissionPeakIndex, len(self.values["s_abs"]) - 1)
        return float(self.values["s_abs"][idx])

    def validate(self, requiredOnly=False):
        """Check required fields and matching wavelength/cross-section lengths."""
        missing = [
            name for name, spec in LASER_PROPERTY_SPECS.items()
            if spec.required and name not in self.values
        ]
        if missing:
            raise ValueError(f"missing required laser properties: {', '.join(missing)}")
        self._validate_pairs()
        if not requiredOnly:
            for name, value in self.values.items():
                self.set(name, value)
        return self

    def _validate_pairs(self):
        if "l_abs" in self.values and "s_abs" in self.values:
            if len(self.values["l_abs"]) != len(self.values["s_abs"]):
                raise ValueError("l_abs and s_abs must have the same length")
        if "l_ems" in self.values and "s_ems" in self.values:
            if len(self.values["l_ems"]) != len(self.values["s_ems"]):
                raise ValueError("l_ems and s_ems must have the same length")
