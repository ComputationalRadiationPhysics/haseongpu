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


@dataclass
class CrossSectionData:
    """Absorption and emission spectra for ASE and pump calculations.

    Wavelength arrays store :math:`\lambda`; matching cross-section arrays
    store :math:`\sigma_a` and :math:`\sigma_e` in ``cm^2``. The wavelength
    unit is kept as supplied, with interpolation helpers handling the common
    ``nm`` table versus ``m`` query mismatch.
    """

    wavelengthsAbsorption: object
    """Wavelength samples for the absorption spectrum."""
    crossSectionAbsorption: object
    """Absorption cross sections :math:`\sigma_a` in ``cm^2``."""
    wavelengthsEmission: object
    """Wavelength samples for the emission spectrum."""
    crossSectionEmission: object
    """Emission cross sections :math:`\sigma_e` in ``cm^2``."""
    resolution: int = 1
    """Spectral interpolation resolution passed to ``calcPhiASE``."""

    def __post_init__(self):
        self.wavelengthsAbsorption = np.asarray(self.wavelengthsAbsorption, dtype=np.float64).reshape(-1)
        self.crossSectionAbsorption = np.asarray(self.crossSectionAbsorption, dtype=np.float64).reshape(-1)
        self.wavelengthsEmission = np.asarray(self.wavelengthsEmission, dtype=np.float64).reshape(-1)
        self.crossSectionEmission = np.asarray(self.crossSectionEmission, dtype=np.float64).reshape(-1)
        self.resolution = int(self.resolution)
        if self.wavelengthsAbsorption.size != self.crossSectionAbsorption.size:
            raise ValueError("wavelengthsAbsorption and crossSectionAbsorption must have the same length")
        if self.wavelengthsEmission.size != self.crossSectionEmission.size:
            raise ValueError("wavelengthsEmission and crossSectionEmission must have the same length")
        if self.resolution < 1:
            raise ValueError("resolution must be positive")

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
        """Interpolate :math:`\sigma_a` at ``wavelength``."""
        return self._interpolate(self.wavelengthsAbsorption, self.crossSectionAbsorption, wavelength)

    def emissionAt(self, wavelength):
        """Interpolate :math:`\sigma_e` at ``wavelength``."""
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


SpectralDecomposition = CrossSectionData


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
        """Maximum absorption cross section :math:`\max(\sigma_a)`."""
        self.validate(requiredOnly=True)
        return float(np.max(self.values["s_abs"]))

    @property
    def maxSigmaE(self):
        """Maximum emission cross section :math:`\max(\sigma_e)`."""
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


@dataclass(init=False)
class PumpProperties:
    """Pump-beam settings used to raise the excited-state fraction ``beta``.

    ``intensity`` is pump intensity :math:`I` in ``W / cm^2``. ``wavelength``
    is the pump wavelength :math:`\lambda`. Built-in Gaussian pumping also
    uses ``radiusX``, ``radiusY``, and ``exponent`` from ``customProperties``.
    A custom ``solver`` may store its own knobs in the same dictionary.
    """

    intensity: float
    """Pump intensity :math:`I` in ``W / cm^2``."""
    wavelength: float | None
    """Pump wavelength :math:`\lambda`; required for monochromatic data."""
    pumpSubsteps: int
    """Number of time samples used by the built-in pump integrator."""
    customProperties: dict
    """Extensible store for beam shape, reflection, spectra, and solver handles."""

    def __init__(self, *, intensity, pumpSubsteps=100, wavelength=None, customProperties=None, **properties):
        """Create pump settings from core fields plus arbitrary custom properties."""
        self.intensity = float(intensity)
        self.wavelength = None if wavelength is None else float(wavelength)
        self.pumpSubsteps = int(pumpSubsteps)
        self.customProperties = dict(customProperties or {})
        self.customProperties.update(properties)
        self._normalizeProperties()

    @classmethod
    def superGaussian(
        cls,
        *,
        intensity,
        wavelength,
        radiusX,
        pumpDuration=None,
        duration=None,
        pumpSubsteps=100,
        temporaryFluorescence=None,
        solver=None,
        crossSections=None,
        spectralProperties=None,
        gainMedium=None,
        crossSectionAbsorption=None,
        crossSectionEmission=None,
        absorption=None,
        emission=None,
        radiusY=None,
        exponent=40.0,
        backReflection=True,
        reflectivity=1.0,
        customProperties=None,
        **extraProperties,
    ):
        """Create settings for the built-in super-Gaussian pump profile.

        The transverse profile is ``intensity * exp(-r ** exponent)`` with
        radii ``radiusX`` and ``radiusY``. ``backReflection`` and
        ``reflectivity`` control the backward pump pass.
        """
        custom = dict(customProperties or {})
        custom.update(extraProperties)
        custom.update(
            {
                "radiusX": radiusX,
                "radiusY": radiusX if radiusY is None else radiusY,
                "exponent": exponent,
                "pumpDuration": pumpDuration if pumpDuration is not None else duration,
                "temporaryFluorescence": temporaryFluorescence,
                "solver": solver,
                "crossSections": crossSections,
                "spectralProperties": spectralProperties,
                "gainMedium": gainMedium,
                "crossSectionAbsorption": crossSectionAbsorption if crossSectionAbsorption is not None else absorption,
                "crossSectionEmission": crossSectionEmission if crossSectionEmission is not None else emission,
                "backReflection": backReflection,
                "reflectivity": reflectivity,
            }
        )
        return cls(
            intensity=intensity,
            wavelength=wavelength,
            pumpSubsteps=pumpSubsteps,
            customProperties=custom,
        )

    def _normalizeProperties(self):
        if "pumpDuration" not in self.customProperties and self.customProperties.get("duration") is not None:
            self.customProperties["pumpDuration"] = self.customProperties["duration"]
        if self.crossSections is None and self.spectralProperties is not None:
            self.customProperties["crossSections"] = self.spectralProperties
        if self.spectralProperties is None and self.crossSections is not None:
            self.customProperties["spectralProperties"] = self.crossSections
        if self.crossSections is None:
            if self.crossSectionAbsorption is None or self.crossSectionEmission is None:
                raise ValueError("PumpProperties requires crossSections or monochromatic cross-section values")
            if self.wavelength is None:
                raise ValueError("PumpProperties requires wavelength when crossSections is not provided")
            self.customProperties["crossSections"] = CrossSectionData(
                wavelengthsAbsorption=[self.wavelength],
                crossSectionAbsorption=[self.crossSectionAbsorption],
                wavelengthsEmission=[self.wavelength],
                crossSectionEmission=[self.crossSectionEmission],
                resolution=1,
            )
            self.customProperties["spectralProperties"] = self.crossSections
        if self.wavelength is None:
            self.wavelength = float(self.crossSections.wavelengthsAbsorption[0])
        if self.solver is None:
            self._validateGaussianPumpParameters()
        if self.pumpDuration is not None and self.pumpDuration <= 0.0:
            raise ValueError("pumpDuration must be positive")
        if self.pumpSubsteps < 2:
            raise ValueError("pumpSubsteps must be at least 2")
        if self.temporaryFluorescence is not None and self.temporaryFluorescence <= 0.0:
            raise ValueError("temporaryFluorescence must be positive")

    def _validateGaussianPumpParameters(self):
        radius_x = self.getProperty("radiusX")
        radius_y = self.radiusY
        if radius_x is None:
            raise ValueError("default Gaussian pump solver requires radiusX")
        if radius_y is None:
            raise ValueError("default Gaussian pump solver requires radiusY or radiusX")
        if float(radius_x) <= 0.0 or float(radius_y) <= 0.0:
            raise ValueError("pump radii must be positive")

    def __getattr__(self, name):
        if "customProperties" in self.__dict__ and name in self.customProperties:
            return self.customProperties[name]
        if name == "crossSections":
            return self.customProperties.get("crossSections")
        if name == "spectralProperties":
            return self.customProperties.get("spectralProperties")
        if name == "radiusY":
            return self.customProperties.get("radiusY", self.customProperties.get("radiusX"))
        if name == "exponent":
            return self.customProperties.get("exponent", 40.0)
        if name == "backReflection":
            return self.customProperties.get("backReflection", True)
        if name == "reflectivity":
            return self.customProperties.get("reflectivity", 1.0)
        if name == "extraction":
            return self.customProperties.get("extraction", False)
        if name in {"pumpDuration", "temporaryFluorescence", "solver", "gainMedium", "duration", "crossSectionAbsorption", "crossSectionEmission"}:
            return self.customProperties.get(name)
        raise AttributeError(name)

    def withProperty(self, name, value):
        """Set one custom pump property and return ``self``."""
        self.customProperties[name] = value
        return self

    def withProperties(self, **properties):
        """Set several custom pump properties and return ``self``."""
        self.customProperties.update(properties)
        return self

    def getProperty(self, name, default=None):
        """Read a custom pump property without raising ``AttributeError``."""
        return self.customProperties.get(name, default)

    def activeDuration(self, timeFrame):
        """Return the pump duration used for one integration frame."""
        if self.pumpDuration is not None:
            return float(self.pumpDuration)
        if timeFrame is None:
            raise ValueError("timeFrame is required when PumpProperties.pumpDuration is not set")
        return float(timeFrame)

    def intensityAt(self, points):
        """Evaluate the super-Gaussian intensity profile at ``(x, y)`` points."""
        self._validateGaussianPumpParameters()
        points = np.asarray(points, dtype=np.float64)
        r = np.sqrt((points[:, 0] ** 2) / (self.radiusY ** 2) + (points[:, 1] ** 2) / (self.radiusX ** 2))
        return self.intensity * np.exp(-(r ** self.exponent))

    def toDict(self, timeFrame=None):
        """Return the low-level pump dictionary consumed by ``pumping.py``."""
        self._validateGaussianPumpParameters()
        sigma_abs = self.crossSections.absorptionAt(self.wavelength)
        sigma_ems = self.crossSections.emissionAt(self.wavelength)
        return {
            "s_abs": np.asarray([sigma_abs], dtype=np.float64),
            "s_ems": np.asarray([sigma_ems], dtype=np.float64),
            "I": float(self.intensity),
            "T": self.activeDuration(timeFrame),
            "wavelength": float(self.wavelength),
            "rx": float(self.radiusX),
            "ry": float(self.radiusY),
            "exp": float(self.exponent),
        }

    def modeDict(self):
        """Return low-level flags for extraction and backward reflection."""
        return {
            "BRM": int(self.backReflection),
            "R": float(self.reflectivity),
            "extr": int(self.extraction),
        }
