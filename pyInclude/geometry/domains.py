# Copyright 2026 Tim Hanel
#
# This file is part of HASEonGPU
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Domain metadata and optics helpers for explicit volume meshes."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class SurfaceOptics:
    """Optical properties for one positive surface-domain id."""

    reflectivity: float = 0.0
    n_inside: float = 1.0
    n_outside: float = 1.0


@dataclass(frozen=True)
class DomainMap:
    """Resolve named domains and gmsh physical tags to HASE domain ids."""

    names: dict[int, str] = field(default_factory=dict)

    @classmethod
    def fromGmsh(cls, gmsh, dimension):
        return cls({int(tag): str(name) for tag, name in gmsh.physical_names.get(int(dimension), {}).items()})

    def resolve(self, domain):
        """Return a positive integer domain id for ``domain``."""
        if isinstance(domain, str):
            matches = [tag for tag, name in self.names.items() if name == domain]
            if not matches:
                raise KeyError(f"unknown domain name '{domain}'")
            if len(matches) > 1:
                raise ValueError(f"domain name '{domain}' is ambiguous")
            return int(matches[0])
        resolved = int(domain)
        if resolved <= 0:
            raise ValueError(f"domain ids must be positive, got {resolved}")
        return resolved

    def arrays(self, opticsByDomain):
        """Build backend surface optics arrays indexed by positive domain id."""
        return surfaceOpticsArrays(self, opticsByDomain)


def surfaceOpticsArrays(domainMap, opticsByDomain, *, minSize=0):
    """Build backend surface optics arrays indexed by positive domain id."""
    resolved = {}
    for domain, optics in opticsByDomain.items():
        resolved[domainMap.resolve(domain)] = _coerceOptics(optics)
    size = max(int(minSize), (max(resolved) + 1) if resolved else 0)
    reflectivity = np.zeros(size, dtype=np.float32)
    inside = np.ones(size, dtype=np.float32)
    outside = np.ones(size, dtype=np.float32)
    for domain, optics in resolved.items():
        reflectivity[domain] = np.float32(optics.reflectivity)
        inside[domain] = np.float32(optics.n_inside)
        outside[domain] = np.float32(optics.n_outside)
    return reflectivity, inside, outside


def _coerceOptics(value):
    if isinstance(value, SurfaceOptics):
        return value
    if isinstance(value, dict):
        data = dict(value)
        if "nInside" in data:
            data["n_inside"] = data.pop("nInside")
        if "nOutside" in data:
            data["n_outside"] = data.pop("nOutside")
        return SurfaceOptics(**data)
    try:
        reflectivity, n_inside, n_outside = value
    except (TypeError, ValueError) as exc:
        raise TypeError("surface optics must be SurfaceOptics, a mapping, or a 3-tuple") from exc
    return SurfaceOptics(reflectivity=reflectivity, n_inside=n_inside, n_outside=n_outside)


SurfaceDomainMap = DomainMap
