from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence


@dataclass
class Result:
    phiAse: Sequence[float] = field(default_factory=list)
    mse: Sequence[float] = field(default_factory=list)
    totalRays: Sequence[int] = field(default_factory=list)
    dndtAse: Sequence[float] = field(default_factory=list)
