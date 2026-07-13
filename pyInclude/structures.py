from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence


@dataclass
class Result:
    phiAse: Sequence[float] = field(default_factory=list)
    mse: Sequence[float] = field(default_factory=list)
    totalRays: Sequence[int] = field(default_factory=list)
    dndtAse: Sequence[float] = field(default_factory=list)
    srmStatus: str = "disabled"
    srmPasses: int = 0
    srmRemainingFraction: float = 0.0
    srmMaxIterations: int = 0
    srmDivergenceStreak: int = 0
