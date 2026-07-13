# laserPumpCladding fixed-legacy-reflection wedge PhiASE reference

This directory stores the deterministic fixed-legacy wedge reference used by
the laserPumpCladding Tet4 integration test. `ptTet4.vtk` represents the same
material geometry with Tet4 cells; the comparison is a retained numerical
contract across topology and transport changes.

The reference was generated from a detached `fix-legacy-reflection-plane`
worktree at:

```text
effd8077edccef93a68d818e8a5eb2f0ebdc03b4 Fix legacy reflection plane height
```

The upstream script was not run through its CLI, because that CLI did not expose
an RNG seed argument.  Instead, the script module was imported and `runExample`
was called directly so the ASE backend seed is part of the fixture contract.

```python
from pathlib import Path
import sys

sys.path.insert(0, "/tmp/hase-laserpump-fixed-legacy/example/python_example")
import laserPumpCladding

state = laserPumpCladding.runExample(
    phiAseConfigPath=Path(
        "/tmp/hase-laserpump-fixed-legacy/example/python_example/config/phiASE.yaml"
    ),
    backend="Host_Cpu_CpuOmpBlocks",
    timeSlices=6,
    pumpSteps=3,
    vtkOutputDir=Path("/tmp/hase-laserpump-fixed-legacy-wedge6"),
    enableAse=True,
    prePump=True,
    rngSeed=5489,
)
```

The build/import environment was:

```text
PYTHONPATH=/tmp/hase-laserpump-fixed-legacy-build/python:/tmp/hase-laserpump-fixed-legacy
cwd=/tmp
```

The fixed-legacy `example/python_example/config/phiASE.yaml` settings were:

```yaml
experiment:
  minRaysPerSample: 2000
  maxRaysPerSample: 1000000
  mseThreshold: 0.01087
  adaptiveSteps: 8
  useReflections: true

compute:
  repetitions: 1
  backend: Cuda_NvidiaGpu_GpuCuda
  parallelMode: single
  numDevices: 4
  nPerNode: 1
```

The generated archive `phiase_reference.npz` contains:

- `metadata`: JSON string with commit, command, config, seed, and checksums.
- `phiASE`: `float64` array with shape `(6, 4210)`, one full legacy wedge
  point buffer per simulation step.
- `points`, `cells`, `cellTypes`: wedge geometry arrays that define the
  reference point order.

The stored `meanPhi` values are only diagnostics. The intended comparison is a
volume integral of the current Tet4 field against the legacy wedge field.

Previous versions use ``mseThreshold: 0.01087``, calibrated here to ``relativeStandardErrorThreshold: 0.05``.

## Reference audit notes

The historical parser interpolates both ASE spectra to 1000 equally spaced
samples over 905--1095 nm before transport; the pump uses the original table
values at 940 nm. The current example reproduces those inputs.

The fixed legacy implementation reflects at the physical top plane,
`(numberOfLevels - 1) * thickness` (0.7 m for the ten sampled planes), the
same boundary used by the Tet4 walker. The Tet4 conversion orders each
wedge's base vertices globally so all shared prism faces are conforming.
