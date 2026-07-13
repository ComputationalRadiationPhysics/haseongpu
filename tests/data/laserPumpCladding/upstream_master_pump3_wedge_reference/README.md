# laserPumpCladding upstream/master wedge PhiASE reference

This directory stores the deterministic legacy wedge reference used by the
laserPumpCladding Tet4 integration test. `ptTet4.vtk` represents the same
material geometry with Tet4 cells; the comparison is a retained numerical
contract across topology and transport changes.

The reference was generated from a detached `upstream/master` worktree at:

```text
469c87770ed13796f2e82385bcf83528e8aeaf1b refactor cuda architecture pass-through
```

The upstream script was not run through its CLI, because that CLI did not expose
an RNG seed argument.  Instead, the script module was imported and `runExample`
was called directly so the ASE backend seed is part of the fixture contract.

```python
from pathlib import Path
import sys

sys.path.insert(0, "/tmp/hase-laserpump-upstream-master/example/python_example")
import laserPumpCladding

state = laserPumpCladding.runExample(
    phiAseConfigPath=Path(
        "/tmp/hase-laserpump-upstream-master/example/python_example/config/phiASE.yaml"
    ),
    backend="Host_Cpu_CpuOmpBlocks",
    timeSlices=6,
    pumpSteps=3,
    vtkOutputDir=Path("/tmp/hase-laserpump-upstream-master-wedge6"),
    enableAse=True,
    prePump=True,
    rngSeed=5489,
)
```

The build/import environment was:

```text
PYTHONPATH=/tmp/hase-laserpump-upstream-master-build/python:/tmp/hase-laserpump-upstream-master
cwd=/tmp
```

The upstream `example/python_example/config/phiASE.yaml` settings were:

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

The historical reflection implementation places its virtual top reflection
plane at `numberOfLevels * thickness` (0.777... m), although the ten sampled
planes end at 0.7 m. The current Tet4 walker reflects at the physical 0.7 m
boundary. For the retained comparison, two physical SRM passes reproduce the
legacy reflection series within 10%; using the current default eight passes
would count the virtual legacy layer repeatedly. The integration contract
therefore pins `reflectionMaxIterations=2`. The Tet4 conversion now orders
each wedge's base vertices globally so all shared prism faces are conforming.
