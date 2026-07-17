# laserPumpCladding true-volume-cladding reference

This fixture is independent of the frozen reflected laserPumpCladding
regression. It tests the legacy and Tet4 material-domain contract with actual
cladding volume cells and with reflections disabled.

The reference is a static PhiASE-only calculation. This deliberately removes
pump interpolation, time integration, beta remapping, and reflected-ray bias
from the material comparison.

## Geometry and material contract

The cladding shell is the set of base triangles that possess an exterior mesh
edge. For `legacy/pt.mat`, that rule selects:

- 24 of 812 base triangle columns;
- 216 of 7,308 wedges over nine z intervals;
- 648 of 21,924 Tet4 cells after mapping every wedge to three children;
- 3.615535851812808 cm³ of 19.566038838161607 cm³, or
  18.4786296384% of the volume.

The shell-triangle centroids lie at radii 2.775--2.805 cm; the largest gain
triangle centroid radius is 2.623 cm. The selection is therefore a contiguous
outer ring rather than an arbitrary set of boundary-touching cells.

`betaVolume` is 0.1 in gain volumes and exactly zero in cladding volumes. This
is important because the legacy importance sampler does not automatically
exclude cladding emitters: source exclusion comes from zero per-volume beta.
`betaCells` remains uniformly 0.1. Shared core/cladding vertices cannot carry
two point values, but they do not control static PhiASE source selection or
segment gain; the authoritative per-wedge/per-Tet `betaVolume` does. A
cladding segment replaces Yb gain with
`exp(-claddingAbsorption * pathLength)`.

Three paired absorption values are stored:

- `0.0`: material-off control;
- `5.5 cm^-1`: the physical laserPumpCladding value used for the cross-model
  comparison;
- `55.0 cm^-1`: a clearly separated high-contrast instrumentation case that
  must fail if cladding types or cladding absorption are ignored.

This contract covers **volumetric cladding absorption only**. The legacy model
does not represent a refractive optical interface between an internal core and
cladding, so this fixture does not claim to test one.

## Legacy provenance and generation

The generator requires the fixed legacy source commit:

```text
effd8077edccef93a68d818e8a5eb2f0ebdc03b4 Fix legacy reflection plane height
```

It imports `example/python_example/laserPumpCladding.py`, verifies that the
frontend and binding resolve inside that worktree, verifies the Git commit,
derives the shell from connectivity, and checks zero cladding source beta
before every PhiASE call.

The reference build and run used isolated Ubuntu 24.04 Docker state. With the
legacy worktree at `/tmp/hase-true-cladding-legacy` and this worktree at
`/tmp/hase-true-cladding-current`, the essential commands were:

```bash
docker run -d --name hase-true-clad-legacy \
  -v /tmp/hase-true-cladding-legacy:/src \
  -v /home/th/workspace/myhaseonpu/.git:/home/th/workspace/myhaseonpu/.git:ro \
  -w /src ubuntu:24.04 sleep infinity

docker exec hase-true-clad-legacy bash -lc \
  'apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
   build-essential cmake ninja-build git python3-dev python3-numpy \
   python3-scipy python3-yaml pybind11-dev libomp-dev'

docker exec hase-true-clad-legacy cmake \
  -S /src -B /src/build-true-clad -G Ninja \
  -DDISABLE_MPI=ON -DHASE_ENABLE_PYTHON=ON \
  -DHASE_SELECT_BACKEND_ALPAKA=ON \
  -Dalpaka_DEP_CUDA=OFF -Dalpaka_DEP_HIP=OFF -Dalpaka_DEP_TBB=OFF \
  -Dalpaka_DEP_OMP=ON -Dalpaka_EXEC_CpuSerial=ON \
  -Dalpaka_EXEC_CpuOmpBlocks=ON -Dalpaka_EXEC_TbbBlocks=OFF \
  -DHASE_NATIVE_OPTIMIZATIONS=OFF -DCMAKE_BUILD_TYPE=Release

docker exec hase-true-clad-legacy \
  cmake --build /src/build-true-clad -j2

docker cp \
  /tmp/hase-true-cladding-current/tests/data/laserPumpCladding/\
true_cladding_no_reflection_phiase_reference/generate_reference.py \
  hase-true-clad-legacy:/tmp/generate_reference.py

docker exec hase-true-clad-legacy bash -lc \
  'cd /tmp && OMP_NUM_THREADS=8 python3 /tmp/generate_reference.py \
   --legacy-repo /src --legacy-build /src/build-true-clad \
   --output /tmp/phiase_reference.npz --rays-per-sample 2000'

docker cp hase-true-clad-legacy:/tmp/phiase_reference.npz \
  /tmp/hase-true-cladding-current/tests/data/laserPumpCladding/\
true_cladding_no_reflection_phiase_reference/phiase_reference.npz
```

The legacy run uses `Host_Cpu_CpuOmpBlocks`, RNG seed 5489, 2,000 rays per
sample, one monochromatic sample at 1030 nm, `sigma_a=0.5e-20 cm²`,
`sigma_e=1.0e-20 cm²`, and no reflections.

## Stored data

`phiase_reference.npz` contains:

- `metadata`: JSON provenance, parameters, geometry, and partition integrals;
- `phiASE`, `mse`, `totalRays`: three full 4,210-point legacy results;
- `points`, `cells`, `cellTypes`: the 7,308-wedge geometry;
- `baseCladdingCellTypes`, `wedgeCladdingCellTypes`: audited material tags;
- `wedgeBetaVolume`: level-major source beta, exactly zero in cladding.

The generated archive SHA-256 is
`9b1823c8eb029afe5ec7970c9c1fb7007f1ca86d81c487f84c69cbb73370d05f`.

The legacy result is projected by the arithmetic mean of the six wedge
vertices, multiplied by wedge volume. The current result is a native
21,924-cell Tet field and is integrated directly with exact Tet volumes. The
cross-model regression gates the native whole-volume integral. Gain and
cladding partitions remain diagnostics because the legacy nodal projection
smears shared core/cladding interface values, whereas the current partition is
cell-local.

The Docker validation produced the following integrals. Bias is
`(current / legacy - 1)`:

| absorption (cm⁻¹) | partition | legacy | current Tet4 | bias |
|---:|:---|---:|---:|---:|
| 0.0 | total | 2.3207021e23 | 2.3820479e23 | +2.64% |
| 0.0 | gain | 2.0891168e23 | 2.1749577e23 | +4.11% |
| 0.0 | cladding | 2.3158534e22 | 2.0709018e22 | -10.58% |
| 5.5 | total | 2.2525996e23 | 2.2255870e23 | -1.20% |
| 5.5 | gain | 2.0632183e23 | 2.1414378e23 | +3.79% |
| 5.5 | cladding | 1.8938130e22 | 8.4149201e21 | -55.57% |
| 55.0 | total | 2.2392766e23 | 2.1334882e23 | -4.72% |
| 55.0 | gain | 2.0538431e23 | 2.1218508e23 | +3.31% |
| 55.0 | cladding | 1.8543345e22 | 1.1637480e21 | -93.72% |

At the physical 5.5 cm⁻¹ value, legacy cladding-partition attenuation is
18.22% and current native cladding-partition attenuation is 59.37%. The total
attenuation is 2.93% and 6.57%, respectively. This difference is retained as
an explicit discretization finding; it is not hidden by substituting a nodal
projection for the native Tet observable. The 55 cm⁻¹ case is an
instrumentation guard, not a physical calibration point.

The maximum legacy MSE was 0.12441 at 2,000 rays per sample. Current four
million-ray maximum relative standard errors were 0.03310, 0.03310, and
0.04264 for absorption 0.0, 5.5, and 55.0 cm⁻¹.

## Current Docker validation

The current source was built separately in an Ubuntu 24.04 container with a
bundled openPMD-api 0.17.0 provider and system HDF5. The focused validation is:

```bash
docker run -d --name hase-true-clad-current \
  -v /tmp/hase-true-cladding-current:/src -w /src \
  ubuntu:24.04 sleep infinity

docker exec hase-true-clad-current bash -lc \
  'apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
   build-essential cmake ninja-build git python3-dev python3-pip \
   python3-numpy python3-scipy python3-yaml python3-pytest \
   libhdf5-dev libomp-dev nlohmann-json3-dev libtoml11-dev pkg-config && \
   python3 -m pip install --break-system-packages pybind11==2.13.6 && \
   mkdir -p /usr/local/lib/cmake && \
   ln -sfn /usr/local/lib/python3.12/dist-packages/pybind11/share/cmake/pybind11 \
     /usr/local/lib/cmake/pybind11'

docker exec hase-true-clad-current cmake \
  -S /src -B /src/build-true-clad -G Ninja \
  -DDISABLE_MPI=ON -DHASE_ENABLE_PYTHON=ON \
  -DHASE_SELECT_BACKEND_ALPAKA=ON \
  -Dalpaka_DEP_CUDA=OFF -Dalpaka_DEP_HIP=OFF -Dalpaka_DEP_TBB=OFF \
  -Dalpaka_DEP_OMP=ON -Dalpaka_EXEC_CpuSerial=ON \
  -Dalpaka_EXEC_CpuOmpBlocks=ON -Dalpaka_EXEC_TbbBlocks=OFF \
  -DHASE_NATIVE_OPTIMIZATIONS=OFF \
  -DHASE_OPENPMD_PROVIDER=bundled \
  -DHASE_OPENPMD_USE_ADIOS2=OFF -DHASE_OPENPMD_USE_HDF5=ON \
  -DHASE_OPENPMD_USE_SST=OFF -DHASE_OPENPMD_FETCH_HDF5=OFF \
  -DHASE_OPENPMD_SUPERBUILD=OFF \
  -DHASE_OPENPMD_BUILD_PYTHON_BINDINGS=ON \
  -DCMAKE_BUILD_TYPE=Release

docker exec hase-true-clad-current \
  cmake --build /src/build-true-clad -j2

docker exec hase-true-clad-current bash -lc \
  'cd /tmp && \
   export HASE_RUNTIME_DIR=/src/build-true-clad && \
   export HASE_OPENPMD_PYTHONPATH=/src/build-true-clad/hase-openpmd-provider/install/lib/python3.12/site-packages && \
   export LD_LIBRARY_PATH=/src/build-true-clad/hase-openpmd-provider/install/lib && \
   export PYTHONPATH=$HASE_OPENPMD_PYTHONPATH:/src && \
   PYTHONNOUSERSITE=1 \
   python3 -m pytest -vv -s \
   /src/tests/python/simulation/test_laserPumpCladdingTrueCladding.py'
```

An import preflight resolved the frontend at `/src/HASEonGPU.py`, openPMD at
the container provider path, and the executable at
`/src/build-true-clad/bin/hase-cpp`. The test asserts that all three Tet
children inherit each wedge tag, classified volumes agree, and cladding cells
have zero source beta.
