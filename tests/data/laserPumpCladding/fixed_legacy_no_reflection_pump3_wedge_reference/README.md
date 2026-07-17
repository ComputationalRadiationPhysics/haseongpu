# laserPumpCladding fixed-legacy no-reflection wedge reference

This directory stores a deterministic six-step legacy wedge reference for the
current Tet4 `laserPumpCladding` regression with surface reflections disabled.
It is separate from, and does not replace, the existing reflected reference.

The source checkout was a detached worktree at exactly:

```text
effd8077edccef93a68d818e8a5eb2f0ebdc03b4 Fix legacy reflection plane height
```

`generate_reference.py` rejects any other `git rev-parse HEAD`. The legacy CLI
did not expose an RNG seed, so `run_legacy.py` imports the example and calls
`runExample` directly with seed 5489. The copied configuration explicitly sets
`useReflections: false`; all other pump and run parameters match the reflected
six-step pump3 fixture.

## Container and legacy build

The reference was generated with Docker image
`hase-no-reflection-toolchain:ubuntu24.04`, image ID
`sha256:8cd6743456e7ca8844f28c84e9aedae35429322eb21ac855bb01a2fdf15407a6`,
built from the checked-in `Dockerfile`:

```bash
fixture_dir="$PWD/tests/data/laserPumpCladding/fixed_legacy_no_reflection_pump3_wedge_reference"
docker build --pull -t hase-no-reflection-toolchain:ubuntu24.04 "$fixture_dir"
git worktree add --detach /tmp/hase-no-reflection-legacy-effd8077 \
    effd8077edccef93a68d818e8a5eb2f0ebdc03b4
mkdir -p /tmp/hase-no-reflection-legacy-build /tmp/hase-no-reflection-home

docker run --rm --user "$(id -u):$(id -g)" \
    -e HOME=/tmp-home \
    -v /tmp/hase-no-reflection-legacy-effd8077:/src:ro \
    -v /tmp/hase-no-reflection-legacy-build:/build \
    -v /tmp/hase-no-reflection-home:/tmp-home \
    hase-no-reflection-toolchain:ubuntu24.04 bash -lc '
      cmake -S /src -B /build -G Ninja \
        -DDISABLE_MPI=ON \
        -DHASE_BUILD_RELEASE=ON \
        -DHASE_NATIVE_OPTIMIZATIONS=OFF \
        -DHASE_ENABLE_PYTHON=ON \
        -DHASE_BUILD_PhiAse=ON \
        -DHASE_TESTING=OFF \
        -Dalpaka_DEP_HWLOC=OFF
      cmake --build /build --parallel 4
    '
```

The recorded toolchain is GCC 13.3.0, CMake 3.28.3, and Python 3.12.3 on
`ubuntu:24.04` base image ID
`sha256:786a8b558f7be160c6c8c4a54f9a57274f3b4fb1491cf65146521ae77ff1dc54`.

## Generation

From the current repository root:

```bash
fixture_dir="$PWD/tests/data/laserPumpCladding/fixed_legacy_no_reflection_pump3_wedge_reference"
legacy_root=/tmp/hase-no-reflection-legacy-effd8077
legacy_git_common_dir="$(git -C "$legacy_root" rev-parse --path-format=absolute --git-common-dir)"
mkdir -p /tmp/hase-no-reflection-legacy-output

docker run --rm --user "$(id -u):$(id -g)" --cpus=12 \
    -e HOME=/tmp-home -e PYTHONNOUSERSITE=1 \
    -e PYTHONPATH=/build/python:/src \
    -v "$legacy_root":/src:ro \
    -v /tmp/hase-no-reflection-legacy-build:/build:ro \
    -v /tmp/hase-no-reflection-home:/tmp-home \
    -v "$fixture_dir":/generation:ro \
    -v /tmp/hase-no-reflection-legacy-output:/output \
    hase-no-reflection-toolchain:ubuntu24.04 \
    python /generation/run_legacy.py

docker run --rm --user "$(id -u):$(id -g)" \
    -e HOME=/tmp-home -e PYTHONNOUSERSITE=1 \
    -v "$PWD":/repo \
    -v "$legacy_root":"$legacy_root":ro \
    -v "$legacy_git_common_dir":"$legacy_git_common_dir":ro \
    -v /tmp/hase-no-reflection-legacy-output:/vtk:ro \
    -v /tmp/hase-no-reflection-home:/tmp-home \
    hase-no-reflection-toolchain:ubuntu24.04 \
    python /repo/tests/data/laserPumpCladding/fixed_legacy_no_reflection_pump3_wedge_reference/generate_reference.py \
      --vtk-dir /vtk \
      --legacy-root "$legacy_root" \
      --output /repo/tests/data/laserPumpCladding/fixed_legacy_no_reflection_pump3_wedge_reference/phiase_reference.npz
```

The container imports the legacy frontend from `/src/HASEonGPU.py`, bindings
from `/build/python/HASEonGPU_Bindings`, and runs the
`Host_Cpu_CpuOmpBlocks` backend. The six wedge-point integrals are:

```text
0
8.353280856963568e21
1.6607832996474322e22
2.476070130202229e22
2.4314895143843087e22
2.3894669267551983e22
```

The archive contains the complete `(6, 4210)` `phiASE` point buffer and the
`(4210, 3)` points, `(7308, 6)` wedge cells, cell types, plus JSON metadata.

## Source audit

```text
legacy laserPumpCladding.py  feb039755932939f3a5455a40616f097166a96083582a077e45e75e0936afc3e
legacy pt.mat                afab3241bb89045a2234006f4c2eca26194bef761e1265e78c5115e6898ed74a
phiASE-no-reflections.yaml   312b269611c464fb2e35a26a2fa09380191596ec11abf2efd52441243a38a9d0
Dockerfile                   3bd9634bdb7aa41870776ca38736065d78d47ff6048d0e294c6292915b6cf61f
phiase_reference.npz         cfc628a32adebea97eecd742e7e04df34f61fd82976c4955fc2769073930c1a2
```

The `.npz` metadata also records SHA-256 checksums for all six generated VTK
snapshots. `pt.mat` is the authoritative legacy geometry input; its checksum is
part of the fixture contract even though the packed archive stores the emitted
wedge geometry directly.

## Current Tet4 container validation and bias audit

The current implementation was independently configured and built in the same
Docker image. It used bundled openPMD 0.17.0 with its HDF5 1.14.6 provider:

```bash
current_build=/tmp/hase-no-reflection-current-build
mkdir -p "$current_build" /tmp/hase-no-reflection-home
docker run --rm --user "$(id -u):$(id -g)" --cpus=12 \
    -e HOME=/tmp-home -e CMAKE_BUILD_PARALLEL_LEVEL=4 \
    -v "$PWD":/src:ro \
    -v "$current_build":/build \
    -v /tmp/hase-no-reflection-home:/tmp-home \
    hase-no-reflection-toolchain:ubuntu24.04 bash -lc '
      cmake -S /src -B /build -G Ninja \
        -DDISABLE_MPI=ON \
        -DHASE_BUILD_RELEASE=ON \
        -DHASE_NATIVE_OPTIMIZATIONS=OFF \
        -DHASE_ENABLE_PYTHON=ON \
        -DHASE_TESTING=ON \
        -Dalpaka_DEP_HWLOC=OFF \
        -DHASE_OPENPMD_PROVIDER=bundled \
        -DHASE_OPENPMD_USE_ADIOS2=OFF \
        -DHASE_OPENPMD_USE_HDF5=ON \
        -DHASE_OPENPMD_FETCH_HDF5=ON \
        -DHASE_OPENPMD_USE_SST=OFF \
        -DHASE_OPENPMD_BUILD_PYTHON_BINDINGS=ON
      cmake --build /build --parallel 4
    '

docker run --rm --user "$(id -u):$(id -g)" --cpus=12 \
    -e HOME=/tmp-home -e PYTHONNOUSERSITE=1 -e PYTHONPATH=/src \
    -v "$PWD":/src:ro \
    -v "$current_build":/build:ro \
    -v /tmp/hase-no-reflection-home:/tmp-home \
    -v "$current_build/python/pyInclude/_native_config.py":/src/pyInclude/_native_config.py:ro \
    hase-no-reflection-toolchain:ubuntu24.04 \
    pytest -q /src/tests/python/simulation/test_laserPumpCladdingNoReflection.py
```

The final test result after the fix was `3 passed in 11.87s`.

Before the fix, increasing the ray count did not reduce the step-2 integral
bias. The fixture-matched adaptive run and a fixed-count run both stayed about
13.5% low, while four times as many rays halved the mean relative sampling
error without changing the integral bias:

| Current run | Step-2 Tet integral | Error vs legacy | Mean RSE | Maximum cell visits |
| --- | ---: | ---: | ---: | ---: |
| adaptive, up to 1M | `7.225260654557545e21` | `-13.5039%` | `0.07078` | `1,000,000` |
| fixed 1M | `7.215956628810387e21` | `-13.6153%` | `0.07097` | `1,000,000` |
| fixed 4M | `7.215963679037562e21` | `-13.6152%` | `0.03540` | `4,000,000` |

The systematic error was an implementation defect, not a changed geometric
domain. `HostMesh` built `betaVolumePrefix`, the source-sampling CDF, when its
initial `betaVolume` was all zero. The compiled time-step runner updated
`betaVolume` after pumping but did not rebuild the CDF. Binary search over the
stale all-zero prefix therefore selected the final Tet for every source ray;
the exact ray-count value in the final cell exposed the failure. Dynamic-only
openPMD updates had the same missing refresh.

A paired PhiASE-only experiment replaced every three Tet beta values with the
corresponding legacy six-vertex wedge mean. This changed the beta-volume
integral by `-0.04049%` and PhiASE by `-0.03956%`, ruling out intra-wedge beta
mapping as the material source of the 13.5% error. With a correctly constructed
CDF, the direct Tet run was already within `+0.08375%` of legacy.

The fix updates `betaVolume` through one setter that immediately rebuilds the
raw, unnormalized CDF. It is used for both time-stepped assignments and for
dynamic openPMD iterations. Post-fix results were:

| Current run | Step-2 Tet integral | Error vs legacy | Maximum cell visits |
| --- | ---: | ---: | ---: |
| fixture-matched adaptive | `8.360277067897611e21` | `+0.08375%` | `2,212` |
| fixed 1M | `8.402970443687753e21` | `+0.59485%` | `2,185` |

The six fixture-matched post-fix errors were `0.0000%`, `+0.0838%`,
`+0.2689%`, `+0.4158%`, `+0.4513%`, and `+0.3963%`, all well inside the
unchanged 5% integral tolerance.
