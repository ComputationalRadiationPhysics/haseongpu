# compareSerial reference data

This directory is reserved for reference data generated from the old C++
`compareSerial` integration test.  That test compared the simple serial physics
implementation against the Alpaka execution path.  It was removed during the
openPMD/Tet4 transition, but its output is still useful as a historical physics
reference.

`phiase_reference.npz` stores the generated reference.  It contains the full
serial point buffers for legacy `cuboid` and `cylindrical` inputs:
`{dataset}_phiASE` and `{dataset}_dndtAse`.  Partial sample ranges are not a
valid physics comparison fixture and should not be added back.

## Generator provenance

Last known commit containing `tests/compSerial_Itest.cpp` and
`tests/data/cfg/compSerial.cfg`:

```text
469c87770ed13796f2e82385bcf83528e8aeaf1b
2026-07-07T16:05:31+02:00
refactor cuda architecture pass-through
```

The old test loaded geometry from:

```text
example/c_example/input
```

The data generator seed is the serial implementation seed:

```text
5489
```

In the old source this is `hase::internal::rng{5489u}` in
`src/core/SerialVersion.cpp`.  The old test also set
`hase::random::SeedGenerator` to `1372085211`, but that is secondary for this
fixture because the recorded data comes from `BaseVersionSerial`.

The old test used this config:

```text
[Experiment]
min-rays=1000000
max-rays=1000000
mse-threshold=0.1
reflection=false
monochromatic=true

[Compute]
parallel-mode=single
repetitions=1
adaptive-steps=1
numDevices=1
min-sample-i=0
max-sample-i=<full point buffer>
```

With ``monochromatic=true``, the legacy parser reduced each stored spectrum to
the first absorption/emission cross-section pair before starting the serial
calculation.  For these inputs that is the pair at 905 nm; the remaining 190
table entries were not sampled by the reference calculation.

At generator commit ``469c8777``, ``BaseVersionSerial`` then consumed
``minRaysPerSample`` (one million), that one cross-section pair, and the legacy
mesh's points, prism connectivity, ``betaVolume``, ``betaCells``, ``nTot``, and
``crystalTFluo``.  Its implementation did not consult the maximum-ray, MSE,
reflection, cladding-optics, repetition, or adaptive-step controls.  A clean
historical build invoking ``BaseVersionSerial`` directly reproduced both
stored arrays bit-for-bit for cylindrical indices 0 through 10; across both
full datasets, every stored ``dndtAse`` value also exactly satisfies the
legacy first-pair ``gainPerDensity * phiASE`` calculation.
The temporary full-buffer dump patch itself is not present in the cited
generator commit, so regenerating the complete archive still requires
reconstructing that small test-only patch.

If new reference data is generated from that old branch, record the exact
generator commit, input geometry, backend, RNG seed, command line/config, and
observable beside the data file.

Previous versions use ``mse-threshold=0.1``, calibrated here to RSE ``0.14`` for ``cuboid`` and ``0.10`` for ``cylindrical``.

## Wedge comparison artifacts

Current Tet4 fixtures can be inverted back to the legacy wedge/prism topology
for inspection:

```sh
python3 utils/compare_serial_wedge_projection.py \
  --output-dir /tmp/hase-compareSerial-wedge-artifacts
```

For each dataset this writes `{name}_original_wedge.vtk` and
`{name}_tet4_roundtrip_wedge.vtk`.  Both files contain `betaCells` and
`betaVolume`.  When `phiase_reference.npz` is present, the original wedge file
also contains full point data `serialPhiASE`.
