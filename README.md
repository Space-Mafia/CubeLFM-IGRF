# CubeLFM-IGRF

**v0.3.0** publishes the 128 × 48 × 48 × 6 Earth magnetosphere event of **23 October 2002**. The model cold-started at 04:00 UT and reached 16:00 UT on 192 physical cores. It uses the native LFM-Cube panel seam, the corrected magnetic-induction reference, and the Yin–Yang GAMERA global ionospheric potential solver with a branch-local return to the MHD inner boundary.

## Build the 128-grid event

The [clean-build guide](event128/CLEAN_BUILD.md) covers compiler, MPI, MKL and HDF5 dependencies, NASA OMNI input preparation, and a 192-core test that advances the first two physical minutes. The two fixed geometry files ship with the code. The clean entry point is `event128/smoke.sbatch`.

## Download and test

```sh
git clone https://github.com/Space-Mafia/CubeLFM-IGRF.git
cd CubeLFM-IGRF
python3 event128/verify.py
```

This offline test uses Python 3.9 or newer and the standard library. It regenerates all 14 C and Fortran inputs compiled for the event, compares them byte for byte by SHA-256 with the Sugon build records, verifies source mappings for 47 inherited objects, and reproduces the two compact-output writers. It takes no cluster time.

To inspect the actual MHD or MI update, start with the [compiled event source](event128/compiled). The [event guide](event128/README.md) describes the scientific configuration and archived 12-hour result. The [clean-build guide](event128/CLEAN_BUILD.md) describes the source build and the short coupled test.

## Reference performance

| Case | Mesh and allocation | Physical interval | Full MHD wall time | Real-time factor |
|---|---|---|---:|---:|
| v0.3.0, dated global-potential event | 128 × 48 × 48 × 6, 192 cores | 04:00–16:00 UT | 87,718.14 s | 0.49249× |
| v0.2.0, idealized native-MI reference | 64 × 24 × 24 × 6, 96 cores | 12 h | 10,789.84 s | 4.00376× |

The factor is simulated duration divided by the complete MHD process wait, including initialization, output and exit. These cases use different meshes, driving and solvers. The 128-grid job completed 495,712 steps and recorded its final MI update at 43,200 s. [128-grid case and result](event128/CASE.json)

The root-level `source/`, `verify_source.py` and `run_reference.sbatch` belong to the earlier 64-grid case. Its original instructions remain available at the [v0.2.0 tag](https://github.com/Space-Mafia/CubeLFM-IGRF/tree/cubelfm-igrf-v0.2.0). The [release notes](RELEASE_NOTES.md) describe the event version.
