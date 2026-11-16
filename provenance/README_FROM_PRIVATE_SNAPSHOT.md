# CubeLFM-IGRF v0.2.0 — accepted 64-grid source snapshot

This Earth research release records the complete **12-hour IGRF + intrinsic tilt + rotating MI** case accepted on 15–16 November 2026. Use the `CubeLFM-IGRF` branch and tag `cubelfm-igrf-v0.2.0`. This is separate from the repository's Jupiter branches.

## Model and reference case

| Setting | Value |
|---|---|
| MHD mesh | Original mapped LFM-Cube, 64 × 24 × 24 × 6; no Yin–Yang or spherical-coordinate MHD mesh substitution |
| Parallelism | 96 MPI ranks / 96 distinct physical cores / one node; OMP=1; layout4×2×2 per face |
| Main field | Pinned Yin–Yang release IGRF coefficients, active degree5; intrinsic tilt retained |
| Rotation | 43200s: one full turn in12h, phase0 |
| Solar wind | Bz +5nT for first7200s, −5nT for next36000s; original convected real32/NINT driver |
| Ionosphere | Native LFM ITMAX100 / ERRMAX1e-3; finite nonconvergence warns and records converged=false |
| Conductance | SigmaP5S, Hall0 |
| Coupling | Natural10s potential solves, held Phi; no MI timestep clipping |
| MHD controls | BDIFFUSE total jump factor0.5; PDM0.75; Boris1100km/s; CFL0.3; density floor0 |
| Archives / restart | MHD+MI every120s; full checkpoint every3600s |

The FAC correction subtracts the complete current-time planetary CT reference before donor CT recovery/curl/projection. It does not impose a low-latitude current mask or threshold. Formal Yin–Yang field/trace/mapping code is reused; native LFM MI remains the solver. The reference solver had finite warned solves; these are not described as converged. Invalid/nonfinite inputs, mapping and communication errors remain fatal.

## Actual qualification

The accepted optimized run advanced0..43200s in246975steps on96cores. Full immediate process wait, including initialization/output/exit, was10789.844918493181s: **4.003764681173281×real time**. It passed the frozen numerical audit and full MHD+MI word comparison with the accepted2.8816891356942693× reference:361archive epochs,34656rank field files,4321natural MI solves and12checkpoint epochs. Different allocations preclude an isolated paired speedup claim.

The **6–7× speed target is still under development**. Later short-test and shared-memory experiments are not part of this release. No global production defaults or Jupiter branches are promoted by this source publication. Yin–Yang full-MHD quantitative benchmarking was not completed; the user accepted the Cube science and stopped further alignment.

## Reading the code

`source/` contains all49 explicitly selected C/Fortran translation units copied byte-for-byte from their executed source paths. Each subdirectory is named for its compiled object. `headers/` retains distinct included header versions by hash. `SOURCE_MANIFEST.json` maps them to exact source/object identities and recorded compilation commands. `BUILD_DEPENDENCIES.json` records the immutable native libraries, flags and module authorities; their inherited source remains under `profiles/deep-tail-256/source`, `cubesphere`, and the pinned Yin–Yang donor repository identified in `provenance/`.

Important entry points: `source/lfm_cube_evolution/`, `source/lfm_cube_ionosphere/`, `source/lfm_cube_yy_mapping/`, `source/lfm_cube_yy_fac_source/`, `source/lfm_cube_yy_fixed_cap/`, `source/lfm_cube_yy_mi_convergence/`, and `source/cube_lfm/`.

```sh
python3 releases/cubelfm-igrf-v0.2.0/verify_source.py
```

This checks source hashes only and does not compile or run science.

## Reproducing on the qualified Sugon installation

`run_reference.py` reuses the accepted executable and exact input files after checking their hashes. It creates a new output directory, checks the96physical-core allocation through the rank launcher, and records the complete immediate process wait. The recorded campaign installation and its permissions must be available. Never run scientific builds, arrays or MHD on a login node.

```sh
git clone --branch CubeLFM-IGRF https://github.com/Space-Mafia/LFM-rescue.git CubeLFM-IGRF
cd CubeLFM-IGRF
git checkout cubelfm-igrf-v0.2.0
export CUBELFM_IGRF_OUTPUT=/your/scratch/new-cubelfm-igrf12h
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION releases/cubelfm-igrf-v0.2.0/run_reference.sbatch
```

For a different machine, use the pinned translation units and exact native build dependencies to establish a new build and numerical qualification. **An independent portable clean rebuild has not yet been qualified.** Repository-root CMake and older v0.1 controllers are historical development entry points, not this release's exact executable recipe. Do not claim a newly rebuilt binary has the reference word identity or4× speed without testing it.

`audit_full12h.py` and `frozen_readers_igrf_v5.py` are the actual frozen readers for newly produced data; `VALIDATED_RESULT.json` and `CASE_CONTRACT.json` preserve the closed reference evidence. Raw simulation fields, checkpoints and ELF files are not committed.

## Access and provenance

The repository remains private. Students with repository access can clone this branch/tag; additional students need collaborator or organization-team access. Existing inherited notices and redistribution rights are unchanged. See `provenance/YIN_YANG_DONOR.json` for the formal donor tag/commit and `provenance/THIRD_PARTY_NOTICES.txt`.
