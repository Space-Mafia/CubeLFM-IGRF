# CubeLFM-IGRF

Public **CubeLFM64 Earth IGRF + tilt + rotating MI** research source snapshot, version **v0.2.0**. This standalone repository contains the accepted scientific version; it has its own Git history.

## Get the code

```sh
git clone https://github.com/Space-Mafia/CubeLFM-IGRF.git
cd CubeLFM-IGRF
git checkout cubelfm-igrf-v0.2.0
python3 verify_source.py
```

The repository is public: students can view and clone it without an invitation. `verify_source.py` checks the pinned source hashes and does not build or run a simulation.

## Accepted case

| Control | Value |
|---|---|
| MHD mesh | Original mapped LFM-Cube, 64 × 24 × 24 × 6 |
| Parallelism | 96 MPI ranks / 96 distinct physical cores / one node; OMP=1 |
| Main field | Pinned Yin–Yang release IGRF coefficients, active degree5; intrinsic tilt retained |
| Rotation | 43200 s: one turn in 12 h, phase0 |
| Solar wind Bz | +5 nT for 2 h, then −5 nT for 10 h |
| MI solver | Native LFM ITMAX=100 / ERRMAX=1e-3; finite nonconvergence warns and records converged=false |
| Conductance | Constant SigmaP=5 S, Hall=0 |
| Coupling | Natural 10 s potential solves, held potential; no MI timestep clipping |
| MHD controls | BDIFFUSE total jump factor0.5; PDM0.75; Boris1100 km/s; CFL0.3; density floor0 |
| Output / checkpoint | 120 s / 3600 s |

The FAC correction subtracts the complete current-time planetary CT reference before donor CT recovery, curl and projection. Field, tracing and mapping code comes from the pinned Yin–Yang implementation; the native LFM MI solver is retained. All invalid/nonfinite input, mapping and communication errors remain fatal.

## Validation and speed

The accepted run advanced **43200 physical seconds** in246975 steps. Its complete immediate MHD process wait, including initialization, output and exit, was10789.844918493181 s: **4.003764681173281× real time on96 cores**. The frozen numerical audit and complete MHD+MI word comparison passed:361 archive epochs,34656 rank field files,4321 natural MI updates and12 checkpoints. Different allocations preclude an isolated paired speedup claim.

The **6–7× speed optimization is ongoing and is not part of this accepted snapshot**. This release records the accepted science; quantitative full-MHD Yin–Yang benchmarking was not completed.

## Source and build status

`source/` contains all49 explicitly selected C/Fortran translation units, copied byte-for-byte from the executed version. `headers/` contains12 distinct included headers. `SOURCE_MANIFEST.json` records exact source/object identities and compilation provenance. `CASE_CONTRACT.json`, `VALIDATED_RESULT.json` and the frozen readers preserve the accepted case and its evidence. Donor identities and notices are in `provenance/`.

Key entry points are `source/lfm_cube_evolution/`, `source/lfm_cube_ionosphere/`, `source/lfm_cube_yy_mapping/`, `source/lfm_cube_yy_fac_source/`, `source/lfm_cube_yy_fixed_cap/`, `source/lfm_cube_yy_mi_convergence/` and `source/cube_lfm/`.

**An independent portable clean rebuild is not yet qualified.** This snapshot contains the selected translation units and headers; it does not bundle the inherited native build tree, immutable libraries, qualified ELF, input data or raw simulation output. `BUILD_DEPENDENCIES.json` identifies the original native-library and module dependencies. The original campaign paths in provenance are historical identities.

On the qualified Sugon installation, `run_reference.py` reuses the accepted ELF and exact input files after checking their hashes. The campaign installation and its permissions must already be available. Submit from an allocation using your own account and partition:

```sh
export CUBELFM_IGRF_OUTPUT=/your/scratch/new-cubelfm-igrf12h
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION run_reference.sbatch
```

All scientific builds, arrays, MPI/MHD and audits belong on allocated compute. A newly built binary requires its own numerical qualification; the recorded4× rate is specific to the accepted executable and case.

## Provenance and notices

This publication exports only the accepted `cubelfm-igrf-v0.2.0` snapshot from source commit14df214a72512a25866661d476b671097b44b286. See `PUBLICATION.json`, `provenance/YIN_YANG_DONOR.json` and `provenance/THIRD_PARTY_NOTICES.txt`. Inherited copyright and redistribution notices remain in force.
