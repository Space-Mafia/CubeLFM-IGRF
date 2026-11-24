# Reproduce the 23 October 2002 Earth event

## Case

The run starts on the fine mesh at 04:00 UT, advances 43,200 s with the dated OMNI solar wind, and ends at 16:00 UT. The initial field uses the degree-five IGRF expansion for the date and its intrinsic tilt. The event rotation rate is 0.000072921150 rad/s. The initial phase and frame conventions are recorded in the event headers and rotation include. The 12 hours in this case describe the simulation interval. [Case and input hashes](CASE.json) · [Rotation constants](reference-source/yy_rotation_profile.inc)

| Quantity | Executed setting |
|---|---|
| MHD mesh | Native mapped LFM-Cube, 128 × 48 × 48 per panel, six panels |
| MPI layout | 8 × 2 × 2 ranks per panel, 16 × 24 × 24 cells per rank, 96 physical cores per node on two nodes, one thread per rank |
| Time stepping | CFL 0.3, Boris speed 1,100 km/s, BDIFFUSE total jump coefficient 0.5 |
| Ionosphere | 192 geographic longitudes × 97 latitudes, global Hall–Pedersen tensor solve with MKL PARDISO |
| Conductance | Event EUV and electron precipitation with diffuse precipitation boundary feedback |
| Output clocks | Natural 10 s MI updates with held potential, frames every 120 s, checkpoints every 3,600 s |

### What changed in the model

The physical initial field retains its radial taper. The induction update constructs a separate, untapered planetary reference on the CT edges. It uses the native curl and the existing shared-face ownership. This change addresses the shell-aligned rings previously seen in solar-wind Bz. The MHD update uses the original cross-panel seam. [Compiled evolution routine](compiled/lfm_cube_evolution.f90) · [Induction source selection](ports/yinyang_igrf/events/20021023/fixes/induction_reference)

The ionosphere uses one geographic potential on the sphere. The solver removes the area-weighted mean FAC for current closure and fixes the potential's area-weighted mean. Both cap views sample that solution. At the MHD wall, adjacent traces can land on different footpoint branches. The wall operator estimates the local voltage increment across such edges from neighboring same-branch differences and reconstructs one scalar voltage for the gas and CT wall updates. The ionospheric solve in this release imposes no north–south conjugate voltage constraint. [Global tensor solve](compiled/global_tensor.c) · [Coupling and output](compiled/cube_event_global.f90) · [Wall reconstruction](compiled/branch_wall.c)

## Offline source test

From the repository root, run `python3 event128/verify.py` with Python 3.9 or newer. The test materializes the source from the dated baseline in a temporary directory, selects the induction correction and branch-local wall, and prepares the main driver and output writers. It checks the resulting files against the executed build. A successful run prints the counts of 14 compiler inputs, 47 inherited source mappings and two writers.

`compiled/` is the byte-exact reading copy of the 14 compiled units. The batch build generates these units from `reference-source/`, `ports/` and `cubesphere/io/`. For development, edit the generating source. The historical comparison in `verify.py` will identify the units affected by an edit.

## Historical site launcher in v0.3.0

For a new installation use the [clean source build and 192-core smoke test](CLEAN_BUILD.md). The commands below describe how the original packaging used a completed site baseline.

Use the completed native two-cap **128-grid baseline** named `event20021023_native128_192_aurora12h_20261120_v2` from the same site installation. Its directory supplies the successful MPI launch layout, shared object/library paths and three inputs. The event starts from the cold initial condition. The launcher compiles the event sources on the allocation, links the 128-grid executable and starts it with the original MPI and input arguments.

The installation must provide:

- Read access to the baseline's `CONFIG.json`, `event12h.receipt.json`, `PHYSICAL_FULL12H_CLOSED.json`, the linked native objects, modules and archives, and the configuration and reference objects named there.
- `shell_geometry.txt`, `ionosphere_reference.txt` and `omni_extended.h5` at the baseline's resolved input paths. Their reference SHA-256 values are in [CASE.json](CASE.json).
- The recorded GCC/GFortran 12.3 and Open MPI environment, HDF5, Intel MKL with PARDISO and DFTI, the OFI/PSM2 transport and Python with NumPy and h5py. The [compiler and link arguments](COMPILE_LINK.json) show the ABI and library order. `@RUN@` denotes the historical build directory and `@SITE_HOME@` the parent of the site's campaign tree.

Run the following from the repository root on Sugon. Set the capitalized placeholders to paths and scheduler values assigned to your account. Use a fresh working directory for each submission.

```sh
repo="$(pwd -P)"
reference=/absolute/path/to/event20021023_native128_192_aurora12h_20261120_v2
run_dir=/your/scratch/cubelfm-event-20021023-v030
mkdir "$run_dir" && cd "$run_dir" || exit 1
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  "$repo/event128/run.sbatch" "$repo" "$reference"
```

The run directory receives the build, process logs, receipts and compact HDF5 outputs under `event12h/`. Each frame's `global_mi` group stores `potential_V`, `outward_fac_A_m2`, `pedersen_S` and `hall_S` on the geographic grid. Use the held MI clock when matching an MHD frame to its ionospheric state. The duplicated residual in the cap-oriented diagnostics refers to one global solve. [Output schema implementation](runtime/archive_global.py)

## Reference result and scope

Sugon job **29545636** exited successfully after 495,712 steps. Its last frame-step record reached 43,200 s and its last MI update index was 4,321. The full MHD process took 87,718.14473083615 s, including initialization, output and exit, giving **0.4924864762× real time**. Compilation time is separate. The SHA-256 of the executed binary is `606770c02d172da36f2f34905aa7d2be926932e0cbf5d250cb5a7f5638181678`. This measurement belongs to the archived Sugon run. The public packaging was tested by regenerating its inputs. The 12-hour measurement applies to the historical build. The scope of the separate clean-build test is described in [the installation guide](CLEAN_BUILD.md).

[BUILD.json](BUILD.json) identifies all 61 explicitly linked object sources and nine native archives. The clean build supplies the 53 archive member sources selected at link time and the two fixed geometry inputs. Intel MKL, MPI and HDF5 are external dependencies. The solar wind is regenerated from NASA OMNI. The repository includes the event-specific IGRF, orientation and geomagnetic lookup headers. The full-sphere solver sources derive from the Yin–Yang GAMERA global tensor implementation; the event uses the sparse direct path and its residual check. [Solver sources](ports/yinyang_igrf/events/20021023/global_mi/yyg) · [Inherited notices](../provenance/THIRD_PARTY_NOTICES.txt)
