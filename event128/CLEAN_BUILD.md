# Build and run the 128-grid event from a fresh checkout

This workflow starts the 23 October 2002 event at 04:00 UT. It builds all 114 C and Fortran translation units in a new directory, including the 53 members selected from nine native archives in the reference link. It uses the original 128 × 48 × 48 × 6 mesh, native seam, corrected CT induction reference and global PARDISO ionospheric solver. All project source, project headers, dated IGRF tables and the two fixed geometry files are in the checkout.

The clean Sugon build compiled all 114 sources and linked the model. The 192-rank short test was cancelled while queued. A runtime result for this clean build is still pending. The archived 12-hour performance figure belongs to the original event run.

## Requirements

| Dependency | What the build or run uses |
|---|---|
| Linux x86-64 with AVX2 | The event includes native AVX2 kernels. Use the same CPU architecture on the build and run nodes. |
| GNU GCC and GFortran, GNU Make, GNU ar | The reference toolchain uses GCC 12.3. Per-file flags preserve the Fortran default real sizes, COMMON layout and floating-point options of the event build. |
| Open MPI with C and Fortran development support | Supply `mpicc`, `mpifort` (including `mpi_f08`) and `mpirun` from the same installation. On Sugon the supplied batch script selects the working Open MPI / OFI / PSM2 configuration. |
| Intel oneMKL development installation | Set `MKLROOT` to the prefix containing `include/mkl.h`, `include/mkl_dfti.h` and the shared `libmkl_rt.so`. The global potential solver uses PARDISO and the wall return uses DFTI. |
| HDF5 C development installation | Set `HDF5_ROOT` to the prefix containing `include/hdf5.h` and `lib/libhdf5.so` or `lib64/libhdf5.so`. |
| Python 3.9+ with NumPy and h5py | Run the input converter, the compact HDF5 writer and the endpoint test with the same Python environment. The compiler, MPI, MKL and native HDF5 libraries are installed separately. |

On a cluster, load the site modules or ask the administrator for their installation prefixes. The Sugon batch script supplies the paths there. On another system, set `MPICC` and `MPIFC` if your GNU-backed MPI wrappers use different names. The two-node run needs 192 MPI ranks. Each rank handles 16 × 24 × 24 cells on one panel.

On a new Python installation, install the two Python packages with `python3 -m pip install numpy h5py`. On Sugon the batch script uses the site Python that already provides them.

## Event inputs

The two text files in `event128/inputs/geometry/` are the exact inputs from the reference run. `shell_geometry.txt` defines the mapped mesh and `ionosphere_reference.txt` supplies the native ionospheric reference. Their SHA-256 values are listed in [CASE.json](CASE.json). The degree-five IGRF coefficients, start phase, spin and coordinate lookup tables are already dated for the event.

`prepare_inputs.py` gets the one-minute OMNI solar wind from NASA CDAWeb. It selects IMF, velocity, proton density and temperature, uses 160 minutes of padding on either side of the 12-hour interval, and transforms GSM IMF and GSE velocity separately into the event frame. The model driver spans 01:20–18:40 UT. The request includes additional samples for interpolation. The converter writes the HDF5 data, source and processing metadata. The archived HDF5 file hash in CASE.json refers to the historical file; a regenerated file can carry different path and retrieval metadata. [NASA OMNI dataset](https://doi.org/10.48322/mj0k-fq60)

For a compute node without Internet access, download the CSV on an Internet-connected machine and copy it into your own scratch directory. This command retrieves the columns in the order expected by the converter:

```sh
curl --fail --location --retry 2 --output nasa-omni.csv \
  'https://cdaweb.gsfc.nasa.gov/hapi/data?id=OMNI_HRO2_1MIN&parameters=BX_GSE%2CBY_GSM%2CBZ_GSM%2CVx%2CVy%2CVz%2Cproton_density%2CT&time.min=2002-10-23T01%3A09%3A00Z&time.max=2002-10-23T18%3A52%3A00Z&format=csv'
```

The reference NASA CSV has SHA-256 `8dade9d352636cd1e694a6fa3c6f9732f858c1a745feec7a851a823df4ee542c`. Check a download with `sha256sum nasa-omni.csv`. A changed digest can indicate a revised OMNI record.

## Submit the short test on Sugon

Start from the root of the public checkout. Replace the scratch path, Slurm account and partition with your own. The checkout and CSV must be readable on both nodes.

```sh
repo="$(pwd -P)"
raw=/absolute/path/to/nasa-omni.csv
work=/your/scratch/cubelfm-128-smoke
mkdir "$work" && cd "$work" || exit 1
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  "$repo/event128/smoke.sbatch" "$repo" "$raw"
```

The job converts the driver, builds and links the model from source, and cold-starts on two nodes with 96 ranks per node and one thread per rank. It advances 120 simulated seconds to **04:02 UT**. The normal 120 s frame interval reaches its first output and the model writes its endpoint checkpoint. A successful test writes `endpoint-result.json` with the actual step count, final time step, MI clock, global residual and full MHD process wall time.

## Build or run the stages separately

Run these commands on allocated compute nodes from a fresh work directory. Set `repo`, `HDF5_ROOT` and `MKLROOT` first. For offline input preparation append `--raw-csv /absolute/path/to/nasa-omni.csv` to the first command.

```sh
python3 "$repo/event128/prepare_inputs.py" inputs
python3 "$repo/event128/build.py" build --jobs 12
python3 "$repo/event128/run_clean.py" build inputs case --seconds 120
python3 "$repo/event128/check_smoke.py" case
```

For a full 12-hour cold start on Sugon, substitute `event128/run.sbatch` for `event128/smoke.sbatch` in the submission command. Both scripts take the checkout and optional NASA CSV as arguments. For a manual run of this interval, request enough allocation time and use `--seconds 43200` with a new output directory. `CUBELFM_MPI_LAUNCH` can supply your site's MPI launch options. The runner appends `-np 192` and the model arguments. The smoke-test batch script sets the qualified Sugon MPI options.

The source inventory and per-unit compile profiles are in [native/units.json](native/units.json). `build.py` creates the Fortran modules and static libraries inside the new build tree. The HDF5 files in `case/` include both the 192 MHD tiles and the global ionosphere fields. The source-only historical check remains available as `python3 event128/verify.py`.
