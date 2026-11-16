#!/bin/bash
#SBATCH --job-name=CubeLFM-IGRF64-12h
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --hint=nomultithread
#SBATCH --time=06:00:00
set -euo pipefail
# Add your site's account/partition options when submitting this template.
module --force purge
module load mpi/openmpi/openmpi-4.1.5-gcc9.3.0
module swap compiler/gcc/9.3.0 compiler/gcc/12.3.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
: "${CUBELFM_IGRF_OUTPUT:?Set a NEW output-directory path before sbatch}"
python3 releases/cubelfm-igrf-v0.2.0/run_reference.py --output "$CUBELFM_IGRF_OUTPUT"
