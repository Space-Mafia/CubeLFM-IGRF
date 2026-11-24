#!/bin/bash
#SBATCH --job-name=cubelfm128-global
#SBATCH --nodes=2
#SBATCH --ntasks=192
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --sockets-per-node=8
#SBATCH --cores-per-socket=12
#SBATCH --ntasks-per-socket=12
#SBATCH --mem=128G
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:cyclic
#SBATCH --time=72:00:00
#SBATCH --no-requeue
set -euo pipefail
repo=${1:?Supply the absolute path to the public checkout}
reference=${2:?Supply the completed two-cap Sugon event directory}
module --force purge
module load mpi/openmpi/openmpi-4.1.5-gcc9.3.0
module swap compiler/gcc/9.3.0 compiler/gcc/12.3.0
export LANG=C LC_ALL=C OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
unset PYTHONOPTIMIZE LFM_TEST_OBLIQUE_B0 LFM_SIDE_EMF_MODE CUBE_DIAGNOSTICS OMPI_MCA_btl
unset OMPI_MCA_mtl_ofi_provider_exclude OMPI_MCA_mtl_ofi_provider_exclude_list
export LD_LIBRARY_PATH=/public/software/mpi/intelmpi/2021.3.0/libfabric/lib:${LD_LIBRARY_PATH:-}
export FI_PROVIDER_PATH=/public/software/mpi/intelmpi/2021.3.0/libfabric/lib/prov
export FI_PROVIDER=psm2 HFI_NO_CPUAFFINITY=1
export OMPI_MCA_pml=cm OMPI_MCA_mtl=ofi OMPI_MCA_mtl_ofi_provider_include=psm2 OMPI_MCA_oob_tcp_if_include=ib0
ulimit -s unlimited
exec /public/software/apps/anaconda3/2023.09/bin/python3 \
  "$repo/event128/ports/yinyang_igrf/events/20021023/event128_comparison/launch.py" \
  "$reference" --source "$repo/event128/reference-source" --branch-local-wall
