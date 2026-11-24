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
repo=${1:?Pass the absolute path of the public checkout}
raw_csv=${2:-}
exec bash "$repo/event128/smoke.sbatch" "$repo" "$raw_csv" 43200
