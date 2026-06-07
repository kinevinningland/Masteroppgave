#!/bin/bash
#SBATCH --partition=CPUQ
#SBATCH --account=share-ie-iel
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=log_base.txt

module load Julia/1.10.4-linux-x86_64

cd ~/Masteroppgave/resddp

julia --threads=$SLURM_CPUS_PER_TASK examples/4area.jl
