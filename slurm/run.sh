#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --partition allgroups
#SBATCH --time 24:00:00

#SBATCH --mem 48G
#SBATCH --cpus-per-task 20

#SBATCH --job-name lds-proj
#SBATCH --output results/slurm/%j-out.txt
#SBATCH --error results/slurm/%j-err.txt
#SBATCH --mail-type ALL

# Use --mem=96G --cpus-per-task 32 for running the full 200 MC runs

srun matlab -batch "$1"
