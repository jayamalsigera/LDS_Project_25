#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --partition allgroups
#SBATCH --time 04:00:00

#SBATCH --mem 32G
#SBATCH --cpus-per-task 15

#SBATCH --job-name lds-nominal
#SBATCH --output results/slurm/%j-out.txt
#SBATCH --error results/slurm/%j-err.txt
#SBATCH --mail-type ALL

srun matlab -batch "estimateSST3d"
