#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --partition allgroups
#SBATCH --time 24:00:00

#SBATCH --mem 8G
#SBATCH --cpus-per-task 10

#SBATCH --job-name lds-proj
#SBATCH --output results/slurm/%j-out.txt
#SBATCH --error results/slurm/%j-err.txt
#SBATCH --mail-type ALL

srun matlab -batch "$1"
