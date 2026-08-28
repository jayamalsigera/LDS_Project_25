#!/bin/bash
#SBATCH --array=0-4
#SBATCH --ntasks=1
#SBATCH --partition=allgroups
#SBATCH --time=00:30:00

#SBATCH --mem=6500M
#SBATCH --cpus-per-task=2

#SBATCH --job-name=lds-smoke
#SBATCH --output=results/slurm/%A_%a-out.txt
#SBATCH --error=results/slurm/%A_%a-err.txt
#SBATCH --mail-type=ALL

# Define your MATLAB commands in a bash array
# Note: Wrapping them in single quotes preserves the internal double quotes
COMMANDS=(
    'estimateSST3d("T=500","totalRuns=2")'
    'estimateSST3dLfm("T=500","totalRuns=2")'
    'estimateSST3dLfm("T=500","totalRuns=2","klTolerance=0")'
    'estimateSST3dLfm("T=500","totalRuns=2","klTolerance=0.005")'
    'estimateSST3dLfm("T=500","totalRuns=2","klTolerance=0.005","lfmKlTolerance=0.01")'
)

# Grab the command corresponding to this specific Slurm task ID
CURRENT_COMMAND="${COMMANDS[$SLURM_ARRAY_TASK_ID]}"

# Execute it
srun matlab -batch "$CURRENT_COMMAND"
