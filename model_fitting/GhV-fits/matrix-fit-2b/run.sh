#!/bin/bash
#SBATCH --job-name=matrixfit2b
#SBATCH --account=ac_brooklab2
#SBATCH --partition=savio3_htc
#SBATCH --array=1001-2000%200
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=48:00:00
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err


# Initialize Lmod module system
module load r/4.4.0-gcc-11.4.0

# Print this sub-job's task ID
echo "Job ID: $SLURM_JOB_ID"
echo "Starting R script with IC list number: $SLURM_ARRAY_TASK_ID"
echo "Running on host: $(hostname)"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

start_time=$(date +%s)

Rscript fit-models.R "$SLURM_ARRAY_TASK_ID"

end_time=$(date +%s)
runtime=$((end_time - start_time))


# Completion message with timing information
echo "Job completed at: $(date)"
echo "Total runtime: $runtime seconds"