#!/bin/bash
#SBATCH --job-name=refitAngVV
#SBATCH --account=ac_brooklab2
#SBATCH --partition=savio3_htc
#SBATCH --array=1-240%200
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=36:00:00
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

IC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ic_list.txt)

echo "Mapped array index ${SLURM_ARRAY_TASK_ID} → IC ${IC}"

Rscript fit-models-AngV.R "$IC"


end_time=$(date +%s)
runtime=$((end_time - start_time))


# Completion message with timing information
echo "Job completed at: $(date)"
echo "Total runtime: $runtime seconds"