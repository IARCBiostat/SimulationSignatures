#!/bin/bash
#SBATCH --job-name=sim_scen2
#SBATCH --array=1-2%100          # cap 15 concurrent jobs
#SBATCH --cpus-per-task=2        # admins halve → ~2 effective workers
#SBATCH --partition=low_p
#SBATCH --output=Slurm_outputs/scenario_2-%A_%a.out
#SBATCH --error=Slurm_outputs/scenario_2-%A_%a.err

# Prevent hidden multithreading
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

/opt/R/4.3.1/bin/Rscript R/run_scenario_2.R