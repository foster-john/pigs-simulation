#!/usr/bin/env bash
#SBATCH --ntasks=1                # number of processes to run
#SBATCH --nodes=1                 # default number of nodes
#SBATCH --partition=cpu_compute   # good enough for what I need
#SBATCH --cpus-per-task=12        # for a multithredded job
#SBATCH --job-name=test_sim       # job name
#SBATCH --output=test_sim%J.txt   # output file 

module load R
Rscript R/workflow.R $SLURM_CPUS_PER_TASK 

