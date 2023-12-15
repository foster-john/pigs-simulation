#!/usr/bin/env bash
#SBATCH --ntasks=1                # number of processes to run
#SBATCH --nodes=1                 # default number of nodes
#SBATCH --partition=cpu_compute   # good enough for what I need
#SBATCH --cpus-per-task=96        # for a multithredded job
#SBATCH --job-name=dens03        # job name
#SBATCH --output=dens03.txt   # output file 

module load R
Rscript R/workflow.R $SLURM_CPUS_PER_TASK 

