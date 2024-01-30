#!/usr/bin/env bash
#SBATCH --ntasks=1                # number of processes to run
#SBATCH --nodes=1                 # default number of nodes
#SBATCH --partition=cpu_compute   # good enough for what I need
#SBATCH --cpus-per-task=1         # for a multithredded job
#SBATCH --job-name=glm                       # job name
#SBATCH --output=outfiles/glm-%a.txt         # output file

module load R
Rscript R/analysis_glm.R $SLURM_ARRAY_TASK_ID

