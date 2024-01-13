#!/usr/bin/env bash
#SBATCH --ntasks=1                # number of processes to run
#SBATCH --nodes=1                 # default number of nodes
#SBATCH --partition=cpu_compute   # good enough for what I need
#SBATCH --cpus-per-task=3         # for a multithredded job
# #SBATCH --mem=20g

# #SBATCH --job-name=collate                       # job name
# #SBATCH --output=outfiles/collate_%J.txt         # output file

#SBATCH --job-name=trapSnare
#SBATCH --output=outfiles/trapSnare-%a.txt   # output file
# #SBATCH --array=1-40

# #SBATCH --dependency=afterany:13589

module load R
Rscript R/workflow.R $SLURM_ARRAY_TASK_ID 
# Rscript R/collate_mcmc_output.R
