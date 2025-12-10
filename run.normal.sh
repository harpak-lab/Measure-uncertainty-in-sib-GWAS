#!/bin/bash
#SBATCH -A OTH24018
#SBATCH -J test.normal
#SBATCH -o test.normal.o%j
#SBATCH -e test.normal.o%j
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH --array=0-10

parameter_list=(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5)
parameter=${parameter_list[$SLURM_ARRAY_TASK_ID]}

echo "Running with parameter = $parameter"

Rscript normal.R $parameter