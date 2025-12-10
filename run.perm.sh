#!/bin/bash
#SBATCH -A OTH24018
#SBATCH -J test.perm
#SBATCH -o test.perm.o%j
#SBATCH -e test.perm.o%j
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 6:00:00
#SBATCH --array=0-5

parameter_list=(0.05 0.1 0.2 0.3 0.4 0.5)
parameter=${parameter_list[$SLURM_ARRAY_TASK_ID]}

echo "Running with parameter = $parameter"

Rscript adaptive_permutation.R $parameter
