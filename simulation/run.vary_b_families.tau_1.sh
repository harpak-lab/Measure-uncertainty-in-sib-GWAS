#!/bin/bash
#SBATCH -J vary_n
#SBATCH -o logs/n_%a.out
#SBATCH -e logs/n_%a.err
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G      
#SBATCH -t 10:00:00    
#SBATCH --array=0-5

# Define sample sizes to test
vals=(1000 2000 3000 4000 5000 10000)
current_n=${vals[$SLURM_ARRAY_TASK_ID]}

echo "STARTING Job for N_families = $current_n"

DIR_NORM="results/vary_n/normal"
DIR_PERM="results/vary_n/permutation"
DIR_BJK="results/vary_n/jackknife"
mkdir -p $DIR_NORM $DIR_PERM $DIR_BJK logs

# 1. Normal
Rscript simulation/normal_simulation.R \
  --n_families $current_n \
  --tau 1.0 \
  --out_dir $DIR_NORM

# 2. Block Jackknife
Rscript simulation/jackknife_simulation.R \
  --n_families $current_n \
  --tau 1.0 \
  --out_dir $DIR_BJK

# 3. Permutation
Rscript simulation/permutation_simulation.R \
  --n_families $current_n \
  --tau 1.0 \
  --out_dir $DIR_PERM
