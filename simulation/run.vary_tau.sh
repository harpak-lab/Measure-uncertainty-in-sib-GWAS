#!/bin/bash
#SBATCH -J vary_tau
#SBATCH -o logs/tau_%a.out
#SBATCH -e logs/tau_%a.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH --array=0-5

# Define Tau values
vals=(0 0.2 0.4 0.6 0.8 1.0)
current_tau=${vals[$SLURM_ARRAY_TASK_ID]}

echo "STARTING Job for Tau = $current_tau"

DIR_NORM="results/vary_tau/normal"
DIR_PERM="results/vary_tau/permutation"
DIR_BJK="results/vary_tau/jackknife"
mkdir -p $DIR_NORM $DIR_PERM $DIR_BJK logs

# 1. Normal
Rscript simulation/normal.R \
  --tau $current_tau \
  --beta 0.0 \
  --out_dir $DIR_NORM

# 2. Block Jackknife
Rscript simulation/jackknife.R \
  --tau $current_tau \
  --beta 0.0 \
  --out_dir $DIR_BJK

# 3. Permutation
Rscript simulation/permutation.R \
  --tau $current_tau \
  --beta 0.0 \
  --out_dir $DIR_PERM
