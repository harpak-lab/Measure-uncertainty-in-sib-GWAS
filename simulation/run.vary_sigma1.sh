#!/bin/bash
#SBATCH -J vary_sigma1
#SBATCH -o logs/sig1_%a.out
#SBATCH -e logs/sig1_%a.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH --array=0-6

# Define Sigma1 values
vals=(-0.3 -0.2 -0.1 0.0 0.1 0.2 0.3)
current_sig1=${vals[$SLURM_ARRAY_TASK_ID]}

echo "STARTING Job for Sigma1 = $current_sig1"

DIR_NORM="results/vary_sigma1/normal"
DIR_PERM="results/vary_sigma1/permutation"
DIR_BJK="results/vary_sigma1/jackknife"
mkdir -p $DIR_NORM $DIR_PERM $DIR_BJK logs

# 1. Normal
Rscript simulation/normal.R \
  --sigma1 $current_sig1 \
  --beta 0.0 \
  --out_dir $DIR_NORM

# 2. Block Jackknife
Rscript simulation/jackknife.R \
  --sigma1 $current_sig1 \
  --beta 0.0 \
  --out_dir $DIR_BJK

# 3. Permutation
Rscript simulation/permutation.R \
  --sigma1 $current_sig1 \
  --beta 0.0 \
  --out_dir $DIR_PERM
