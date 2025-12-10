#!/bin/bash
#SBATCH -J vary_maf
#SBATCH -o logs/maf_%a.out
#SBATCH -e logs/maf_%a.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH --array=0-5

# Define MAF values (0.05 to 0.5)
vals=(0.05 0.1 0.2 0.3 0.4 0.5)
current_maf=${vals[$SLURM_ARRAY_TASK_ID]}

echo "STARTING Job for MAF = $current_maf"

DIR_NORM="results/vary_maf/normal"
DIR_PERM="results/vary_maf/permutation"
DIR_BJK="results/vary_maf/jackknife"
mkdir -p $DIR_NORM $DIR_PERM $DIR_BJK logs

# 1. Normal
Rscript simulation/normal.R \
  --maf $current_maf \
  --beta 1.0 \
  --out_dir $DIR_NORM

# 2. Block Jackknife
Rscript simulation/jackknife.R \
  --maf $current_maf \
  --beta 1.0 \
  --out_dir $DIR_BJK

# 3. Permutation
Rscript simulation/permutation.R \
  --maf $current_maf \
  --beta 1.0 \
  --out_dir $DIR_PERM
