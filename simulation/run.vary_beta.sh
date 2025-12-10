#!/bin/bash
#SBATCH -J vary_beta_all
#SBATCH -o logs/beta_%a.out
#SBATCH -e logs/beta_%a.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G       
#SBATCH -t 08:00:00   
#SBATCH --array=0-10  

# 1. Define the parameters to vary
vals=(-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5)
current_beta=${vals[$SLURM_ARRAY_TASK_ID]}

echo "STARTING Job for Beta = $current_beta"

# 2. Define Output Directories 
DIR_NORM="results/vary_beta/normal"
DIR_PERM="results/vary_beta/permutation"
DIR_BJK="results/vary_beta/jackknife"

# Create folders if they don't exist
mkdir -p $DIR_NORM $DIR_PERM $DIR_BJK
mkdir -p logs

# ---------------------------------------------------------
# 3. RUN METHOD A: Normal Theory 
# ---------------------------------------------------------
echo "Running Normal Theory..."
Rscript simulation/normal.R \
  --beta $current_beta \
  --n_families 5000 \
  --out_dir $DIR_NORM

# ---------------------------------------------------------
# 4. RUN METHOD B: Block Jackknife
# ---------------------------------------------------------
echo "Running Block Jackknife..."
Rscript simulation/jackknife.R \
  --beta $current_beta \
  --n_families 5000 \
  --m_iter 5000 \
  --out_dir $DIR_BJK

# ---------------------------------------------------------
# 5. RUN METHOD C: Permutation 
# ---------------------------------------------------------
echo "Running Permutation..."
Rscript simulation/permutation.R \
  --beta $current_beta \
  --n_families 5000 \
  --n_sims 5000 \
  --out_dir $DIR_PERM

echo "ALL DONE for Beta = $current_beta"
