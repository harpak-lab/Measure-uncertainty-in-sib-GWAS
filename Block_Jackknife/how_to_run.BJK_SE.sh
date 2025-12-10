# 1. Generate Phenotypes
python3 scripts/01_generate_phenos.py \
  --pheno_name height \
  --input_file data/raw/height.txt \
  --out_dir data/processed

# 2. Generate and Run PLINK Jobs
bash scripts/02_make_plink_jobs.sh height data/genotypes data/processed/height.phenos.txt output/plink_results
# (Then run the generated .sh files using parallel or SLURM)

# 3. Construct Beta Matrix
Rscript scripts/03_beta_matrix.R output/plink_results output/beta_matrices

# 4. Compute Final SEs
Rscript scripts/04_compute_se.R \
  --trait height \
  --n_families 17178 \ #manually set according to the available sibling information
  --input_dir output/beta_matrices \
  --out_dir results
