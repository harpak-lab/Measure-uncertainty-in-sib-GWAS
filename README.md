# Block Jackknife Standard Error Estimation

The Block_Jackknife/ directory implements the **block jackknife (BJK) procedure** for estimating standard errors of sib-GWAS effect sizes. The pipeline repeatedly omits blocks of families, re-runs association testing, and derives jackknife-based SE estimates across all SNPs.

## Pipeline Overview

```
01_generate_phenos.py ──► 02_make_plink_jobs.sh ──► 03_beta_matrix.R ──► 04_compute_se.R

   Create 500                Generate PLINK           Collect betas          Compute final
   permuted pheno            scripts for each         across iterations      jackknife SEs
   files (omitting           chr × iteration          into one matrix        per SNP
   500 families each)                                 per chromosome
```

## Requirements

### Software
- Python >= 3.8 (`numpy`, `pandas`)
- R >= 4.0 (`data.table`, `matrixStats`, `optparse`)
- [PLINK 1.9](https://www.cog-genomics.org/plink/)

### Input Data
- **Residual phenotype file**: tab-separated, 3 columns (FID, IID, phenotype), no header
- **Sibling genotype files**: PLINK binary format (`.bed/.bim/.fam`), one per chromosome

## Step-by-Step Usage

### Step 1: Generate Jackknife Phenotype Files

Creates 500 permuted phenotype columns, each with a different block of 500 families set to NA.

```bash
python3 scripts/01_generate_phenos.py \
  --pheno_name height \
  --input_file data/raw/height.txt \
  --out_dir data/processed
```

**Output**: `height.blockjackknife.phenos.txt` — one file with columns `FID, IID, pheno1, ..., pheno500`

### Step 2: Generate PLINK Job Scripts

Creates 22 shell scripts (one per chromosome), each containing 500 PLINK `--qfam` commands.

```bash
bash scripts/02_make_plink_jobs.sh \
  height \
  data/genotypes \
  data/processed/height.blockjackknife.phenos.txt \
  output/plink_results
```

**Output**: `run.bjk_se.height.chr1.sh` ... `run.bjk_se.height.chr22.sh`

> **Note**: Submit these scripts to your cluster (e.g., via SLURM or GNU parallel). Each script runs 500 PLINK calls, so this is the most compute-intensive step.

### Step 3: Construct Beta Matrices

Collects PLINK `.qfam.within` output files and assembles a SNP × iteration beta matrix for each chromosome.

```bash
Rscript scripts/03_beta_matrix.R \
  output/plink_results \
  output/beta_matrices
```

**Output**: `beta_matrix.chr1.txt` ... `beta_matrix.chr22.txt`

### Step 4: Compute Jackknife Standard Errors

Applies the block jackknife variance formula to the beta matrices.

```bash
Rscript scripts/04_compute_se.R \
  --trait height \
  --n_families 17178 \
  --input_dir output/beta_matrices \
  --out_dir results
```

Set `--n_families` to the total number of sibling families available for your trait.

**Output**:
- Per-chromosome files: `height_chr1.bjk_se.txt` ... `height_chr22.bjk_se.txt`
- Combined file: `height.bjk_se.all_snps.csv`

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `m` | 500 | Number of jackknife replicates |
| `d` | 500 | Families omitted per replicate |
| `n_families` | trait-specific | Total number of sibling families (set manually) |

The jackknife variance for each SNP is calculated as:

```
Var(β) = (r / (d × (m − 1))) × Σ(β_i − β̄)²
```

where `r = n_families − d`.
