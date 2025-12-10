#!/bin/bash
TRAIT=$1
GENO_DIR=$2   # e.g., sibgeneticfiles
PHENO_FILE=$3 # e.g., pheno_files/trait.blockjackknife.phenos.txt
OUT_DIR=$4

if [ -z "$TRAIT" ]; then
  echo "Usage: ./02_make_plink_jobs.sh <TRAIT> <GENO_DIR> <PHENO_FILE> <OUT_DIR>"
  exit 1
fi

mkdir -p "$OUT_DIR"

for chr in {1..22}; do
  script_name="run.bjk_se.${TRAIT}.chr${chr}.sh"
  
  # Create/Overwrite file
  echo "#!/bin/bash" > "$script_name"
  
  for iter in $(seq 1 500); do
    echo "plink1.9 \\
      --bfile ${GENO_DIR}/sibs.relabeled.chr${chr} \\
      --pheno ${PHENO_FILE} \\
      --pheno-name pheno${iter} \\
      --qfam mperm=1 --memory 32000 --silent --threads 64 \\
      --out ${OUT_DIR}/blockjackknife.se.chr${chr}.iter${iter}" >> "$script_name"
  done
  
  chmod +x "$script_name"
  echo "Generated $script_name"
done
