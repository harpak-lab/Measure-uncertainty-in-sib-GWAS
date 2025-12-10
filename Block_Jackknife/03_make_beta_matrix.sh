args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript make_beta_matrix.R <TRAIT_DIR> <OUTPUT_DIR>")

trait_dir <- args[1]
output_dir <- args[2]

library(data.table)

dir.create(output_dir, showWarnings = FALSE)

# Iterate through chromosomes
for(chr in 1:22) {
  beta_mat <- NULL
  
  # Note: Ideally, check if files exist before loop
  for(i in 1:500) {
    fn <- file.path(trait_dir, paste0("blockjackknife.se.chr", chr, ".iter", i, ".qfam.within"))
    
    if(!file.exists(fn)) next # Skip if missing
    
    dt <- fread(fn, select = c("SNP", "BETA"))
    setkey(dt, SNP)

    if (is.null(beta_mat)) {
      beta_mat <- data.table(SNP = dt$SNP)
    }
    
    # Dynamic column assignment
    beta_mat[[paste0("iter", i)]] <- dt$BETA
  }
  
  if(!is.null(beta_mat)) {
    outfn <- file.path(output_dir, paste0("beta_matrix.chr", chr, ".txt"))
    fwrite(beta_mat, outfn, sep = " ")
    message("Wrote: ", outfn)
  }
}
