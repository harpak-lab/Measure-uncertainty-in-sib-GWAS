library(data.table)
library(matrixStats)
library(optparse)

option_list <- list(
  make_option(c("-t", "--trait"), type="character", help="Trait name"),
  make_option(c("-n", "--n_families"), type="integer", help="Number of families"),
  make_option(c("-i", "--input_dir"), type="character", help="Directory containing beta matrices"),
  make_option(c("-o", "--out_dir"), type="character", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

m <- 500 # The number of resampling replicates
d <- 500 # The Block Size: the number of families deleted in one Block Jackknife resampling replicate
r <- opt$n_families - d # The number of remaining families after deleting d of them 

dir.create(opt$out_dir, showWarnings = FALSE)
bjk_list <- list()

for(chr in 1:22){
  fname <- file.path(opt$input_dir, paste0("beta_matrix.chr", chr, ".txt"))
  
  if(!file.exists(fname)) {
    warning(paste("Missing file:", fname))
    next
  }

  beta_matrix <- fread(fname) # fread is faster than read.table
  
  # Convert to matrix (excluding SNP col)
  beta_vals <- as.matrix(beta_matrix[, -1, with=FALSE])
  
  beta_means <- rowMeans(beta_vals, na.rm = TRUE)
  ssq <- (beta_vals - beta_means)^2
  ssq_sum <- rowSums(ssq, na.rm = TRUE)

  # Jackknife Variance Formula
  jack_var <- (r / (d * (m - 1))) * ssq_sum
  jack_se  <- sqrt(jack_var)

  bjk_df <- data.frame(SNP = beta_matrix$SNP, bjk_se = jack_se)
  
  # Save individual chromosome
  write.table(bjk_df, file.path(opt$out_dir, paste0(opt$trait, "_chr", chr, ".bjk_se.txt")), 
              row.names=FALSE, quote=FALSE)
  
  bjk_list[[chr]] <- bjk_df
  message(paste("Finished CHR", chr))
}

# Bind and save all
bjk_all <- do.call(rbind, bjk_list)
write.csv(bjk_all, file.path(opt$out_dir, paste0(opt$trait, ".bjk_se.all_snps.csv")), 
          row.names=FALSE, quote=FALSE)
