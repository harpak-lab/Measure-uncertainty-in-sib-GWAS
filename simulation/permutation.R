#!/usr/bin/env Rscript
library(MASS)
library(optparse)

# ── 1) Define Command Line Arguments 
option_list <- list(
  make_option(c("--maf"), type = "numeric", default = 0.1, 
              help = "Minor Allele Frequency [default: %default]"),
  make_option(c("--beta"), type = "numeric", default = 0.0, 
              help = "Beta 0 value [default: %default]"),
  make_option(c("--sigma1"), type = "numeric", default = 0.0, 
              help = "Sigma 1 value [default: %default]"),
  make_option(c("--tau"), type = "numeric", default = 0.0, 
              help = "Between-family effect SD (tau) [default: %default]"),
  make_option(c("--n_families"), type = "integer", default = 5000, 
              help = "Number of families per simulation [default: %default]"),
  make_option(c("--n_sims"), type = "integer", default = 5000, 
              help = "Number of simulations to run [default: %default]"),
  make_option(c("--seed_start"), type = "integer", default = 1, 
              help = "Starting seed number [default: %default]"),
  make_option(c("--out_dir"), type = "character", default = "permutation", 
              help = "Output directory [default: %default]")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Print configuration for log files
print("Running Simulation with parameters:")
print(opt)

# ── 2) Setup Variables 
maf_0       <- opt$maf
n_families  <- opt$n_families
beta0       <- opt$beta
sigma1      <- opt$sigma1
tau         <- opt$tau
n_sims      <- opt$n_sims
seed_start  <- opt$seed_start

# Create output directory if it doesn't exist
if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE)
}

# ── 3) Helper Functions
simulate_parents <- function(maf, n_families) {
  p0 <- (1 - maf)^2
  p1 <- 2 * maf * (1 - maf)
  p2 <- maf^2
  replicate(length(maf),
            sample(0:2, n_families, replace = TRUE, prob = c(p0, p1, p2)))
}

transmit_allele <- function(parent_geno) {
  ifelse(parent_geno == 1,
         rbinom(n_families, 1, 0.5),
         parent_geno / 2)
}

run_one_sign_flip_permutation <- function(delta_y, delta_X) {
  n <- length(delta_y)
  s <- sample(c(-1, 1), size = n, replace = TRUE)
  delta_y_perm <- s * delta_y
  
  # Fit model
  model <- lm(delta_y_perm ~ 0 + delta_X)
  betas <- coef(model)[1]
  names(betas) <- "SNP1"
  return(betas)
}

# ── 4) Main Simulation Loop 
df_list <- list()
seed_sequence <- seed_start:(seed_start + n_sims - 1)

# Use a counter 'i' for list indexing to avoid massive sparse lists
for (i in seq_along(seed_sequence)) {
  current_seed <- seed_sequence[i]
  set.seed(current_seed)
  
  if (i %% 10 == 0) print(paste("Processing seed:", current_seed)) # Progress tracker

  n_snps <- 1         
  maf    <- rep(maf_0, n_snps)
  
  # A) Genotypes
  father_geno <- simulate_parents(maf, n_families)
  mother_geno <- simulate_parents(maf, n_families)
  
  # B) Transmit
  sibling1_geno <- transmit_allele(father_geno) + transmit_allele(mother_geno)
  sibling2_geno <- transmit_allele(father_geno) + transmit_allele(mother_geno)
  delta_X <- sibling1_geno - sibling2_geno 
  
  # C) Phenotypes
  sigma0 <- 4
  rho_e  <- 0.06 #Within-family background correlation
  
  var_delta_eps <- sigma0 + sigma1 * (delta_X^2)
  var_delta_eps <- pmax(var_delta_eps, 1e-8)
  
  scale_i <- sqrt(var_delta_eps / (2 * (1 - rho_e)))
  
  Sigma_Z <- matrix(c(1, rho_e, rho_e, 1), nrow = 2)
  Z <- mvrnorm(n = n_families, mu = c(0, 0), Sigma = Sigma_Z)
  
  epsilon1 <- scale_i * Z[,1]
  epsilon2 <- scale_i * Z[,2]
  
  u <- rnorm(n_families, 0, tau)
  beta_i_vec <- beta0 + u
  
  y1 <- beta_i_vec * sibling1_geno + epsilon1
  y2 <- beta_i_vec * sibling2_geno + epsilon2
  delta_y <- y1 - y2        
  
  # D) Analysis (Original Model)
  model <- lm(delta_y ~ 0 + delta_X)
  beta_hat <- coef(model)[1]
  
  # E) Permutation
  B <- 1000
  perm_beta_matrix <- replicate(B, run_one_sign_flip_permutation(delta_y, delta_X))
  perm_var <- var(perm_beta_matrix)
  
  # Save result for this seed
  df_list[[i]] <- data.frame(
    seed = current_seed, # Good practice to save the seed in the data
    SNP = "SNP1",
    perm_var = perm_var,
    beta = beta_hat
  )
}

# ── 5) Save Output ────────────────────────────────────────────────
df_all <- do.call(rbind, df_list)

# Generate filename automatically based on parameters
filename <- paste0("perm_n", n_families, 
                   "_beta", beta0, 
                   "_tau", tau, 
                   "_maf", maf_0, 
                   "_sig1_", sigma1, 
                   "_start", seed_start,
                   ".txt")

full_path <- file.path(opt$out_dir, filename)
write.table(df_all, full_path, row.names = FALSE, quote = FALSE)

print(paste("Done! Results saved to:", full_path))
