#!/usr/bin/env Rscript
library(MASS)
library(optparse)

# ── 1) Define Command Line Arguments
option_list <- list(
  make_option(c("--maf"), type = "numeric", default = 0.1, 
              help = "Minor Allele Frequency [default: %default]"),
  make_option(c("--beta"), type = "numeric", default = 0.0, 
              help = "Beta 0 value (Mean effect size) [default: %default]"),
  make_option(c("--sigma1"), type = "numeric", default = 0.0, 
              help = "Sigma 1 value (Heteroskedasticity) [default: %default]"),
  make_option(c("--tau"), type = "numeric", default = 0.5, 
              help = "Between-family effect SD (tau) [default: %default]"),
  make_option(c("--n_families"), type = "integer", default = 5000, 
              help = "Number of families per simulation [default: %default]"),
  make_option(c("--n_sims"), type = "integer", default = 100, 
              help = "Number of simulations to run [default: %default]"),
  make_option(c("--seed_start"), type = "integer", default = 40001, 
              help = "Starting seed number [default: %default]"),
  make_option(c("--m_iter"), type = "integer", default = 1000, 
              help = "Number of jackknife iterations (m) [default: %default]"),
  make_option(c("--out_dir"), type = "character", default = "jackknife", 
              help = "Output directory [default: %default]")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Print configuration
print("Running Block Jackknife Simulation with parameters:")
print(opt)

# ── 2) Setup Variables 
maf_0       <- opt$maf
n_families  <- opt$n_families
beta0       <- opt$beta
sigma1      <- opt$sigma1
tau         <- opt$tau
n_sims      <- opt$n_sims
seed_start  <- opt$seed_start
m_iter      <- opt$m_iter

# Constants
n_snps <- 1
sigma0 <- 4
rho_e  <- 0.06

# Create output directory
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
         rbinom(length(parent_geno), 1, 0.5),
         parent_geno / 2)
}

# Updated to robustly handle vector or matrix input for delta_X
block_jackknife_random <- function(delta_X, delta_y, d, m) {

  if (is.null(dim(delta_X))) {
    delta_X <- as.matrix(delta_X)
  }
  
  n_families <- nrow(delta_X)
  n_snps <- ncol(delta_X)
  
  block_jackknife_results <- data.frame(
    SNP = character(),
    beta = numeric(),
    block_jackknife_var = numeric(),
    block_jackknife_bias = numeric(),
    block_jackknife_mean = numeric(),
    stringsAsFactors = FALSE
  )
  
  r <- n_families - d  # size of retained sample
  
  for (j in 1:n_snps) {
    beta_block_jack <- numeric(m)
    
    # Calculate full model beta once
    model_full <- lm(delta_y ~ 0 + delta_X[, j])
    beta_full  <- coef(model_full)[1]
    
    for (b in 1:m) {
      omit <- sample(1:n_families, d, replace = FALSE)
      keep <- setdiff(1:n_families, omit)
      
      # Safety check
      if (length(keep) < 2) next 
      
      # Fit model on subset
      model_sub <- lm(delta_y[keep] ~ 0 + delta_X[keep, j])
      beta_block_jack[b] <- coef(model_sub)[1]
    }
    
    beta_mean <- mean(beta_block_jack)
    
    # Variance formula for delete-d jackknife
    v_jd <- (r / (d * m)) * sum((beta_block_jack - beta_full)^2)
    
    bias <- (m - 1) * (beta_mean - beta_full)
    
    block_jackknife_results <- rbind(block_jackknife_results, data.frame(
      SNP = paste0("SNP", j),
      beta = beta_full,
      block_jackknife_var = v_jd,
      block_jackknife_bias = bias,
      block_jackknife_mean = beta_mean
    ))
  }
  
  return(block_jackknife_results)
}

# ── 4) Main Simulation Loop
df_list <- list()
seed_sequence <- seed_start:(seed_start + n_sims - 1)

for (i in seq_along(seed_sequence)) {
  current_seed <- seed_sequence[i]
  set.seed(current_seed)
  
  if (i %% 10 == 0) print(paste("Processing seed:", current_seed))
  
  maf <- rep(maf_0, n_snps)
  
  # A) Genotypes
  father_geno <- simulate_parents(maf, n_families)
  mother_geno <- simulate_parents(maf, n_families)
  
  sibling1_geno <- transmit_allele(father_geno) + transmit_allele(mother_geno)
  sibling2_geno <- transmit_allele(father_geno) + transmit_allele(mother_geno)
  
  delta_X <- sibling1_geno - sibling2_geno
  
  # B) Phenotypes
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
  
  # C) Normal Theory Analysis (as baseline)
  model <- lm(delta_y ~ 0 + delta_X)
  beta_hat <- coef(model)[1]
  se_hat   <- summary(model)$coefficients[1, "Std. Error"]
  
  # D) Block Jackknife Analysis
  d_size <- min(40, floor(n_families / 5)) # 'floor' adds safety for odd numbers
  
  block_df <- block_jackknife_random(delta_X, delta_y, d = d_size, m = m_iter)
  
  # E) Store Results
  df_list[[i]] <- data.frame(
    seed = current_seed,
    SNP = "SNP1",
    beta = beta_hat,
    bjk_var = block_df$block_jackknife_var,
    normal_var = se_hat^2
  )
}

# ── 5) Save Output ────────────────────────────────────────────────
df_all <- do.call(rbind, df_list)

filename <- paste0("jackknife",
                   ".n_", n_families, 
                   ".beta_", beta0, 
                   ".tau_", tau, 
                   ".maf_", maf_0, 
                   ".sigma1_", sigma1, 
                   ".start_", seed_start,
                   ".ite_", m_iter,
                   ".txt")

full_path <- file.path(opt$out_dir, filename)
write.table(df_all, full_path, row.names = FALSE, quote = FALSE)

print(paste("Done! Results saved to:", full_path))
