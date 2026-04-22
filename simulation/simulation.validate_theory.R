#!/usr/bin/env Rscript
library(MASS)
library(optparse)
library(sandwich) # For Cluster-Robust SE
library(lmtest) 
# ── 1) Define Command Line Arguments
option_list <- list(
  make_option(c("--maf"), type = "numeric", default = 0.1, 
              help = "Minor Allele Frequency [default: %default]"),
  make_option(c("--beta"), type = "numeric", default = 0.0, 
              help = "Beta 0 value (Mean effect size) [default: %default]"),
  make_option(c("--sigma1"), type = "numeric", default = 0.0, 
              help = "Sigma 1 value (Heteroskedasticity) [default: %default]"),
  make_option(c("--tau"), type = "numeric", default = 1, 
              help = "Between-family effect SD (tau) [default: %default]"),
  make_option(c("--n_families"), type = "integer", default = 5000, 
              help = "Number of families per simulation [default: %default]"),
  make_option(c("--n_sims"), type = "integer", default = 5000, 
              help = "Number of simulations to run [default: %default]"),
  make_option(c("--m_iters"), type = "integer", default = 500, 
              help = "Number of iterations for perm and bjk [default: %default]"),
  make_option(c("--seed_start"), type = "integer", default = 1, 
              help = "Starting seed number [default: %default]"),
  make_option(c("--out_dir"), type = "character", default = "results/", 
              help = "Output directory [default: %default]")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Print configuration
print("Running Normal Theory Simulation with parameters:")
print(opt)

# ── 2) Setup Variables ────────────────────────────────────────────
# Assign CLI arguments to variables
maf_0       <- opt$maf
n_families  <- opt$n_families
beta0       <- opt$beta
sigma1      <- opt$sigma1
tau         <- opt$tau
n_sims      <- opt$n_sims
seed_start  <- opt$seed_start
m_iters = opt$m_iters
out_dir = opt$out_dir

# Constants (Fixed parameters)
n_snps <- 1
sigma0 <- 4
rho_e  <- 0.06 # Within-family shared background
library(MASS) # For mvrnorm

# 1. Helper: Transmit allele 
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


# --- MAIN SIMULATION ---

# Parameters
maf <- rep(maf_0, n_snps)

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
  }
  return(beta_block_jack)
}

# ── 4) Main Simulation Loop
df_list <- list()
seed_sequence <- seed_start:(seed_start + n_sims - 1)

# Use index 'i' for list storage, but 'current_seed' for RNG
for (i in seq_along(seed_sequence)) {
  current_seed <- seed_sequence[i]
  set.seed(current_seed)
  
  if (i %% 100 == 0) print(paste("Processing simulation:", i, "/", n_sims)) 
  
  # A) Genotypes: The Pedigree Approach
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
  # D) Analysis (Fit Linear Model)
  model <- lm(delta_y ~ 0 + delta_X)
  # Robust Standard Errors 
  robust_vcov  <- vcovHC(model, type = "HC1") 
  robust_stats <- coeftest(model, vcov = robust_vcov)
  se_robust  <- robust_stats["delta_X", "Std. Error"]
  # Extract results
  beta_hat <- coef(model)[1]
  se_hat   <- summary(model)$coefficients[1, "Std. Error"]
  p_val    <- summary(model)$coefficients[1, "Pr(>|t|)"]
  
  B <- m_iters
  perm_beta_matrix <- replicate(B, run_one_sign_flip_permutation(delta_y, delta_X))
  perm_var <- var(perm_beta_matrix)
  
  d_size <- min(40, floor(n_families / 5)) # 'floor' adds safety for odd numbers
  block_beta <- block_jackknife_random(delta_X, delta_y, d = d_size, m = m_iters)
  beta_mean <- mean(block_beta)
  r = n_families - d_size
  bjk_var <- (r / (d_size * (m_iters - 1)  )) * sum((block_beta - beta_mean)^2)
  # Store results
  df_list[[i]] <- data.frame(
    beta = beta_hat,
    cluster_var   = se_robust^2,
    normal_var = se_hat^2,
    perm_var = perm_var,
    bjk_var = bjk_var
  )
}

# ── 5) Save Output 
df_all <- do.call(rbind, df_list)
true_var <- function(p, n, tau, sigma1, sigma0) {
  q <- 1 - p
  term_1 <- 1 / (n * 2 * p * q)
  term_2 <- (1 + 3 * p * q) * (tau^2 + sigma1)
  result <- term_1 * (term_2 + sigma0)
  
  return(result)
}
variance = true_var(p = maf_0, n = n_families, tau = tau, sigma1 = sigma1, sigma0 = sigma0)

df_all$cluster_bias <- df_all$cluster_var - variance
df_all$normal_bias  <- df_all$normal_var - variance
df_all$perm_bias    <- df_all$perm_var - variance
df_all$bjk_bias     <- df_all$bjk_var - variance

filename <- paste0(out_dir,
                   "compare.n_", n_families, 
                   ".beta_", beta0, 
                   ".indep_tau_", tau, 
                   ".maf_", maf_0, 
                   ".indep_sigma1_", sigma1,
                   ".m_iter.",  m_iters,
                   ".txt")

write.table(df_all, filename, row.names = FALSE, quote = FALSE, sep = "\t")

print(paste("Done! Results saved to:", filename))
