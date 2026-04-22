#!/usr/bin/env Rscript
library(MASS)
library(optparse)
library(sandwich) # For vcovCL (Cluster-Robust SE)
library(lmtest) 

# ── 1) Define Command Line Arguments ──────────────────────────────
option_list <- list(
  make_option(c("--maf1"), type = "numeric", default = 0.1, 
              help = "Minor Allele Frequency for Pop 1 [default: %default]"),
  make_option(c("--maf2"), type = "numeric", default = 0.4, 
              help = "Minor Allele Frequency for Pop 2 [default: %default]"),
  make_option(c("--beta"), type = "numeric", default = 0.0, 
              help = "True allelic effect size [default: %default]"),
  make_option(c("--pop_effect"), type = "numeric", default = 1.0, 
              help = "Baseline phenotypic difference between Pop 1 and Pop 2 [default: %default]"),
  make_option(c("--sigma1"), type = "numeric", default = 0.0, 
              help = "Heteroskedasticity scale [default: %default]"),
  make_option(c("--tau"), type = "numeric", default = 0.0, 
              help = "Between-family effect SD (family cluster variance) [default: %default]"),
  make_option(c("--n_families"), type = "integer", default = 5000, 
              help = "Total number of families per simulation [default: %default]"),
  make_option(c("--n_sims"), type = "integer", default = 5000, 
              help = "Number of simulations to run [default: %default]"),
  make_option(c("--m_iters"), type = "integer", default = 500, 
              help = "Number of iterations for perm and bjk [default: %default]"),
  make_option(c("--seed_start"), type = "integer", default = 1, 
              help = "Starting seed number [default: %default]"), 
  make_option(c("--out_dir"), type = "character", default = "results/", 
              help = "Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── 2) Setup Variables ────────────────────────────────────────────
n_families  <- opt$n_families
n_sims      <- opt$n_sims
m_iters     <- opt$m_iters
sigma0      <- 4
maf_1       <- opt$maf1
maf_2 = opt$maf2
beta0       <- opt$beta
sigma1      <- opt$sigma1
tau         <- opt$tau
seed_start  <- opt$seed_start
out_dir = opt$out_dir
rho_e = 0.06
if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE)
}

# ── 3) Helper Functions ───────────────────────────────────────────
simulate_parents <- function(maf_vec) {
  # maf_vec is now a vector of length n_families containing the specific MAF for that family's population
  p0 <- (1 - maf_vec)^2
  p1 <- 2 * maf_vec * (1 - maf_vec)
  p2 <- maf_vec^2
  
  # Sample genotypes row-by-row based on the specific MAF probabilities
  sapply(1:length(maf_vec), function(i) {
    sample(0:2, 1, prob = c(p0[i], p1[i], p2[i]))
  })
}

transmit_allele <- function(parent_geno) {
  ifelse(parent_geno == 1,
         rbinom(length(parent_geno), 1, 0.5),
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

block_jackknife_random <- function(delta_X, delta_y, d, m) {
  
  if (is.null(dim(delta_X))) {
    delta_X <- as.matrix(delta_X)
  }
  
  n_families <- nrow(delta_X)
  n_snps <- ncol(delta_X)
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
# ── 4) Main Simulation Loop ───────────────────────────────────────
df_list <- list()

for (i in 1:n_sims) {
  set.seed(opt$seed_start + i - 1)
  if (i %% 100 == 0) print(paste("Processing simulation:", i, "/", n_sims)) 
  
  # A) Define Populations (Half Pop 1, Half Pop 2)
  pop_id <- rep(c(1, 2), each = n_families / 2)
  maf_vec <- ifelse(pop_id == 1, opt$maf1, opt$maf2)
  
  # B) Genotypes
  father_geno <- simulate_parents(maf_vec)
  mother_geno <- simulate_parents(maf_vec)
  
  sib1_geno <- transmit_allele(father_geno) + transmit_allele(mother_geno)
  sib2_geno <- transmit_allele(father_geno) + transmit_allele(mother_geno)
  delta_X <- sib1_geno - sib2_geno
  
  var_delta_eps <- sigma0 + sigma1 * (delta_X^2)
  var_delta_eps <- pmax(var_delta_eps, 1e-8)
  
  scale_i <- sqrt(var_delta_eps / (2 * (1 - rho_e)))
  
  Sigma_Z <- matrix(c(1, rho_e, rho_e, 1), nrow = 2)
  Z <- mvrnorm(n = n_families, mu = c(0, 0), Sigma = Sigma_Z)
  
  epsilon1 <- scale_i * Z[,1]
  epsilon2 <- scale_i * Z[,2]
  u <- rnorm(n_families, 0, tau)
  beta_i_vec <- beta0 + u
  
  y1 <- beta_i_vec * sib1_geno + epsilon1
  y2 <- beta_i_vec * sib2_geno + epsilon2
  delta_y <- y1 - y2        
  
  # E) Analysis (Fit Linear Model controlling for Population)
  model <- lm(delta_y ~ 0 + delta_X)
  # Robust Standard Errors 
  robust_vcov  <- vcovHC(model, type = "HC1") 
  robust_stats <- coeftest(model, vcov = robust_vcov)
  se_robust  <- robust_stats["delta_X", "Std. Error"]
  # Extract results
  beta_hat <- coef(model)[1]
  se_hat   <- summary(model)$coefficients[1, "Std. Error"]

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

# ── 5) Save Output ────────────────────────────────────────────────
df_all <- do.call(rbind, df_list)
emp_var <- var(df_all$beta, na.rm = TRUE)

df_all$cluster_bias <- df_all$cluster_var - emp_var
df_all$normal_bias  <- df_all$normal_var - emp_var
df_all$perm_bias    <- df_all$perm_var - emp_var
df_all$bjk_bias     <- df_all$bjk_var - emp_var

print("Simulation complete.")
filename <- paste0(opt$out_dir,
                   "pop_structure_parents.compare.n_", n_families, 
                   ".beta_", beta0, 
                   ".tau_", tau, 
                   ".maf1_", maf_1,
                   ".maf2_", maf_2,
                   ".sigma1_", sigma1,
                   ".txt")

# full_path <- file.path(opt$out_dir, filename)
write.table(df_all, filename, row.names = FALSE, quote = FALSE)

print(paste("Done! Results saved to:", filename))
