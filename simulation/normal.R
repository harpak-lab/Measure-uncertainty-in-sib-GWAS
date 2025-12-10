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
  make_option(c("--n_sims"), type = "integer", default = 5000, 
              help = "Number of simulations to run [default: %default]"),
  make_option(c("--seed_start"), type = "integer", default = 1, 
              help = "Starting seed number [default: %default]"),
  make_option(c("--out_dir"), type = "character", default = "normal_theory", 
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

# Constants (Fixed parameters)
n_snps <- 1
sigma0 <- 4
rho_e  <- 0.06 # Within-family shared background

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
  # If genotype==1 (het), pick 0/1 with prob .5; else deterministic
  ifelse(parent_geno == 1,
         rbinom(length(parent_geno), 1, 0.5),
         parent_geno / 2)
}

# ── 4) Main Simulation Loop
df_list <- list()
seed_sequence <- seed_start:(seed_start + n_sims - 1)

# Use index 'i' for list storage, but 'current_seed' for RNG
for (i in seq_along(seed_sequence)) {
  current_seed <- seed_sequence[i]
  set.seed(current_seed)
  
  if (i %% 100 == 0) print(paste("Processing simulation:", i, "/", n_sims)) 
  
  maf <- rep(maf_0, n_snps)
  
  # A) Genotypes
  father_geno <- simulate_parents(maf, n_families)
  mother_geno <- simulate_parents(maf, n_families)
  
  # B) Transmit
  sibling1_geno <- transmit_allele(father_geno) + transmit_allele(mother_geno)
  sibling2_geno <- transmit_allele(father_geno) + transmit_allele(mother_geno)
  
  delta_X <- sibling1_geno - sibling2_geno  # within‐family genotype difference
  
  # C) Simulate Phenotypes
  # 1) Desired var(Δε_i)
  var_delta_eps <- sigma0 + sigma1 * (delta_X^2)
  var_delta_eps <- pmax(var_delta_eps, 1e-8)
  
  # 2) Scale factor
  scale_i <- sqrt(var_delta_eps / (2 * (1 - rho_e)))
  
  # 3) Correlated standard normals
  Sigma_Z <- matrix(c(1, rho_e, rho_e, 1), nrow = 2)
  Z <- mvrnorm(n = n_families, mu = c(0, 0), Sigma = Sigma_Z)
  
  # 4) Scale residuals
  epsilon1 <- scale_i * Z[,1]
  epsilon2 <- scale_i * Z[,2]
  
  # 5) Family-specific allelic effect
  u <- rnorm(n_families, 0, tau)
  beta_i_vec <- beta0 + u
  
  # 6) Final Phenotypes
  y1 <- beta_i_vec * sibling1_geno + epsilon1
  y2 <- beta_i_vec * sibling2_geno + epsilon2
  
  delta_y <- y1 - y2
  
  # D) Analysis (Fit Linear Model)
  model <- lm(delta_y ~ 0 + delta_X)
  
  # Extract results
  beta_hat <- coef(model)[1]
  se_hat   <- summary(model)$coefficients[1, "Std. Error"]
  p_val    <- summary(model)$coefficients[1, "Pr(>|t|)"]
  
  # Store results
  df_list[[i]] <- data.frame(
    seed = current_seed,
    SNP = "SNP1", 
    beta = beta_hat, 
    se = se_hat, 
    p = p_val
  )
}

# ── 5) Save Output 
df_all <- do.call(rbind, df_list)

# Construct dynamic filename
filename <- paste0("normal_method",
                   ".n_", n_families, 
                   ".beta_", beta0, 
                   ".tau_", tau, 
                   ".maf_", maf_0, 
                   ".sigma1_", sigma1, 
                   ".start_", seed_start,
                   ".txt")

full_path <- file.path(opt$out_dir, filename)
write.table(df_all, full_path, row.names = FALSE, quote = FALSE)

print(paste("Done! Results saved to:", full_path))
