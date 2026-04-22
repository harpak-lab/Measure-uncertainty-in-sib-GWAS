##############################################################################
# Within-Family GWAS Simulation: Variance Estimator Comparison
#
# Model:
#   Y_ij = (beta0 + u_i) * G_ij + epsilon_ij
#   u_i  ~ N(0, tau^2)           [family-level random slope]
#   epsilon ~ MVN(0, Sigma_eps)  [correlated within-family noise]
#   Sigma_eps[j,j]  = sigma0 + sigma1 * V_i   [heteroskedastic variance]
#   Sigma_eps[j,k]  = rho_e * sqrt(diag_j * diag_k)  [within-family correlation]
#
# Estimand:
#   beta_within = sum(Y_ij_w * G_ij_w) / sum(G_ij_w^2)
#   where _w denotes within-family demeaning (family FE removed)
#
# Variance estimators compared:
#   1. Naive OLS SE
#   2. Cluster-robust (HC) SE  [clustered by family]
#   3. Within-family permutation
#   4. Block jackknife
##############################################################################

suppressPackageStartupMessages({
  library(MASS)      # mvrnorm
  library(sandwich)  # vcovCL
  library(lmtest)    # coeftest
  library(optparse)
})

# ── 1) Command-Line Arguments ─────────────────────────────────────────────────
option_list <- list(
  make_option("--maf",        type="numeric",   default=0.1),
  make_option("--n_sibs",     type="numeric",   default=3),
  make_option("--beta",       type="numeric",   default=0.0),
  make_option("--sigma1",     type="numeric",   default=0.0),
  make_option("--tau",        type="numeric",   default=0.0),
  make_option("--n_families", type="integer",   default=5000L),
  make_option("--n_sims",     type="integer",   default=5000L),
  make_option("--seed_start", type="integer",   default=1L),
  make_option("--out_dir",    type="character", default="results/")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ── 2) Fixed Parameters ───────────────────────────────────────────────────────
maf_0      <- opt$maf
n_sibs     <- opt$n_sibs         # siblings per family  (columns)
n_families <- opt$n_families     # families             (rows)
beta0      <- opt$beta
sigma1     <- opt$sigma1
tau        <- opt$tau
n_sims     <- opt$n_sims
out_dir    <- opt$out_dir
seed_start <- opt$seed_start

sigma0 <- 4      # baseline within-family noise variance
rho_e  <- 0.06   # within-family noise correlation
m_iter <- 500    # permutation / jackknife iterations

# ── 3) Helper: genotype transmission ─────────────────────────────────────────
# transmit() samples one allele from a parent:
#   homozygous (0 or 2) -> deterministic, heterozyous (1) -> Bernoulli(0.5)
transmit <- function(parent_geno) {
  ifelse(parent_geno == 1,
         rbinom(length(parent_geno), 1L, 0.5),
         parent_geno / 2L)
}

# ── 4) Within-Family Permutation ──────────────────────────────────────────────
# Permute genotype labels WITHIN each family, keeping phenotypes fixed.
#
# Arguments
#   Y_w_vec    : numeric vector length (n_families * n_sibs), column-major
#                [already within-family demeaned phenotype]
#   G_w_matrix : numeric matrix [n_families x n_sibs], within-family demeaned G
#   denom      : scalar = sum(G_w_vec^2), precomputed for speed
#   B          : number of permutations
#
# Column-major alignment guarantee:
#   as.vector(G_w_matrix) == G_w_vec  (both column-major)
#   After row-permutation: t(apply(G_w_matrix, 1, sample)) stays [n_fam x n_sib]
#   as.vector() on that result is still column-major -> aligns with Y_w_vec ✓
#
# Note: permuting within-family preserves sum(G_w_perm^2) = denom exactly,
# so reusing the original denom is correct.

run_within_family_permutations <- function(Y_w_vec, G_w_matrix, denom, B=m_iter) {
  if (denom < 1e-9) return(rep(NA_real_, B))

  n_fam <- nrow(G_w_matrix)
  n_sib <- ncol(G_w_matrix)
  stopifnot(length(Y_w_vec) == n_fam * n_sib)

  replicate(B, {
    # Shuffle sibling order within each family (row-wise)
    # Result is [n_fam x n_sib]; as.vector() flattens column-major
    G_w_perm_vec <- as.vector(t(apply(G_w_matrix, 1L, sample)))
    sum(Y_w_vec * G_w_perm_vec) / denom
  })
}

# ── 5) Block Jackknife ────────────────────────────────────────────────────────
# Randomly delete d families per draw, recompute beta, aggregate variance.
#
# Variance formula (Kunsch 1989 random-delete-d jackknife):
#   V_jk = (r / (d * m)) * sum( (beta_b - beta_full)^2 )
#   where r = n_families - d
#
# d choice: sqrt(n_families) is a standard bias/variance trade-off.
# Using floor(n_families * 0.2) with small n_families gives r/d >> 1
# and inflates the variance estimate; sqrt() is safer.

block_jackknife_random <- function(Y_w_vec, G_w_vec, fam_ids, d, m=m_iter) {
  stopifnot(length(Y_w_vec) == length(G_w_vec),
            length(G_w_vec) == length(fam_ids))

  unique_fams    <- unique(fam_ids)
  n_fams_total   <- length(unique_fams)
  r              <- n_fams_total - d

  denom_full <- sum(G_w_vec^2)
  if (denom_full < 1e-9) return(NA_real_)
  beta_full <- sum(Y_w_vec * G_w_vec) / denom_full

  beta_jack <- vapply(seq_len(m), function(b) {
    omit      <- sample(unique_fams, d, replace=FALSE)
    keep_idx  <- which(!(fam_ids %in% omit))
    G_sub     <- G_w_vec[keep_idx]
    Y_sub     <- Y_w_vec[keep_idx]
    denom_sub <- sum(G_sub^2)
    if (denom_sub < 1e-9) return(NA_real_)
    sum(Y_sub * G_sub) / denom_sub
  }, numeric(1L))

  beta_jack <- beta_jack[!is.na(beta_jack)]
  m_valid   <- length(beta_jack)
  if (m_valid == 0L) return(NA_real_)

  # Random-delete-d jackknife variance
  v_jd <- (r / (d * m_valid)) * sum((beta_jack - beta_full)^2)
  return(v_jd)
}

# ── 6) Main Simulation Loop ───────────────────────────────────────────────────

# Pre-build the within-family noise covariance template (correlation only;
# the diagonal scaling by fam_sds is applied per-family below).
Sigma_corr <- matrix(rho_e, nrow=n_sibs, ncol=n_sibs)
diag(Sigma_corr) <- 1.0

# Jackknife block size: sqrt(n_families) balances bias and variance
d_size <- max(1L, floor(sqrt(n_families)))

results_list <- vector("list", n_sims)

for (i in seq_len(n_sims)) {
  set.seed(seed_start + i - 1L)
  if (i %% 100L == 0L) message("Simulation ", i, " / ", n_sims)

  # ── A) Generate parental and offspring genotypes ──────────────────────────
  # Hardy-Weinberg probabilities
  p0 <- (1 - maf_0)^2
  p1 <- 2 * maf_0 * (1 - maf_0)
  p2 <- maf_0^2

  father_geno <- sample(0:2, n_families, replace=TRUE, prob=c(p0, p1, p2))
  mother_geno <- sample(0:2, n_families, replace=TRUE, prob=c(p0, p1, p2))

  # sib_geno_matrix : [n_families x n_sibs]  (row = family, col = sibling)
  sib_geno_matrix <- vapply(seq_len(n_sibs), function(s) {
    transmit(father_geno) + transmit(mother_geno)
  }, numeric(n_families))
  # vapply with FUN.VALUE=numeric(n_families) returns [n_families x n_sibs] ✓

  # ── B) Family-level random slopes ────────────────────────────────────────
  u <- rnorm(n_families, mean=0, sd=tau)  # length n_families

  # beta_i_matrix : [n_families x n_sibs]
  # Each row gets the same family-specific slope; expand to matrix for broadcast
  beta_i_matrix <- matrix(beta0 + u, nrow=n_families, ncol=n_sibs)

  # ── C) Simulate phenotypes ────────────────────────────────────────────────
  # Family-specific genotype variance: V_i = sum_j (G_ij - G_i_bar)^2
  fam_means_G <- rowMeans(sib_geno_matrix)   # length n_families
  V_i         <- rowSums((sib_geno_matrix - fam_means_G)^2)  # length n_families

  # Phenotypic noise SD per family (heteroskedastic)
  fam_sds <- sqrt(sigma0 + sigma1 * V_i)     # length n_families

  # Draw correlated noise: mvrnorm returns [n_families x n_sibs]
  raw_noise     <- mvrnorm(n=n_families, mu=rep(0, n_sibs), Sigma=Sigma_corr)
  # Scale each family's noise by its SD (row-wise multiplication)
  epsilon_matrix <- raw_noise * fam_sds      # fam_sds recycled column-wise

  # Phenotype matrix: [n_families x n_sibs]
  pheno_matrix <- beta_i_matrix * sib_geno_matrix + epsilon_matrix

  # ── D) Flatten to vectors (COLUMN-MAJOR throughout) ──────────────────────
  #
  # as.vector() on an [n_families x n_sibs] matrix uses column-major order:
  #   indices 1:n_families   -> sib 1 for all families
  #   indices (n_fam+1):(2*n_fam) -> sib 2 for all families, etc.
  #
  # fam_ids must use rep(..., times=n_sibs) to match this layout:
  #   rep(1:n_fam, times=3) = [1,2,...,n_fam, 1,2,...,n_fam, 1,2,...,n_fam]
  #   NOT rep(1:n_fam, each=n_sibs) which would be row-major
  #
  N_obs   <- n_families * n_sibs
  fam_ids <- rep(seq_len(n_families), times=n_sibs)   # length N_obs, col-major

  G_vec <- as.vector(sib_geno_matrix)   # length N_obs, col-major
  Y_vec <- as.vector(pheno_matrix)      # length N_obs, col-major

  # G_bar_vec: family mean of G repeated for each sib, column-major
  # rep(..., times=n_sibs) matches the same col-major layout as G_vec
  G_bar_vec <- rep(fam_means_G, times=n_sibs)         # length N_obs

  # Within-family demeaned genotype (col-major, aligned with G_vec)
  G_w_vec <- G_vec - G_bar_vec                         # length N_obs

  # Verify: G - fam_means on the matrix uses R's column-wise recycling,
  # so as.vector(sib_geno_matrix - fam_means_G) == G_w_vec ✓
  # (fam_means_G has length n_families = nrow, recycled across columns)
  G_w_matrix <- sib_geno_matrix - fam_means_G          # [n_families x n_sibs]

  # ── E) Within-family demeaned phenotype (FE residual) ────────────────────
  #
  # We want Y_w such that beta_within = sum(Y_w * G_w) / sum(G_w^2).
  # This is the Frisch-Waugh-Lovell (FWL) residual of Y after projecting
  # out family fixed effects.
  #
  # With balanced panels, FE residual = row-mean demeaning:
  #   Y_w_ij = Y_ij - Y_i_bar
  # This is EQUIVALENT to lm(Y ~ factor(fam_id)) residuals when n_sibs
  # is the same for every family (balanced), which is our case here.
  #
  # Note: the original code used rowMeans(pheno_matrix) which IS correct
  # for the balanced case.  We make it explicit below.

  fam_means_Y <- rowMeans(pheno_matrix)                # length n_families
  Y_w_matrix  <- pheno_matrix - fam_means_Y            # [n_families x n_sibs]
  Y_w_vec     <- as.vector(Y_w_matrix)                 # length N_obs, col-major

  # ── F) OLS with cluster-robust SE ────────────────────────────────────────
  # Model: Y ~ G_within + G_bar
  # G_bar is the between-family component (Mundlak correction).
  # Clustering by family accounts for within-family error correlation.

  df <- data.frame(
    fam_id   = factor(fam_ids),
    Y        = Y_vec,
    G_within = G_w_vec,
    G_bar    = G_bar_vec
  )

  model      <- lm(Y ~ G_within + G_bar, data=df)
  se_naive   <- summary(model)$coefficients["G_within", "Std. Error"]

  robust_vcov  <- vcovCL(model, cluster=~fam_id, data=df)
  robust_stats <- coeftest(model, vcov=robust_vcov)
  se_cluster   <- robust_stats["G_within", "Std. Error"]

  # ── G) Permutation variance ──────────────────────────────────────────────
  # Numerator identity: sum(Y_w * G_w) = sum(Y_w_mat * G_w_mat) element-wise
  # because both are column-major flattened from [n_fam x n_sib].
  # Denom = sum(G_w^2) is preserved under within-family permutation. ✓

  denom_SS   <- sum(G_w_vec^2)
  perm_betas <- run_within_family_permutations(Y_w_vec, G_w_matrix,
                                               denom=denom_SS, B=m_iter)
  perm_var   <- var(perm_betas, na.rm=TRUE)

  # ── H) Block jackknife variance ──────────────────────────────────────────
  bjk_var <- block_jackknife_random(Y_w_vec, G_w_vec,
                                    fam_ids=fam_ids,
                                    d=d_size, m=m_iter)

  results_list[[i]] <- data.frame(
    sim        = i,
    beta       = coef(model)["G_within"],
    normal_var  = se_naive^2,
    cluster_var = se_cluster^2,
    perm_var    = perm_var,
    bjk_var     = bjk_var
  )
}

# ── 7) Aggregate & Output ─────────────────────────────────────────────────────
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

final_results <- do.call(rbind, results_list)
emp_var <- var(final_results$beta)

# Bias = mean(estimated variance) - empirical variance
final_results$normal_bias  <- final_results$normal_var  - emp_var
final_results$cluster_bias <- final_results$cluster_var - emp_var
final_results$perm_bias    <- final_results$perm_var    - emp_var
final_results$bjk_bias     <- final_results$bjk_var     - emp_var

fname <- paste0(out_dir, "compare.beta_",beta0, ".n_fam_", n_families, ".n_sibs_", n_sibs   ,".maf_", maf_0, ".tau_", tau, ".sigma1_", sigma1, ".txt" )
write.table(final_results, fname, row.names=FALSE, quote=FALSE)
