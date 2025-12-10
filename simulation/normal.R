args = commandArgs(trailingOnly=TRUE) 
# beta_i <- as.numeric(args[1])
maf_0 <- as.numeric(args[1])
library(MASS)
seeds=1
df_list <- list()
for(seed in seeds:(seeds+4999)){
  set.seed(seed)
  # print(paste("Start for seed", seed))
  n_snps      <- 1        
  n_families  <- 5000
  # maf_0= 0.1
  maf         <- rep(maf_0, n_snps)
  
  # ── 1) Simulate parental genotypes under HWE ─────────────────────
  # returns an n_families × n_snps matrix
  simulate_parents <- function(maf, n_families) {
    p0 <- (1 - maf)^2
    p1 <- 2 * maf * (1 - maf)
    p2 <- maf^2
    replicate(length(maf),
              sample(0:2, n_families, replace = TRUE, prob = c(p0, p1, p2)))
  }
  father_geno <- simulate_parents(maf, n_families)
  mother_geno <- simulate_parents(maf, n_families)
  
  # ── 2) Transmit one allele from each parent ────────────────────────
  # If genotype==1 (het), pick 0/1 with prob .5; else deterministic
  transmit_allele <- function(parent_geno) {
    ifelse(parent_geno == 1,
           rbinom(n_families, 1, 0.5),
           parent_geno / 2)
  }
  
  # Sibling 1
  paternal1 <- transmit_allele(father_geno)
  maternal1 <- transmit_allele(mother_geno)
  sibling1_geno <- paternal1 + maternal1
  
  # Sibling 2 (independent draw of alleles)
  paternal2 <- transmit_allele(father_geno)
  maternal2 <- transmit_allele(mother_geno)
  sibling2_geno <- paternal2 + maternal2
  
  delta_X <- sibling1_geno - sibling2_geno  # within‐family genotype difference
  
  # ── 3) Simulate phenotypes ─────────────────────────────────────────
  library(MASS)
  
  # parameters
  beta0  <- 0
  tau    <- 0.5           # between-family effect SD
  sigma0 <- 4
  sigma1 <- 0     # your chosen value
  rho_e  <- 0.06
  
  # 1) Desired var(Δε_i) for each family
  var_delta_eps <- sigma0 + sigma1 * (delta_X^2)      # length n_families
  var_delta_eps <- pmax(var_delta_eps, 1e-8)          # safety
  
  # 2) Scale factor so that Var(ε1 - ε2) = var_delta_eps
  #    because Var(Z1 - Z2) = 2(1 - rho_e)
  scale_i <- sqrt( var_delta_eps / (2 * (1 - rho_e)) )
  
  # 3) Draw correlated standard normals
  Sigma_Z <- matrix(c(1, rho_e,
                      rho_e, 1), nrow = 2)
  Z <- mvrnorm(n = n_families, mu = c(0, 0), Sigma = Sigma_Z)
  
  # 4) Scale by family-specific factor
  epsilon1 <- scale_i * Z[,1]
  epsilon2 <- scale_i * Z[,2]
  
  # 5) Family-specific allelic effect
  u      <- rnorm(n_families, 0, tau)
  beta_i <- beta0 + u
  
  # 6) Phenotypes
  y1 <- beta_i * sibling1_geno + epsilon1
  y2 <- beta_i * sibling2_geno + epsilon2
  
  delta_y <- y1 - y2             # sibling‐GWAS outcome
  
  results <- data.frame(SNP = character(), beta = numeric(), se = numeric(), p = numeric())
    model <- lm(delta_y ~ 0 + delta_X)
    summary_model <- summary(model)
    beta_hat <- coef(model)[1]
    se_hat <-  summary(model)$coefficients[1, "Std. Error"]
    p_val <- summary(model)$coefficients[1, "Pr(>|t|)"]
    results <- rbind(results, data.frame(SNP = paste0("SNP", 1), beta = beta_hat, se = se_hat, p = p_val))
  rownames(results) <- NULL
  df_list[[seed]] <- results
}
df_all <- do.call(rbind, df_list)
write.table(df_all, paste0("normal_theory/normal_method.n_", n_families, ".beta_0_", beta0, ".tau_", tau, ".maf_", maf_0, ".sigma_1_", sigma1,  ".sim_5k.txt"))
