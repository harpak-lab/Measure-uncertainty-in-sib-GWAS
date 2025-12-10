#VARY_SIGMA1
# parameters
n_families <- 5000
p <- 0.1
q <- 1 - p
sigma0 <- 4
sigma1_list <- seq(-0.5, 0.5, by = 0.1)
tau <- 0   # given

# constants from sib transmission
E_dX2  <- 2 * p * q
E_dX4  <- 2 * p * q * (1 + 3 * p * q)   # for reference
fac <- (1 + 3 * p * q)                  # appears in formula

# true variance formula for beta_hat
# true_var_list = list()
# true_var <- ( (fac * (tau^2 + sigma1_list)) + sigma0 ) / (n_families * 2 * p * q)
df_list = list()
for(sigma in sigma1_list){
  df_list[[as.character(sigma)]] = read.table(
    paste0("permutation/permutation.n_",
           n_families,
           ".beta_", 0,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1.", sigma,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(sigma)]]$sigma_1 = sigma
  df_list[[as.character(sigma)]]$perm_bias = df_list[[as.character(sigma)]]$perm_var  - true_var
}
df <- do.call(rbind, df_list)
rownames(df) = NULL
write.table(df,  paste0("vary_sigma/permutation.n_",
                               n_families,
                               ".beta_", 0,
                               ".tau_", tau,
                               ".maf_", p,
                               ".sigma_1_-0.5_0.5.sim_5k.txt"))


# VARY_N_FAMILIES
n_families = 5000
p_lsit <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
q_list <- 1 - maf_list
sigma0 <- 4
sigma1 = 0
beta = 1
tau = 0

# constants from sib transmission
E_dX2  <- 2 * p * q
E_dX4  <- 2 * p * q * (1 + 3 * p * q)   # for reference
fac <- (1 + 3 * p * q)                  # appears in formula

# true_var <- ( (fac * (tau^2 + sigma1_list)) + sigma0 ) / (n_families * 2 * p * q)
df_list = list()
for(n_families in n_families_list){
  df_list[[as.character(n_families)]] = read.table(
    paste0("normal_theory/normal_method.n_",
           n_families,
           ".beta_0_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1_", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(n_families)]]$n_families = n_families
  df_list[[as.character(n_families)]]$normal_bias = df_list[[as.character(n_families)]]$se^2  - true_var
}

df_list = list()
for(n_families in n_families_list){
  df_list[[as.character(n_families)]] = read.table(
    paste0("permutation/permutation.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1.", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(n_families)]]$n_families = n_families
  df_list[[as.character(n_families)]]$perm_bias = df_list[[as.character(n_families)]]$perm_var  - true_var
}

df_list = list()
for(n_families in n_families_list){
  df_list[[as.character(n_families)]] = read.table(
    paste0("jackknife/jackknife.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(n_families)]]$n_families = n_families
  df_list[[as.character(n_families)]]$bjk_bias = df_list[[as.character(n_families)]]$bjk_var  - true_var
}
df <- do.call(rbind, df_list)
rownames(df) = NULL
write.table(df,  paste0("vary_n_families/normal.n_",
                        "1k_10k",
                        ".beta_", beta,
                        ".tau_", tau ,
                        ".maf_", p,
                        ".sigma_1_", sigma1, "sim_5k.txt"))

#VARY N_FAMILIES
n_families_list <- c(1000, 2000, 3000, 4000, 5000, 10000)
p <- 0.1
q <- 1 - p
sigma0 <- 4
sigma1 = 0
beta = 1
tau = 0

# constants from sib transmission
E_dX2  <- 2 * p * q
E_dX4  <- 2 * p * q * (1 + 3 * p * q)   # for reference
fac <- (1 + 3 * p * q)                  # appears in formula

# true_var <- ( (fac * (tau^2 + sigma1_list)) + sigma0 ) / (n_families * 2 * p * q)
df_list = list()
for(n_families in n_families_list){
  df_list[[as.character(n_families)]] = read.table(
    paste0("normal_theory/normal_method.n_",
           n_families,
           ".beta_0_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1_", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(n_families)]]$n_families = n_families
  df_list[[as.character(n_families)]]$normal_bias = df_list[[as.character(n_families)]]$se^2  - true_var
}

df_list = list()
for(n_families in n_families_list){
  df_list[[as.character(n_families)]] = read.table(
    paste0("permutation/permutation.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1.", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(n_families)]]$n_families = n_families
  df_list[[as.character(n_families)]]$perm_bias = df_list[[as.character(n_families)]]$perm_var  - true_var
}

df_list = list()
for(n_families in n_families_list){
  df_list[[as.character(n_families)]] = read.table(
    paste0("jackknife/jackknife.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(n_families)]]$n_families = n_families
  df_list[[as.character(n_families)]]$bjk_bias = df_list[[as.character(n_families)]]$bjk_var  - true_var
}
df <- do.call(rbind, df_list)
rownames(df) = NULL
write.table(df,  paste0("jackknife.n_",
                        "1k_10k",
                        ".beta_", beta,
                        ".tau_", tau ,
                        ".maf_", p,
                        ".sigma_1_", sigma1, "sim_5k.txt"))

n_families_list <- c(1000, 2000, 3000, 4000, 5000, 10000)
p <- 0.1
q <- 1 - p
sigma0 <- 4
sigma1 = 0
beta = 1
tau = 0

# constants from sib transmission
E_dX2  <- 2 * p * q
E_dX4  <- 2 * p * q * (1 + 3 * p * q)   # for reference
fac <- (1 + 3 * p * q)                  # appears in formula

# true_var <- ( (fac * (tau^2 + sigma1_list)) + sigma0 ) / (n_families * 2 * p * q)
df_list = list()
for(n_families in n_families_list){
  df_list[[as.character(n_families)]] = read.table(
    paste0("normal_theory/normal_method.n_",
           n_families,
           ".beta_0_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1_", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(n_families)]]$n_families = n_families
  df_list[[as.character(n_families)]]$normal_bias = df_list[[as.character(n_families)]]$se^2  - true_var
}

df_list = list()
for(n_families in n_families_list){
  df_list[[as.character(n_families)]] = read.table(
    paste0("permutation/permutation.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1.", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(n_families)]]$n_families = n_families
  df_list[[as.character(n_families)]]$perm_bias = df_list[[as.character(n_families)]]$perm_var  - true_var
}

df_list = list()
for(n_families in n_families_list){
  df_list[[as.character(n_families)]] = read.table(
    paste0("jackknife/jackknife.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(n_families)]]$n_families = n_families
  df_list[[as.character(n_families)]]$bjk_bias = df_list[[as.character(n_families)]]$bjk_var  - true_var
}
df <- do.call(rbind, df_list)
rownames(df) = NULL
write.table(df,  paste0("vary_n_families/normal.n_",
                        "1k_10k",
                        ".beta_", beta,
                        ".tau_", tau ,
                        ".maf_", p,
                        ".sigma_1_", sigma1, "sim_5k.txt"))
#VARY BETA
n_families <- 5000
p <- 0.1
q <- 1 - p
sigma0 <- 4
sigma1 = 0
beta_list = seq(-0.2, 0.5, by = 0.1)
tau <- 0

# constants from sib transmission
E_dX2  <- 2 * p * q
E_dX4  <- 2 * p * q * (1 + 3 * p * q)   # for reference
fac <- (1 + 3 * p * q)                  # appears in formula

# true variance formula for beta_hat
# true_var_list = list()
# true_var <- ( (fac * (tau^2 + sigma1_list)) + sigma0 ) / (n_families * 2 * p * q)
df_list = list()
for(beta in beta_list){
  df_list[[as.character(beta)]] = read.table(
    paste0("normal_theory/normal_method.n_",
           n_families,
           ".beta_0_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1_", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(beta)]]$beta = beta
  df_list[[as.character(beta)]]$normal_bias = df_list[[as.character(beta)]]$se^2  - true_var
}

df_list = list()
for(beta in beta_list){
  df_list[[as.character(beta)]] = read.table(
    paste0("permutation/permutation.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1.", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(beta)]]$beta = beta
  df_list[[as.character(beta)]]$perm_bias = df_list[[as.character(beta)]]$perm_var  - true_var
}
df_list = list()
for(beta in beta_list){
  df_list[[as.character(beta)]] = read.table(
    paste0("jackknife/jackknife.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(beta)]]$beta = beta
  df_list[[as.character(beta)]]$bjk_bias = df_list[[as.character(beta)]]$bjk_var  - true_var
}
df <- do.call(rbind, df_list)
rownames(df) = NULL
write.table(df,  paste0("vary_tau/normal.n_",
                        n_families,
                        ".beta_", beta,
                        ".tau_", "0_1" ,
                        ".maf_", p,
                        ".sigma_1_0.sim_5k.txt"))
#VARY TAU
n_families <- 5000
p <- 0.1
q <- 1 - p
sigma0 <- 4
sigma1 = 0
beta = 0
tau_list <- seq(0, 0.7, by = 0.1)   # given

# constants from sib transmission
E_dX2  <- 2 * p * q
E_dX4  <- 2 * p * q * (1 + 3 * p * q)   # for reference
fac <- (1 + 3 * p * q)                  # appears in formula

# true variance formula for beta_hat
# true_var_list = list()
# true_var <- ( (fac * (tau^2 + sigma1_list)) + sigma0 ) / (n_families * 2 * p * q)
df_list = list()
for(tau in tau_list){
  df_list[[as.character(tau)]] = read.table(
    paste0("normal_theory/normal_method.n_",
           n_families,
           ".beta_0_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1_", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(tau)]]$tau = tau
  df_list[[as.character(tau)]]$normal_bias = df_list[[as.character(beta)]]$se^2  - true_var
}
df <- do.call(rbind, df_list)

df_list = list()
for(tau in tau_list){
  df_list[[as.character(tau)]] = read.table(
    paste0("permutation/permutation.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           ".sigma_1.", sigma1,   # remove the dot before sigma_1
           ".sim_5k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(tau)]]$tau = tau
  df_list[[as.character(tau)]]$perm_bias = df_list[[as.character(tau)]]$perm_var  - true_var
}
df <- do.call(rbind, df_list)
rownames(df) = NULL
write.table(df,  paste0("vary_beta/normal.n_",
                        n_families,
                        ".beta_", "0_0.5",
                        ".tau_", tau,
                        ".maf_", p,
                        ".sigma_1_0.sim_5k.txt"))

df_list = list()
for(tau in tau_list){
  df_list[[as.character(tau)]] = read.table(
    paste0("jackknife/jackknife.n_",
           n_families,
           ".beta_", beta,
           ".tau_", tau,
           ".maf_", p,
           "sigma_1", sigma1,   # remove the dot before sigma_1
           ".sim_5k.ite_1k.txt"),
    header = TRUE
  )
  true_var =  ( (fac * (tau^2 + sigma1)) + sigma0 ) / (n_families * 2 * p * q)
  df_list[[as.character(tau)]]$tau = tau
  df_list[[as.character(tau)]]$bjk_bias = df_list[[as.character(tau)]]$bjk_var  - true_var
}
df <- do.call(rbind, df_list)
rownames(df) = NULL
write.table(df,  paste0("vary_tau/jackknife.n_",
                        n_families,
                        ".beta_", beta,
                        ".tau_", "0_0.7" ,
                        ".maf_", p,
                        ".sigma_1_0.sim_5k.txt"))


# fixed parameters
n_families <- 5000
maf_list   <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
sigma0     <- 4
sigma1     <- 0
beta       <- 1
tau        <- 0

df_list <- list() 

for (p in maf_list) {
  q   <- 1 - p
  fac <- 1 + 3 * p * q    # (1 + 3pq) term

  fname <- paste0(
    "normal_theory/normal_method.n_",
    n_families,
    ".beta_0_", beta,
    ".tau_", tau,
    ".maf_", p,
    ".sigma_1_", sigma1,
    ".sim_5k.txt"
  )
  
  df_p <- read.table(fname, header = TRUE)
  
  # true variance for this MAF
  true_var <- ((fac * (tau^2 + sigma1)) + sigma0) / (n_families * 2 * p * q)
  
  df_p$maf         <- p
  df_p$true_var    <- true_var
  df_p$normal_bias <- df_p$se^2 - true_var
  
  df_list[[as.character(p)]] <- df_p
}

df_list = list()
for (p in maf_list) {
  q   <- 1 - p
  fac <- 1 + 3 * p * q    # (1 + 3pq) term
  
  fname <-  paste0("permutation/permutation.n_",
                   n_families,
                   ".beta_", beta,
                   ".tau_", tau,
                   ".maf_", p,
                   ".sigma_1.", sigma1,   # remove the dot before sigma_1
                   ".sim_5k.txt")
  
  df_p <- read.table(fname, header = TRUE)
  
  # true variance for this MAF
  true_var <- ((fac * (tau^2 + sigma1)) + sigma0) / (n_families * 2 * p * q)
  
  df_p$maf         <- p
  df_p$true_var    <- true_var
  df_p$perm_bias <- df_p$perm_var - true_var
  
  df_list[[as.character(p)]] <- df_p
}

df_list = list()

for (p in maf_list) {
  q   <- 1 - p
  fac <- 1 + 3 * p * q    # (1 + 3pq) term
  
  fname <-  paste0("jackknife/jackknife.n_",
                   n_families,
                   ".beta_", beta,
                   ".tau_", tau,
                   ".maf_", p,
                   ".sigma_1", sigma1,   # remove the dot before sigma_1
                   ".sim_5k.txt")
  
  df_p <- read.table(fname, header = TRUE)
  
  # true variance for this MAF
  true_var <- ((fac * (tau^2 + sigma1)) + sigma0) / (n_families * 2 * p * q)
  
  df_p$maf         <- p
  df_p$true_var    <- true_var
  df_p$bjk_bias <- df_p$bjk_var - true_var
  
  df_list[[as.character(p)]] <- df_p
}
df <- do.call(rbind, df_list)
rownames(df) <- NULL
