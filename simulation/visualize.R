#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(optparse)
})

# ==============================================================================
# 1. SETUP: Command Line Arguments
# ==============================================================================
library(optparse)

option_list <- list(
  make_option(c("-v", "--vary"), type="character", default="maf", 
              help="Parameter varied"),
  make_option(c("--n_fixed"), type="numeric", default=5000),
  make_option(c("--beta_fixed"), type="numeric", default=0), 
  make_option(c("--maf_fixed"), type="numeric", default=0.1),
  make_option(c("--tau_fixed"), type="numeric", default=0),
  make_option(c("--sigma1_fixed"), type="numeric", default=0),
  make_option(c("--out"), type="character", default="plot.png")
)

opt <- parse_args(OptionParser(option_list=option_list))
results_dir <- "results" 

# 3. READ: Use the SAME Standardized Format
input_filename <- sprintf("comparison_df.%s_vary.beta_%s.n_%d.tau_%s.maf_%s.sigma1_%s.txt",
                          opt$vary,                        # 1. Key Param
                          as.character(opt$beta_fixed),    # 2. Beta
                          opt$n_fixed,                     # 3. N
                          as.character(opt$tau_fixed),     # 4. Tau
                          as.character(opt$maf_fixed),     # 5. MAF
                          as.character(opt$sigma1_fixed))  # 6. Sigma1

full_path <- file.path(results_dir, input_filename)
# 4. Read the file
if (!file.exists(full_path)) {
  stop(paste("File not found:", full_path, "\nCheck that your fixed parameters match the filename exactly."))
}

df <- read.table(full_path, header = TRUE)
print(paste("Successfully read:", full_path))

# ==============================================================================
# 2. HELPER: Calculate Theoretical Bias
# ==============================================================================
get_theoretical_bias <- function(vary_param, min_val, max_val, fixed_opts) {
  
  # Create a sequence for smooth lines
  # Handle special case where min=max (single point) to avoid errors
  if(min_val == max_val) {
    x_seq <- min_val 
  } else {
    x_seq <- seq(min_val, max_val, length.out = 300)
  }
  
  # Expand grid for methods
  base_grid <- expand.grid(x = x_seq, method = c("normal", "perm", "bjk"))
  
  # ----------------------------------------------------------------------------
  # FORMULAS
  # Normal: - [1 / (n-1)] * [(1 + pq) / (2pq)] * (tau^2 + sigma1)
  # Perm:   [1 / n] * [(1 + 3pq) / (2pq)] * beta^2
  # BJK:    o(1/dr) + o(1/n^2) -> Treated as 0 for theoretical plotting
  # ----------------------------------------------------------------------------
  
  # We define a localized helper to compute bias given current parameters
  calc_bias <- function(curr_n, curr_p, curr_beta, curr_tau, curr_sigma1, method_name) {
    q <- 1 - curr_p
    
    if (method_name == "normal") {
      term1 <- -1 / (curr_n - 1)
      term2 <- (1 + curr_p * q) / (2 * curr_p * q)
      term3 <- (curr_tau^2 + curr_sigma1)
      return(term1 * term2 * term3)
    } 
    else if (method_name == "perm") {
      term1 <- 1 / curr_n
      term2 <- (1 + 3 * curr_p * q) / (2 * curr_p * q)
      term3 <- curr_beta^2
      return(term1 * term2 * term3)
    } 
    else if (method_name == "bjk") {
      # The bias is negligible order terms o(...)
      return(0)
    }
    return(0)
  }
  
  # Apply logic based on what is varying
  # Extract fixed defaults first
  N_def <- fixed_opts$n_fixed
  B_def <- fixed_opts$beta_fixed
  P_def <- fixed_opts$maf_fixed
  T_def <- fixed_opts$tau_fixed
  S_def <- fixed_opts$sigma1_fixed
  
  out <- base_grid %>% rowwise() %>% mutate(
    theo = case_when(
      
      # --- SCENARIO A: VARY MAF ---
      vary_param == "maf" ~ calc_bias(N_def, x, B_def, T_def, S_def, method),
      
      # --- SCENARIO B: VARY BETA ---
      vary_param == "beta" ~ calc_bias(N_def, P_def, x, T_def, S_def, method),
      
      # --- SCENARIO C: VARY N_FAMILIES ---
      vary_param == "n_families" ~ calc_bias(x, P_def, B_def, T_def, S_def, method),
      
      # --- SCENARIO D: VARY SIGMA1 ---
      vary_param == "sigma1" ~ calc_bias(N_def, P_def, B_def, T_def, x, method),
      
      # --- SCENARIO E: VARY TAU ---
      vary_param == "tau" ~ calc_bias(N_def, P_def, B_def, x, S_def, method),
      
      TRUE ~ 0
    )
  ) %>% ungroup()
  
  # Align names for plotting
  out$method <- recode(out$method, 
                       "normal" = "Normal Theory", 
                       "perm"   = "Permutation", 
                       "bjk"    = "Block Jackknife")
  return(out)
}

# ==============================================================================
# 3. PREPARE DATA
# ==============================================================================

# A. Summarize Simulation Data
#    Dynamically group by the input argument
bias_summary <- df %>%
  pivot_longer(cols = ends_with("bias"), names_to = "method", values_to = "bias") %>%
  mutate(method = recode(method, 
                         normal_bias = "Normal Theory", 
                         perm_bias = "Permutation", 
                         bjk_bias = "Block Jackknife")) %>%
  group_by(method, !!sym(opt$vary)) %>%
  summarise(
    mean_bias = mean(bias, na.rm=TRUE),
    se_bias   = sd(bias, na.rm=TRUE) / sqrt(n()),
    ci_lower  = mean_bias - 1.96 * se_bias,
    ci_upper  = mean_bias + 1.96 * se_bias,
    .groups   = "drop"
  )

# B. Generate Theoretical Lines
min_x <- min(bias_summary[[opt$vary]])
max_x <- max(bias_summary[[opt$vary]])
theory_df <- get_theoretical_bias(opt$vary, min_x, max_x, opt)

# ==============================================================================
# 4. PLOTTING
# ==============================================================================
method_cols <- c("Block Jackknife"="#D81B60", "Normal Theory"="#1E88E5", "Permutation"="#f8971f")

labels_map <- list(
  maf = "Minor Allele Frequency (MAF)",
  beta = "Allelic effect size",
  n_families = "Sample Size (Families)",
  sigma1 = "Heteroskedasticity (Sigma 1)",
  tau = "Between-Family SD (Tau)"
)

# Use labels_map if key exists, else default to raw string
x_label <- if(!is.null(labels_map[[opt$vary]])) labels_map[[opt$vary]] else opt$vary

p <- ggplot(bias_summary, aes(x = .data[[opt$vary]], color = method)) +
  
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  # Theoretical lines
  geom_line(data = theory_df, aes(x = x, y = theo), linewidth = 1, alpha = 0.8) +
  
  # Simulation points
  geom_pointrange(aes(y = mean_bias, ymin = ci_lower, ymax = ci_upper), 
                  size = 0.5, fatten = 2) +
  
  # Styling
  scale_color_manual(values = method_cols) +
  labs(
    x = x_label,
    y = "Bias in Variance Estimate",
    title = paste("Bias vs.", x_label)
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_text(face="bold", hjust=0.5)
  )

# Save
output_file_name = sprintf("plots/comparison_plot.%s_vary.beta_%s.n_%d.tau_%s.sigma1_%s.pdf",
                  opt$vary,                # e.g., "maf"
                  as.character(opt$beta_fixed),   # e.g., "1"
                  opt$n_fixed,             # e.g., 5000
                  as.character(opt$tau_fixed),    # e.g., "0"
                  as.character(opt$sigma1_fixed)) # e.g., "0"
ggsave(output_file_name, plot = p, width = 7, height = 5)
print(paste("Plot saved to",output_file_name))
