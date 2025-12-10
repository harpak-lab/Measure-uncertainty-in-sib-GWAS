normal = read.table("vary_beta/normal.n_5000.beta_-0.2_0.5.tau_0.maf_0.1.sigma_1_0.sim_5k.txt", header = TRUE)
perm = read.table("vary_beta/permutation.n_5000.beta_-0.2_0.2.tau_0.maf_0.1.sigma_1_0.sim_5k.txt", header = TRUE)
bjk = read.table("vary_beta/jackknife.n_5000.beta_-0.2_0.5.tau_0.maf_0.1.sigma_1_0.sim_5k.txt", header = TRUE)

comparison_df = cbind(normal[, c("se", "beta", "normal_bias")], perm[, c("perm_var", "beta", "perm_bias")])
comparison_df = cbind(comparison_df, bjk[, c("bjk_var", "beta", "bjk_bias")])
comparison_df <- comparison_df[, !duplicated(names(comparison_df))]
library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Long format: one row per (replicate, sigma_1, method)
bias_long <- comparison_df %>%
  select(beta,
         normal_bias,
         perm_bias,
         bjk_bias) %>%
  pivot_longer(
    cols = ends_with("bias"),
    names_to  = "method",
    values_to = "bias"
  ) %>%
  mutate(
    method = recode(method,
                    normal_bias = "normal",
                    perm_bias   = "perm",
                    bjk_bias    = "bjk")
  )

# 2) Summarise: mean bias and 95% CI by method & sigma_1
bias_summary <- bias_long %>%
  group_by(method, beta) %>%
  summarise(
    mean_bias = mean(bias),
    se_bias   = sd(bias) / sqrt(n()),
    ci_lower  = mean_bias - 1.96 * se_bias,
    ci_upper  = mean_bias + 1.96 * se_bias,
    .groups   = "drop"
  )

# 3) Colors
method_cols <- c(
  bjk    = "#D81B60",
  normal = "#1E88E5",
  perm   = "#F89620"
)
library(ggplot2)

p = ggplot(bias_summary,
           aes(x = beta, y = mean_bias, color = method)) +
  
  # zero bias line
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, color = "gray40") +
  
  # CI
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                  size = 0.5, fatten = 2) +
  
  # mean trend line
  geom_line(size = 1.2) +
  
  scale_color_manual(values = method_cols, name = "Method") +
  
  labs(
    x = "Allelic effect",
    y = "Bias in Variance Estimate",
    title = "Bias vs. allelic effect"
  ) +
  
  # Scientific-style theme
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    
    # Remove grid (classic already removes major/minor grid)
    panel.grid = element_blank(),
    
    # Axes styling
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 0.9),
    
    # spacing
    plot.margin = margin(10, 15, 10, 15)
  )
library(dplyr)
library(tidyr)
library(ggplot2)
normal = read.table("vary_sigma/normal_method.n_5000.beta_0.tau_0.maf_0.1.sigma_1_-0.5_0.5.sim_5k.txt", header = TRUE)
perm = read.table("vary_sigma/permutation.n_5000.beta_0.tau_0.maf_0.1.sigma_1_-0.5_0.5.sim_5k.txt", header = TRUE)
bjk = read.table("vary_sigma/jackknife.n_5000.beta_0.tau_0.maf_0.1.sigma_1_-0.5_0.5.sim_5k.txt", header = TRUE)

comparison_df = cbind(normal[, c("se", "sigma_1", "normal_bias")], perm[, c("perm_var", "sigma_1", "perm_bias")])
comparison_df = cbind(comparison_df, bjk[, c("bjk_var", "sigma_1", "bjk_bias")])
comparison_df <- comparison_df[, !duplicated(names(comparison_df))]

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Long format: one row per (replicate, sigma_1, method)
bias_long <- comparison_df %>%
  select(sigma_1,
         normal_bias,
         perm_bias,
         bjk_bias) %>%
  pivot_longer(
    cols = ends_with("bias"),
    names_to  = "method",
    values_to = "bias"
  ) %>%
  mutate(
    method = recode(method,
                    normal_bias = "normal",
                    perm_bias   = "perm",
                    bjk_bias    = "bjk")
  )

# 2) Summarise: mean bias and 95% CI by method & sigma_1
bias_summary <- bias_long %>%
  group_by(method, sigma_1) %>%
  summarise(
    mean_bias = mean(bias),
    se_bias   = sd(bias) / sqrt(n()),
    ci_lower  = mean_bias - 1.96 * se_bias,
    ci_upper  = mean_bias + 1.96 * se_bias,
    .groups   = "drop"
  )

# 3) Colors
method_cols <- c(
  bjk    = "#D81B60",
  normal = "#1E88E5",
  perm   = "#F89620"
)

# 4) Plot
library(ggplot2)

p = ggplot(bias_summary,
       aes(x = sigma_1, y = mean_bias, color = method)) +
  
  # zero bias line
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, color = "gray40") +
  
  # CI
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                  size = 0.5, fatten = 2) +
  
  # mean trend line
  geom_line(size = 1.2) +
  
  scale_color_manual(values = method_cols, name = "Method") +
  
  labs(
    x = expression(sigma[1]),
    y = "Bias in Variance Estimate",
    title = "Bias vs. Non–focal Variance Parameter"
  ) +
  
  # Scientific-style theme
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    
    # Remove grid (classic already removes major/minor grid)
    panel.grid = element_blank(),
    
    # Axes styling
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 0.9),
    
    # spacing
    plot.margin = margin(10, 15, 10, 15)
  )


normal = read.table("vary_n_families/normal.n_1k_10k.beta_0.tau_0.maf_0.1.sigma_1_-0.5sim_5k.txt", header = TRUE)
perm = read.table("vary_n_families/perm.n_1k_10k.beta_0.tau_0.maf_0.1.sigma_1_0.5sim_5k.txt", header = TRUE)
bjk = read.table("vary_n_families/jackknife.n_1k_10k.beta_0.tau_0.maf_0.1.sigma_1_0.5sim_5k.txt", header = TRUE)

comparison_df = cbind(normal[, c("se", "n_families", "normal_bias")], perm[, c("perm_var", "n_families", "perm_bias")])
comparison_df = cbind(comparison_df, bjk[, c("bjk_var", "n_families", "bjk_bias")])
comparison_df <- comparison_df[, !duplicated(names(comparison_df))]

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Long format: one row per (replicate, sigma_1, method)
bias_long <- comparison_df %>%
  select(n_families,
         normal_bias,
         perm_bias,
         bjk_bias) %>%
  pivot_longer(
    cols = ends_with("bias"),
    names_to  = "method",
    values_to = "bias"
  ) %>%
  mutate(
    method = recode(method,
                    normal_bias = "normal",
                    perm_bias   = "perm",
                    bjk_bias    = "bjk")
  )

# 2) Summarise: mean bias and 95% CI by method & sigma_1
bias_summary <- bias_long %>%
  group_by(method, n_families) %>%
  summarise(
    mean_bias = mean(bias),
    se_bias   = sd(bias) / sqrt(n()),
    ci_lower  = mean_bias - 1.96 * se_bias,
    ci_upper  = mean_bias + 1.96 * se_bias,
    .groups   = "drop"
  )

# 3) Colors
method_cols <- c(
  bjk    = "#D81B60",
  normal = "#1E88E5",
  perm   = "#F89620"
)

# 4) Plot
library(ggplot2)

p = ggplot(bias_summary,
           aes(x = n_families, y = mean_bias, color = method)) +
  
  # zero bias line
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, color = "gray40") +
  
  # CI
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                  size = 0.5, fatten = 2) +
  
  # mean trend line
  geom_line(size = 1.2) +
  
  scale_color_manual(values = method_cols, name = "Method") +
  
  labs(
    x = expression(sigma[1]),
    y = "Bias in Variance Estimate",
    title = "Bias vs.number of sibling paris"
  ) +
  
  # Scientific-style theme
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    
    # Remove grid (classic already removes major/minor grid)
    panel.grid = element_blank(),
    
    # Axes styling
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 0.9),
    
    # spacing
    plot.margin = margin(10, 15, 10, 15)
  )

normal = read.table("vary_n_families/normal.n_1k_10k.beta_1.tau_0.maf_0.1.sigma_1_0sim_5k.txt", header = TRUE)
perm = read.table("vary_n_families/permutation.n_1k_10k.beta_1.tau_0.maf_0.1.sigma_1_0sim_5k.txt", header = TRUE)
bjk = read.table("vary_n_families/jackknife.n_1k_10k.beta_1.tau_0.maf_0.1.sigma_1_0.sim_5k.txt", header = TRUE)

comparison_df = cbind(normal[, c("se", "n_families", "normal_bias")], perm[, c("perm_var", "n_families", "perm_bias")])
comparison_df = cbind(comparison_df, bjk[, c("bjk_var", "n_families", "bjk_bias")])
comparison_df <- comparison_df[, !duplicated(names(comparison_df))]

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Long format: one row per (replicate, sigma_1, method)
bias_long <- comparison_df %>%
  select(n_families,
         normal_bias,
         perm_bias,
         bjk_bias) %>%
  pivot_longer(
    cols = ends_with("bias"),
    names_to  = "method",
    values_to = "bias"
  ) %>%
  mutate(
    method = recode(method,
                    normal_bias = "normal",
                    perm_bias   = "perm",
                    bjk_bias    = "bjk")
  )

# 2) Summarise: mean bias and 95% CI by method & sigma_1
bias_summary <- bias_long %>%
  group_by(method, n_families) %>%
  summarise(
    mean_bias = mean(bias),
    se_bias   = sd(bias) / sqrt(n()),
    ci_lower  = mean_bias - 1.96 * se_bias,
    ci_upper  = mean_bias + 1.96 * se_bias,
    .groups   = "drop"
  )

# 3) Colors
method_cols <- c(
  bjk    = "#D81B60",
  normal = "#1E88E5",
  perm   = "#F89620"
)

# 4) Plot
library(ggplot2)

p = ggplot(bias_summary,
           aes(x = n_families, y = mean_bias, color = method)) +
  
  # zero bias line
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, color = "gray40") +
  
  # CI
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                  size = 0.5, fatten = 2) +
  
  # mean trend line
  geom_line(size = 1.2) +
  
  scale_color_manual(values = method_cols, name = "Method") +
  
  labs(
    x = expression(sigma[1]),
    y = "Bias in Variance Estimate",
    title = "Bias vs.number of sibling paris"
  ) +
  
  # Scientific-style theme
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    
    # Remove grid (classic already removes major/minor grid)
    panel.grid = element_blank(),
    
    # Axes styling
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 0.9),
    
    # spacing
    plot.margin = margin(10, 15, 10, 15)
  )

normal = read.table("vary_maf/normal.n_5k.beta_0.tau_0.maf_0.05_0.5.sigma_1_0.5.txt", header = TRUE)
perm = read.table("vary_maf/permutation.n_5k.beta_0.tau_0.maf_0.05_0.5.sigma_1_0.5.sim_5k.txt", header = TRUE) 
bjk = read.table("vary_maf/jackknife.n_5000.beta_0.tau_0.maf_0.05_0.5.sigma_1_0.5.sim_5k.txt", header = TRUE)


comparison_df = cbind(normal[, c("se", "maf", "normal_bias")], perm[, c("perm_var", "maf", "perm_bias")])
comparison_df = cbind(comparison_df, bjk[, c("bjk_var", "maf", "bjk_bias")])
comparison_df <- comparison_df[, !duplicated(names(comparison_df))]

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Long format: one row per (replicate, sigma_1, method)
bias_long <- comparison_df %>%
  select(maf,
         normal_bias,
         perm_bias,
         bjk_bias) %>%
  pivot_longer(
    cols = ends_with("bias"),
    names_to  = "method",
    values_to = "bias"
  ) %>%
  mutate(
    method = recode(method,
                    normal_bias = "normal",
                    perm_bias   = "perm",
                    bjk_bias    = "bjk")
  )

# 2) Summarise: mean bias and 95% CI by method & sigma_1
bias_summary <- bias_long %>%
  group_by(method, maf) %>%
  summarise(
    mean_bias = mean(bias),
    se_bias   = sd(bias) / sqrt(n()),
    ci_lower  = mean_bias - 1.96 * se_bias,
    ci_upper  = mean_bias + 1.96 * se_bias,
    .groups   = "drop"
  )

# 3) Colors
method_cols <- c(
  bjk    = "#D81B60",
  normal = "#1E88E5",
  perm   = "#F89620"
)

# 4) Plot
library(ggplot2)

p = ggplot(bias_summary,
           aes(x = maf, y = mean_bias, color = method)) +
  
  # zero bias line
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, color = "gray40") +
  
  # CI
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                  size = 0.5, fatten = 2) +
  
  # mean trend line
  geom_line(size = 1.2) +
  
  scale_color_manual(values = method_cols, name = "Method") +
  
  labs(
    x = "minor allele frequency",
    y = "Bias in Variance Estimate",
    title = "Bias vs.minor allele frequency"
  ) +
  
  # Scientific-style theme
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    
    # Remove grid (classic already removes major/minor grid)
    panel.grid = element_blank(),
    
    # Axes styling
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 0.9),
    
    # spacing
    plot.margin = margin(10, 15, 10, 15)
  )

comparison_df = cbind(normal[, c("se", "tau", "normal_bias")], perm[, c("perm_var", "tau", "perm_bias")])
comparison_df = cbind(comparison_df, bjk[, c("bjk_var", "tau", "bjk_bias")])
comparison_df <- comparison_df[, !duplicated(names(comparison_df))]

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Long format: one row per (replicate, sigma_1, method)
bias_long <- comparison_df %>%
  select(tau,
         normal_bias,
         perm_bias,
         bjk_bias) %>%
  pivot_longer(
    cols = ends_with("bias"),
    names_to  = "method",
    values_to = "bias"
  ) %>%
  mutate(
    method = recode(method,
                    normal_bias = "normal",
                    perm_bias   = "perm",
                    bjk_bias    = "bjk")
  )

# 2) Summarise: mean bias and 95% CI by method & sigma_1
bias_summary <- bias_long %>%
  group_by(method, tau) %>%
  summarise(
    mean_bias = mean(bias),
    se_bias   = sd(bias) / sqrt(n()),
    ci_lower  = mean_bias - 1.96 * se_bias,
    ci_upper  = mean_bias + 1.96 * se_bias,
    .groups   = "drop"
  )

# 3) Colors
method_cols <- c(
  bjk    = "#D81B60",
  normal = "#1E88E5",
  perm   = "#F89620"
)

# 4) Plot
library(ggplot2)

p = ggplot(bias_summary,
           aes(x = sigma_1, y = mean_bias, color = method)) +
  
  # zero bias line
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.6, color = "gray40") +
  
  # CI
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                  size = 0.5, fatten = 2) +
  
  # mean trend line
  geom_line(size = 1.2) +
  
  scale_color_manual(values = method_cols, name = "Method") +
  
  labs(
    x = expression(sigma[1]),
    y = "Bias in Variance Estimate",
    title = "Bias vs. Non–focal Variance Parameter"
  ) +
  
  # Scientific-style theme
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    
    # Remove grid (classic already removes major/minor grid)
    panel.grid = element_blank(),
    
    # Axes styling
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 0.9),
    
    # spacing
    plot.margin = margin(10, 15, 10, 15)
  )
