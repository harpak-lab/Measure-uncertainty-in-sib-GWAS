trait="height" # example trait
comparison_df = read.table(paste0("results/", trait, ".normal.bjk.perm_1k.txt"), header = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)
make_ratio_df <- function(comparison_df) {
  comparison_df %>%
    transmute(
      SNP, EMP1, bjk_var,
      perm_var , MAF,normal_var
    ) %>%
    mutate(
      ratio_perm   = perm_var   / bjk_var,
      ratio_normal = normal_var / bjk_var
    ) %>%
    pivot_longer(
      cols = starts_with("ratio_"),
      names_to  = "method",
      values_to = "ratio"
    ) %>%
    mutate(
      method = recode(method,
                      ratio_perm   = "Permutation",
                      ratio_normal = "Normal"
                      )
    ) %>%
    filter(is.finite(ratio), is.finite(EMP1))
}

ratio_df <- make_ratio_df(comparison_df)

# --- parameters ---
n_per_bin <- 300L
col_perm  <- "#FDB813"
col_norm  <- "#1E88E5"
line_perm <- "#B34700" # dark orange-brown 
line_norm <- "#0D47A1" # dark navy blue

# # --- prepare data ---
d0 <- ratio_df %>%
  filter(method %in% c("Permutation","Normal"),
         is.finite(EMP1), EMP1 > 0,
         is.finite(ratio), ratio > 0) %>%
  mutate(
    neglog10_p = -log10(EMP1),
    log2_ratio = log2(ratio),
    log10_maf = log10(MAF)
  )

# --- equal-count bins per method ---
# precompute on d0
d1 <- d0 %>%
  group_by(method) %>%
  arrange(neglog10_p, .by_group = TRUE) %>%        # sort by -log10(p)
  mutate(bin_id = (row_number() - 1L) %/% n_per_bin + 1L) %>%
  ungroup()

bin_means <- d1 %>%
  group_by(method, bin_id) %>%
  summarise(
    p_mean   = mean(neglog10_p, na.rm = TRUE),     # x: average -log10(p)
    log2_mean = mean(log2_ratio, na.rm = TRUE),    # y: mean log2 ratio (geometric mean)
    n_bin    = n(),
    .groups  = "drop"
  )

# Shared limits from ALL binned means 
get_shared_limits_from_bins <- function(df_bins, pad_frac = 0.03, keep_zero = TRUE) {
  yb <- df_bins$log2_mean[is.finite(df_bins$log2_mean)]
  if (!length(yb)) return(c(-0.5, 0.5))
  r   <- range(yb)
  pad <- pad_frac * diff(r)
  lim <- c(r[1] - pad, r[2] + pad)
  if (keep_zero) lim <- range(c(lim, 0))  # ensure ratio=1 is visible
  lim
}
make_panel_binned_pts_raw_spline <- function(
    df_raw, df_bins, method,
    pt_col, line_col,
    k_spline = 10, sample_n = NULL,
    y_limits_shared = NULL,
    fit_on = c("bins", "raw"),
    x_pad = 0.25
) {
  fit_on <- match.arg(fit_on)
  library(dplyr); library(ggplot2); library(mgcv); library(scales); library(grid)
  
  paper_theme <- theme_classic(base_size = 13) +
    theme(
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6))
    )
  
  # --- data subsets ---
  raw_m  <- df_raw %>%
    filter(method == !!method,
           is.finite(neglog10_p), is.finite(log2_ratio))
  if (!is.null(sample_n) && nrow(raw_m) > sample_n) {
    raw_m <- slice_sample(raw_m, n = sample_n)
  }
  bins_m <- df_bins %>%
    filter(method == !!method,
           is.finite(p_mean), is.finite(log2_mean))
  
  # --- shared Y limits ---
  if (is.null(y_limits_shared)) {
    yb <- df_bins$log2_mean[is.finite(df_bins$log2_mean)]
    r   <- range(yb)
    pad <- 0.03 * diff(r)
    y_limits_shared <- range(c(r[1] - pad, r[2] + pad, 0))
  }
  
  # Smallest p-value across raw & binned data
  p_min_raw  <- min(10^(-df_raw$neglog10_p), na.rm = TRUE)
  p_min_bins <- min(10^(-df_bins$p_mean),     na.rm = TRUE)
  p_min      <- min(p_min_raw, p_min_bins)
  
  # Convert to -log10 scale for plotting
  x_max_auto <- ceiling(-log10(p_min))
  
  # Final x-axis limits
  x_limits <- c(0, x_max_auto)
  
  # --- automatic x-axis limit from data ---
  p_min_raw  <- min(10^(-df_raw$neglog10_p), na.rm = TRUE)
  p_min_bins <- min(10^(-df_bins$p_mean),     na.rm = TRUE)
  p_min      <- min(p_min_raw, p_min_bins)
  
  # convert to -log10 scale
  x_max_auto <- ceiling(-log10(p_min))
  x_limits <- c(0, x_max_auto)
  
  # --- pretty labels on linear ratio scale (y) ---
  linear_limits <- 2 ^ y_limits_shared
  linear_breaks <- pretty_breaks(6)(linear_limits)
  linear_breaks <- linear_breaks[linear_breaks >= min(linear_limits) &
                                   linear_breaks <= max(linear_limits)]
  y_breaks <- log2(linear_breaks)
  y_labels <- formatC(linear_breaks, format = "f", digits = 2)
  
  # --- choose data for spline ---
  smooth_df <- if (fit_on == "bins") bins_m else raw_m
  aes_x <- if (fit_on == "bins") aes(x = p_mean, y = log2_mean) else aes(x = neglog10_p, y = log2_ratio)
  
  # --- log ticks positions for x ---
  # Map true p-values 10^(-x) to evenly spaced ticks in -log10(p)
  p_breaks <- 10^(0:-6)
  x_tick_positions <- -log10(p_breaks)     # e.g. 0, 1, 2, 3, 4, 5, 6
  x_tick_minor     <- seq(0, max(x_tick_positions), by = 0.5)  # short minor ticks
  
  ggplot() +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.6, color = "grey45") +
    geom_point(data = bins_m,
               aes(x = p_mean, y = log2_mean),
               shape = 21, size = 2.6, stroke = 0.5,
               color = pt_col, fill = alpha(pt_col, 0.35)) +
    geom_smooth(data = smooth_df, mapping = aes_x,
                method = "gam", formula = y ~ s(x, bs = "cs", k = k_spline),
                se = FALSE, linewidth = 1.2, color = line_col) +
    annotation_logticks(sides = "l", short = unit(0.15, "cm"), mid = unit(0.25, "cm"), long = unit(0.3, "cm")) +
    annotation_logticks(
      sides = "b", 
      short = unit(0.10, "cm"),
      mid   = unit(0.18, "cm"),
      long  = unit(0.25, "cm"),
      scaled = TRUE
    ) +
    
    scale_x_continuous(
      name   = expression("Permutation " * p * "-value"),
      limits = x_limits,
      breaks = x_tick_positions,
      labels = function(b) paste0("1e-", b),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      name   = bquote(widehat(Var)[.(method)] / widehat(Var)[BJK]),
      breaks = y_breaks,
      labels = y_labels
    ) +
    coord_cartesian(ylim = y_limits_shared) +
    paper_theme
}

col_perm  <- "#FDB813"   # point color (Permutation)
col_norm  <- "#1E88E5"   # point color (Normal)
line_perm <- "#B34700"   # line color (Permutation)
line_norm <- "#0D47A1"   # line color (Normal)

ylims <- get_shared_limits_from_bins(bin_means, pad_frac = 0.03, keep_zero = TRUE)

p_norm <- make_panel_binned_pts_raw_spline(
  d0, bin_means, "Normal",
  pt_col = col_norm, line_col = line_norm,
  y_limits_shared = ylims, fit_on = "raw", x_pad = 0.4
)
p_perm <- make_panel_binned_pts_raw_spline(
  d0, bin_means, "Permutation",
  pt_col = col_perm, line_col = line_perm,
  y_limits_shared = ylims, fit_on = "raw", x_pad = 0.4
)

panels  <- plot_grid(p_norm, p_perm, ncol = 2, align = "hv", rel_widths = c(1, 1))
title_g <- ggdraw() + draw_label(trait, x = 0, hjust = 0, fontface = "bold", size = 16)
p_final <- plot_grid(title_g, panels, ncol = 1, rel_heights = c(0.10, 1))
ggsave(paste0(trait , ".binned.ratio_perm_norm_bjk.pdf"), p_final, width = 10, height = 4, dpi = 300)


d0p <- d0 %>%
  mutate(
    maf_lin   = ifelse(is.finite(MAF), MAF, 10^(log10_maf)),
    maf_log10 = ifelse(is.finite(log10_maf), log10_maf, log10(MAF)),
    p_stratum = ifelse(EMP1 < 0.05, "p < 0.05", "p >= 0.05")
  ) %>%
  filter(is.finite(maf_lin), maf_lin > 0, is.finite(ratio), ratio > 0) %>%
  mutate(p_stratum = forcats::fct_relevel(p_stratum, "p >= 0.05", "p < 0.05"))

# bin means for ALL methods (so limits are shared)
bin_means_all <- d0p %>%
  group_by(method, p_stratum) %>%
  arrange(log10_maf, .by_group = TRUE) %>%
  mutate(bin_id = (row_number() - 1L) %/% 300L + 1L) %>%
  group_by(method, p_stratum, bin_id) %>%
  summarise(
    log10_maf_mean = mean(log10_maf, na.rm = TRUE),
    log2_mean      = mean(log2(ratio), na.rm = TRUE),
    .groups = "drop"
  )

get_shared_limits_from_bins <- function(df_bins, pad_frac = 0.03, keep_zero = TRUE) {
  yb <- df_bins$log2_mean[is.finite(df_bins$log2_mean)]
  if (!length(yb)) return(c(-0.5, 0.5))
  r   <- range(yb)
  pad <- pad_frac * diff(r)
  lim <- c(r[1] - pad, r[2] + pad)
  if (keep_zero) lim <- range(c(lim, 0))  # keep ratio=1 visible (log2=0)
  lim
}

ylims_log2 <- get_shared_limits_from_bins(bin_means_all, pad_frac = 0.03, keep_zero = TRUE)
make_maf_scatter_rawspline_fixedbins <- function(df, method_name, point_col, line_col,
                                                 bin_size = 300L,
                                                 share_y_limits_log2 = NULL,
                                                 k_spline = 12,
                                                 smooth_binned = FALSE) {
  # --- raw data (for spline) ---
  raw_m <- df %>%
    dplyr::filter(method == method_name,
                  is.finite(maf_lin), maf_lin > 0,
                  is.finite(log10_maf),
                  is.finite(ratio), ratio > 0) %>%
    dplyr::mutate(
      log2_ratio = if (!"log2_ratio" %in% names(.)) log2(ratio) else log2_ratio
    )
  
  # --- bins (still defined in log10 space) ---
  bins <- df %>%
    dplyr::filter(method == method_name,
                  is.finite(log10_maf),
                  is.finite(ratio), ratio > 0) %>%
    dplyr::group_by(p_stratum) %>%
    dplyr::arrange(log10_maf, .by_group = TRUE) %>%
    dplyr::mutate(bin_id = (dplyr::row_number() - 1L) %/% bin_size + 1L) %>%
    dplyr::group_by(p_stratum, bin_id) %>%
    dplyr::summarise(
      n_bin          = dplyr::n(),
      log10_maf_mean = mean(log10_maf, na.rm = TRUE),
      maf_mid        = 10^log10_maf_mean,              # back to linear MAF for plotting
      log2_mean      = mean(log2(ratio),  na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- y-limits: use shared if provided; otherwise derive from bins ---
  if (is.null(share_y_limits_log2)) {
    yb <- bins$log2_mean[is.finite(bins$log2_mean)]
    if (!length(yb)) {
      y_limits_log2 <- c(-0.5, 0.5)
    } else {
      r   <- range(yb)
      pad <- 0.03 * diff(r)
      y_limits_log2 <- range(c(r[1] - pad, r[2] + pad, 0))  # keep 0 (ratio=1) visible
    }
  } else {
    y_limits_log2 <- share_y_limits_log2
  }
  
  # --- derive nice ratio ticks from those limits ---
  linear_limits <- 2 ^ y_limits_log2
  linear_breaks <- scales::pretty_breaks(5)(linear_limits)
  # keep only within range
  linear_breaks <- linear_breaks[linear_breaks >= min(linear_limits) &
                                   linear_breaks <= max(linear_limits)]
  y_breaks <- log2(linear_breaks)
  y_labels <- formatC(linear_breaks, format = "f", digits = 2)
  
  # --- x ticks (MAF) ---
  maf_ticks_lin <- c(1e-4, 1e-3, 1e-2, 1e-1, 0.5)
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.6, color = "grey45") +
    
    # binned means
    geom_point(
      data  = bins,
      aes(x = maf_mid, y = log2_mean),
      shape  = 21,
      size   = 2.0,
      stroke = 0.45,
      fill   = alpha(point_col, 0.35),
      color  = point_col
    )
  
  # spline (on bins or raw)
  if (smooth_binned) {
    p <- p + geom_smooth(
      data   = bins,
      aes(x = maf_mid, y = log2_mean),
      method = "gam",
      formula = y ~ s(x, bs = "cs", k = k_spline),
      se = FALSE, linewidth = 1.2, color = line_col
    )
  } else {
    p <- p + geom_smooth(
      data   = raw_m,
      aes(x = maf_lin, y = log2_ratio, group = p_stratum),
      method = "gam",
      formula = y ~ s(x, bs = "cs", k = k_spline),
      se = FALSE, linewidth = 1.2, color = line_col
    )
  }
  
  p +
    # X: log10 MAF with labels in scientific notation + log ticks
    scale_x_log10(
      name   = "Minor allele frequency",
      breaks = maf_ticks_lin,
      labels = scales::label_scientific(digits = 1)
    ) +
    annotation_logticks(
      sides = "b",
      short = grid::unit(0.05, "lines"),
      mid   = grid::unit(0.08, "lines"),
      long  = grid::unit(0.12, "lines")
    ) +
    
    # Y: log2 scale with labels in linear ratio space
    scale_y_continuous(
      name   = bquote(widehat(Var)[.(method_name)]/widehat(Var)[BJK]),
      breaks = y_breaks,
      labels = y_labels
    ) +
    annotation_logticks(sides = "l", short = unit(0.15, "cm"), mid = unit(0.25, "cm"), long = unit(0.3, "cm")) +
    facet_wrap(
      ~ p_stratum, nrow = 1, scales = "fixed",
      labeller = labeller(p_stratum = label_parsed)
    ) +
    coord_cartesian(ylim = y_limits_log2) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(face = "bold"),
      axis.title.y     = element_text(margin = margin(r = 6)),
      axis.title.x     = element_text(margin = margin(t = 6)),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    labs(title = method_name)
}


# ---- calls (now shared_lin_limits exists in the caller) ----
ylims_log2 <- get_shared_limits_from_bins(bin_means_all, pad_frac = 0.03, keep_zero = TRUE)

p_norm <- make_maf_scatter_rawspline_fixedbins(
  d0p, "Normal", col_norm, line_norm,
  bin_size = 300L,
  share_y_limits_log2 = ylims_log2,
  k_spline = 12,
  smooth_binned = TRUE
)

p_perm <- make_maf_scatter_rawspline_fixedbins(
  d0p, "Permutation", col_perm, line_perm,
  bin_size = 300L,
  share_y_limits_log2 = ylims_log2,
  k_spline = 12,
  smooth_binned = TRUE
)

# combine & save
row_panels <- cowplot::plot_grid(p_norm, p_perm, ncol = 2, align = "hv", rel_widths = c(1, 1))
title_g    <- cowplot::ggdraw() + cowplot::draw_label(trait, x = 0, hjust = 0, fontface = "bold", size = 16)
p_final    <- cowplot::plot_grid(title_g, row_panels, ncol = 1, rel_heights = c(0.10, 1))

ggsave(paste0(trait, ".binned_maf.ratio_perm_norm_bjk.pdf"),
       p_final, width = 10, height = 4, dpi = 300)
