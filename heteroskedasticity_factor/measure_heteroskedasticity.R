## Estimate within-family residual (non-focal) variance as a function of the focal genotype contrast DeltaX^2, for each SNP.
##
## This script produces two outputs:
##   1. Per-SNP residual variance grouped by DeltaX^2 in {0, 1, 4}
##   2. Per-SNP linear fit: Var(residual) = sigma0 + sigma1 * DeltaX^2
##
## The estimated slope sigma1 quantifies heteroskedasticity in the
## non-focal variance (see Eq. 4 in the main text).
##
## Usage:
##   Rscript estimate_residual_variance.R <trait>
##
## NOTE: This script requires access to individual-level UK Biobank
## data, which cannot be shared publicly. Researchers can apply for
## access at https://www.ukbiobank.ac.uk/.
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(data.table)
  library(broom)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript estimate_residual_variance.R <trait>")
trait <- args[1]

## ---- Configuration: update these paths for your environment ----
BASE_DIR        <- getwd()
COMPARISON_DIR  <- file.path(BASE_DIR, "data", "comparison_df")
CHUNKS_BASE_DIR <- BASE_DIR

## ============================================================
## PART 1: Compute per-SNP residual variance by DeltaX^2 group
## ============================================================

## Helper: compute residual variance grouped by DeltaX^2 for all SNPs
all_snp_variances_with_counts <- function(df, df_snp,
                                          use_marker_id = FALSE) {
  key_col <- if (use_marker_id) "marker.ID" else "SNP"

  stopifnot("beta" %in% names(df))
  if (!key_col %in% names(df))
    stop(sprintf("df missing '%s'", key_col))
  if (!all(c(key_col, "DeltaX", "DeltaY") %in% names(df_snp)))
    stop(sprintf("df_snp must have '%s', 'DeltaX', 'DeltaY'", key_col))

  df     <- df     %>% mutate(across(all_of(key_col),
                                     ~ trimws(as.character(.))))
  df_snp <- df_snp %>% mutate(across(all_of(key_col),
                                     ~ trimws(as.character(.))))

  snp_index <- df %>%
    distinct(.data[[key_col]]) %>%
    rename(SNP_key = !!key_col)

  ## Join beta to per-family data and compute residuals
  dat <- df_snp %>%
    inner_join(
      df %>% select(all_of(key_col), beta) %>% distinct(),
      by = key_col
    ) %>%
    mutate(
      deltaX2   = DeltaX^2,
      residual  = DeltaY - DeltaX * beta,
      deltaX2_i = as.integer(round(deltaX2))
    )

  ## Summarise by SNP x DeltaX^2 group
  agg <- dat %>%
    group_by(.data[[key_col]], deltaX2_i) %>%
    summarise(
      var_resid = var(residual, na.rm = TRUE),
      n         = dplyr::n(),
      .groups   = "drop"
    )

  if (nrow(agg) == 0) {
    agg <- tibble(!!key_col := character(), deltaX2_i = integer(),
                  var_resid = numeric(), n = integer())
  }

  ## Ensure every SNP has rows for DeltaX^2 in {0, 1, 4}
  agg_completed <- agg %>%
    group_by(.data[[key_col]]) %>%
    tidyr::complete(deltaX2_i = c(0L, 1L, 4L)) %>%
    ungroup()

  ## Pivot wide
  wide <- agg_completed %>%
    tidyr::pivot_wider(
      names_from  = deltaX2_i,
      values_from = c(var_resid, n),
      names_glue  = "{.value}_{deltaX2_i}"
    )

  ## Attach to SNP index so ALL SNPs appear
  out <- snp_index %>%
    left_join(wide %>% rename(SNP_key = !!key_col), by = "SNP_key") %>%
    rename(!!key_col := SNP_key) %>%
    select(all_of(key_col),
           var_resid_0, var_resid_1, var_resid_4,
           n_0, n_1, n_4)

  if (key_col != "SNP") out <- out %>% rename(SNP = all_of(key_col))
  out
}

## Process all chunks for one trait
run_by_chunks <- function(trait,
                          chunk_start = 0L,
                          chunk_end   = 252L,
                          out_dir     = NULL,
                          verbose     = TRUE) {
  if (is.null(out_dir))
    out_dir <- file.path(CHUNKS_BASE_DIR, paste0("chunks_", trait), "out")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  result_list <- vector("list", length = chunk_end - chunk_start + 1L)
  idx <- 0L

  for (chunk_id in seq.int(chunk_start, chunk_end)) {
    idx <- idx + 1L
    pad <- sprintf("%05d", chunk_id)
    chunk_dir <- file.path(
      CHUNKS_BASE_DIR,
      paste0("chunks_", trait), "out", trait,
      paste0(trait, "_chunk_", pad)
    )
    fn_sum <- file.path(chunk_dir,
                        paste0(basename(chunk_dir), "_perSNP_summary.tsv.gz"))
    fn_snp <- file.path(chunk_dir,
                        paste0(basename(chunk_dir), "_perFamily_perSNP.tsv.gz"))

    dt_sum <- tryCatch(
      fread(fn_sum, sep = ",", showProgress = FALSE),
      error = function(e) {
        if (verbose) message("  read summary failed: ", e$message)
        NULL
      }
    )
    dt_snp <- tryCatch(
      fread(fn_snp, sep = ",", showProgress = FALSE),
      error = function(e) {
        if (verbose) message("  read perFamily failed: ", e$message)
        NULL
      }
    )

    if (is.null(dt_sum) || is.null(dt_snp)) {
      result_list[[idx]] <- NULL
      next
    }

    colnames(dt_snp) <- c("CHUNK", "FID", "IID1", "IID2", "pair_id",
                           "marker_ID", "DeltaX", "DeltaY", "SNP")

    need_sum <- c("SNP", "beta")
    need_snp <- c("SNP", "DeltaX", "DeltaY")
    if (!all(need_sum %in% names(dt_sum)) ||
        !all(need_snp %in% names(dt_snp))) {
      if (verbose) message("  Missing required columns, skipping.")
      result_list[[idx]] <- NULL
      next
    }

    chunk_res <- tryCatch(
      all_snp_variances_with_counts(as.data.frame(dt_sum),
                                    as.data.frame(dt_snp),
                                    use_marker_id = FALSE),
      error = function(e) {
        if (verbose) message("  compute failed: ", e$message)
        NULL
      }
    )

    result_list[[idx]] <- chunk_res
    rm(dt_sum, dt_snp, chunk_res); gc()
  }

  out <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
  setkey(out, SNP)

  out_file <- file.path(out_dir,
                        paste0(trait, "_residual_var_by_deltaX2.tsv.gz"))
  fwrite(out, out_file, sep = "\t")
  if (verbose) message("Done. Wrote: ", out_file)

  invisible(out)
}

res_all <- run_by_chunks(trait = trait, chunk_start = 0L, chunk_end = 252L)

## ============================================================
## PART 2: Fit linear model Var(residual) ~ DeltaX^2 per SNP
## ============================================================

## Load the comparison data (for ratio computation)
comparison_df <- read.table(
  file.path(COMPARISON_DIR,
            paste0(trait, ".normal_bjk_perm.maf.adapt_p.txt")),
  header = TRUE
)
comparison_df$normal_var <- comparison_df$se_norm^2
comparison_df$bjk_var    <- comparison_df$bjk_se^2

## Load residual variance data
df_res <- read.table(
  file.path(CHUNKS_BASE_DIR, paste0("chunks_", trait), "out",
            paste0(trait, "_residual_var_by_deltaX2.tsv.gz")),
  header = TRUE, sep = "\t"
)

## Fit Var(residual) = sigma0 + sigma1 * DeltaX^2 for each SNP
## using weighted least squares (weights = number of families in
## each DeltaX^2 group).
fit_var_model <- function(row) {
  df_long <- tibble(
    deltaX2 = c(0, 1, 4),
    var     = c(row$var_resid_0, row$var_resid_1, row$var_resid_4),
    n       = c(row$n_0, row$n_1, row$n_4)
  ) %>%
    filter(!is.na(var) & !is.na(n) & n > 0)

  if (nrow(df_long) < 2) return(NULL)

  fit <- lm(var ~ deltaX2, data = df_long, weights = n)

  tibble(
    SNP        = row$SNP,
    sigma0_hat = coef(fit)[["(Intercept)"]],
    sigma1_hat = coef(fit)[["deltaX2"]],
    r2         = summary(fit)$r.squared,
    p_slope    = summary(fit)$coefficients["deltaX2", "Pr(>|t|)"]
  )
}

df_fit <- df_res %>%
  split(.$SNP) %>%
  map_dfr(~ fit_var_model(.x))

out_fit_file <- file.path(
  CHUNKS_BASE_DIR, paste0("chunks_", trait), "out",
  paste0(trait, ".fit.res_delta_X2.txt")
)
write.table(df_fit, out_fit_file, row.names = FALSE)
message("Wrote fit results: ", out_fit_file)
