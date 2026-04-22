## ============================================================
## compute_normal_var.R
##
## Compute OLS variance estimates for sib-GWAS
## allelic effect estimates, along with per-family per-SNP
## genotype and phenotype differences.
##
##
## Inputs (per trait):
##   - PLINK binary files (.bed/.bim/.fam) split into chunks
##     under <CHUNK_DIR>/<trait>/bfile/
##   - Residualized phenotype file with columns: FID, IID, pheno
##
## Outputs (per trait):
##   - <OUT_DIR>/<trait>_perSNP_summary.tsv.gz
##       Per-SNP: beta, se_norm, sum(DeltaX^2), n_pairs, n_eff
##   - <OUT_DIR>/<trait>_perFamily_perSNP.tsv.gz
##       Per-family per-SNP: DeltaX, DeltaY (for downstream
##       heteroskedasticity analysis)
##
## NOTE: This script requires access to individual-level UK Biobank
## genotype and phenotype data, which cannot be shared publicly.
## Researchers can apply for access at https://www.ukbiobank.ac.uk/.
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(bigsnpr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript compute_normal_se.R <trait>")
trait <- args[1]

## ---- Configuration: update these paths for your environment ----
BASE_DIR   <- getwd()
CHUNK_DIR  <- file.path(BASE_DIR, paste0("chunks_", trait), "bfile")
PHENO_FILE <- file.path(BASE_DIR, "data", "phenotype_residual",
                        paste0("residual_pheno_", trait, ".txt"))
TRAIT_COL  <- "pheno"

OUT_DIR <- file.path(BASE_DIR, paste0("chunks_", trait), "out", trait)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_SUMMARY <- file.path(OUT_DIR, paste0(trait, "_perSNP_summary.tsv.gz"))
OUT_LONG    <- file.path(OUT_DIR, paste0(trait, "_perFamily_perSNP.tsv.gz"))

## Write global headers if new
if (!file.exists(OUT_SUMMARY)) {
  fwrite(data.table(
    CHUNK = character(), marker.ID = character(), SNP = character(),
    beta = double(), se_norm = double(), Sx2 = double(),
    n_pairs = integer(), n_eff = integer()
  ), OUT_SUMMARY)
}
if (!file.exists(OUT_LONG)) {
  fwrite(data.table(
    CHUNK = character(), FID = character(), IID1 = character(),
    IID2 = character(), pair_id = integer(), marker.ID = character(),
    SNP = character(), DeltaX = double(), DeltaY = double()
  ), OUT_LONG)
}

## ---- Helper: build sibling pairs within one family ----
## Prefers PID/MID-based pairing if parental IDs are present;
## otherwise pairs all individuals within a family (by FID).
.make_pairs_one_family <- function(df) {
  if (nrow(df) < 2) return(NULL)
  has_par <- !is.na(df$PID) & df$PID != "0" &
             !is.na(df$MID) & df$MID != "0"
  if (any(has_par)) {
    d <- df[has_par, , drop = FALSE]
    if (nrow(d) < 2) return(NULL)
    d %>%
      dplyr::group_by(PID, MID) %>%
      dplyr::group_map(~ {
        if (nrow(.x) < 2) return(NULL)
        idx <- t(combn(seq_len(nrow(.x)), 2))
        tibble::tibble(
          FID  = .x$FID[idx[, 1]],
          IID1 = .x$IID[idx[, 1]],
          IID2 = .x$IID[idx[, 2]]
        )
      }) %>%
      dplyr::bind_rows()
  } else {
    idx <- t(combn(seq_len(nrow(df)), 2))
    tibble::tibble(
      FID  = df$FID[idx[, 1]],
      IID1 = df$IID[idx[, 1]],
      IID2 = df$IID[idx[, 2]]
    )
  }
}

## ---- Helper: OLS (Normal-theory) slope + SE for a DeltaX block ----
## Computes beta = sum(DX * dY) / sum(DX^2) and the OLS standard error
## for each SNP column in the DeltaX matrix. NA-safe: only uses pairs
## where DeltaX is observed.
normal_se_block <- function(dY, DXb) {
  W   <- !is.na(DXb)
  DX0 <- DXb
  DX0[!W] <- 0

  n_eff <- colSums(W)
  Sx2   <- colSums(DX0 * DX0)
  num   <- as.numeric(crossprod(DX0, dY))
  dy2   <- as.numeric(crossprod(W, dY^2))

  beta <- rep(NA_real_, length(Sx2))
  okS  <- (Sx2 > 0)
  beta[okS] <- num[okS] / Sx2[okS]

  RSS <- rep(NA_real_, length(Sx2))
  RSS[okS] <- dy2[okS] - 2 * beta[okS] * num[okS] +
              (beta[okS]^2) * Sx2[okS]

  se_norm <- rep(NA_real_, length(Sx2))
  okSE    <- okS & (n_eff > 1L)
  se_norm[okSE] <- sqrt((RSS[okSE] / (n_eff[okSE] - 1L)) / Sx2[okSE])

  list(beta = beta, se_norm = se_norm, Sx2 = Sx2,
       n_eff = n_eff, RSS = RSS)
}

## ---- Helper: stream over FBM columns in blocks ----
compute_normal_by_blocks <- function(G, idx1, idx2, dY,
                                     block_size = 4096L) {
  p <- ncol(G)
  beta    <- rep(NA_real_, p)
  se_norm <- rep(NA_real_, p)
  Sx2     <- rep(NA_real_, p)
  n_eff   <- rep(0L, p)

  for (start in seq.int(1L, p, by = block_size)) {
    stop_ <- min(start + block_size - 1L, p)
    j_sub <- start:stop_

    X1  <- as.matrix(G[idx1, j_sub, drop = FALSE])
    X2  <- as.matrix(G[idx2, j_sub, drop = FALSE])
    DXb <- X1 - X2
    rm(X1, X2)

    stats <- normal_se_block(dY, DXb)
    beta[j_sub]    <- stats$beta
    se_norm[j_sub] <- stats$se_norm
    Sx2[j_sub]     <- stats$Sx2
    n_eff[j_sub]   <- stats$n_eff

    rm(DXb, stats); gc(FALSE)
  }

  list(beta = beta, se_norm = se_norm, Sx2 = Sx2, n_eff = n_eff)
}

## ============================================================
## Main: read phenotypes, process each chunk
## ============================================================

## 1) Read phenotypes
pheno <- fread(PHENO_FILE, header = FALSE)
setnames(pheno, c("FID", "IID", TRAIT_COL))
pheno <- pheno %>% mutate(across(c(FID, IID), as.character))

## 2) Locate chunk prefixes
bed_files <- list.files(CHUNK_DIR, pattern = "\\.bed$", full.names = TRUE)
if (length(bed_files) == 0) stop("No .bed files found under: ", CHUNK_DIR)
chunk_prefixes <- sub("\\.bed$", "", bed_files)
message("Found ", length(chunk_prefixes), " chunk(s).")

## 3) Process each chunk
for (cc in seq_along(chunk_prefixes)) {
  pref <- chunk_prefixes[cc]
  fam_file <- paste0(pref, ".fam")
  bed_file <- paste0(pref, ".bed")

  chunk_name <- basename(pref)
  chunk_dir  <- file.path(OUT_DIR, chunk_name)
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

  OUT_SUMMARY_CHUNK <- file.path(chunk_dir,
                                 paste0(chunk_name, "_perSNP_summary.tsv.gz"))
  OUT_LONG_CHUNK    <- file.path(chunk_dir,
                                 paste0(chunk_name, "_perFamily_perSNP.tsv.gz"))

  if (file.exists(OUT_SUMMARY_CHUNK)) file.remove(OUT_SUMMARY_CHUNK)
  if (file.exists(OUT_LONG_CHUNK))    file.remove(OUT_LONG_CHUNK)

  fwrite(data.table(
    CHUNK = character(), marker.ID = character(), SNP = character(),
    beta = double(), se_norm = double(), Sx2 = double(),
    n_pairs = integer(), n_eff = integer()
  ), OUT_SUMMARY_CHUNK)

  fwrite(data.table(
    CHUNK = character(), FID = character(), IID1 = character(),
    IID2 = character(), pair_id = integer(), marker.ID = character(),
    SNP = character(), DeltaX = double(), DeltaY = double()
  ), OUT_LONG_CHUNK)

  message(sprintf("\n[Chunk %d/%d] %s", cc, length(chunk_prefixes),
                  chunk_name))

  ## 3a) Read FAM, join phenotype
  fam <- fread(
    fam_file,
    col.names = c("FID", "IID", "PID", "MID", "SEX", "PHENO"),
    na.strings = c("NA", "-9")
  ) %>% mutate(across(c(FID, IID, PID, MID), as.character))

  ind <- fam %>%
    select(FID, IID, PID, MID) %>%
    inner_join(pheno %>% select(FID, IID, all_of(TRAIT_COL)),
               by = c("FID", "IID"))

  message("  fam N=", nrow(fam), "  pheno N=", nrow(pheno),
          "  joined N=", nrow(ind))

  ## 3b) Build sibling pairs per family
  pairs <- ind %>%
    group_split(FID) %>%
    lapply(.make_pairs_one_family) %>%
    bind_rows()

  if (nrow(pairs) == 0) {
    message("  (No sibling pairs in this chunk; skipping.)")
    next
  }
  pairs <- pairs %>% distinct(FID, IID1, IID2)
  pairs$pair_id <- seq_len(nrow(pairs))

  ## DeltaY for available pairs
  y_map <- ind %>% select(IID, all_of(TRAIT_COL))
  dY <- y_map[[TRAIT_COL]][match(pairs$IID1, y_map$IID)] -
        y_map[[TRAIT_COL]][match(pairs$IID2, y_map$IID)]
  keep_pairs <- !is.na(dY)
  if (!any(keep_pairs)) {
    message("  (No usable DeltaY after NA filter; skipping.)")
    next
  }
  pairs <- pairs[keep_pairs, ]
  dY    <- dY[keep_pairs]
  n_pairs <- length(dY)

  ## 3c) Attach genotypes
  rds_path <- paste0(pref, ".rds")
  if (!file.exists(rds_path)) {
    snp_readBed(bed_file, backingfile = pref)
  }
  obj <- snp_attach(rds_path)
  G   <- obj$genotypes
  map <- obj$map
  if (!("marker.ID" %in% names(map))) stop("map missing 'marker.ID'")
  if (!("rsid" %in% names(map))) map$rsid <- map$marker.ID

  fam_chunk <- fam$IID
  idx1 <- match(pairs$IID1, fam_chunk)
  idx2 <- match(pairs$IID2, fam_chunk)
  keep <- !is.na(idx1) & !is.na(idx2)
  if (!any(keep)) {
    message("  (Pairs not in chunk genotype rows; skipping.)")
    rm(obj, G, map); gc(FALSE); next
  }
  idx1 <- idx1[keep]; idx2 <- idx2[keep]
  pairs <- pairs[keep, ]
  dY    <- dY[keep]
  n_pairs <- length(dY)

  message(sprintf("  Usable sibling pairs: %d   SNPs: %d",
                  n_pairs, ncol(G)))

  ## 3d) Compute OLS stats in column blocks
  stats <- compute_normal_by_blocks(G, idx1, idx2, dY,
                                    block_size = 4096L)

  too_small <- (stats$n_eff <= 1L) | !(stats$Sx2 > 0)
  if (any(too_small)) stats$se_norm[too_small] <- NA_real_

  ## 3e) Write per-SNP summary
  sum_dt <- data.table(
    CHUNK     = chunk_name,
    marker.ID = as.character(map$marker.ID),
    SNP       = as.character(map$rsid),
    beta      = as.numeric(stats$beta),
    se_norm   = as.numeric(stats$se_norm),
    Sx2       = as.numeric(stats$Sx2),
    n_pairs   = as.integer(n_pairs),
    n_eff     = as.integer(stats$n_eff)
  )
  cols <- c("CHUNK", "marker.ID", "SNP", "beta", "se_norm",
            "Sx2", "n_pairs", "n_eff")
  setcolorder(sum_dt, cols)

  fwrite(sum_dt, OUT_SUMMARY_CHUNK, append = TRUE)
  fwrite(sum_dt, OUT_SUMMARY,       append = TRUE)

  ## 3f) Write per-family per-SNP long format
  col_chunk <- 500L
  starts <- seq(1L, ncol(G), by = col_chunk)
  for (st in starts) {
    en  <- min(st + col_chunk - 1L, ncol(G))
    X1  <- as.matrix(G[idx1, st:en, drop = FALSE])
    X2  <- as.matrix(G[idx2, st:en, drop = FALSE])
    DXb <- X1 - X2

    long_dt <- as.data.table(DXb)
    setnames(long_dt, as.character(map$rsid[st:en]))
    long_dt[, pair_id := pairs$pair_id]
    long_dt <- melt(long_dt, id.vars = "pair_id",
                    variable.name = "SNP", value.name = "DeltaX",
                    variable.factor = FALSE)

    annot <- data.table(pair_id = pairs$pair_id, FID = pairs$FID,
                        IID1 = pairs$IID1, IID2 = pairs$IID2,
                        DeltaY = dY)
    long_dt <- annot[long_dt, on = "pair_id"]

    key_map <- data.table(
      SNP       = as.character(map$rsid[st:en]),
      marker.ID = as.character(map$marker.ID[st:en])
    )
    long_dt <- key_map[long_dt, on = "SNP"]
    long_dt[, CHUNK := chunk_name]
    setcolorder(long_dt, c("CHUNK", "FID", "IID1", "IID2", "pair_id",
                           "SNP", "DeltaX", "DeltaY"))

    fwrite(long_dt, OUT_LONG_CHUNK, append = TRUE)
    fwrite(long_dt, OUT_LONG,       append = TRUE)

    rm(X1, X2, DXb, long_dt, annot, key_map); gc(FALSE)
  }

  rm(sum_dt, stats, obj, G, map, pairs, dY, idx1, idx2); gc(FALSE)
}

message("\nDone.")
message("Summary per-SNP:  ", OUT_SUMMARY)
message("Long per-family:  ", OUT_LONG)
