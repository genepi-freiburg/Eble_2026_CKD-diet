#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# MR direct run — dietary components → kidney function (without TwoSampleMR)
# Uses ieugwasr directly for data retrieval + a manual IVW/WM/Egger
# implementation. Comparable to the MR-Agent output (same methods, same
# GWAS IDs) for the validation reported in the Supplement.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ieugwasr)
  library(dplyr)
  library(tidyr)
})

TOKEN <- Sys.getenv("OPENGWAS_JWT")  # Set via: Sys.setenv(OPENGWAS_JWT = "your_token")
Sys.setenv(OPENGWAS_JWT = TOKEN)

set.seed(42)
OUT_DIR <- "analysis/mragent_validation/results"
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# ─── GWAS IDs ────────────────────────────────────────────────────────────────
EXPOSURES <- list(
  "Tea intake"                  = "ukb-b-6066",
  "Salad/raw vegetable intake"  = "ukb-b-1996",
  "Processed meat intake"       = "ukb-b-6324",
  "Fresh fruit intake"          = "ukb-b-3881",
  "Cooked vegetable intake"     = "ukb-b-8089",
  "Dried fruit intake"          = "ukb-b-16576",
  "Oily fish intake"            = "ukb-b-2209",
  "Beef intake"                 = "ukb-b-2862",
  "Pork intake"                 = "ukb-b-5640",
  "Lamb/mutton intake"          = "ukb-b-14179"
)

# Outcomes (best available proxies on OpenGWAS for the validation):
# ebi-a-GCST90103634 = eGFR (creatinine), CKDGen 2019, EUR, n=1,004,040
# ebi-a-GCST90018822 = Chronic renal failure, Sakaue 2021, EUR, n=482,858, 24M SNPs
# Note: the binary CKDGen 2019 CKD statistics used in the primary analysis are
# not fully available on OpenGWAS (ieu-b-7440 = 404; ieu-a-1102 = 2015 array,
# 2.1M SNPs), hence the Sakaue 2021 proxy is used here for the API-based check.
OUTCOMES <- list(
  "eGFR (CKDGen 2019, EUR)"        = "ebi-a-GCST90103634",
  "CKD/renal failure (Sakaue 2021)" = "ebi-a-GCST90018822"
)

# Reference IVW results from the primary analysis (this repository)
PRIMARY_IVW <- data.frame(
  exposure = c("Tea intake","Tea intake","Salad/raw vegetable intake",
               "Processed meat intake","Fresh fruit intake","Dried fruit intake"),
  outcome  = c("eGFR (log-transformed)","CKD (binary)","eGFR (log-transformed)",
               "eGFR (log-transformed)","eGFR (log-transformed)","eGFR (log-transformed)"),
  beta_primary = c(0.01621, -0.4026, 0.02515, -0.01252, -0.004843, 0.007920),
  se_primary   = c(0.00382,  0.1022, 0.01064,  0.00675,  0.006415, 0.005273),
  pval_primary = c(2.17e-5,  8.15e-5, 0.0181,  0.0636,   0.4502,   0.1331),
  nsnp_primary = c(37, 37, 15, 18, 43, 38),
  stringsAsFactors = FALSE
)

# ─── MR functions (IVW, Weighted Median, MR-Egger) ───────────────────────────

mr_ivw <- function(b_exp, b_out, se_out) {
  # Random-effects IVW (multiplicative)
  W     <- 1 / se_out^2
  beta  <- sum(W * b_out * b_exp) / sum(W * b_exp^2)
  se    <- sqrt(1 / sum(W * b_exp^2))
  # Cochran's Q
  Q     <- sum(W * (b_out - beta * b_exp)^2)
  Q_df  <- length(b_exp) - 1
  Q_p   <- pchisq(Q, Q_df, lower.tail=FALSE)
  # RE se (Bowden)
  if (Q / Q_df > 1) se <- se * sqrt(Q / Q_df)
  pval  <- 2 * pnorm(abs(beta / se), lower.tail=FALSE)
  list(method="IVW", beta=beta, se=se, pval=pval,
       Q=Q, Q_df=Q_df, Q_pval=Q_p, nsnp=length(b_exp))
}

mr_weighted_median <- function(b_exp, b_out, se_out, nboot=1000) {
  # Weighted median estimator (Bowden 2016)
  b_iv  <- b_out / b_exp
  w     <- (b_exp^2 / se_out^2)
  w     <- w / sum(w)
  ord   <- order(b_iv)
  b_iv_s<- b_iv[ord]; w_s <- w[ord]
  cumw  <- cumsum(w_s)
  # Median
  med_i <- min(which(cumw >= 0.5))
  beta  <- b_iv_s[med_i]
  # Bootstrap SE
  bs <- replicate(nboot, {
    idx  <- sample(length(b_iv), replace=TRUE)
    b_s  <- b_iv[idx]; w_bs <- w[idx]; w_bs <- w_bs/sum(w_bs)
    ord2 <- order(b_s)
    cm   <- cumsum(w_bs[ord2])
    b_s[ord2][min(which(cm >= 0.5))]
  })
  se   <- sd(bs)
  pval <- 2 * pnorm(abs(beta/se), lower.tail=FALSE)
  list(method="Weighted Median", beta=beta, se=se, pval=pval, nsnp=length(b_exp))
}

mr_egger <- function(b_exp, b_out, se_out) {
  # MR-Egger (Bowden 2015) — orient to positive b_exp
  sign_e <- sign(b_exp)
  b_exp2 <- abs(b_exp)
  b_out2 <- b_out * sign_e
  W      <- 1 / se_out^2
  # WLS
  Sx2    <- sum(W * b_exp2^2)
  Sx     <- sum(W * b_exp2)
  Sy     <- sum(W * b_out2)
  Sxy    <- sum(W * b_exp2 * b_out2)
  Sw     <- sum(W)
  denom  <- Sw * Sx2 - Sx^2
  slope  <- (Sw * Sxy - Sx * Sy) / denom
  intercept <- (Sy * Sx2 - Sx * Sxy) / denom
  # Residuals → SE (simulation extrapolation not applied here)
  res    <- b_out2 - intercept - slope * b_exp2
  sigma2 <- sum(W * res^2) / (length(b_exp) - 2)
  se_slope <- sqrt(sigma2 * Sw / denom)
  se_int   <- sqrt(sigma2 * Sx2 / denom)
  pval_slope <- 2*pnorm(abs(slope/se_slope), lower.tail=FALSE)
  pval_int   <- 2*pnorm(abs(intercept/se_int), lower.tail=FALSE)
  list(method="MR-Egger", beta=slope, se=se_slope, pval=pval_slope,
       intercept=intercept, intercept_se=se_int, intercept_pval=pval_int,
       nsnp=length(b_exp))
}

mr_weighted_mode <- function(b_exp, b_out, se_out, phi=1) {
  # Weighted Mode (Hartwig 2017)
  b_iv   <- b_out / b_exp
  w      <- (b_exp^2 / se_out^2); w <- w / sum(w)
  bw     <- phi * (0.9 * min(sd(b_iv), IQR(b_iv)/1.34)) * length(b_iv)^(-1/5)
  # Density at each candidate
  dens <- density(b_iv, weights=w, bw=bw, n=512)
  beta <- dens$x[which.max(dens$y)]
  # Bootstrap SE
  bs <- replicate(500, {
    idx  <- sample(length(b_iv), replace=TRUE)
    b_s  <- b_iv[idx]; w_bs <- w[idx]; w_bs <- w_bs/sum(w_bs)
    bw2  <- phi * (0.9 * min(sd(b_s), IQR(b_s)/1.34)) * length(b_s)^(-1/5)
    d2   <- density(b_s, weights=w_bs, bw=bw2, n=512)
    d2$x[which.max(d2$y)]
  })
  se   <- sd(bs)
  pval <- 2*pnorm(abs(beta/se), lower.tail=FALSE)
  list(method="Weighted Mode", beta=beta, se=se, pval=pval, nsnp=length(b_exp))
}

# ─── Harmonisierung ──────────────────────────────────────────────────────────

harmonise <- function(exp_dat, out_dat) {
  merged <- inner_join(exp_dat, out_dat, by="rsid", suffix=c("_exp","_out"))
  # Palindromische SNPs identifizieren
  is_palindrome <- function(ea, oa) {
    pairs <- paste0(toupper(ea), toupper(oa))
    pairs %in% c("AT","TA","CG","GC")
  }
  # Allele-Alignment
  merged <- merged %>%
    mutate(
      ea_e = toupper(effect_allele_exp),
      oa_e = toupper(other_allele_exp),
      ea_o = toupper(effect_allele_out),
      oa_o = toupper(other_allele_out),
      palindrome = is_palindrome(ea_e, oa_e),
      # Direct match
      aligned = ea_e == ea_o & oa_e == oa_o,
      # Flipped
      flipped = ea_e == oa_o & oa_e == ea_o,
      beta_out_aligned = case_when(
        aligned ~ beta_out,
        flipped ~ -beta_out,
        TRUE    ~ NA_real_
      ),
      eaf_aligned = case_when(
        aligned ~ eaf_out,
        flipped ~ 1 - eaf_out,
        TRUE    ~ eaf_out
      )
    ) %>%
    # Palindromische SNPs mit EAF ~0.5 entfernen
    filter(!(palindrome & abs(eaf_aligned - 0.5) < 0.08)) %>%
    filter(!is.na(beta_out_aligned))
  merged
}

# ─── Hauptanalyse ────────────────────────────────────────────────────────────

all_results <- list()
snp_counts  <- list()

cat("=== MRAgent-kompatibler MR-Run ===\n")
cat("Exposures:", length(EXPOSURES), "| Outcomes:", length(OUTCOMES), "\n\n")

for (exp_name in names(EXPOSURES)) {
  exp_id <- EXPOSURES[[exp_name]]
  cat("─── Exposure:", exp_name, "(", exp_id, ") ───\n")

  # Instrumente extrahieren
  exp_dat <- tryCatch({
    hits <- tophits(id=exp_id, pval=5e-8, clump=1, r2=0.001, kb=10000,
                    pop="EUR", opengwas_jwt=TOKEN)
    if (is.null(hits) || nrow(hits)==0) {
      cat("  Kein GWS-Signal bei p<5e-8, versuche p<5e-6...\n")
      hits <- tophits(id=exp_id, pval=5e-6, clump=1, r2=0.001, kb=10000,
                      pop="EUR", opengwas_jwt=TOKEN)
    }
    hits %>%
      rename(rsid = rsid, effect_allele_exp = ea, other_allele_exp = nea,
             beta_exp = beta, se_exp = se, pval_exp = p, eaf_exp = eaf) %>%
      select(rsid, effect_allele_exp, other_allele_exp, beta_exp, se_exp, pval_exp, eaf_exp) %>%
      filter(!is.na(rsid), !is.na(beta_exp))
  }, error = function(e) { cat("  ERROR Exposure:", e$message, "\n"); NULL })

  if (is.null(exp_dat) || nrow(exp_dat) < 3) {
    cat("  Zu wenige Instrumente — überspringe\n\n"); next
  }
  cat("  Instrumente:", nrow(exp_dat), "SNPs\n")

  # F-Statistik
  exp_dat <- exp_dat %>%
    mutate(F_stat = (beta_exp / se_exp)^2)
  cat("  F-Statistik (Median):", round(median(exp_dat$F_stat), 1),
      "| Schwache Instrumente (<10):", sum(exp_dat$F_stat < 10), "\n")

  for (out_name in names(OUTCOMES)) {
    out_id <- OUTCOMES[[out_name]]
    cat("  → Outcome:", out_name, "\n")

    # Outcome-Daten abrufen
    out_dat <- tryCatch({
      assoc <- associations(variants=exp_dat$rsid, id=out_id,
                            proxies=0, opengwas_jwt=TOKEN)
      assoc %>%
        rename(rsid = rsid, effect_allele_out = ea, other_allele_out = nea,
               beta_out = beta, se_out = se, pval_out = p, eaf_out = eaf) %>%
        select(rsid, effect_allele_out, other_allele_out, beta_out, se_out, pval_out, eaf_out) %>%
        filter(!is.na(beta_out))
    }, error = function(e) { cat("    ERROR Outcome:", e$message, "\n"); NULL })

    if (is.null(out_dat) || nrow(out_dat) < 3) {
      cat("    Zu wenige Outcome-Daten\n"); next
    }

    # Harmonisierung
    harm <- tryCatch(harmonise(exp_dat, out_dat), error=function(e){ NULL })
    if (is.null(harm) || nrow(harm) < 3) {
      cat("    Zu wenige harmonisierte SNPs:", if(is.null(harm)) 0 else nrow(harm), "\n"); next
    }
    cat("    Harmonisierte SNPs:", nrow(harm), "\n")

    b_e <- harm$beta_exp
    b_o <- harm$beta_out_aligned
    s_o <- harm$se_out

    # MR-Methoden
    ivw  <- tryCatch(mr_ivw(b_e, b_o, s_o),  error=function(e) NULL)
    wm   <- tryCatch(mr_weighted_median(b_e, b_o, s_o), error=function(e) NULL)
    egg  <- tryCatch(mr_egger(b_e, b_o, s_o), error=function(e) NULL)
    wmod <- tryCatch(mr_weighted_mode(b_e, b_o, s_o),  error=function(e) NULL)

    # Assemble result
    res_list <- list(ivw, wm, egg, wmod)
    for (res in res_list) {
      if (is.null(res)) next
      row <- data.frame(
        exposure = exp_name, exposure_id = exp_id,
        outcome  = out_name,  outcome_id  = out_id,
        method   = res$method,
        beta     = res$beta,  se = res$se,
        pval     = res$pval,  nsnp = res$nsnp,
        Q = if(!is.null(res$Q)) res$Q else NA,
        Q_pval = if(!is.null(res$Q_pval)) res$Q_pval else NA,
        intercept     = if(!is.null(res$intercept))     res$intercept     else NA,
        intercept_se  = if(!is.null(res$intercept_se))  res$intercept_se  else NA,
        intercept_pval= if(!is.null(res$intercept_pval))res$intercept_pval else NA,
        stringsAsFactors=FALSE
      )
      # OR für binäre Outcomes
      row$or    <- if(grepl("CKD", out_name)) exp(row$beta) else NA
      row$or_lo <- if(grepl("CKD", out_name)) exp(row$beta - 1.96*row$se) else NA
      row$or_hi <- if(grepl("CKD", out_name)) exp(row$beta + 1.96*row$se) else NA
      row$ci_lo <- row$beta - 1.96*row$se
      row$ci_hi <- row$beta + 1.96*row$se

      all_results[[length(all_results)+1]] <- row
    }
    snp_counts[[paste(exp_name, out_name)]] <- nrow(harm)

    # Print IVW comparison against the primary analysis
    if (!is.null(ivw)) {
      cat(sprintf("    IVW: β=%.4f (SE=%.4f), p=%.3e, nSNP=%d",
                  ivw$beta, ivw$se, ivw$pval, ivw$nsnp))
      # Match against the primary analysis
      primary_out_match <- if(grepl("eGFR", out_name)) "eGFR (log-transformed)" else "CKD (binary)"
      prim <- PRIMARY_IVW %>% filter(exposure==exp_name, outcome==primary_out_match)
      if(nrow(prim)>0) {
        cat(sprintf("\n    Primary:   β=%.4f (SE=%.4f), p=%.3e, nSNP=%d",
                    prim$beta_primary[1], prim$se_primary[1], prim$pval_primary[1], prim$nsnp_primary[1]))
        dir_match <- sign(ivw$beta)==sign(prim$beta_primary[1])
        cat(sprintf("  [direction: %s]", if(dir_match) "✓ same" else "✗ different"))
      }
      cat("\n")
    }
  }
  cat("\n")
}

# ─── Save results ────────────────────────────────────────────────────

if (length(all_results) > 0) {
  results_df <- bind_rows(all_results)
  write.csv(results_df, file.path(OUT_DIR, "mragent_results_all.csv"), row.names=FALSE)

  # IVW results separately
  ivw_df <- results_df %>% filter(method=="IVW") %>%
    mutate(pval_fdr = p.adjust(pval, method="BH")) %>%
    arrange(pval)
  write.csv(ivw_df, file.path(OUT_DIR, "mragent_ivw_results.csv"), row.names=FALSE)

  cat("\n=== RESULTS SUMMARY (IVW) ===\n")
  cat(sprintf("%-30s %-22s %8s %8s %10s %5s %8s\n",
              "Exposure","Outcome","beta","SE","p","nSNP","p_FDR"))
  cat(strrep("-", 95), "\n")
  for(i in seq_len(nrow(ivw_df))) {
    r <- ivw_df[i,]
    sig <- if(r$pval_fdr<0.05) "**FDR" else if(r$pval<0.05) "*nom" else ""
    cat(sprintf("%-30s %-22s %8.4f %8.4f %10.3e %5d %8.3e %s\n",
                substr(r$exposure,1,30), substr(r$outcome,1,22),
                r$beta, r$se, r$pval, r$nsnp, r$pval_fdr, sig))
  }

  # Head-to-head comparison with the primary analysis
  cat("\n=== DIRECT COMPARISON WITH THE PRIMARY ANALYSIS ===\n")
  cat(sprintf("%-30s %-25s | %12s %12s | %12s %12s | %s\n",
              "Exposure","Outcome","β MR-Agent","β primary","p MR-Agent","p primary","direction"))
  cat(strrep("-",110),"\n")
  for(i in seq_len(nrow(PRIMARY_IVW))) {
    prow <- PRIMARY_IVW[i,]
    out_match <- if(grepl("eGFR", prow$outcome)) "eGFR (CKDGen 2019)" else "CKD (CKDGen 2019)"
    mra <- results_df %>%
      filter(exposure==prow$exposure, outcome==out_match, method=="IVW")
    if(nrow(mra)==0) { next }
    dir_ok <- sign(mra$beta[1])==sign(prow$beta_primary)
    cat(sprintf("%-30s %-25s | %12.4f %12.4f | %12.3e %12.3e | %s\n",
                substr(prow$exposure,1,30), substr(prow$outcome,1,25),
                mra$beta[1], prow$beta_primary,
                mra$pval[1], prow$pval_primary,
                if(dir_ok) "✓ same" else "✗ different"))
  }

  cat("\nResults saved:\n")
  cat(" ", file.path(OUT_DIR, "mragent_results_all.csv"), "\n")
  cat(" ", file.path(OUT_DIR, "mragent_ivw_results.csv"), "\n")
} else {
  cat("[WARNING] No results — check token and API connection\n")
}
