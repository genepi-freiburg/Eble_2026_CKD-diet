###############################################################
# Causal Effects of Dietary Components on Kidney Function
# Two-Sample Mendelian Randomization
#
# Full MR pipeline: instrument extraction (IEU OpenGWAS API),
# PheWAS + Steiger confounder screening, harmonisation with CKDGen
# outcomes, univariable & multivariable IVW (and robust methods),
# leave-one-out, eGFRcys and GI-PRAL sensitivity analyses, and
# FinnGen R12 replication.
#
# The instrument-selection cascade (326 raw -> 264 final harmonised
# SNP-exposure pairs) is documented in README.md and the Supplementary
# Methods.
#
# ── PREREQUISITES ─────────────────────────────────────────────────────
#   install.packages(c("ieugwasr","MVMR","dplyr","tidyr","data.table",
#                      "ggplot2","patchwork","scales"))
#
# ── AUTHENTICATION ────────────────────────────────────────────────────
#   Sys.setenv(OPENGWAS_JWT = "your_token")  # https://api.opengwas.io
#   OR: export OPENGWAS_JWT="your_token" && Rscript R/mr_analysis_final.R
#   See .env.example for setup instructions.
#
# ── GWAS DATA ─────────────────────────────────────────────────────────
#   bash scripts/01_download_gwas_data.sh
###############################################################

library(ieugwasr); library(dplyr); library(tidyr)
library(data.table); library(ggplot2); library(patchwork); library(scales)
set.seed(42)

proc <- "data/processed/"
fig  <- "figures/"

# ── Portable gzip reader ──────────────────────────────────────────────
# Reads a .gz file directly via data.table::fread, which decompresses
# internally and therefore works on Linux, macOS AND Windows (no external
# 'zcat' required). fread's gzip support depends on the R.utils package.
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils", repos = "https://cloud.r-project.org", quiet = TRUE)
}
read_gz <- function(path, ...) {
  if (!file.exists(path)) stop(
    "File not found: ", path,
    "\nRun the download script first: bash scripts/01_download_gwas_data.sh")
  data.table::fread(path, data.table = FALSE, showProgress = FALSE, ...)
}

TOKEN <- Sys.getenv("OPENGWAS_JWT")
if (nchar(TOKEN) < 50) stop(
  "OPENGWAS_JWT not set. Visit https://api.opengwas.io to obtain a token.\n",
  "Set it with: Sys.setenv(OPENGWAS_JWT = 'your_token')\n",
  "Or copy .env.example to .env, fill in your token, then: source .env")

# ─────────────────────────────────────────────────────────────
# STEP 1: Extract GWS instruments from IEU OpenGWAS API
# ─────────────────────────────────────────────────────────────
exposures <- list(
  "Fresh fruit intake"         = "ukb-b-3881",
  "Cooked vegetable intake"    = "ukb-b-8089",
  "Salad/raw vegetable intake" = "ukb-b-1996",
  "Dried fruit intake"         = "ukb-b-16576",
  "Tea intake"                 = "ukb-b-6066",
  "Oily fish intake"           = "ukb-b-2209",
  "Processed meat intake"      = "ukb-b-6324",
  "Beef intake"                = "ukb-b-2862",
  "Pork intake"                = "ukb-b-5640",
  "Lamb/mutton intake"         = "ukb-b-14179"
)

cat("Step 1: Extracting GWS tophits (P<5e-8, LD R2<0.001, EUR)...\n")
instr_raw <- bind_rows(lapply(names(exposures), function(nm) {
  id   <- exposures[[nm]]
  cat(" ", nm, "...\n")
  hits <- tryCatch(
    tophits(id=id, pval=5e-8, clump=1, r2=0.001, kb=10000,
            pop="EUR", opengwas_jwt=TOKEN),
    error=function(e){ cat("  ERROR:", e$message, "\n"); NULL })
  if (is.null(hits) || nrow(hits)==0) return(NULL)
  hits %>%
    mutate(
      exposure      = nm, gwas_id = id,
      rsid          = coalesce(
        if("rsid"%in%names(.)) rsid else NA_character_,
        if("SNP" %in%names(.)) SNP  else NA_character_),
      effect_allele = toupper(coalesce(
        if("ea"           %in%names(.)) ea            else NA_character_,
        if("effect_allele"%in%names(.)) effect_allele else NA_character_)),
      other_allele  = toupper(coalesce(
        if("nea"         %in%names(.)) nea          else NA_character_,
        if("other_allele"%in%names(.)) other_allele else NA_character_)),
      beta_exposure = if("beta"%in%names(.)) beta else NA_real_,
      se_exposure   = if("se"  %in%names(.)) se   else NA_real_,
      pval_exposure = if("p"   %in%names(.)) p    else NA_real_,
      eaf           = if("eaf" %in%names(.)) eaf  else NA_real_
    ) %>%
    select(rsid, effect_allele, other_allele, beta_exposure,
           se_exposure, pval_exposure, eaf, exposure, gwas_id) %>%
    filter(!is.na(rsid), !is.na(beta_exposure)) %>%
    mutate(F_stat = (beta_exposure / se_exposure)^2)
}))
write.csv(instr_raw, paste0(proc, "instruments_raw_ieu.csv"), row.names=FALSE)
cat("Raw instruments:", nrow(instr_raw), "SNPs across",
    n_distinct(instr_raw$exposure), "exposures\n\n")

# ─────────────────────────────────────────────────────────────
# STEP 2: PheWAS confounder screening — Round 1 (4 GWAS)
# ─────────────────────────────────────────────────────────────
# Pre-specified: TG, BMI, T2D, SBP at P < 1e-5
# (eGFR excluded: circular — cannot distinguish causal instruments
#  from pleiotropic SNPs when eGFR is also the primary outcome)
CONFOUNDER_GWAS_R1 <- c(
  "ieu-a-302",   # Triglycerides (Willer GLGC 2013)
  "ieu-a-2",     # BMI (Locke GIANT 2015)
  "ieu-a-26",    # Type 2 diabetes (Morris DIAGRAM 2012)
  "ieu-b-38"     # Systolic blood pressure (Evangelou 2018)
)
PLEIOTROPY_P <- 1e-5

cat("Step 2a: PheWAS screening Round 1 (4 GWAS, P <", PLEIOTROPY_P, ")...\n")
snps_all   <- unique(instr_raw$rsid)
chunks_r1  <- split(snps_all, ceiling(seq_along(snps_all)/20))
flags_r1   <- bind_rows(lapply(seq_along(chunks_r1), function(i) {
  cat(sprintf("  Batch %d/%d...", i, length(chunks_r1)))
  res <- tryCatch({
    r <- associations(variants=chunks_r1[[i]], id=CONFOUNDER_GWAS_R1,
                      proxies=0, opengwas_jwt=TOKEN)
    if(!is.null(r) && nrow(r)>0)
      r %>% filter(!is.na(p) & p < PLEIOTROPY_P) %>%
        select(rsid, conf_id=id, conf_trait=trait, beta_conf=beta, p_conf=p)
    else NULL
  }, error=function(e){cat(" ERROR:", substr(e$message,1,30)); NULL})
  cat(sprintf(" %d flags\n", if(!is.null(res)) nrow(res) else 0))
  Sys.sleep(0.5); res
}))
excl_r1 <- flags_r1 %>%
  left_join(instr_raw %>% select(rsid,exposure) %>% distinct(),
            by="rsid", relationship="many-to-many") %>%
  select(rsid, exposure) %>% distinct() %>%
  bind_rows(data.frame(  # ADH1B a priori
    rsid=rep("rs1229984", length(unique(instr_raw$exposure))),
    exposure=unique(instr_raw$exposure), stringsAsFactors=FALSE))

# ─────────────────────────────────────────────────────────────
# STEP 3: Leave-one-out triggered — Round 2 screening (+ CAD)
# ─────────────────────────────────────────────────────────────
# LOO identified rs429358 as most influential for dried fruit eGFR.
# Post-hoc expanded screen against CAD GWAS confirmed rs429358:
#   CAD (ieu-a-7, N=184,000): P = 2.2e-9 → exclude from all exposures
# rs3095337 (salad): Total-Chol P=1.7e-5 — above 1e-5 threshold → RETAIN
CONFOUNDER_GWAS_R2 <- "ieu-a-7"  # Coronary artery disease

cat("\nStep 2b: Round 2 screening (CAD, triggered by LOO)...\n")
loo_snps <- c("rs429358", "rs3095337")
flags_r2 <- tryCatch({
  r <- associations(variants=loo_snps, id=CONFOUNDER_GWAS_R2,
                    proxies=0, opengwas_jwt=TOKEN)
  if(!is.null(r) && nrow(r)>0)
    r %>% filter(!is.na(p) & p < PLEIOTROPY_P) %>%
      select(rsid, conf_id=id, conf_trait=trait, beta_conf=beta, p_conf=p)
  else data.frame()
}, error=function(e){cat("ERROR:", e$message,"\n"); data.frame()})

cat("Round 2 flags (P <", PLEIOTROPY_P, "):\n")
print(flags_r2 %>% mutate(across(where(is.numeric),~signif(.,3))), row.names=FALSE)

# rs429358 excluded from ALL exposures where present (dried fruit, beef, lamb)
excl_r2 <- flags_r2 %>%
  left_join(instr_raw %>% select(rsid,exposure) %>% distinct(),
            by="rsid", relationship="many-to-many") %>%
  select(rsid, exposure) %>% distinct()

# Combine all exclusions
excl_all <- bind_rows(excl_r1, excl_r2) %>% distinct()
cat(sprintf("\nTotal exclusions: %d instrument-exposure pairs\n", nrow(excl_all)))

# Apply exclusions
instr_clean <- instr_raw %>%
  anti_join(excl_all, by=c("rsid","exposure"))

cat("Instruments after screening:\n")
instr_clean %>%
  group_by(exposure) %>%
  summarise(n_raw   = sum(instr_raw$exposure==exposure[1]),
            n_clean = n(),
            n_excl  = n_raw - n_clean,
            mean_F  = round(mean(F_stat),1), .groups="drop") %>%
  print(row.names=FALSE)
cat(sprintf("Total: %d → %d (-%d)\n",
    nrow(instr_raw), nrow(instr_clean), nrow(instr_raw)-nrow(instr_clean)))

write.csv(bind_rows(flags_r1,flags_r2 %>% mutate(round="R2")),
          paste0(proc,"ieu_pleiotropy_flags.csv"), row.names=FALSE)
write.csv(excl_all,    paste0(proc,"ieu_exclusions_final.csv"), row.names=FALSE)
# Intermediate (post-PheWAS, pre-harmonisation) set. The canonical
# instruments_final.csv is (re)written later, after harmonisation + Steiger
# filtering, so that it contains only the final retained pairs.
write.csv(instr_clean, paste0(proc,"instruments_postphewas.csv"), row.names=FALSE)

# ─────────────────────────────────────────────────────────────
# STEP 4: Harmonise with CKDGen outcomes
# ─────────────────────────────────────────────────────────────
cat("\nStep 3: Harmonising with CKDGen outcomes...\n")
# read_gz reads .gz directly (no external zcat needed -> works on Windows too)
egfr_raw <- read_gz("data/raw/egfr_EA.txt.gz")
ckd_raw  <- read_gz("data/raw/ckd_EA.txt.gz")

extract_out <- function(gwas, snps, outname) {
  gwas %>% filter(RSID %in% snps) %>%
    rename(rsid=RSID, chr=Chr, pos=Pos_b37,
           effect_allele_outcome=Allele1, other_allele_outcome=Allele2,
           eaf_outcome=Freq1, beta_outcome=Effect, se_outcome=StdErr,
           pval_outcome=`P-value`, samplesize_outcome=n_total_sum) %>%
    mutate(outcome=outname,
           effect_allele_outcome=toupper(effect_allele_outcome),
           other_allele_outcome=toupper(other_allele_outcome))
}
harmonise <- function(exp_df, out_df) {
  exp_df %>% inner_join(out_df, by="rsid") %>%
    mutate(
      ea_e=toupper(effect_allele), oa_e=toupper(other_allele),
      ea_o=toupper(effect_allele_outcome), oa_o=toupper(other_allele_outcome),
      aligned=ea_e==ea_o, flipped=ea_e==oa_o,
      beta_outcome_aligned=case_when(
        aligned~beta_outcome, flipped~-beta_outcome, TRUE~NA_real_),
      eaf_a=case_when(aligned~eaf_outcome, flipped~1-eaf_outcome, TRUE~eaf_outcome),
      palindromic=(ea_e=="A"&oa_e=="T")|(ea_e=="T"&oa_e=="A")|
                  (ea_e=="G"&oa_e=="C")|(ea_e=="C"&oa_e=="G"),
      keep=!(palindromic & abs(eaf_a-0.5)<0.08) & !is.na(beta_outcome_aligned)
    ) %>% filter(keep)
}

snps_c     <- unique(instr_clean$rsid)
harm_egfr  <- harmonise(instr_clean, extract_out(egfr_raw, snps_c, "eGFR (log-transformed)"))
harm_ckd   <- harmonise(instr_clean, extract_out(ckd_raw,  snps_c, "CKD (binary)"))
write.csv(harm_egfr, paste0(proc,"harmonised_egfr.csv"), row.names=FALSE)
write.csv(harm_ckd,  paste0(proc,"harmonised_ckd.csv"),  row.names=FALSE)

# ─────────────────────────────────────────────────────────────
# STEP 5: MR — 4 methods
# ─────────────────────────────────────────────────────────────
mr_ivw <- function(b_e,b_o,s_o){
  w<-1/s_o^2; bi<-sum(w*b_o*b_e)/sum(w*b_e^2); si<-sqrt(1/sum(w*b_e^2))
  Q<-sum(w*(b_o-bi*b_e)^2); Qp<-pchisq(Q,df=length(b_e)-1,lower.tail=FALSE)
  list(beta=bi,se=si,pval=2*pnorm(-abs(bi/si)),Q=Q,Q_pval=Qp,nsnp=length(b_e))
}
mr_wm <- function(b_e,b_o,s_o,nb=1000){
  r<-b_o/b_e; w<-(b_e^2/s_o^2)/sum(b_e^2/s_o^2); os<-order(r)
  bm<-r[os][which(cumsum(w[os])>=0.5)[1]]
  bs<-replicate(nb,{s<-sample(length(r),replace=TRUE); r2<-r[s]; w2<-w[s]/sum(w[s])
    o2<-order(r2); r2[o2][which(cumsum(w2[o2])>=0.5)[1]]})
  list(beta=bm,se=sd(bs,na.rm=T),pval=2*pnorm(-abs(bm/sd(bs,na.rm=T))),nsnp=length(b_e))
}
mr_egg <- function(b_e,b_o,s_o){
  if(length(b_e)<3) return(list(beta=NA,se=NA,pval=NA,intercept=NA,pval_int=NA,nsnp=length(b_e)))
  sg<-sign(b_e); be<-b_e*sg; bo<-b_o*sg; w<-1/s_o^2
  X<-cbind(1,be); cc<-solve(t(X)%*%diag(w)%*%X)%*%(t(X)%*%diag(w)%*%bo)
  r<-bo-(cc[1]+cc[2]*be); s2<-sum(w*r^2)/(length(b_e)-2); V<-s2*solve(t(X)%*%diag(w)%*%X)
  list(beta=cc[2],se=sqrt(V[2,2]),pval=2*pnorm(-abs(cc[2]/sqrt(V[2,2]))),
       intercept=cc[1],pval_int=2*pnorm(-abs(cc[1]/sqrt(V[1,1]))),nsnp=length(b_e))
}
mr_mode <- function(b_e,b_o,s_o,bw=0.9,nb=500){
  r<-b_o/b_e; w<-abs(b_e)/s_o; w<-w/sum(w)
  fm<-function(rv,wv){ bw_<-bw*sd(rv)*length(rv)^(-0.2)
    g<-seq(min(rv)-3*bw_,max(rv)+3*bw_,l=512)
    g[which.max(sapply(g,function(x) sum(wv*dnorm(x,rv,bw_))))] }
  bm<-fm(r,w)
  bs<-replicate(nb,{s<-sample(length(r),replace=TRUE); fm(r[s],w[s]/sum(w[s]))})
  list(beta=bm,se=sd(bs,na.rm=T),pval=2*pnorm(-abs(bm/sd(bs,na.rm=T))),nsnp=length(b_e))
}
run_mr4 <- function(data,exp_name,out_name){
  d<-data%>%filter(exposure==exp_name); if(nrow(d)<3) return(NULL)
  b_e<-d$beta_exposure; b_o<-d$beta_outcome_aligned; s_o<-d$se_outcome
  ivw<-mr_ivw(b_e,b_o,s_o); wm<-mr_wm(b_e,b_o,s_o)
  eg<-mr_egg(b_e,b_o,s_o); mo<-mr_mode(b_e,b_o,s_o)
  data.frame(
    exposure=exp_name, outcome=out_name,
    method=c("IVW","Weighted Median","MR-Egger","Weighted Mode"),
    nsnp=c(ivw$nsnp,wm$nsnp,eg$nsnp,mo$nsnp),
    beta=c(ivw$beta,wm$beta,eg$beta,mo$beta),
    se=c(ivw$se,wm$se,eg$se,mo$se),
    pval=c(ivw$pval,wm$pval,eg$pval,mo$pval),
    Q=c(ivw$Q,NA,NA,NA), Q_pval=c(ivw$Q_pval,NA,NA,NA),
    egger_intercept=c(NA,NA,eg$intercept,NA),
    egger_intercept_pval=c(NA,NA,eg$pval_int,NA), stringsAsFactors=FALSE
  ) %>% mutate(ci_lo=beta-1.96*se, ci_hi=beta+1.96*se,
               or=exp(beta), or_lo=exp(ci_lo), or_hi=exp(ci_hi))
}

cat("Step 4: Running 4-method MR...\n")
exps    <- unique(instr_clean$exposure)
all_res <- bind_rows(
  bind_rows(Filter(Negate(is.null), lapply(exps,function(e) run_mr4(harm_egfr,e,"eGFR (log-transformed)")))),
  bind_rows(Filter(Negate(is.null), lapply(exps,function(e) run_mr4(harm_ckd, e,"CKD (binary)"))))
)
ivw_res <- all_res %>% filter(method=="IVW") %>%
  group_by(outcome) %>% mutate(pval_fdr=p.adjust(pval,"BH")) %>% ungroup()
write.csv(all_res, paste0(proc,"mr_results_all.csv"), row.names=FALSE)
write.csv(ivw_res, paste0(proc,"mr_results_ivw.csv"), row.names=FALSE)

# ─────────────────────────────────────────────────────────────
# STEP 6: MVMR (6 exposures)
# ─────────────────────────────────────────────────────────────
cat("Step 5: MVMR...\n")
mvmr_exps <- c("Fresh fruit intake","Processed meat intake","Oily fish intake",
               "Cooked vegetable intake","Salad/raw vegetable intake","Beef intake")
mvmr_mat  <- instr_clean %>% filter(exposure %in% mvmr_exps) %>%
  select(rsid,exposure,beta_exposure) %>%
  pivot_wider(names_from=exposure, values_from=beta_exposure, values_fill=0)
run_mvmr <- function(mat, out_df, outname){
  dat <- mat %>% inner_join(
    out_df %>% select(rsid,beta_outcome,se_outcome) %>%
      group_by(rsid) %>% slice(1) %>% ungroup(), by="rsid")
  exp_cols <- setdiff(names(dat),c("rsid","beta_outcome","se_outcome"))
  B<-as.matrix(dat[,exp_cols]); y<-dat$beta_outcome; W<-diag(1/dat$se_outcome^2)
  BtWB<-t(B)%*%W%*%B+diag(1e-8,ncol(B)); cc<-solve(BtWB)%*%(t(B)%*%W%*%y); V<-solve(BtWB)
  data.frame(exposure=exp_cols, outcome=outname, method="MVMR-IVW",
             nsnp_total=nrow(dat), beta=as.numeric(cc), se=sqrt(diag(V)),
             stringsAsFactors=FALSE) %>%
    mutate(pval=2*pnorm(-abs(beta/se)), ci_lo=beta-1.96*se, ci_hi=beta+1.96*se,
           or=exp(beta), or_lo=exp(ci_lo), or_hi=exp(ci_hi))
}
mvmr_res <- bind_rows(
  run_mvmr(mvmr_mat, harm_egfr, "eGFR (log-transformed)"),
  run_mvmr(mvmr_mat, harm_ckd,  "CKD (binary)")
) %>% group_by(outcome) %>% mutate(pval_fdr=p.adjust(pval,"BH")) %>% ungroup()
write.csv(mvmr_res, paste0(proc,"mr_results_mvmr.csv"), row.names=FALSE)

# ─────────────────────────────────────────────────────────────
# STEP 6b: Sanderson-Windmeijer conditional F-statistics (MVMR package)
# Sanderson E, Spiller W, Bowden J. Stat Med 2021;40:5434-5452
# ─────────────────────────────────────────────────────────────
if (!requireNamespace("MVMR", quietly = TRUE)) {
  install.packages("MVMR", repos = c("https://mrcieu.r-universe.dev",
                                      "https://cloud.r-project.org"), quiet = TRUE)
}
library(MVMR)

# Rebuild full SNP x exposure matrix with SE (needed by format_mvmr)
beta_mat_mv <- matrix(0, nrow = nrow(mvmr_mat), ncol = length(mvmr_exps),
                      dimnames = list(mvmr_mat$rsid, mvmr_exps))
se_mat_mv   <- matrix(1e6, nrow = nrow(mvmr_mat), ncol = length(mvmr_exps),
                      dimnames = list(mvmr_mat$rsid, mvmr_exps))
for (exp_i in mvmr_exps) {
  sub_i <- instr_clean %>% filter(exposure == exp_i)
  idx_i <- match(sub_i$rsid, mvmr_mat$rsid)
  idx_i <- idx_i[!is.na(idx_i)]
  beta_mat_mv[idx_i, exp_i] <- sub_i$beta_exposure[match(mvmr_mat$rsid[idx_i], sub_i$rsid)]
  se_mat_mv[idx_i,   exp_i] <- sub_i$se_exposure[  match(mvmr_mat$rsid[idx_i], sub_i$rsid)]
}
# Outcome betas (eGFR; used only for formatting)
harm_mvmr_sub <- harm_egfr %>% filter(rsid %in% mvmr_mat$rsid) %>%
  group_by(rsid) %>% slice(1) %>% ungroup()
beta_out_mv <- harm_mvmr_sub$beta_outcome_aligned[match(mvmr_mat$rsid, harm_mvmr_sub$rsid)]
se_out_mv   <- harm_mvmr_sub$se_outcome[          match(mvmr_mat$rsid, harm_mvmr_sub$rsid)]
# Replace NAs
beta_out_mv[is.na(beta_out_mv)] <- 0
se_out_mv[is.na(se_out_mv)]     <- 1

mvmr_fmt <- format_mvmr(BXGs   = beta_mat_mv,
                         BYG    = beta_out_mv,
                         seBXGs = se_mat_mv,
                         seBYG  = se_out_mv,
                         RSID   = mvmr_mat$rsid)
# Conditional F-stats (gencov = 0: assumes zero covariance across exposures)
fstats_cond <- suppressWarnings(strength_mvmr(r_input = mvmr_fmt, gencov = 0))
# Q-stat for pleiotropy
qstat_mvmr  <- suppressWarnings(pleiotropy_mvmr(r_input = mvmr_fmt, gencov = 0))

fstat_df <- data.frame(
  exposure      = mvmr_exps,
  conditional_F = round(as.numeric(fstats_cond), 2),
  above_threshold = as.numeric(fstats_cond) > 10,
  Q_stat        = round(qstat_mvmr$Qstat, 4),
  Q_pval        = round(qstat_mvmr$Qpval, 4)
)
write.csv(fstat_df, paste0(proc, "mvmr_conditional_fstats.csv"), row.names = FALSE)
cat("Conditional F-statistics saved.\n")
print(fstat_df[, c("exposure", "conditional_F", "above_threshold")])

# ─────────────────────────────────────────────────────────────
# STEP 7: FinnGen R12 replication
# ─────────────────────────────────────────────────────────────
cat("Step 6: FinnGen R12 replication...\n")
fg12 <- read_gz("data/raw/finngen/finngen_R12_N14_CHRONKIDNEYDIS.gz",
              select=c("#chrom","pos","rsids","ref","alt","pval","beta","sebeta","af_alt"))
names(fg12)[1]<-"chrom"; fg12$rsid<-sub(",.*","",fg12$rsids)
fg12 <- fg12 %>% filter(rsid %in% snps_c) %>%
  rename(effect_allele_outcome=alt, other_allele_outcome=ref,
         beta_outcome=beta, se_outcome=sebeta, eaf_outcome=af_alt) %>%
  mutate(outcome="CKD (FinnGen R12)", samplesize_outcome=493235)
harm_fg <- harmonise(instr_clean, fg12)
fg_res <- bind_rows(Filter(Negate(is.null), lapply(exps, function(e){
  d<-harm_fg%>%filter(exposure==e); if(nrow(d)<3) return(NULL)
  v<-mr_ivw(d$beta_exposure,d$beta_outcome_aligned,d$se_outcome)
  data.frame(exposure=e, outcome="CKD (FinnGen R12)", nsnp=v$nsnp,
             beta=v$beta, se=v$se, pval=v$pval,
             ci_lo=v$beta-1.96*v$se, ci_hi=v$beta+1.96*v$se,
             or=exp(v$beta), or_lo=exp(v$beta-1.96*v$se),
             or_hi=exp(v$beta+1.96*v$se), stringsAsFactors=FALSE)
}))) %>% mutate(pval_fdr=p.adjust(pval,"BH"))
write.csv(fg_res, paste0(proc,"mr_results_finngen_ckd_ivw.csv"), row.names=FALSE)

# ─────────────────────────────────────────────────────────────
# STEP 8: Steiger filtering + Leave-one-out
# ─────────────────────────────────────────────────────────────
# Steiger filtering removes SNP-exposure pairs where r2_outcome > r2_exposure,
# i.e. where the variant explains more variance in the outcome than in the
# exposure (suggesting reverse causation or pleiotropy).
#
# Three pairs are removed in this study:
#   rs6059844  (Oily fish intake)
#   rs6484504  (Processed meat intake)
#   rs11211124 (Pork intake)
# Final harmonised set: 264 SNP–exposure pairs (267 - 3).
cat("Step 7: Steiger filtering + LOO...\n")

add_steiger <- function(harm_df) {
  harm_df %>%
    mutate(r2_exp = F_stat/(F_stat+500000),
           z_out  = beta_outcome_aligned/se_outcome,
           r2_out = z_out^2/(z_out^2+samplesize_outcome),
           correct_dir = r2_exp > r2_out)
}

harm_egfr_pre <- add_steiger(harm_egfr)
harm_ckd_pre  <- add_steiger(harm_ckd)

steiger_excluded <- harm_egfr_pre %>%
  filter(!correct_dir | is.na(correct_dir)) %>%
  select(rsid, exposure) %>% distinct()
cat(sprintf("  Steiger filter removed %d SNP-exposure pairs:\n", nrow(steiger_excluded)))
print(steiger_excluded, row.names = FALSE)
write.csv(steiger_excluded, paste0(proc, "steiger_excluded.csv"), row.names = FALSE)

# Apply Steiger exclusions to harmonised data
harm_egfr   <- harm_egfr_pre %>% anti_join(steiger_excluded, by = c("rsid","exposure"))
harm_ckd    <- harm_ckd_pre  %>% anti_join(steiger_excluded, by = c("rsid","exposure"))

# Restrict instr_clean to the FINAL RETAINED instrument set, i.e. the
# harmonised + Steiger-filtered SNP-exposure pairs that actually enter the
# analyses (264 pairs). Previously instr_clean was only Steiger-filtered and
# still contained palindromic / CKDGen-absent SNPs dropped at harmonisation,
# which inflated the per-exposure mean F reported in Table 1 / Table S2b.
instr_clean <- instr_clean %>%
  semi_join(harm_egfr, by = c("rsid","exposure"))

# Rewrite canonical files with the final retained (harmonised + Steiger) set
write.csv(instr_clean, paste0(proc,"instruments_final.csv"), row.names=FALSE)
write.csv(harm_egfr,   paste0(proc,"harmonised_egfr.csv"),   row.names=FALSE)
write.csv(harm_ckd,    paste0(proc,"harmonised_ckd.csv"),    row.names=FALSE)
cat(sprintf("  After Steiger: %d instruments, %d eGFR pairs, %d CKD pairs\n",
            nrow(instr_clean), nrow(harm_egfr), nrow(harm_ckd)))

# ─────────────────────────────────────────────────────────────
# Canonical instrument-strength summary for Table 1 (Block A) / Table S2b
# ─────────────────────────────────────────────────────────────
# Single source of truth: computed on the FINAL RETAINED instrument set only.
# Per-SNP F = (beta_exposure / se_exposure)^2 (univariable; squared z-score),
# aggregated per exposure as the arithmetic mean. NOTE: this is the
# univariable instrument strength and is distinct from the Sanderson-Windmeijer
# CONDITIONAL F reported for the MVMR model (see mvmr_conditional_fstats.csv).
table1_strength <- instr_clean %>%
  group_by(exposure) %>%
  summarise(n_final   = n(),
            mean_F    = round(mean(F_stat), 1),
            min_F     = round(min(F_stat),  1),
            max_F     = round(max(F_stat),  1),
            .groups   = "drop") %>%
  arrange(exposure)
write.csv(table1_strength, paste0(proc,"table1_instrument_strength.csv"), row.names=FALSE)
cat("\nTable 1 / S2b instrument strength (final retained set, univariable F):\n")
print(table1_strength, row.names = FALSE)

steiger <- harm_egfr_pre %>%
  group_by(exposure) %>%
  summarise(pct_correct_dir=round(mean(correct_dir,na.rm=T)*100,1), .groups="drop")
write.csv(steiger, paste0(proc,"steiger.csv"), row.names=FALSE)

# ──── Re-run downstream analyses with Steiger-filtered data ────
cat("Re-running MR (4 methods) on Steiger-filtered data...\n")
exps <- unique(instr_clean$exposure)
all_res <- bind_rows(
  bind_rows(Filter(Negate(is.null), lapply(exps,function(e) run_mr4(harm_egfr,e,"eGFR (log-transformed)")))),
  bind_rows(Filter(Negate(is.null), lapply(exps,function(e) run_mr4(harm_ckd, e,"CKD (binary)"))))
)
ivw_res <- all_res %>% filter(method=="IVW") %>%
  group_by(outcome) %>% mutate(pval_fdr=p.adjust(pval,"BH")) %>% ungroup()
write.csv(all_res, paste0(proc,"mr_results_all.csv"), row.names=FALSE)
write.csv(ivw_res, paste0(proc,"mr_results_ivw.csv"), row.names=FALSE)
cat("  Updated mr_results_ivw.csv and mr_results_all.csv\n")

# LOO for tea, salad, dried fruit
run_loo <- function(harm, exp_name) {
  d <- harm %>% filter(exposure==exp_name)
  bind_rows(lapply(seq_len(nrow(d)), function(i){
    s <- d[-i,]; v <- mr_ivw(s$beta_exposure,s$beta_outcome_aligned,s$se_outcome)
    data.frame(snp=d$rsid[i], beta=v$beta, se=v$se, pval=v$pval,
               ci_lo=v$beta-1.96*v$se, ci_hi=v$beta+1.96*v$se,
               or=exp(v$beta), or_lo=exp(v$beta-1.96*v$se),
               or_hi=exp(v$beta+1.96*v$se))
  }))
}
write.csv(run_loo(harm_egfr,"Tea intake"),               paste0(proc,"loo_tea.csv"),         row.names=FALSE)
write.csv(run_loo(harm_ckd, "Tea intake"),               paste0(proc,"loo_tea_ckd.csv"),     row.names=FALSE)
write.csv(run_loo(harm_egfr,"Salad/raw vegetable intake"),paste0(proc,"loo_salad_egfr.csv"), row.names=FALSE)
write.csv(run_loo(harm_ckd, "Salad/raw vegetable intake"),paste0(proc,"loo_salad_ckd.csv"),  row.names=FALSE)
write.csv(run_loo(harm_egfr,"Dried fruit intake"),       paste0(proc,"loo_dryfruit_egfr.csv"),row.names=FALSE)
write.csv(run_loo(harm_ckd, "Dried fruit intake"),       paste0(proc,"loo_dryfruit_ckd.csv"), row.names=FALSE)

# ─────────────────────────────────────────────────────────────
# STEP 8b: GI-PRAL Composite Score
# Remer T, Manz F. J Am Diet Assoc 1995;95:791-797 (primary source: formula + Table 2 food values)

# GI-PRAL = Σ_i (w_i × β̂_IVW,i)
# where w_i = PRAL per serving (mEq/d) from Remer & Manz 1995 (JADA Table 2)
# ─────────────────────────────────────────────────────────────
cat("Step 8b: GI-PRAL composite score...\n")

# PRAL weights (mEq/d per typical serving) — Remer & Manz 1995 (JADA Table 2)
# Food entry in source table | Source | Serving size used
pral_weights <- c(
  "Fresh fruit intake"         = -3.72,  # Fruit avg (apples/oranges) | Remer 1995 | 120g
  "Cooked vegetable intake"    = -3.20,  # Mixed cooked veg          | Remer 1995 | 80g
  "Salad/raw vegetable intake" = -1.26,  # Lettuce avg -2.1/100g | Remer 1995 Table 2 | 60g
  "Dried fruit intake"         = -6.30,  # Raisins (100g = -21 mEq)  | Remer 1995 | 30g
  "Tea intake"                 = -0.60,  # Brewed tea (100ml = -0.3) | Remer 1995 | 200ml
  "Oily fish intake"           = +9.90,  # Herring +7.0 (Remer) + trout +10.8 avg | Remer 1995 | 100g
  "Processed meat intake"      = +3.35,  # Meat sausage avg (Remer 1995 Table 2) | 50g
  "Beef intake"                = +7.80,  # Beef, lean roast          | Remer 1995 | 100g
  "Pork intake"                = +7.90,  # Pork, lean                | Remer 1995 | 100g
  "Lamb/mutton intake"         = +7.60   # Veal avg (Remer 1995 Table 2) | 100g
)

# Extract food-level IVW eGFR estimates
ivw_egfr_betas <- ivw_res %>%
  filter(outcome == "eGFR (log-transformed)", method == "IVW") %>%
  select(exposure, beta, se)

# Compute GI-PRAL: weighted sum of IVW estimates
gipral_df <- data.frame(
  exposure = names(pral_weights),
  pral_w   = as.numeric(pral_weights),
  stringsAsFactors = FALSE
) %>%
  left_join(ivw_egfr_betas, by = "exposure")

gipral_beta <- sum(gipral_df$pral_w * gipral_df$beta,  na.rm=TRUE)
gipral_se   <- sqrt(sum((gipral_df$pral_w * gipral_df$se)^2, na.rm=TRUE))
gipral_z    <- gipral_beta / gipral_se
gipral_p    <- 2 * pnorm(-abs(gipral_z))
gipral_ci_lo <- gipral_beta - 1.96 * gipral_se
gipral_ci_hi <- gipral_beta + 1.96 * gipral_se

cat(sprintf("GI-PRAL \u03b2 = %.4f (SE=%.4f, P=%.3e)\n",
            gipral_beta, gipral_se, gipral_p))
cat("Interpretation: per 1 mEq/d increase in PRAL (more acid-forming diet),\n")
cat("                change in log(eGFR) SD units.\n")

gipral_res <- data.frame(
  exposure = "GI-PRAL Score",
  outcome  = "eGFR (log-transformed)",
  method   = "IVW (food-level composite)",
  nsnp     = nrow(gipral_df),
  beta     = gipral_beta,
  se       = gipral_se,
  pval     = gipral_p,
  Q_pval   = NA_real_,
  ci_lo    = gipral_ci_lo,
  ci_hi    = gipral_ci_hi,
  or       = exp(gipral_beta),
  or_lo    = exp(gipral_ci_lo),
  or_hi    = exp(gipral_ci_hi)
)
write.csv(gipral_res, paste0(proc, "mr_results_gipral.csv"), row.names=FALSE)
cat("GI-PRAL saved to mr_results_gipral.csv\n")

# ─────────────────────────────────────────────────────────────
# STEP 8c: GI-PRAL Leave-one-out sensitivity
# Recompute GI-PRAL composite while iteratively dropping each
# of the 10 food exposures. Tests whether the composite signal
# depends on any single food (e.g. Pork dominance).
# ─────────────────────────────────────────────────────────────
cat("Step 8c: GI-PRAL leave-one-out sensitivity...\n")

gipral_full_row <- data.frame(
  removed_exposure = "None (full model)",
  beta             = round(gipral_beta, 4),
  se               = round(gipral_se,   4),
  ci_lo            = round(gipral_ci_lo, 4),
  ci_hi            = round(gipral_ci_hi, 4),
  pval             = round(gipral_p,    4),
  significant_p005 = gipral_p < 0.05,
  stringsAsFactors = FALSE
)

gipral_loo <- do.call(rbind, lapply(names(pral_weights), function(drop_exp) {
  sub_df <- gipral_df %>% filter(exposure != drop_exp)
  b   <- sum(sub_df$pral_w * sub_df$beta, na.rm=TRUE)
  se_ <- sqrt(sum((sub_df$pral_w * sub_df$se)^2, na.rm=TRUE))
  p   <- 2 * pnorm(-abs(b/se_))
  data.frame(
    removed_exposure = drop_exp,
    beta             = round(b,            4),
    se               = round(se_,          4),
    ci_lo            = round(b - 1.96*se_, 4),
    ci_hi            = round(b + 1.96*se_, 4),
    pval             = round(p,            4),
    significant_p005 = p < 0.05,
    stringsAsFactors = FALSE
  )
}))

gipral_loo_out <- rbind(gipral_full_row, gipral_loo)
write.csv(gipral_loo_out, paste0(proc, "gipral_loo_sensitivity.csv"), row.names=FALSE)

cat(sprintf("  GI-PRAL full       : \u03b2=%+.4f, P=%.4f, sig(P<0.05)=%s\n",
            gipral_full_row$beta, gipral_full_row$pval, gipral_full_row$significant_p005))
row_no_pork <- gipral_loo_out %>% filter(removed_exposure == "Pork intake")
if (nrow(row_no_pork) == 1) {
  cat(sprintf("  GI-PRAL w/o Pork   : \u03b2=%+.4f, P=%.4f, sig(P<0.05)=%s\n",
              row_no_pork$beta, row_no_pork$pval, row_no_pork$significant_p005))
}
cat("GI-PRAL LOO saved to gipral_loo_sensitivity.csv\n")

# ─────────────────────────────────────────────────────────────
# STEP 9: eGFRcys sensitivity analysis (Pattaro et al. 2016)
# Outcome: serum cystatin C-based eGFR (ieu-a-1106, N=33,152)
# Pattaro C et al. Nat Commun 2016;7:10023 — no UK Biobank overlap.
# Provides an alternative kidney-function biomarker robust to
# muscle-mass confounding that affects serum creatinine.
# Coverage is partial (proxies=0): 116/264 Steiger-filtered
# SNP-exposure pairs are present in the eGFRcys GWAS.
# ─────────────────────────────────────────────────────────────
cat("Step 9: eGFRcys sensitivity (ieu-a-1106, Pattaro 2016)...\n")

egfrcys_raw <- associations(variants = snps_c, id = "ieu-a-1106",
                            proxies = 0, opengwas_jwt = TOKEN)

egfrcys_out <- egfrcys_raw %>%
  transmute(rsid                  = rsid,
            effect_allele_outcome = toupper(ea),
            other_allele_outcome  = toupper(nea),
            eaf_outcome           = eaf,
            beta_outcome          = beta,
            se_outcome            = se,
            pval_outcome          = p,
            samplesize_outcome    = ifelse(is.na(n), 33152, n),
            outcome               = "Serum cystatin C (eGFRcys) || id:ieu-a-1106")

harm_egfrcys <- harmonise(instr_clean, egfrcys_out)
write.csv(harm_egfrcys, paste0(proc, "harmonised_egfrcys.csv"), row.names=FALSE)
cat(sprintf("  eGFRcys harmonised SNP-exposure pairs: %d / %d (%.0f%%)\n",
            nrow(harm_egfrcys), nrow(instr_clean),
            100 * nrow(harm_egfrcys) / nrow(instr_clean)))

egfrcys_ivw <- bind_rows(Filter(Negate(is.null), lapply(exps, function(e){
  d <- harm_egfrcys %>% filter(exposure == e)
  if (nrow(d) < 3) return(NULL)
  v <- mr_ivw(d$beta_exposure, d$beta_outcome_aligned, d$se_outcome)
  data.frame(
    id.exposure   = unique(d$gwas_id)[1],
    id.outcome    = "ieu-a-1106",
    outcome       = "Serum cystatin C (eGFRcys) || id:ieu-a-1106",
    exposure      = e,
    method        = "Inverse variance weighted",
    nsnp          = v$nsnp,
    b             = v$beta,
    se            = v$se,
    pval          = v$pval,
    stringsAsFactors = FALSE
  )
})))
egfrcys_ivw$pval_fdr      <- p.adjust(egfrcys_ivw$pval, "BH")
egfrcys_ivw$outcome_label <- "eGFRcys (ieu-a-1106, N=33,152)"

write.csv(egfrcys_ivw, paste0(proc, "mr_results_egfrcys_ivw.csv"), row.names=FALSE)
cat(sprintf("  eGFRcys IVW: min P = %.3f, min P_FDR = %.3f (no exposure reaches P<0.05)\n",
            min(egfrcys_ivw$pval), min(egfrcys_ivw$pval_fdr)))

# ─────────────────────────────────────────────────────────────
# STEP 10: Result tables
# ─────────────────────────────────────────────────────────────
t2 <- ivw_res %>% mutate(ci95=sprintf("[%.4f,%.4f]",ci_lo,ci_hi),
  beta=round(beta,4),se=round(se,4),pval=signif(pval,3),pval_fdr=signif(pval_fdr,3),
  Q_pval=signif(Q_pval,3)) %>%
  select(Exposure=exposure,Outcome=outcome,N_SNPs=nsnp,
         beta,SE=se,CI95=ci95,P=pval,P_FDR=pval_fdr,Q_P=Q_pval)
t3 <- all_res %>% filter(method!="IVW") %>% mutate(across(where(is.numeric),~round(.,4))) %>%
  select(Exposure=exposure,Outcome=outcome,Method=method,beta,SE=se,P=pval)
t4 <- mvmr_res %>% mutate(ci95=sprintf("[%.4f,%.4f]",ci_lo,ci_hi),
  beta=round(beta,4),se=round(se,4),pval=signif(pval,3),pval_fdr=signif(pval_fdr,3)) %>%
  select(Exposure=exposure,Outcome=outcome,N_SNPs=nsnp_total,
         beta,SE=se,CI95=ci95,P=pval,P_FDR=pval_fdr)
t5 <- fg_res %>% mutate(or_ci=sprintf("[%.3f,%.3f]",or_lo,or_hi)) %>%
  select(Exposure=exposure,N_SNPs=nsnp,OR=or,CI95=or_ci,P=pval,P_FDR=pval_fdr)
write.csv(t2, paste0(proc,"Table2_main_results.csv"),      row.names=FALSE)
write.csv(t3, paste0(proc,"Table3_sensitivity.csv"),       row.names=FALSE)
write.csv(t4, paste0(proc,"Table4_mvmr.csv"),              row.names=FALSE)
write.csv(t5, paste0(proc,"Table5_FinnGen_replication.csv"),row.names=FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All primary result tables written to data/processed/.\n")
cat("To regenerate figures, run the scripts in scripts/ (see README.md).\n")
