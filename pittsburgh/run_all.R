#!/usr/bin/env Rscript
###############################################################################
#  run_all.R  --  runs Figure 5A, 5B and 5C and writes a summary
#
#    <outdir>/
#      Fig5A_* Fig5B_* Fig5C_* T1_cutoff_and_stage_counts.csv
#      SUMMARY.txt          headline numbers
#      SUMMARY.csv          same, machine-readable
#      run_log.txt          full console output
#
#  Usage
#    Rscript run_all.R --input plasma.csv --clinical mmse.csv --outdir Fig5_output
#
#  Optional
#    --skip-mmse                skip Figure 5C
#    any --*-col override accepted by the individual scripts
###############################################################################
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "validation_common.R"))

if (has_flag("--help") || !length(ARGS)) {
  cat("Run Figure 5A/B/C.\n",
      "  --input FILE  --clinical FILE  --outdir DIR  [--skip-mmse]\n"); quit(status = 0) }

input    <- get_arg("--input");    if (is.null(input)) stop("--input is required.")
clinical <- get_arg("--clinical")
outdir   <- get_arg("--outdir", "Fig5_output")
skip_mmse<- has_flag("--skip-mmse") || is.null(clinical)
passthru <- character(0)
for (k in c(names(COLUMN_MAP), "sheet", "clinical-sheet", "mmse-date-col", "mmse-col", "clinical-id-col")) {
  f <- if (k %in% names(COLUMN_MAP)) paste0("--", gsub("_", "-", k), "-col") else paste0("--", k)
  v <- get_arg(f); if (!is.null(v)) passthru <- c(passthru, f, v)
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(outdir, "run_log.txt"); if (file.exists(logf)) file.remove(logf)

run_one <- function(script, args) {
  cmd <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(cmd, c(shQuote(file.path(here, script)), shQuote(args)), stdout = TRUE, stderr = TRUE))
  cat(sprintf("\n$ Rscript %s %s\n", script, paste(args, collapse = " ")), out, sep = "\n", file = logf, append = TRUE)
  status <- attr(out, "status"); if (is.null(status)) status <- 0
  list(ok = status == 0, out = out)
}
grab <- function(path, ...) if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE, ...) else NULL

banner("Figure 5 external validation")
say("input    : ", input); if (!skip_mmse) say("clinical : ", clinical); say("outdir   : ", normalizePath(outdir))

base <- c("--input", input, "--outdir", outdir, passthru)
a <- run_one("01_ADNC_validation.R", base);       say(sprintf("  5A ADNC       %s", if (a$ok) "done" else "FAILED (see run_log.txt)"))
b <- run_one("02_34protein_replication.R", base); say(sprintf("  5B proteins   %s", if (b$ok) "done" else "FAILED (see run_log.txt)"))
cres <- NULL
if (!skip_mmse) { cres <- run_one("03_MMSE_validation.R", c(base, "--clinical", clinical)); say(sprintf("  5C MMSE       %s", if (cres$ok) "done" else "FAILED (see run_log.txt)")) }

t1 <- grab(file.path(outdir, "T1_cutoff_and_stage_counts.csv")); tv <- function(f) if (is.null(t1)) NA else t1$value[t1$field == f]
am <- grab(file.path(outdir, "Fig5A_ADNC_model_results.csv"));  ap <- if (!is.null(am)) am[am$model == "Primary", ] else NULL
bs <- grab(file.path(outdir, "Fig5B_replication_summary.csv")); bx <- if (!is.null(bs)) bs[bs$set == "all analysed proteins", ] else NULL
cg <- grab(file.path(outdir, "Fig5C_MMSE_global_test.csv"));    cp <- if (!is.null(cg)) cg[cg$set == "Primary", ] else NULL
S <- data.frame(
  t1_cutoff = tv("t1_cutoff"), n_index = tv("n_index"),
  n_stage0 = tv("n_Stage 0"), n_stage1 = tv("n_Stage 1"), n_stage2 = tv("n_Stage 2"), n_stage3 = tv("n_Stage 3"), n_stage4 = tv("n_Stage 4"), n_unstageable = tv("n_Unstageable"),
  ADNC_n = if (!is.null(ap)) ap$n else NA, ADNC_OR_per_stage = if (!is.null(ap)) ap$OR_per_stage else NA,
  ADNC_OR_low = if (!is.null(ap)) ap$OR_95CI_low else NA, ADNC_OR_high = if (!is.null(ap)) ap$OR_95CI_high else NA, ADNC_P = if (!is.null(ap)) ap$P else NA,
  proteins_n = if (!is.null(bx)) bx$n_proteins else NA, proteins_concordant_pct = if (!is.null(bx)) bx$direction_concordance_percent else NA,
  proteins_spearman_rho = if (!is.null(bx)) bx$Spearman_rho else NA, proteins_spearman_P = if (!is.null(bx)) bx$Spearman_P else NA,
  proteins_replicated_q05 = if (!is.null(bx)) bx$n_replicated_q_lt_0.05_same_direction else NA,
  MMSE_n = if (!is.null(cp)) cp$n_participants else NA, MMSE_LR_chi2 = if (!is.null(cp)) cp$LR_chi2 else NA, MMSE_P = if (!is.null(cp)) cp$P else NA,
  stringsAsFactors = FALSE)
write_csv_safe(S, file.path(outdir, "SUMMARY.csv"))

f <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = d))
txt <- c(
  "Figure 5 external validation  --  summary",
  sprintf("input: %s   run: %s", basename(input), format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  sprintf("%-32s %12s", "T1 cutoff (log2 NPQ)",      f(S$t1_cutoff, 9)),
  sprintf("%-32s %12s", "Stage 0 / 1 / 2 / 3 / 4 (n)", paste(S$n_stage0, S$n_stage1, S$n_stage2, S$n_stage3, S$n_stage4, sep = "/")),
  sprintf("%-32s %12s", "Unstageable (n)",             S$n_unstageable),
  "",
  "5A  ADNC severity ~ stage (ordinal logistic, adjusted)",
  sprintf("%-32s %12s", "  OR per stage",              f(S$ADNC_OR_per_stage, 2)),
  sprintf("%-32s %12s", "  95% CI low",                f(S$ADNC_OR_low, 2)),
  sprintf("%-32s %12s", "  95% CI high",               f(S$ADNC_OR_high, 2)),
  sprintf("%-32s %12s", "  P",                         fp(as.numeric(S$ADNC_P))),
  "",
  "5B  34-protein replication",
  sprintf("%-32s %12s", "  proteins analysed",         S$proteins_n),
  sprintf("%-32s %12s", "  direction concordance (%)", f(S$proteins_concordant_pct, 1)),
  sprintf("%-32s %12s", "  Spearman rho (betas)",      f(S$proteins_spearman_rho, 2)),
  sprintf("%-32s %12s", "  replicated (q<0.05, same dir)", S$proteins_replicated_q05),
  "",
  "5C  MMSE ~ time x stage (mixed model)",
  sprintf("%-32s %12s", "  participants",              S$MMSE_n),
  sprintf("%-32s %12s", "  stage x time LR chi2(4)",   f(S$MMSE_LR_chi2, 1)),
  sprintf("%-32s %12s", "  P",                         fp(as.numeric(S$MMSE_P)))
)
writeLines(txt, file.path(outdir, "SUMMARY.txt"))
cat("\n", paste(txt, collapse = "\n"), "\n\n")
say("Done. Output: ", normalizePath(outdir))
