#!/usr/bin/env Rscript
###############################################################################
#  03_MMSE_validation.R  --  Figure 5C
#  Longitudinal MMSE trajectories by baseline plasma stage.
#
#  The SAME closest-to-death index plasma sample and stage used in Figure 5A/B
#  is used here.
#
#  MMSE observations retained
#    within 1 year before the index plasma, and all observations after it
#  Participant eligibility
#    at least 2 retained MMSE assessments, at least 1 on/after the index plasma
#
#  Primary model (nlme::lme, ML)
#    MMSE ~ time x categorical stage + age + sex
#    participant-specific random intercept and slope
#  Global stage-by-time interaction: likelihood-ratio test (4 df) against the
#  model without the interaction.
#
#  Sensitivity: participants whose first retained MMSE is >= 10 (floor effect)
#
#  Usage
#    Rscript 03_MMSE_validation.R --input plasma.csv --clinical mmse.csv --outdir Fig5C
###############################################################################
suppressPackageStartupMessages(library(nlme))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "validation_common.R"))

if (has_flag("--help") || !length(ARGS)) {
  cat("Figure 5C: longitudinal MMSE by plasma stage\n", COMMON_HELP,
      "  --clinical FILE         longitudinal MMSE file (participant id, date, score)\n",
      "  --clinical-sheet NAME   worksheet, for .xlsx\n",
      "  --mmse-date-col / --mmse-col   override clinical column detection\n"); quit(status = 0) }
opt <- common_opts()
opt$clinical <- get_arg("--clinical"); opt$clin_sheet <- get_arg("--clinical-sheet")
if (is.null(opt$input) || is.null(opt$clinical)) stop("--input and --clinical are both required.")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
PRE_INDEX_YEARS <- 1.0; PLOT_YEARS <- 3.0

banner("Figure 5C  Longitudinal MMSE by plasma stage")
raw  <- read_any(opt$input, opt$sheet)
prep <- prepare_index(raw)
idx  <- prep$idx[prep$idx$stage %in% STAGES, ]
say(sprintf("Fixed cutoffs: BD-p-tau217 %.9f; eMTBR-tau %s", CUT_P217, paste(c(CUT_EMTBR_1, CUT_EMTBR_2, CUT_EMTBR_3), collapse = " / ")))
write_t1_info(prep, opt$outdir)

## --- clinical file ---------------------------------------------------------
cl <- read_any(opt$clinical, opt$clin_sheet)
cid <- get_arg("--clinical-id-col"); if (is.null(cid)) cid <- resolve_col(cl, "id", TRUE)
find_col <- function(df, override, cands, what) {
  if (!is.null(override)) { if (!override %in% names(df)) stop(sprintf("%s column '%s' not found.", what, override)); return(override) }
  nm <- setNames(names(df), norm_name(names(df)))
  for (a in cands) if (norm_name(a) %in% names(nm)) return(unname(nm[norm_name(a)]))
  stop(sprintf("%s column not found in the clinical file. Columns: %s\nPass --%s-col.", what, paste(names(df), collapse = ", "), tolower(gsub(" ", "-", what))))
}
cdate <- find_col(cl, get_arg("--mmse-date-col"), c("MMSE_date","mmse_date","assessment_date","visit_date","exam_date","date"), "MMSE date")
cscore<- find_col(cl, get_arg("--mmse-col"),      c("MMSE","mmse","MMSE_score","mmse_score","mmse_total"), "MMSE")
c <- data.frame(id = norm_id(cl[[cid]]), mmse_date = parse_date(cl[[cdate]]),
                mmse = suppressWarnings(as.numeric(cl[[cscore]])), stringsAsFactors = FALSE)
c <- c[complete.cases(c), ]

info <- idx[, c("id","sample_date","stage","stage_num","age","female")]
d <- merge(c, info, by = "id")
d$time <- as.numeric(d$mmse_date - d$sample_date) / 365.25
d <- d[d$time >= -PRE_INDEX_YEARS, ]
n_all  <- table(d$id); n_post <- table(d$id[d$time >= 0])
keep <- names(n_all)[n_all >= 2 & names(n_all) %in% names(n_post)]
d <- d[d$id %in% keep, ]
n_cov_na <- sum(is.na(d$age) | is.na(d$female))
if (n_cov_na) { message(sprintf("NOTE: %d assessment rows dropped for missing age or sex.", n_cov_na)); d <- d[!is.na(d$age) & !is.na(d$female), ] }
if (length(unique(d$id)) < 30) stop(sprintf("Only %d participants meet the MMSE criteria.", length(unique(d$id))))
## Stages with no eligible participants are dropped from the factor; the
## interaction test then has (number of stages present - 1) degrees of freedom.
present <- STAGES[STAGES %in% unique(d$stage)]
if (length(present) < 2) stop("Fewer than two plasma stages have eligible MMSE data.")
if (length(present) < 5) message("NOTE: stage(s) with no eligible MMSE participants: ", paste(setdiff(STAGES, present), collapse = ", "))
d$stage_cat <- factor(d$stage, levels = present)
d$age_c <- d$age - mean(d$age); d$female_c <- d$female - mean(d$female)
d <- d[order(d$id, d$time), ]
IDF <- length(present) - 1

## --- models ----------------------------------------------------------------
fit_lme <- function(f, dd) {
  ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200, returnObject = TRUE)
  for (rf in list(~ time | id, ~ 1 | id)) {
    fit <- tryCatch(lme(f, data = dd, random = rf, method = "ML", control = ctrl), error = function(e) NULL)
    if (!is.null(fit)) { attr(fit, "random_used") <- deparse(rf); return(fit) }
  }
  stop("lme did not converge with random slope or random intercept.")
}
run_mmse <- function(dd, label) {
  dd$stage_cat <- droplevels(dd$stage_cat)
  full <- fit_lme(mmse ~ time * stage_cat + age_c + female_c, dd)
  red  <- fit_lme(mmse ~ time + stage_cat + age_c + female_c, dd)
  lv <- levels(droplevels(dd$stage_cat)); kdf <- length(lv) - 1
  lr <- max(0, 2 * (logLik(full) - logLik(red))); p_lr <- pchisq(lr, kdf, lower.tail = FALSE)
  fe <- fixef(full); V <- vcov(full)
  slopes <- do.call(rbind, lapply(lv, function(s) {
    w <- setNames(rep(0, length(fe)), names(fe)); w["time"] <- 1
    if (s != lv[1]) { term <- paste0("time:stage_cat", s); if (!term %in% names(fe)) stop("missing ", term); w[term] <- 1 }
    est <- sum(w * fe); se <- sqrt(as.numeric(t(w) %*% V %*% w))
    data.frame(set = label, stage = s, n_participants = length(unique(dd$id[dd$stage == s])),
               annual_MMSE_slope = est, se = se, ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
               P = 2 * pnorm(-abs(est / se)))
  }))
  list(full = full, slopes = slopes,
       glob = data.frame(set = label, n_participants = length(unique(dd$id)), n_assessments = nrow(dd),
                         stages_present = paste(lv, collapse = ";"),
                         LR_chi2 = lr, df = kdf, P = p_lr, random_effects = attr(full, "random_used")))
}

prim <- run_mmse(d, "Primary")
sl <- prim$slopes; gl <- prim$glob
first <- d[!duplicated(d$id), c("id","mmse")]
keep10 <- first$id[first$mmse >= 10]
if (length(keep10) >= 30 && length(keep10) < length(unique(d$id))) {
  s10 <- tryCatch(run_mmse(d[d$id %in% keep10, ], "Sensitivity: first MMSE >= 10"), error = function(e) NULL)
  if (!is.null(s10)) { sl <- rbind(sl, s10$slopes); gl <- rbind(gl, s10$glob) }
}
write_csv_safe(sl, file.path(opt$outdir, "Fig5C_MMSE_stage_slopes.csv"))
write_csv_safe(gl, file.path(opt$outdir, "Fig5C_MMSE_global_test.csv"))
cnt <- data.frame(stage = STAGES,
                  participants = sapply(STAGES, function(s) length(unique(d$id[d$stage == s]))),
                  assessments  = sapply(STAGES, function(s) sum(d$stage == s)))
write_csv_safe(cnt, file.path(opt$outdir, "Fig5C_MMSE_counts_by_stage.csv"))
sink(file.path(opt$outdir, "Fig5C_MMSE_model_summary.txt")); print(summary(prim$full)); cat(sprintf("\nGlobal stage-by-time LRT: chi2(%d) = %.3f, P = %s\n", gl$df[1], gl$LR_chi2[1], fp(gl$P[1]))); sink()

## --- predicted trajectories ------------------------------------------------
fe <- fixef(prim$full); V <- vcov(prim$full)
grid <- seq(-PRE_INDEX_YEARS, PLOT_YEARS, length.out = 81)
pred <- do.call(rbind, lapply(present, function(s) {
  X <- t(sapply(grid, function(t) {
    w <- setNames(rep(0, length(fe)), names(fe)); w["(Intercept)"] <- 1; w["time"] <- t
    if (s != present[1]) { w[paste0("stage_cat", s)] <- 1; w[paste0("time:stage_cat", s)] <- t }
    w }))
  mu <- as.numeric(X %*% fe); se <- sqrt(rowSums((X %*% V) * X))
  data.frame(stage = s, time_years = grid, predicted_MMSE = mu, ci_low = mu - 1.96 * se, ci_high = mu + 1.96 * se)
}))
write_csv_safe(pred, file.path(opt$outdir, "Fig5C_MMSE_predicted_trajectories.csv"))

save_fig(file.path(opt$outdir, "Fig5C_MMSE_validation"), function() {
  par(mar = c(4.5, 4.5, 1.5, 6.5), xpd = FALSE)
  plot(NA, xlim = c(-PRE_INDEX_YEARS, PLOT_YEARS), ylim = range(c(pred$ci_low, pred$ci_high)),
       xlab = "Years from index plasma", ylab = "MMSE", las = 1)
  abline(v = 0, lty = 2, col = grey(0.55))
  lty <- c(1, 2, 4, 3, 5)
  for (i in match(present, STAGES)) {
    q <- pred[pred$stage == STAGES[i], ]
    polygon(c(q$time_years, rev(q$time_years)), c(q$ci_low, rev(q$ci_high)), col = adjustcolor(grey_col(i), 0.08), border = NA)
    lines(q$time_years, q$predicted_MMSE, col = grey_col(i), lwd = 1.7, lty = lty[i])
  }
  par(xpd = TRUE)
  ii <- match(present, STAGES)
  legend(PLOT_YEARS + 0.15, max(pred$ci_high), legend = present, col = sapply(ii, grey_col), lwd = 1.7, lty = lty[ii], bty = "n", cex = 0.8)
  par(xpd = FALSE)
  legend("topleft", bty = "n", cex = 0.75, legend = sprintf("Stage x time: chi2(%d) = %.1f, P = %s", gl$df[1], gl$LR_chi2[1], fp(gl$P[1])))
}, width = 6.2, height = 4.2)

say(""); print(sl[sl$set == "Primary", c("stage","n_participants","annual_MMSE_slope","ci_low","ci_high","P")], row.names = FALSE)
say(sprintf("\nGlobal stage-by-time LRT: chi2(%d) = %.2f, P = %s   [%s]", gl$df[1], gl$LR_chi2[1], fp(gl$P[1]), gl$random_effects[1]))
say(sprintf("Figure 5C complete: %s", normalizePath(opt$outdir)))
