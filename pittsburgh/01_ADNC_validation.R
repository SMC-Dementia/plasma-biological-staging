#!/usr/bin/env Rscript
###############################################################################
#  01_ADNC_validation.R  --  Figure 5A
#  Plasma Stage 0-4 versus ADNC severity in the Pittsburgh autopsy cohort.
#
#  Primary model
#    ordered ADNC (Not/Low < Intermediate < High)
#      ~ ordinal plasma stage (0-4, one OR per stage) + age + sex
#        + blood-to-death interval
#    proportional-odds logistic regression (MASS::polr)
#
#  Secondary
#    same model with stage as a categorical factor: OR for each stage vs Stage 0
#    Spearman correlation between stage and ADNC
#
#  Sensitivity
#    + plate fixed effects
#    restricted to blood-to-death interval <= 10 years
#
#  Source ADNC labels are harmonised for reporting:
#    "Low/Not ADNC" -> Not/Low,  "Int ADNC" -> Intermediate,  "Sev ADNC" -> High
#
#  Usage
#    Rscript 01_ADNC_validation.R --input plasma.csv --outdir Fig5A
###############################################################################
suppressPackageStartupMessages(library(MASS))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "validation_common.R"))

if (has_flag("--help") || !length(ARGS)) { cat("Figure 5A: plasma stage versus ADNC severity\n", COMMON_HELP); quit(status = 0) }
opt <- common_opts()
if (is.null(opt$input)) stop("--input is required.")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

banner("Figure 5A  Plasma stage versus ADNC severity")
raw  <- read_any(opt$input, opt$sheet)
prep <- prepare_index(raw)
idx  <- prep$idx
say(sprintf("Fixed cutoffs: BD-p-tau217 %.9f; eMTBR-tau %s", CUT_P217, paste(c(CUT_EMTBR_1, CUT_EMTBR_2, CUT_EMTBR_3), collapse = " / ")))
write_t1_info(prep, opt$outdir)

## --- ADNC ------------------------------------------------------------------
adnc_col <- resolve_col(raw, "adnc", required = TRUE, contains_fallback = "adnc")
harmonise_adnc <- function(v) {
  ## Text labels
  z <- norm_name(v)
  out <- rep(NA_character_, length(z))
  out[grepl("^(lownot|notlow|not|none|no|low|lowadnc|notadnc|noneadnc|notad|lowad)", z)] <- "Not/Low"
  out[grepl("^(int|intermediate)", z)] <- "Intermediate"
  out[grepl("^(sev|severe|high)", z)] <- "High"
  ## Numeric codes: decide by how many distinct levels are present.
  ##   4 levels 0-3  (NIA-AA: Not / Low / Intermediate / High)  -> 0,1 = Not/Low, 2 = Int, 3 = High
  ##   3 levels 0-2                                            -> 0 = Not/Low, 1 = Int, 2 = High
  num <- suppressWarnings(as.numeric(as.character(v)))
  isnum <- !is.na(num) & is.na(out)
  if (any(isnum)) {
    lv <- sort(unique(num[isnum]))
    if (all(lv %in% 0:3) && max(lv) == 3) {
      out[isnum & num <= 1] <- "Not/Low"; out[isnum & num == 2] <- "Intermediate"; out[isnum & num == 3] <- "High"
      message("NOTE: ADNC is coded 0-3; interpreted as Not(0)/Low(1)/Intermediate(2)/High(3). Check Fig5A_ADNC_label_audit.csv.")
    } else if (all(lv %in% 0:2)) {
      out[isnum & num == 0] <- "Not/Low"; out[isnum & num == 1] <- "Intermediate"; out[isnum & num == 2] <- "High"
      message("NOTE: ADNC is coded 0-2; interpreted as Not/Low(0)/Intermediate(1)/High(2). Check Fig5A_ADNC_label_audit.csv.")
    } else if (all(lv %in% 1:4)) {
      out[isnum & num <= 2] <- "Not/Low"; out[isnum & num == 3] <- "Intermediate"; out[isnum & num == 4] <- "High"
      message("NOTE: ADNC is coded 1-4; interpreted as Not(1)/Low(2)/Intermediate(3)/High(4). Check Fig5A_ADNC_label_audit.csv.")
    }
  }
  out
}
idx$adnc_src <- as.character(raw[[adnc_col]][idx$raw_row])
idx$adnc     <- factor(harmonise_adnc(idx$adnc_src), levels = c("Not/Low", "Intermediate", "High"), ordered = TRUE)
audit <- unique(data.frame(source_label = idx$adnc_src, harmonised = as.character(idx$adnc), stringsAsFactors = FALSE))
write_csv_safe(audit, file.path(opt$outdir, "Fig5A_ADNC_label_audit.csv"))
if (any(is.na(audit$harmonised) & !is.na(audit$source_label)))
  warning("Some ADNC labels were not recognised and are excluded: ",
          paste(audit$source_label[is.na(audit$harmonised)], collapse = ", "))

ana <- idx[idx$stage %in% STAGES & !is.na(idx$adnc), ]
if (!nrow(ana)) stop("No stageable participants with a recognised ADNC label.")
say(sprintf("Analysed: %d participants (Stage 0-4 with ADNC)", nrow(ana)))

## --- proportions table -----------------------------------------------------
tab <- table(factor(ana$stage, levels = STAGES), ana$adnc)
prop <- prop.table(tab, 1) * 100
prop_df <- do.call(rbind, lapply(STAGES, function(s) data.frame(
  stage = s, stage_n_with_ADNC = sum(tab[s, ]), ADNC = colnames(tab),
  n = as.integer(tab[s, ]), percent_within_stage = round(as.numeric(prop[s, ]), 1))))
write_csv_safe(prop_df, file.path(opt$outdir, "Fig5A_ADNC_proportions.csv"))

## --- models ----------------------------------------------------------------
ana$age_c      <- ana$age - mean(ana$age, na.rm = TRUE)
ana$female_c   <- ana$female - mean(ana$female, na.rm = TRUE)
ana$interval_c <- ana$b2d_years - mean(ana$b2d_years, na.rm = TRUE)

fit_polr <- function(d, extra = NULL, label) {
  f <- as.formula(paste("adnc ~ stage_num + age_c + female_c + interval_c", if (!is.null(extra)) paste("+", extra) else ""))
  d <- d[complete.cases(d[, c("adnc","stage_num","age_c","female_c","interval_c")]), ]
  fit <- tryCatch(polr(f, data = d, Hess = TRUE, method = "logistic"), error = function(e) NULL)
  if (is.null(fit)) return(list(fit = NULL, row = data.frame(model = label, n = nrow(d), note = "did not converge")))
  cf <- polr_coef(fit, "stage_num")
  list(fit = fit, row = data.frame(model = label, n = nrow(d), beta_log_odds_per_stage = cf["beta"], se = cf["se"],
                                   OR_per_stage = cf["OR"], OR_95CI_low = cf["OR_low"], OR_95CI_high = cf["OR_high"],
                                   P = cf["P"], row.names = NULL))
}

prim <- fit_polr(ana, NULL, "Primary")
rows <- list(prim$row)
if (!all(is.na(ana$plate)) && length(unique(na.omit(ana$plate))) > 1) {
  ana$plate_f <- factor(ana$plate)
  rows[[length(rows) + 1]] <- fit_polr(ana, "plate_f", "Sensitivity: + plate")$row
}
sub10 <- ana[!is.na(ana$b2d_years) & ana$b2d_years <= 10, ]
if (nrow(sub10) >= 50) rows[[length(rows) + 1]] <- fit_polr(sub10, NULL, "Sensitivity: interval <= 10 y")$row
model_tbl <- do.call(rbind, lapply(rows, function(r) { r[setdiff(c("model","n","beta_log_odds_per_stage","se","OR_per_stage","OR_95CI_low","OR_95CI_high","P","note"), names(r))] <- NA; r }))
write_csv_safe(model_tbl, file.path(opt$outdir, "Fig5A_ADNC_model_results.csv"))

## secondary: categorical stage
ana$stage_cat <- droplevels(ana$stage_cat)
cat_fit <- tryCatch(polr(adnc ~ stage_cat + age_c + female_c + interval_c, data = ana, Hess = TRUE), error = function(e) NULL)
if (!is.null(cat_fit)) {
  cat_rows <- do.call(rbind, lapply(levels(ana$stage_cat)[-1], function(s) {
    cf <- polr_coef(cat_fit, paste0("stage_cat", s))
    data.frame(contrast = paste(s, "vs Stage 0"), OR = cf["OR"], OR_95CI_low = cf["OR_low"],
               OR_95CI_high = cf["OR_high"], P = cf["P"], row.names = NULL)
  }))
  write_csv_safe(cat_rows, file.path(opt$outdir, "Fig5A_ADNC_categorical_stage.csv"))
}
sp <- suppressWarnings(cor.test(ana$stage_num, as.integer(ana$adnc), method = "spearman"))
write_csv_safe(data.frame(Spearman_rho = unname(sp$estimate), P = sp$p.value, n = nrow(ana)),
               file.path(opt$outdir, "Fig5A_ADNC_spearman.csv"))

if (!is.null(prim$fit)) {
  sink(file.path(opt$outdir, "Fig5A_ADNC_model_summary.txt")); print(summary(prim$fit)); sink()
}

## --- figure ----------------------------------------------------------------
pr <- prim$row
save_fig(file.path(opt$outdir, "Fig5A_ADNC_validation"), function() {
  par(mar = c(4.5, 4.5, 1.5, 9), xpd = TRUE)
  m <- t(as.matrix(prop))[, STAGES]
  bp <- barplot(m, col = c("white", grey(0.68), grey(0.18)), border = "black", lwd = 0.8,
                names.arg = sprintf("%s\n(n=%d)", STAGES, rowSums(tab)[STAGES]),
                ylim = c(0, 100), ylab = "Participants (%)", xlab = "Plasma stage", cex.names = 0.8, las = 1)
  legend("topright", inset = c(-0.48, 0), legend = rev(c("Not/Low ADNC", "Intermediate ADNC", "High ADNC")),
         fill = rev(c("white", grey(0.68), grey(0.18))), bty = "n", cex = 0.85)
  if (!is.na(pr$OR_per_stage))
    text(max(bp) + 0.85, 55, adj = c(0, 1), cex = 0.78,
         sprintf("Adjusted OR per stage\n  %.2f\n95%% CI %.2f to %.2f\nP = %s", pr$OR_per_stage, pr$OR_95CI_low, pr$OR_95CI_high, fp(pr$P)))
}, width = 6.4, height = 4.2)

write_csv_safe(data.frame(index_participants = nrow(idx), stageable_0_to_4 = sum(idx$stage %in% STAGES),
                          stageable_with_ADNC = nrow(ana), median_blood_to_death_years = round(median(ana$b2d_years, na.rm = TRUE), 2),
                          ADNC_source_column = adnc_col, t1_cutoff = CUT_P217),
               file.path(opt$outdir, "Fig5A_cohort_summary.csv"))

say(""); print(model_tbl, row.names = FALSE)
say(sprintf("\nFigure 5A complete: %s", normalizePath(opt$outdir)))
