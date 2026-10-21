#!/usr/bin/env Rscript
###############################################################################
#  02_34protein_replication.R  --  Figure 5B
#  Replication of the 34 prespecified K-ROAD stage-associated proteins.
#
#  For each protein
#    standardised protein (z within Pittsburgh)
#      ~ ordinal plasma stage (0-4) + age + sex
#    ordinary least squares with HC3 heteroscedasticity-robust standard errors
#
#  This matches the K-ROAD secondary ordinal-stage model from which the
#  reference betas below were taken. The Pittsburgh beta (SD per stage) is
#  compared with the K-ROAD beta for direction and magnitude.
#
#  Sensitivity: + plate fixed effects
#
#  Five of the 34 are rare-case targets on the Neuro 220 Panel (mHTT-exon1,
#  pQ-ATXN3, pRAB29-T71, pRAB12-S106, pPRKN-S65). They sit near the detection
#  limit in most people, so non-replication for these is expected. They are
#  flagged in the results and the summary is reported with and without them.
#
#  Usage
#    Rscript 02_34protein_replication.R --input plasma.csv --outdir Fig5B
###############################################################################
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "validation_common.R"))

if (has_flag("--help") || !length(ARGS)) { cat("Figure 5B: 34-protein stage-effect replication\n", COMMON_HELP); quit(status = 0) }
opt <- common_opts()
if (is.null(opt$input)) stop("--input is required.")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## --- K-ROAD reference betas (SD per stage, ordinal-stage model, HC3) --------
KROAD_BETA <- c(
  GFAP=0.441, SPP1=0.134, IL33=0.108, MMP9=0.079, MOG=0.091, TREM1=0.064, IL1B=0.063,
  NfL=0.221, STMN2=0.097, MAG=0.071, ACHE=0.297, GAP43=0.224, SNAP25=0.113, BASP1=0.096,
  APOA2=-0.099, NTRK3=0.085, PAFAH1B3=0.035, NPTX2=0.081, NCAM1=0.084, VGF=-0.051,
  PPY=0.073, APOA1=0.062, EFEMP1=0.044, BACE1=0.195, `mHTT-exon1`=0.092, `pQ-ATXN3`=0.039,
  `Ab42`=-0.055, `pRAB29-T71`=0.109, `pRAB12-S106`=0.088, PSME1=0.063, GBA=0.076,
  `pPRKN-S65`=0.114, GLP1R=0.073, VSNL1=0.124
)
RARE_CASE <- c("mHTT-exon1", "pQ-ATXN3", "pRAB29-T71", "pRAB12-S106", "pPRKN-S65")

## If a protein column is not detected automatically, put the column name here.
PROTEIN_COLUMN_MAP <- list()

protein_aliases <- function(nm) {
  a <- c(nm, gsub("-", "_", nm), paste0("neuro220_", nm), paste0("neuro220_", gsub("-", "_", nm)), paste0("NPQ_", nm))
  if (nm == "Ab42") a <- c(a, "Abeta42", "AB42", "A_beta_42", "amyloid_beta_42", "Aβ42")
  if (nm == "NfL")  a <- c(a, "NEFL", "NFL")
  unique(a)
}
resolve_protein <- function(df, nm) {
  if (!is.null(PROTEIN_COLUMN_MAP[[nm]])) return(PROTEIN_COLUMN_MAP[[nm]])
  nmap <- setNames(names(df), norm_name(names(df)))
  for (a in protein_aliases(nm)) if (norm_name(a) %in% names(nmap)) return(unname(nmap[norm_name(a)]))
  pats <- norm_name(protein_aliases(nm))
  hits <- names(df)[vapply(norm_name(names(df)), function(c) any(endsWith(c, pats)), logical(1))]
  if (length(hits) == 1) hits else NA_character_
}

banner("Figure 5B  34-protein stage-effect replication")
raw  <- read_any(opt$input, opt$sheet)
prep <- prepare_index(raw)
idx  <- prep$idx
say(sprintf("Fixed cutoffs: BD-p-tau217 %.9f; eMTBR-tau %s", CUT_P217, paste(c(CUT_EMTBR_1, CUT_EMTBR_2, CUT_EMTBR_3), collapse = " / ")))
write_t1_info(prep, opt$outdir)

d <- idx[idx$stage %in% STAGES, ]
d$age_c <- d$age - mean(d$age, na.rm = TRUE); d$female_c <- d$female - mean(d$female, na.rm = TRUE)
has_plate <- !all(is.na(d$plate)) && length(unique(na.omit(d$plate))) > 1
if (has_plate) d$plate_f <- factor(d$plate)

fit_beta <- function(d, y, plate = FALSE) {
  d$y <- suppressWarnings(as.numeric(y))
  d <- d[complete.cases(d[, c("y","stage_num","age_c","female_c")]), ]
  if (nrow(d) < 30) return(NULL)
  s <- sd(d$y); if (!is.finite(s) || s <= 0) return(NULL)
  d$yz <- (d$y - mean(d$y)) / s
  f <- if (plate) yz ~ stage_num + age_c + female_c + plate_f else yz ~ stage_num + age_c + female_c
  fit <- lm(f, data = d)
  V <- tryCatch(vcov_hc3(fit), error = function(e) vcov(fit))
  b <- coef(fit)["stage_num"]; se <- sqrt(V["stage_num", "stage_num"])
  data.frame(n = nrow(d), pittsburgh_beta = unname(b), se = unname(se),
             ci_low = unname(b - 1.96 * se), ci_high = unname(b + 1.96 * se),
             P = unname(2 * pnorm(-abs(b / se))))
}

res <- list(); sens <- list(); resolution <- list()
for (nm in names(KROAD_BETA)) {
  col <- resolve_protein(raw, nm)
  resolution[[nm]] <- data.frame(protein = nm, resolved_column = ifelse(is.na(col), "", col), resolved = !is.na(col))
  if (is.na(col)) next
  y <- raw[[col]][d$raw_row]
  r <- fit_beta(d, y, FALSE); if (is.null(r)) next
  r$protein <- nm; r$resolved_column <- col; r$KROAD_beta <- KROAD_BETA[[nm]]
  r$rare_case_target <- nm %in% RARE_CASE
  r$same_direction <- sign(r$pittsburgh_beta) == sign(r$KROAD_beta)
  res[[nm]] <- r
  if (has_plate) { rs <- fit_beta(d, y, TRUE); if (!is.null(rs)) { rs$protein <- nm; rs$KROAD_beta <- KROAD_BETA[[nm]]; rs$same_direction <- sign(rs$pittsburgh_beta) == sign(rs$KROAD_beta); sens[[nm]] <- rs } }
}
resolution <- do.call(rbind, resolution)
write_csv_safe(resolution, file.path(opt$outdir, "Fig5B_protein_column_resolution.csv"))
unres <- resolution$protein[!resolution$resolved]
if (length(unres)) say(sprintf("  NOTE: %d protein(s) not detected: %s\n        Add them to PROTEIN_COLUMN_MAP in 02_34protein_replication.R", length(unres), paste(unres, collapse = ", ")))
if (length(res) < 10) stop(sprintf("Only %d of 34 proteins could be analysed. See Fig5B_protein_column_resolution.csv.", length(res)))

res <- do.call(rbind, res); rownames(res) <- NULL
res$q <- p.adjust(res$P, "BH")
res <- res[order(-res$KROAD_beta), c("protein","resolved_column","rare_case_target","n","KROAD_beta","pittsburgh_beta","se","ci_low","ci_high","P","q","same_direction")]
write_csv_safe(res, file.path(opt$outdir, "Fig5B_34protein_results.csv"))
if (length(sens)) { s <- do.call(rbind, sens); rownames(s) <- NULL; s$q <- p.adjust(s$P, "BH"); write_csv_safe(s, file.path(opt$outdir, "Fig5B_34protein_plate_adjusted_sensitivity.csv")) }

summarise <- function(r, label) {
  nc <- sum(r$same_direction); nt <- nrow(r)
  pe <- suppressWarnings(cor.test(r$KROAD_beta, r$pittsburgh_beta, method = "pearson"))
  sp <- suppressWarnings(cor.test(r$KROAD_beta, r$pittsburgh_beta, method = "spearman"))
  data.frame(set = label, n_proteins = nt, n_same_direction = nc, direction_concordance_percent = round(100 * nc / nt, 1),
             binomial_P_vs_50pct = binom.test(nc, nt, 0.5, alternative = "greater")$p.value,
             n_replicated_q_lt_0.05_same_direction = sum(r$same_direction & r$q < 0.05),
             Pearson_r = unname(pe$estimate), Pearson_P = pe$p.value,
             Spearman_rho = unname(sp$estimate), Spearman_P = sp$p.value)
}
summ <- rbind(summarise(res, "all analysed proteins"),
              summarise(res[!res$rare_case_target, ], "excluding rare-case targets"))
write_csv_safe(summ, file.path(opt$outdir, "Fig5B_replication_summary.csv"))

## --- figure ----------------------------------------------------------------
save_fig(file.path(opt$outdir, "Fig5B_34protein_replication"), function() {
  par(mar = c(4.5, 4.5, 1.5, 1))
  x <- res$KROAD_beta; y <- res$pittsburgh_beta
  lim <- range(c(x, y, 0)); pad <- 0.12 * diff(lim); lim <- lim + c(-pad, pad)
  plot(NA, xlim = lim, ylim = lim, xlab = "K-ROAD stage effect (SD per stage)", ylab = "Pittsburgh stage effect (SD per stage)", las = 1)
  abline(0, 1, lty = 2, col = grey(0.45)); abline(h = 0, v = 0, col = grey(0.8), lwd = 0.6)
  con <- res$same_direction & !res$rare_case_target
  points(x[con], y[con], pch = 21, bg = grey(0.15), col = "black", cex = 1.1)
  points(x[!res$same_direction & !res$rare_case_target], y[!res$same_direction & !res$rare_case_target], pch = 21, bg = "white", col = "black", cex = 1.2)
  points(x[res$rare_case_target], y[res$rare_case_target], pch = 24, bg = grey(0.75), col = "black", cex = 1.0)
  lab <- res$protein %in% c("GFAP","ACHE","NfL","GAP43","BACE1","SPP1") | !res$same_direction
  text(x[lab], y[lab], res$protein[lab], pos = 4, offset = 0.3, cex = 0.65)
  s1 <- summ[2, ]
  legend("topleft", bty = "n", cex = 0.75,
         legend = c(sprintf("Spearman rho = %.2f, P = %s (excl. rare-case)", s1$Spearman_rho, fp(s1$Spearman_P)),
                    sprintf("Direction concordance = %d/%d (%.0f%%)", s1$n_same_direction, s1$n_proteins, s1$direction_concordance_percent)))
  legend("bottomright", bty = "n", cex = 0.75, pch = c(21, 21, 24), pt.bg = c(grey(0.15), "white", grey(0.75)),
         legend = c("Concordant direction", "Opposite direction", "Rare-case target"))
}, width = 5.2, height = 4.8)

say(""); print(summ[, c("set","n_proteins","n_same_direction","direction_concordance_percent","Spearman_rho","Spearman_P")], row.names = FALSE)
say(sprintf("\nFigure 5B complete: %d/34 proteins analysed. %s", nrow(res), normalizePath(opt$outdir)))
