###############################################################################
#  validation_common.R
#  Shared code for the Figure 5 external-validation scripts
#  (00_check_inputs.R, 01_ADNC_validation.R, 02_34protein_replication.R,
#  03_MMSE_validation.R, run_all.R). Not meant to be run on its own.
#
#  Requirements: base R plus MASS and nlme.
###############################################################################

## ===========================================================================
## 1. Fixed K-ROAD cutoffs (log2 NPQ)
## ===========================================================================
CUT_P217    <- 9.562758948
CUT_EMTBR_1 <- 9.735217224
CUT_EMTBR_2 <- 10.35845619
CUT_EMTBR_3 <- 11.60847036
STAGES      <- paste("Stage", 0:4)

## ===========================================================================
## 2. Column detection
## ===========================================================================
## Set a name here only if automatic detection picks the wrong column.
COLUMN_MAP <- list(
  id = NULL, sample_id = NULL, sample_date = NULL, death_date = NULL,
  p217 = NULL, emtbr = NULL, age = NULL, sex = NULL, plate = NULL, adnc = NULL
)
ALIASES <- list(
  id          = c("PTID","participant_id","subject_id","case_id","BBLID","id"),
  sample_id   = c("sample_id","specimen_id","plasma_id","aliquot_id"),
  sample_date = c("draw_date","sample_date","blood_date","plasma_date","collection_date","blood_draw_date",
                  "date_of_draw","date_drawn","plasma_collection_date","blood_collection_date","specimen_date","visit_date"),
  death_date  = c("DateDeath","death_date","date_of_death","dod","autopsy_date","date_death","deathdate","dateofautopsy"),
  p217        = c("neuro220_BD_pTau_217","BD_pTau_217","BD-p-tau217","BD_p_tau217","bd_ptau217"),
  emtbr       = c("multiplex_MTBR","eMTBR_Tau","eMTBR-tau","eMTBR_tau","emtbr_tau","emtbr_tau243"),
  age         = c("blood_draw_age","age_at_draw","age","age_at_blood","age_at_sample","age_draw","draw_age","age_at_visit"),
  sex         = c("sex","gender"),
  plate       = c("neuro220_plateID","plate","plate_id","plateid"),
  adnc        = c("ADNC","adnc","ADNC_severity","adnc_severity","ADNC_category","ADNC_level","adnc_stage")
)

norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))

resolve_col <- function(df, slot, required = FALSE, contains_fallback = NULL) {
  if (!is.null(COLUMN_MAP[[slot]])) {
    if (!COLUMN_MAP[[slot]] %in% names(df))
      stop(sprintf("COLUMN_MAP$%s = '%s' is not a column in the input.", slot, COLUMN_MAP[[slot]]))
    return(COLUMN_MAP[[slot]])
  }
  nm <- setNames(names(df), norm_name(names(df)))
  for (a in ALIASES[[slot]]) if (norm_name(a) %in% names(nm)) return(unname(nm[norm_name(a)]))
  if (!is.null(contains_fallback)) {
    hits <- names(df)[grepl(contains_fallback, norm_name(names(df)))]
    if (length(hits) == 1) return(hits)
  }
  if (required)
    stop(sprintf(paste0("Required column '%s' could not be found.\nAvailable columns:\n  %s\n",
                        "Set COLUMN_MAP$%s in validation_common.R or pass --%s-col."),
                 slot, paste(names(df), collapse = ", "), slot, gsub("_", "-", slot)))
  NULL
}

read_any <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Excel input needs the 'readxl' package. Either install.packages(\"readxl\") or export to CSV.")
    df <- as.data.frame(readxl::read_excel(path, sheet = if (is.null(sheet)) 1 else sheet))
  } else {
    ## Read as raw UTF-8 text first so that non-ASCII header characters
    ## (for example a Greek beta in a protein name), a byte-order mark, or
    ## Windows line endings never cause a silent zero-row read.
    txt <- readLines(path, encoding = "UTF-8", warn = FALSE)
    txt <- sub("^\ufeff", "", txt); txt <- sub("\r$", "", txt)
    txt <- txt[nzchar(trimws(txt))]
    sep <- if (ext == "tsv") "\t" else ","
    df <- read.csv(text = txt, sep = sep, check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
    if (!nrow(df)) stop("The file was read but contains no data rows: ", path)
  }
  names(df) <- gsub("^\ufeff", "", trimws(names(df)))
  df
}

## Participant ids: trim, and turn "10023.0" (a numeric id that went through
## a float column) back into "10023" so the two files merge.
norm_id <- function(x) {
  z <- trimws(as.character(x))
  z <- sub("\\.0+$", "", z)
  z[z %in% c("", "NA", "NaN")] <- NA
  z
}

parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  x <- trimws(as.character(x)); x[x %in% c("", "NA", "NaN", ".", "-", "NULL", "n/a", "N/A")] <- NA
  if (all(is.na(x) | grepl("^[0-9]{5}(\\.[0-9]+)?$", x))) return(as.Date(as.numeric(x), origin = "1899-12-30"))
  ## Candidate formats. A parse only counts as a success if the resulting year
  ## is plausible: R will happily read "07/03/2014" with %Y/%m/%d as the year 7.
  fm <- c("%Y-%m-%d","%Y/%m/%d","%m/%d/%Y","%d/%m/%Y","%m-%d-%Y","%d-%m-%Y","%Y%m%d",
          "%d%b%Y","%d-%b-%Y","%d %b %Y","%b %d %Y","%b %d, %Y","%Y-%m-%d %H:%M:%S","%m/%d/%Y %H:%M","%m/%d/%y","%d/%m/%y")
  plausible <- function(d) !is.na(d) & as.integer(format(d, "%Y")) >= 1900 & as.integer(format(d, "%Y")) <= 2100
  best <- NULL; bok <- -1
  for (f in fm) {
    d <- suppressWarnings(as.Date(x, format = f)); d[!plausible(d)] <- NA
    ok <- sum(!is.na(d))
    if (ok > bok) { bok <- ok; best <- d }
  }
  best
}

encode_female <- function(v) {
  if (is.numeric(v) && all(v[!is.na(v)] %in% c(0, 1))) return(as.numeric(v))
  z <- tolower(trimws(as.character(v)))
  out <- rep(NA_real_, length(z))
  out[z %in% c("f","female","woman","women","2")] <- 1
  out[z %in% c("m","male","man","men","1")]       <- 0
  out
}

## ===========================================================================
## 3. Stage assignment and index-sample selection
## ===========================================================================
assign_stage <- function(p, e) {
  s <- rep(NA_character_, length(p))
  ok <- !is.na(p) & !is.na(e)
  t1 <- ok & p >= CUT_P217
  s[ok & !t1 & e <  CUT_EMTBR_1] <- "Stage 0"
  s[ok &  t1 & e <  CUT_EMTBR_1] <- "Stage 1"
  s[ok &  t1 & e >= CUT_EMTBR_1 & e < CUT_EMTBR_2] <- "Stage 2"
  s[ok &  t1 & e >= CUT_EMTBR_2 & e < CUT_EMTBR_3] <- "Stage 3"
  s[ok &  t1 & e >= CUT_EMTBR_3] <- "Stage 4"
  s[ok & !t1 & e >= CUT_EMTBR_1] <- "Unstageable"
  s[!ok] <- "Missing biomarker"
  s
}

## Build the base frame, then keep one pre-mortem sample per participant
## (closest to death; single-sample participants kept; post-mortem-dated
## samples excluded; ties broken by later draw date then sample id).
prepare_index <- function(raw) {
  col <- list(
    id = resolve_col(raw, "id", TRUE), sample_id = resolve_col(raw, "sample_id"),
    sample_date = resolve_col(raw, "sample_date", TRUE, contains_fallback = "draw|collect"),
    death_date = resolve_col(raw, "death_date", TRUE, contains_fallback = "death|autops"),
    p217 = resolve_col(raw, "p217", TRUE), emtbr = resolve_col(raw, "emtbr", TRUE),
    age = resolve_col(raw, "age", TRUE), sex = resolve_col(raw, "sex", TRUE),
    plate = resolve_col(raw, "plate")
  )
  d <- data.frame(raw_row = seq_len(nrow(raw)), stringsAsFactors = FALSE)
  d$id          <- norm_id(raw[[col$id]])
  d$sample_id   <- if (!is.null(col$sample_id)) norm_id(raw[[col$sample_id]]) else paste0("row_", d$raw_row)
  d$sample_date <- parse_date(raw[[col$sample_date]])
  d$death_date  <- parse_date(raw[[col$death_date]])
  d$p217        <- suppressWarnings(as.numeric(raw[[col$p217]]))
  d$emtbr       <- suppressWarnings(as.numeric(raw[[col$emtbr]]))
  d$age         <- suppressWarnings(as.numeric(raw[[col$age]]))
  d$female      <- encode_female(raw[[col$sex]])
  d$plate       <- if (!is.null(col$plate)) as.character(raw[[col$plate]]) else NA_character_
  d$b2d_days    <- as.numeric(d$death_date - d$sample_date)
  d$b2d_years   <- d$b2d_days / 365.25

  el <- d[!is.na(d$sample_date) & !is.na(d$death_date) & d$b2d_days >= 0, ]
  el <- el[order(el$id, el$b2d_days, -as.numeric(el$sample_date), el$sample_id), ]
  idx <- el[!duplicated(el$id), ]
  rownames(idx) <- NULL

  idx$stage     <- assign_stage(idx$p217, idx$emtbr)
  idx$stage_num <- match(idx$stage, STAGES) - 1
  idx$stage_cat <- factor(idx$stage, levels = STAGES)
  list(idx = idx, col = col,
       n_rows = nrow(raw), n_participants = length(unique(d$id)),
       n_dropped_postmortem = sum(!is.na(d$b2d_days) & d$b2d_days < 0),
       n_missing_dates = sum(is.na(d$sample_date) | is.na(d$death_date)))
}

## ===========================================================================
## 4. Command line helpers
## ===========================================================================
ARGS <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(ARGS == flag); if (!length(i)) return(default)
  if (i[1] == length(ARGS)) stop(sprintf("%s requires a value.", flag)); ARGS[i[1] + 1]
}
has_flag <- function(flag) flag %in% ARGS

## Apply --*-col overrides to COLUMN_MAP
apply_col_overrides <- function() {
  for (k in names(COLUMN_MAP)) {
    v <- get_arg(paste0("--", gsub("_", "-", k), "-col"))
    if (!is.null(v)) COLUMN_MAP[[k]] <<- v
  }
}

common_opts <- function() {
  apply_col_overrides()
  list(
    input    = get_arg("--input"),
    sheet    = get_arg("--sheet"),
    outdir   = get_arg("--outdir", "Fig5_output")
  )
}

COMMON_HELP <- "
Shared options
  --input FILE            plasma marker file (.csv / .tsv / .xlsx)
  --sheet NAME            worksheet, for .xlsx
  --outdir DIR            output directory (default Fig5_output)
  --<field>-col NAME      override column detection: id, sample-id, sample-date,
                          death-date, p217, emtbr, age, sex, plate, adnc
"

## ===========================================================================
## 5. Statistics helpers
## ===========================================================================
## HC3 heteroscedasticity-robust covariance for an lm fit (base R only).
vcov_hc3 <- function(fit) {
  X <- model.matrix(fit); e <- residuals(fit); h <- hatvalues(fit)
  XtXi <- solve(crossprod(X))
  meat <- crossprod(X * (e / (1 - h)))
  XtXi %*% meat %*% XtXi
}

## OR, CI and p for one coefficient of a MASS::polr fit
polr_coef <- function(fit, term) {
  b <- coef(fit)[term]; se <- sqrt(diag(vcov(fit)))[term]
  z <- b / se
  c(beta = unname(b), se = unname(se), OR = exp(unname(b)),
    OR_low = exp(unname(b - 1.96 * se)), OR_high = exp(unname(b + 1.96 * se)),
    P = unname(2 * pnorm(-abs(z))))
}

fp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 2), formatC(p, format = "f", digits = 3)))

## ===========================================================================
## 6. Output helpers
## ===========================================================================
write_csv_safe <- function(df, path) { write.csv(df, path, row.names = FALSE, na = ""); invisible(path) }

write_t1_info <- function(prep, out) {
  sc <- table(factor(prep$idx$stage, levels = c(STAGES, "Unstageable", "Missing biomarker")))
  df <- data.frame(field = c("t1_cutoff","emtbr_cutoffs",
                             "n_rows_in_file","n_participants","n_index","n_dropped_postmortem",
                             "n_missing_dates", paste0("n_", names(sc))),
                   value = c(sprintf("%.9f", CUT_P217), paste(c(CUT_EMTBR_1, CUT_EMTBR_2, CUT_EMTBR_3), collapse = " / "),
                             prep$n_rows, prep$n_participants, nrow(prep$idx),
                             prep$n_dropped_postmortem, prep$n_missing_dates, as.integer(sc)),
                   stringsAsFactors = FALSE)
  write_csv_safe(df, file.path(out, "T1_cutoff_and_stage_counts.csv"))
  invisible(df)
}

## Save a figure as vector PDF and 600 dpi PNG. `draw` is a function.
save_fig <- function(path_stem, draw, width = 5, height = 4.2) {
  pdf(paste0(path_stem, ".pdf"), width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  draw(); invisible(dev.off())
  png(paste0(path_stem, ".png"), width = width, height = height, units = "in", res = 600)
  draw(); invisible(dev.off())
}

GREYS <- c("0.10", "0.30", "0.48", "0.66", "0.80")
grey_col <- function(i) grey(as.numeric(GREYS[i]))

banner <- function(txt) { cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "") }
say <- function(...) cat(paste0(..., "\n"))
