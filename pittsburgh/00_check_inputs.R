#!/usr/bin/env Rscript
###############################################################################
#  00_check_inputs.R  --  dry run, no modelling
#
#  Reads the input files, reports which columns were detected for every
#  field the analyses need, counts the index samples and stages, and lists
#  any protein whose column could not be found.
#
#  Usage
#    Rscript 00_check_inputs.R --input plasma.csv --clinical mmse.csv
###############################################################################
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "validation_common.R"))
if (has_flag("--help") || !length(ARGS)) { cat("Dry run.\n  --input FILE  [--clinical FILE]  plus any --*-col override\n"); quit(status = 0) }
opt <- common_opts(); clinical <- get_arg("--clinical")
if (is.null(opt$input)) stop("--input is required.")

banner("Input check (no analysis is run)")
say("R version : ", R.version$major, ".", R.version$minor)
say("MASS      : ", if (requireNamespace("MASS", quietly = TRUE)) "available" else "MISSING (needed for 5A)")
say("nlme      : ", if (requireNamespace("nlme", quietly = TRUE)) "available" else "MISSING (needed for 5C)")

raw <- read_any(opt$input, opt$sheet)
say(sprintf("\nPlasma file: %d rows, %d columns", nrow(raw), ncol(raw)))
say("\nColumn detection")
for (k in c("id","sample_id","sample_date","death_date","p217","emtbr","age","sex","plate","adnc")) {
  v <- tryCatch(resolve_col(raw, k, FALSE, if (k == "adnc") k else NULL), error = function(e) NULL)
  say(sprintf("  %-12s -> %s", k, if (is.null(v)) "(not found)" else v))
}

say("\nIndex sample selection and stage counts")
p <- prepare_index(raw)
sc <- table(factor(p$idx$stage, levels = c(STAGES, "Unstageable", "Missing biomarker")))
say(sprintf("  index n=%d  stages 0/1/2/3/4 = %s  unstageable=%d  missing biomarker=%d",
            nrow(p$idx), paste(sc[STAGES], collapse = "/"), sc["Unstageable"], sc["Missing biomarker"]))
say(sprintf("\n  participants in file: %d   post-mortem-dated samples dropped: %d   rows with missing dates: %d",
            p$n_participants, p$n_dropped_postmortem, p$n_missing_dates))
b2d <- p$idx$b2d_years
say(sprintf("  blood-to-death (years): median %.2f, IQR %.2f to %.2f, max %.2f", median(b2d, na.rm = TRUE), quantile(b2d, .25, na.rm = TRUE), quantile(b2d, .75, na.rm = TRUE), max(b2d, na.rm = TRUE)))
ad <- tryCatch(resolve_col(raw, "adnc", FALSE, "adnc"), error = function(e) NULL)
if (!is.null(ad)) say("  ADNC values seen : ", paste(names(table(raw[[ad]])), collapse = " | "))

## proteins
say("\nProtein columns for Figure 5B")
src <- readLines(file.path(here, "02_34protein_replication.R"))
env <- new.env(); eval(parse(text = src[grep("^KROAD_BETA <- c\\(", src):grep("^\\)$", src)[1]]), env)
eval(parse(text = src[grep("^protein_aliases <- function", src):(grep("^resolve_protein <- function", src) - 1)]), env)
eval(parse(text = src[grep("^resolve_protein <- function", src):(grep("^banner\\(", src)[1] - 1)]), env)
env$PROTEIN_COLUMN_MAP <- list(); env$norm_name <- norm_name
found <- sapply(names(env$KROAD_BETA), function(nm) env$resolve_protein(raw, nm))
say(sprintf("  detected: %d / 34", sum(!is.na(found))))
if (any(is.na(found))) say("  NOT detected: ", paste(names(found)[is.na(found)], collapse = ", "),
                           "\n  -> add these to PROTEIN_COLUMN_MAP in 02_34protein_replication.R, e.g.\n",
                           "     PROTEIN_COLUMN_MAP <- list(`NfL` = \"neuro220_NEFL\")")

## clinical
if (!is.null(clinical)) {
  cl <- read_any(clinical, get_arg("--clinical-sheet"))
  say(sprintf("\nClinical file: %d rows, %d columns. Columns: %s", nrow(cl), ncol(cl), paste(names(cl), collapse = ", ")))
  cid <- tryCatch(resolve_col(cl, "id", TRUE), error = function(e) NULL)
  ov <- length(intersect(unique(norm_id(cl[[cid]])), p$idx$id))
  say(sprintf("  id column: %s   participants overlapping with plasma index set: %d", if (is.null(cid)) "(not found)" else cid, ov))
} else say("\nNo --clinical file given; Figure 5C will be skipped.")

say("\nIf the above looks right:\n  Rscript run_all.R --input <plasma> --clinical <mmse> --outdir Fig5_output")
