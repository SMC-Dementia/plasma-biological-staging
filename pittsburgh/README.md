# PITT-ADRC external validation (Figure 5)

R scripts that apply the fixed K-ROAD staging thresholds to the PITT-ADRC cohort and run the three Figure 5 analyses. They were run by the PITT-ADRC team on participant-level data; the returned aggregate results are in `../figure5/source_data/`.

Requirements: R 4.x with `MASS` and `nlme` (both ship with R). No other packages.

| Panel | Analysis | Script |
|---|---|---|
| 5a | Ordered ADNC severity vs plasma stage 0 to 4 (proportional-odds regression, adjusted for age, sex and blood-to-death interval; categorical-stage model, Spearman correlation, plate-adjusted and interval-restricted sensitivity analyses) | `01_ADNC_validation.R` |
| 5b | Replication of the 34 stage-associated proteins (z-scored protein ~ ordinal stage + age + sex, HC3 standard errors; comparison with the K-ROAD coefficients hard-coded in `KROAD_BETA`; BH-FDR across the 34 tests; plate-adjusted sensitivity) | `02_34protein_replication.R` |
| 5c | MMSE trajectories by baseline stage (ML linear mixed model with time x stage, age and sex, random intercept and slope per participant; 4-df likelihood-ratio test) | `03_MMSE_validation.R` |

`validation_common.R` holds the fixed cutoffs, index-sample selection (one pre-mortem sample per participant, closest to death) and stage assignment. `00_check_inputs.R` is a dry run that reports column detection and stage counts without fitting any model. `run_all.R` runs all three analyses and writes a summary.

## Running

Demo (synthetic data):

```sh
python run_pittsburgh.py --input demo/synthetic_plasma.csv --clinical demo/synthetic_mmse.csv --outdir pitt_demo_out
```

`run_pittsburgh.py` calls the three scripts and writes `run_summary.json`. Expected output is in `demo/expected_output/`. Use `--rscript <path>` if `Rscript` is not on the PATH. Column and sheet overrides accepted by the R scripts can be appended.

Scripts can also be run directly:

```sh
Rscript 00_check_inputs.R --input plasma.csv --clinical mmse.csv
Rscript run_all.R --input plasma.csv --clinical mmse.csv --outdir Fig5_output
```

## Input files

Columns are detected by name (alias lists in `validation_common.R`; `--<field>-col` overrides, or `COLUMN_MAP`).

Plasma CSV, one row per sample:

| Field | Content |
|---|---|
| `PTID` | participant id |
| `sample_id` | sample id (tie breaking) |
| `draw_date` | plasma collection date, ISO format |
| `DateDeath` | death date, ISO format |
| `neuro220_BD_pTau_217` | BD-p-tau217, log2 NPQ |
| `multiplex_MTBR` | eMTBR-tau, log2 NPQ |
| `blood_draw_age` | age at draw |
| `sex` | Male / Female |
| `ADNC` | Not/Low, Intermediate, High (source labels are harmonised; see `Fig5A_ADNC_label_audit.csv`) |
| `neuro220_plateID` | plate id (sensitivity analyses) |
| 34 protein columns | log2 NPQ; names resolved as listed in `Fig5B_protein_column_resolution.csv`, overrides in `PROTEIN_COLUMN_MAP` |

Longitudinal MMSE CSV, one row per assessment: `PTID`, `MMSE_date`, `MMSE`. Assessments from 1 year before the index plasma sample onwards are retained; participants need at least two retained assessments including one on or after the index date.

## Output

Per run: `T1_cutoff_and_stage_counts.csv`, `Fig5A_*` (proportions, model results, categorical-stage model, Spearman, label audit, cohort summary, figure), `Fig5B_*` (per-protein results, replication summary, plate-adjusted sensitivity, column resolution, figure), `Fig5C_*` (stage slopes, global test, counts, predicted trajectories, figure) and model summaries. Aggregate values only; no participant-level rows are written.

`tests/check_staging.R` checks the stage assignment at the exact cutoffs:

```sh
Rscript --vanilla tests/check_staging.R
```
