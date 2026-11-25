# GenRomica: genetic ICA component scores and child outcomes in Generation R

This pipeline tests whether genetic scores are associated with child
outcomes in the Generation R Study (GenR). The scores are ICA components of
genome-wide SNP effects, derived in an independent study, plus other polygenic
scores. Outcomes are behaviour (CBCL), brain imaging (MRI, DTI), SES, diet and
online activity.

It is a restructured version of the original analysis code
(`Code_Hailey_Components.docx`). **The statistical logic is unchanged, with one
deliberate exception:** all mixed-model tests now use the same Satterthwaite
F-test (see [Change from the original](#change-from-the-original)). Apart from
that, on simulated data the new code reproduces the original code's output
exactly (see [Verification](#verification)). Concerns about the code and the science are
listed in [Concerns](#concerns); none of them has been acted on yet.

## Pipeline

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `R/01_export_genr_data.R` | GenR RDS file | one `csv.gz` per data type (run once) |
| 2 | `R/02_merge_components.R` | variant catalogue, ICA Z-scores and p-values | `component<i>.txt` |
| – | *(not in this pipeline)* | `component<i>.txt` | `component<i>_corrected.txt` (see concern A1) |
| 3 | `R/03_clean_components.R` | `component<i>_corrected.txt` | `component_<i>_QC.txt` |
| 4 | `R/04_compute_scores.R` | GenR genotypes (bigSNP), QC'd weights | subject-level scores `csv.gz` |
| 5 | `R/05_run_regressions.R` | step 1 + step 4 output | results tables |

Every step takes a config file, and is run from this folder:

```bash
Rscript R/05_run_regressions.R config/project.R
```

- `config/project.R` holds all paths and analysis settings: waves, outcomes,
  score sets and reference age. The paths are the original ones; edit them for
  your machine.
- `config/toy.R` points the same settings at simulated data.
- `R/utils.R` holds the shared helpers.

Required R packages: `data.table`, `R.utils` (reads `.gz`), `lme4`, `lmerTest`,
`performance`, `bigsnpr`, `dplyr` (used by `bigsnpr::snp_match` for strand flips).

### Data conventions

- **Effect allele.** In the variant catalogue and all summary statistics, `A2`
  is the effect allele: the sign of the Z-score refers to `A2`. `A1` is the
  other allele.
- **`BETA` holds Z-scores.** The ICA component files already contain Z-scores.
  The column is called `BETA` only so that standard tools recognise it.
- **Positional matching.** The ICA Z-score and p-value files have no SNP
  identifier. Row *k* belongs to the *k*-th SNP of the variant catalogue with
  `include == 1`. Step 2 stops if the row counts differ.
- **GenR data.** Only GenR variables are exported in step 1. Genetic scores
  are never read from the RDS. They are computed from the genotypes in step 4,
  because their weights come from another study.

### Checks that stop the pipeline

Every step checks its inputs and stops with a message saying what is wrong and
which R commands to run to investigate. It stops when:

- an input file is missing;
- a required column is missing (the message gives `grep()` commands to find similarly named columns);
- `idc` values are duplicated;
- the row counts of the variant catalogue and the component files differ;
- a variable that enters a model as a number is a factor or text;
- a column shared by two GenR files holds different values.

In the original code, a measure missing at a wave was skipped silently. Here
every `<measure>_<wave>` listed in `cfg$longitudinal` must exist. If a measure
was not collected at a wave, remove that row from the config (there is an
example there). For example, CBCL/1½-5 at wave f05 has no rule-breaking scale.
This gives the same results as the original code, but the omission is now explicit.

### Answers to review questions

**Is the cleaning done in place?** No, in both versions. Each file is read from
an input folder and a cleaned copy is written to a separate output folder; the
input files are never changed. In memory, `df <- df[...]` creates a filtered
copy. The `:=` assignments change the freshly read table in place, which has
no effect outside the script.

**Why not use `bigsnpr::snp_match`?** It is now used in step 4, where it matters
most. `snp_match` aligns the weights with the genotype data:

- it joins on chromosome and position;
- it flips the sign of a weight when the genotype data count the other allele;
- it matches alleles reported on the opposite strand;
- it removes duplicates.

Without it, any SNP whose counted allele differs from the effect allele would
enter the score with the wrong sign.

`snp_match` does not replace the QC in step 3, which also removes non-finite
Z-scores, p-values outside (0, 1], non-autosomal SNPs and non-ACGT alleles.
Palindromic SNPs are still removed in step 3 to replicate the original.

Strand flipping follows the usual procedure. Matching first runs with
`strand_flip = TRUE`. If fewer than `cfg$min_strand_flips` SNPs (default 10)
needed a strand flip, it reruns with `strand_flip = FALSE`, which keeps
strand-ambiguous SNPs. Because step 3 already removes palindromic SNPs, the
rerun can currently only drop those few flipped SNPs; see concern B2. The
choice for each score is logged in `*_snp_matching.csv`.

## The statistical models

The full description is in the comments of `R/05_run_regressions.R`. In
summary:

**Repeated measures, one score at a time (linear mixed model, ML).** For
observation *j* of child *i*:

```
M0: y_ij = b0 + b1 a_ij + b2 a_ij² + g'X_i                                + u_i + e_ij
M1: M0 + b3 S_i
M2: M1 + b4 S_i·a_ij + b5 S_i·a_ij²
```

- *y* is the outcome, z-scored over all observations.
- *a* is age − 17 years.
- *S* is the score, z-scored over children.
- *X* holds sex, genotyping cohort and gPC1-4.
- *u_i* is a random intercept per child.

**Fixed effects:**
- the age terms describe the mean developmental trajectory;
- sex, genotyping batch and ancestry (gPCs) are confounder adjustments;
- *S* is the score being tested;
- *S*×age terms let the association change with age.

**Random effect:** only a child intercept. It accounts for the correlation
between repeated observations of the same child. There are no random slopes
and no family effect.

**Tests** (all `lmerTest::contestMD`, Satterthwaite denominator df):
- **H1** is the *F*-test of b3 = 0 in M1 (1 df; the same p-value as the
  Satterthwaite *t*-test).
- **H2** is the joint *F*-test of b4 = b5 = 0 in M2.
- The incremental Nakagawa marginal R² is reported for each.

**Repeated measures, a whole score set jointly.** Same model with all *K*
scores as fixed effects. H1 is the joint *F*-test that all *K* score
coefficients are zero in M1.

**Single-time outcomes (OLS).** `y_i = b0 + g'X_i [+ age at measurement] + Σ c_k S_ik + e_i`.
H1 is the *F*-test of the nested models; for a single score, the estimate is c_1.

### Change from the original

The original code used three different tests for the mixed models:

- a Satterthwaite *t*-test for single-score H1;
- a Satterthwaite *F*-test (`contestMD`) for H2;
- a likelihood-ratio χ² test, `anova(m0, m1)`, for H1 of the combined score sets.

All three now use `contestMD` through one helper, `joint_f_test()`. This
changes only the combined-set p-values, and the FDR values within those sets.
For single-score H1, the 1-df *F*-test gives the same p-value as the *t*-test.

To reproduce the original results, uncomment one line in `fit_lmm_score_set()`
in `R/05_run_regressions.R`:

```r
# h1_p <- anova(m0, m1)[2, "Pr(>Chisq)"]
```

With that line uncommented, the toy run matches the original exactly. A
commented line also restores the old single-score *t*-test p-value, which is
numerically identical.

**Multiple testing.** Benjamini-Hochberg FDR, applied separately within each
score set × hypothesis.

## Verification

```bash
bash toy/run_toy.sh
```

This takes about one minute on a laptop. It does three things:

1. It simulates toy data (`toy/make_toy_data.R`): 500 GenR children, 480
   genotyped children and 3,000 SNPs. The data include SNPs that QC must remove,
   SNPs with swapped alleles or reported on the opposite strand, missing
   genotype calls, the `sum_wid_f05` naming quirk, stale score columns in the
   RDS, and a non-alphabetical factor reference level.
2. It runs steps 1-5 with `config/toy.R`.
3. It runs `toy/compare_with_original.R`, which:
   - runs the original code, transcribed in `reference/` with only the paths
     changed, and checks that the merged and cleaned files are identical;
   - recomputes one score by hand and checks it against bigsnpr;
   - runs the original regression code and checks that every result table
     matches the new code to 1e-10. The exception is the combined-set
     mixed-model p-values (and FDR values), which changed deliberately. For
     those, it reports how the new F-test and the old likelihood-ratio
     p-values compare. On the toy data the correlation of −log10 p is 1.000.

All checks pass.

**Reference transcription.** In `reference/original_regression.R`, two parts
could not be read from the screenshots and had to be assumed (marked `# EDIT`):

- line 134, the body of `scale_gpcs`, assumed to be `z_by_id` per gPC;
- everything after line 476 (the end of the FDR merge and writing the output).

The assumption about line 134 only affects the gPC rows of
`covariate_results.csv`. Rescaling a covariate does not change any other
estimate, p-value or R². Check both parts against the real script.

**Precision.** Step 1 writes decimal numbers with 17 significant digits.
`fwrite`'s default of 15 digits changed the Nakagawa R² at the 1e-6 level.

## Concerns

These are recorded for discussion and have **not** been acted on, so that
the original results can be replicated first.

### A. Pipeline and data handling

1. **Unknown "corrected" step.** The QC script reads `component<i>_corrected.txt`,
   but step 2 writes `component<i>.txt`. The correction in between is not part of
   the code we received, so this pipeline cannot rerun it. It should be added
   or documented.
2. **Positional matching of ICA output.** Row order is checked only through the
   row count. If the ICA software kept SNP identifiers, joining on them would be safer.
3. **Meaning of the catalogue's 7th column (`include`).** It is undocumented.
4. **Genotype IDs.** Step 4 assumes the `.fam` `sample.ID` equals GenR `idc`
   (`cfg$genotype_id_column`). If it does not, a linking file is needed.
5. **Genome build.** `snp_match` joins on position, so the ICA study and the
   GenR genotypes must use the same build (GRCh37 vs GRCh38). This is not checked.
   `snp_match` errors when fewer than 20% of SNPs match, which catches gross
   mismatches only.
6. **Imputed dosages vs hard calls.** Step 4 mean-imputes missing calls
   (`mean2`, rounded to 2 decimals). If the scores in the original RDS were made
   differently (imputed dosages, other software), they will not match exactly.
   Compare step-4 scores with the RDS scores once before replacing them.
7. **UK_PCS and PGS sources.** It is unknown which summary statistics produced
   `UK_PCS_*` and `PGS_*`. Add them to `cfg$score_weights`, or list a file of
   precomputed scores in `cfg$score_files`.
8. **Duplicated `idc`.** The original code would silently use the first row
   (`match()`); the pipeline now stops.
9. **Warnings were invisible.** The original code captured errors but not
   warnings, such as non-convergence. They are now logged in
   `model_warnings.csv` without changing any result. Singular fits are flagged
   in `results.csv` (`singular`), but their p-values are still reported.
10. **Silent skips.** A model with too little variation returns no row and is
    not logged in `failed_models.csv`. Compare the number of result rows with
    the number expected.

### B. Summary statistics and scores

1. **p = 0 is removed.** QC keeps `P > 0`, which drops SNPs whose p-value
   underflowed to 0. These are the strongest signals. Consider setting them to
   the smallest representable number, or recomputing P from Z.
2. **Palindromic SNPs are removed by default.** This costs roughly 15% of SNPs.
   With well-defined strands (same reference panel), `snp_match(strand_flip = FALSE)`
   could keep them. The strand-flip switch in step 4 only helps if they are not
   already removed in step 3.
3. **Lowercase alleles are dropped** rather than upper-cased.
4. **Score construction.** The score is Σ dosage × Z over all SNPs, with no LD
   clumping or shrinkage. Correlated SNPs therefore count several times, and
   Z-scores scale with the discovery study's sample size. That is fine for an
   ICA-component projection, but the score is not an optimised PGS.
5. **Ancestry.** GenR is multi-ethnic. If the ICA study used European-ancestry
   data, score performance will differ by ancestry. Four gPCs may not fully
   control stratification; consider restricting to, or stratifying by,
   genetic ancestry.

### C. Statistical models

1. **Age centred at 17, the oldest wave.** `age_c` and `age_c²` are then highly
   correlated, and in M2 the score main effect is the effect at age 17. In M1
   (H1) centring does not change b3, but centring at the mean age, or using
   orthogonal polynomials, would make M2 more stable and easier to interpret.
2. **Different CBCL instruments.** Wave f05 used the CBCL/1½-5, and later waves
   used the CBCL/6-18. The scales differ in items and range, but are pooled and
   z-scored together. Standardising within wave, or modelling the instrument,
   would be more defensible.
3. **Outcome standardisation is pooled over waves.** It does not change the
   p-values, but estimates are in pooled-SD units.
4. **Skewed outcomes.** CBCL sums are counts with many zeros, and the models
   assume normal residuals. Consider a transform, or a count model, for robustness.
5. **Siblings and twins.** GenR includes related children. The random intercept
   is per child only, so related children are treated as independent. Keep one
   child per family or add a family random effect.
6. **ML vs REML.** Models are fitted with ML, which is needed for the
   likelihood-ratio test. Satterthwaite tests of fixed effects are usually done
   on REML fits; the differences are small at GenR sample sizes.
7. **Different tests for mixed and OLS models.** All mixed-model tests now use
   a Satterthwaite *F*-test (see [Change from the original](#change-from-the-original)).
   The single-time OLS models use the exact nested-model *F*-test. Both are
   *F*-tests, but the mixed-model denominator df is approximated.
8. **Constant covariates are dropped silently.** If a covariate is constant
   (e.g. one genotyping cohort in a subsample), models for different outcomes
   use different covariate sets.
9. **Complete-case samples differ** between outcomes and scores, so `n_id` varies.
10. **FDR grouping.** Within each score set, longitudinal and single-time outcomes
    share one FDR family, while the jointly tested sets are separate families.
    This is defensible, but should be pre-specified. Whether to correct across
    score sets as well should be decided too.
11. **Incremental R².** Nakagawa R² is NA for singular fits, and differences
    can be negative. For OLS, the raw (not adjusted) R² difference is reported.
12. **SES outcomes are parental.** Maternal and paternal education and household
    income are parent characteristics. An association with the child's score
    reflects shared parental genetics (genetic nurture, or confounding),
    not an effect on the child.
13. **Ordinal outcomes treated as continuous** (`INCOME5`, `EDUCM`, `EDUCP`).
    Step 1 stops if they are stored as factors, so their coding has to be
    chosen explicitly.
14. **Timing.** `DietScore_8y` is measured at about age 8 and the online measure
    at about 17. Single-time outcomes are adjusted for age only where an age
    variable is given.
