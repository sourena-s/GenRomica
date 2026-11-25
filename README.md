# GenomICA polygenic scores (GPS) of the environment and childhood behavioral outcomes in the Generation R cohort

This project uses data-driven, multivariate genetic dispositions to a large range
of heritable traits (up to 691) to estimate a broad range of environmental,
lifestyle and psychiatric traits in individuals, and to identify gene–environment
and gene–behaviour correlations. Specifically, it tests whether these genetic
predispositions predict the following variables in the Generation R Study
(GenR; N = **[N]** individuals, up to four study waves, ages **[x]**–**[y]** years):

- behavioural and emotional problems (Child Behavior Checklist, CBCL);
- brain structure and white-matter microstructure (MRI, DTI);
- family socioeconomic status (SES);
- diet;
- online activity.

The genetic dispositions are multivariate components extracted by genomic
independent component analysis (GenomICA)
([Soheili-Nezhad et al., 2021](https://doi.org/10.1101/2021.01.06.425535);
[Oblong, Soheili-Nezhad et al., 2024](https://doi.org/10.1111/gbb.12876)). They
are extracted in a data-driven way with the noisy (probabilistic) ICA algorithm
MELODIC ([Beckmann & Smith, 2004](https://doi.org/10.1109/TMI.2003.822821)) and
FastICA ([Hyvärinen, 1999](https://doi.org/10.1109/72.761722)), applied to up to
691 genome-wide association studies (GWAS) of UK Biobank traits from the Pan-UK
Biobank project ([Karczewski et al., 2025](https://doi.org/10.1038/s41588-025-02335-7)).
Each component is a population-wide set of genome-wide variant weights, and
every genotyped GenR individual receives one score per component.

### References

- Soheili-Nezhad S, Beckmann CF, Sprooten E. Independent genomic sources of
  brain structure and function. *bioRxiv* 2021. doi:[10.1101/2021.01.06.425535](https://doi.org/10.1101/2021.01.06.425535)
- Oblong LM, Soheili-Nezhad S, *et al.* Principal and independent genomic
  components of brain structure and function. *Genes, Brain and Behavior*
  2024;23(1):e12876. doi:[10.1111/gbb.12876](https://doi.org/10.1111/gbb.12876)
- Beckmann CF, Smith SM. Probabilistic independent component analysis for
  functional magnetic resonance imaging. *IEEE Transactions on Medical Imaging*
  2004;23(2):137–152. doi:[10.1109/TMI.2003.822821](https://doi.org/10.1109/TMI.2003.822821)
- Hyvärinen A. Fast and robust fixed-point algorithms for independent component
  analysis. *IEEE Transactions on Neural Networks* 1999;10(3):626–634.
  doi:[10.1109/72.761722](https://doi.org/10.1109/72.761722)
- Karczewski KJ, *et al.* (Pan-UK Biobank). Pan-UK Biobank genome-wide
  association analyses enhance discovery and resolution of ancestry-enriched
  effects. *Nature Genetics* 2025;57:2408–2417. doi:[10.1038/s41588-025-02335-7](https://doi.org/10.1038/s41588-025-02335-7)

## Pipeline

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `R/01_export_genr_data.R` | GenR RDS file | one `csv.gz` per data type (run once) |
| 2 | `R/02_merge_components.R` | variant catalogue, ICA Z-scores | `component<i>.txt`, Z-score QC plots |
| – | *(not in this pipeline)* | `component<i>.txt` | `component<i>_corrected.txt` (see concern A1) |
| 3 | `R/03_clean_components.R` | `component<i>_corrected.txt` | `component_<i>_QC.txt` |
| 4 | `R/04_compute_scores.R` | GenR genotypes (bigSNP), QC'd weights | subject-level scores `csv.gz` |
| 5 | `R/05_run_regressions.R` | step 1 + step 4 output | results tables |
| 6 | `R/06_plot_results.R` | step 5 output | PheWAS plot and heatmaps |

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
- **`BETA` holds the component weights.** These are the `melodic_oIC` values
  exported by genomica, treated as Z-scores (open issue: see concern B6). The
  column is called `BETA` only so that standard tools recognise it.
- **P is computed from Z.** Only the Z-score files are read. Step 2 computes
  two-sided p-values as `P = 2·Φ(−|Z|)`. Old `p_val_component<i>.txt` files are
  not used; the toy data writes them only for the original code in `reference/`.
- **Positional matching.** The ICA Z-score files have no SNP identifier.
  Row *k* belongs to the *k*-th SNP of the variant catalogue with
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
Results match the original code, and the config shows which waves each measure uses.

### Answers to review questions

**Is the cleaning done in place?** No. Each file is read from
an input folder and a cleaned copy is written to a separate output folder; the
input files are never changed. In memory, `df <- df[...]` creates a filtered
copy. The `:=` assignments change the freshly read table in place, which has
no effect outside the script.

**Why not use `bigsnpr::snp_match`?** It is used in step 4, where it aligns the
weights with the genotype data:

- it joins on chromosome and position;
- it flips the sign of a weight when the genotype data count the other allele;
- it matches alleles reported on the opposite strand;
- it removes duplicates.

Without it, any SNP whose counted allele differs from the effect allele would
enter the score with the wrong sign.

`snp_match` does not replace the QC in step 3, which also removes non-finite
Z-scores, p-values outside (0, 1], non-autosomal SNPs and non-ACGT alleles.
Palindromic SNPs are still removed in step 3 to replicate the original.

Matching first runs with
`strand_flip = TRUE`. If fewer than `cfg$min_strand_flips` SNPs (default 10)
needed a strand flip, it reruns with `strand_flip = FALSE`, which keeps
strand-ambiguous SNPs. Because step 3 already removes palindromic SNPs, the
rerun can currently only drop those few flipped SNPs; see concern B2. The
choice for each score is logged in `*_snp_matching.csv`.

## GenR variables used

Each variable is tagged with its data type: *continuous*, *count*, *ordinal*,
*binary* or *categorical*. Candidate variables not used yet are listed at the
end, in [Candidate GenR variables not used yet](#candidate-genr-variables-not-used-yet).

- **Behaviour:** CBCL syndrome scales at waves f05, f09, f13 and f17 (`sum_<scale>_<wave>`): *count* sum scores.
  - anxious/depressed;
  - somatic complaints;
  - withdrawn;
  - social problems;
  - thought problems;
  - attention problems;
  - rule-breaking behaviour;
  - aggressive behaviour.
- **Brain imaging:**
  - total grey-matter volume and mean cortical thickness, by wave: *continuous*;
  - mean fractional anisotropy (FA) and mean diffusivity (MD) from DTI, by wave: *continuous*.
- **Lifestyle:**
  - diet quality score at 8 years (`DietScore_8y`): *continuous*;
  - weekly online activity at about 17 years (`onlineweek_child_F17`): *continuous*.
- **Socioeconomic status:**
  - household income (`INCOME5`): *ordinal* (5 categories), modelled as continuous;
  - maternal and paternal education (`EDUCM`, `EDUCP`): *ordinal*, modelled as continuous.
- **Covariates:**
  - sex: *binary*;
  - genotyping cohort: *categorical*;
  - gPC1–4: *continuous*;
  - age at each measurement: *continuous*.

## Quality control and filtering

Everything below happens before any model is fitted, in pipeline order.

### GenR phenotype data (step 1)

- **Required variables:** every variable listed in the config must exist; otherwise the pipeline stops.
- **One row per individual:** duplicated `idc` values stop the pipeline.
- **Numeric types:** outcomes, ages and gPCs must be numeric; a factor or text column stops the pipeline.
- **Naming fix:** `sum_wid_f05` is renamed to `sum_wit_f05`, the CBCL withdrawn scale at wave f05.
- **Exported variables:** only the variables used by the pipeline are exported. Score columns in the RDS are ignored.
- **No row filtering:** all individuals are kept at this step.

### GenomICA component Z-scores (step 2)

- **Included SNPs:** only SNPs flagged `include == 1` in the variant catalogue are kept.
- **Input format:** each Z-score file must be one numeric column, without a header.
- **Row counts:** the number of Z-scores must equal the number of included SNPs, because rows are matched by position.
- **p-values:** P is computed from Z as $`2\,\Phi(-|Z|)`$.

### Component summary statistics (step 3)

SNPs are removed, in this order, when they have:

1. a missing chromosome or position, a non-finite Z-score, or a p-value outside (0, 1]. This also removes |Z| > ~37.5, where P underflows to 0;
2. a chromosome other than 1–22;
3. an allele other than A, C, G or T, or two identical alleles.
   **TODO:** this also removes indels. Check whether indels can be aligned with
   the GenR genotypes and kept;
4. strand-ambiguous (palindromic) alleles, A/T or C/G;
5. a position (CHR:BP) that occurs more than once. All copies are removed.

The number of SNPs removed at each filter is printed for every component.

### Genotypes and scoring (step 4)

- **Missing genotype calls:** imputed with the mean dosage (`snp_fastImputeSimple`, method `mean2`).
- **Allele matching:** with `bigsnpr::snp_match()`, joining on chromosome and position.
  - A SNP is flipped in sign when the effect allele is the genotype's other allele.
  - Strand flips are allowed only if at least 10 SNPs need one (`cfg$min_strand_flips`); otherwise matching is redone without them.
  - SNPs duplicated after matching are removed.
  - The pipeline stops if fewer than 20% of SNPs match (`snp_match` default).
- **Logging:** the number of SNPs used, sign-flipped and strand-flipped is logged per score.
- **Individuals:** genotyped individuals not in GenR are dropped. GenR individuals without genotypes get missing scores.

### Model data (step 5)

- **Repeated measures:** an observation is used only if both the outcome and the age at that wave are present.
- **Complete cases, per model:** individuals or observations with any missing value in the outcome, age, sex, genotyping cohort, gPC1–4 or the score(s) are excluded.
- **Too little variation:** a model is skipped if, after exclusions,
  - the outcome has fewer than 2 distinct values;
  - any score has fewer than 2 distinct values;
  - or, for mixed models, age has fewer than 3 distinct values.
- **Constant covariates:** covariates that are constant in the analysis sample are dropped from that model.
- **Standardisation:** outcomes, scores and gPCs are z-transformed after the exclusions, within each model's sample.
- **Problem fits:** singular fits are flagged (`singular`) and their Nakagawa $`R^2`$ is set to NA. Errors are logged in `failed_models.csv` and warnings in `model_warnings.csv`.

### Not applied

These are left to upstream processing, or listed as open points under [Concerns](#concerns):

- **SNP-level filters:** no minor-allele-frequency, imputation-quality (INFO), call-rate or Hardy–Weinberg filtering of SNPs. Genotype QC is assumed done upstream.
- **Sample-level filters:**
  - no exclusion of related individuals (concern C5);
  - no restriction by genetic ancestry (concern B5);
  - no sample call-rate or sex-check filtering.
- **Outcomes:** no outlier removal or transformation (concern C4).

## The statistical models

The comments in `R/05_run_regressions.R` describe the same models next to the
code.

**Terminology.** Subject-level scores built from the ICA component weights are
called **GPS (GenomICA Polygenic Scores)**. For SNP $`k`$, let $`Z_k`$ be its
component Z-score and $`G_{ik} \in [0, 2]`$ the number of effect alleles individual
$`i`$ carries. The GPS of individual $`i`$ is

```math
\mathrm{GPS}_i = \sum_{k} G_{ik}\, Z_k .
```

It is computed in step 4. Score sets `C` and `UK_PCS` are GPS; set `PGS` holds
conventional polygenic scores. All are analysed in the same way.

### 1. Repeated measures, one score at a time — equation [1]

Linear mixed model, fitted by maximum likelihood. For measurement $`j`$
(wave f05, f09, f13 or f17) of individual $`i`$:

Full model, with every term labelled:

```math
\begin{aligned}
\underbrace{y_{ij}}_{\text{outcome}}
={}& \underbrace{\beta_0}_{\text{intercept}}
 + \underbrace{\beta_1\,a_{ij}}_{\text{age}}
 + \underbrace{\beta_2\,a_{ij}^2}_{\text{age}^2}
 + \underbrace{\beta_3\,S_i}_{\text{GPS}}
 + \underbrace{\beta_4\,S_i\,a_{ij}}_{\text{GPS}\,\times\,\text{age}}
 + \underbrace{\beta_5\,S_i\,a_{ij}^2}_{\text{GPS}\,\times\,\text{age}^2}\\[10pt]
&+ \underbrace{\gamma_{\mathrm{sex}}\,\mathrm{sex}_i}_{\text{sex}}
 + \underbrace{\boldsymbol{\gamma}_{\mathrm{coh}}^{\top}\mathbf{c}_i}_{\text{genotyping cohort}}
 + \underbrace{\textstyle\sum_{p=1}^{4}\gamma_{p}\,\mathrm{gPC}_{p,i}}_{\text{ancestry: gPC 1 to 4}}
 + {\color{red}\underbrace{u_i}_{\text{individual random intercept}}}
 + \underbrace{\varepsilon_{ij}}_{\text{residual}}
\end{aligned}
\qquad\qquad \textbf{[1]}
```

The three nested models differ only in which GPS terms are included:

```math
\mathcal{M}_0:\ \beta_3=\beta_4=\beta_5=0 \ \ (\text{no GPS})
\qquad
\mathcal{M}_1:\ \beta_4=\beta_5=0 \ \ (\text{GPS main effect})
\qquad
\mathcal{M}_2:\ \text{all terms} \ \ (\text{GPS}\times\text{age})
```

Fitted in R (`R/05_run_regressions.R`; the formula is built with `reformulate()`, and covariates that are constant in the sample are left out):

```math
\begin{aligned}
&\mathcal{M}_0:\ \mathtt{lmer(outcome }\sim\mathtt{ age\_c + age2} \\
&\qquad\quad \mathtt{+ sex + genotype\_cohort + gPC1 + gPC2 + gPC3 + gPC4 + }{\color{red}\mathtt{(1 | idc)}}\mathtt{, REML = FALSE)} \\[2pt]&\mathcal{M}_1:\ \mathtt{lmer(outcome }\sim\mathtt{ age\_c + age2 + score} \\
&\qquad\quad \mathtt{+ sex + genotype\_cohort + gPC1 + gPC2 + gPC3 + gPC4 + }{\color{red}\mathtt{(1 | idc)}}\mathtt{, REML = FALSE)} \\[2pt]&\mathcal{M}_2:\ \mathtt{lmer(outcome }\sim\mathtt{ age\_c + age2 + score + score:age\_c + score:age2} \\
&\qquad\quad \mathtt{+ sex + genotype\_cohort + gPC1 + gPC2 + gPC3 + gPC4 + }{\color{red}\mathtt{(1 | idc)}}\mathtt{, REML = FALSE)}
\end{aligned}
```

```math
{\color{red}u_i \sim \mathcal{N}(0, \sigma_u^2)}, \qquad \varepsilon_{ij} \sim \mathcal{N}(0, \sigma_\varepsilon^2), \qquad u_i \perp \varepsilon_{ij}
```

| Symbol | Meaning |
|---|---|
| $`y_{ij}`$ | **Outcome** of individual $`i`$ at wave $`j`$. Examples: the CBCL anxious/depressed sum score (`sum_anx_<wave>`), total grey-matter volume (`TotalGrayVol_<wave>`), mean fractional anisotropy (`mean_FA_genr_<wave>`). It is z-transformed once, over all observations in the model's analysis sample, so across individuals *and* waves together. |
| $`a_{ij}`$ | **Age** at that measurement minus 17 years. $`a_{ij}^2`$ is its square; together they give a quadratic developmental trajectory centred at age 17. |
| $`S_i`$ | **Score** of individual $`i`$ (a GPS or PGS). It is constant over waves and z-transformed across individuals: the mean and SD are taken over individuals, each counted once however many measurements they have. |
| $`\mathrm{sex}_i`$, $`\mathbf{c}_i`$, $`\mathrm{gPC}_{p,i}`$ | **Covariates:** sex (indicator); genotyping cohort (genotyping batch/array, as dummy variables $`\mathbf{c}_i`$); gPC1–4, the first four principal components of the genome-wide genotypes, which capture genetic ancestry (z-transformed across individuals). A covariate that is constant in the analysis sample is dropped. |
| $`u_i`$ | **Random intercept** of individual $`i`$: how far that individual's average outcome level lies from the population intercept $`\beta_0`$, the same at every wave. It models the correlation between repeated measurements of the same individual. |
| $`\varepsilon_{ij}`$ | Residual error of measurement $`j`$ of individual $`i`$. |

**Why these terms.**
- **Fixed effects:**
  - $`\beta_1, \beta_2`$ describe the average outcome trajectory with age.
  - $`\gamma_{\mathrm{sex}}`$, $`\boldsymbol{\gamma}_{\mathrm{coh}}`$ and $`\gamma_1,\dots,\gamma_4`$ adjust for sex differences, genotyping batch effects and population stratification (ancestry differences correlated with both score and outcome).
  - $`\beta_3`$ is the score association.
  - $`\beta_4, \beta_5`$ let that association change with age.
- **Random effects:** only $`u_i`$. There are no random slopes, so all individuals share one trajectory shape apart from the score terms. There is no family-level effect (see concern C5).
- **Maximum likelihood** is used rather than REML, so that models with different fixed effects can be compared.

**Hypotheses** (all tested with `lmerTest::contestMD`, Satterthwaite denominator df):

```math
\textbf{H1:}\;\; \beta_3 = 0 \;\text{ in } \mathcal{M}_1
\qquad\qquad
\textbf{H2:}\;\; \beta_4 = \beta_5 = 0 \;\text{ in } \mathcal{M}_2
```

- **H1** is a 1-df *F*-test; its p-value equals that of the Satterthwaite *t*-test. Because $`\mathcal{M}_1`$ has no interaction terms, $`\beta_3`$ is the association averaged over ages, in SD of outcome per SD of score.
- **H2** is a 2-df *F*-test of whether the association changes with age. In $`\mathcal{M}_2`$, $`\beta_3`$ is the association at age 17.
- **Effect size:** the incremental Nakagawa marginal $`R^2`$ (variance explained by fixed effects): $`R^2_{\mathcal{M}_1} - R^2_{\mathcal{M}_0}`$ for H1 and $`R^2_{\mathcal{M}_2} - R^2_{\mathcal{M}_1}`$ for H2.

### 2. Repeated measures, a whole score set jointly — equation [2]

All $`K`$ scores of a set (e.g. the 20 component GPS) enter one model, each
z-transformed across individuals:

```math
\begin{aligned}
\underbrace{y_{ij}}_{\text{outcome}}
={}& \underbrace{\beta_0}_{\text{intercept}}
 + \underbrace{\beta_1\,a_{ij}}_{\text{age}}
 + \underbrace{\beta_2\,a_{ij}^2}_{\text{age}^2}
 + \underbrace{\textstyle\sum_{k=1}^{K}\theta_k\,S_{ik}}_{K\ \text{GPS of one set}}\\[10pt]
&+ \underbrace{\gamma_{\mathrm{sex}}\,\mathrm{sex}_i}_{\text{sex}}
 + \underbrace{\boldsymbol{\gamma}_{\mathrm{coh}}^{\top}\mathbf{c}_i}_{\text{genotyping cohort}}
 + \underbrace{\textstyle\sum_{p=1}^{4}\gamma_{p}\,\mathrm{gPC}_{p,i}}_{\text{ancestry: gPC 1 to 4}}
 + {\color{red}\underbrace{u_i}_{\text{individual random intercept}}}
 + \underbrace{\varepsilon_{ij}}_{\text{residual}}
\end{aligned}
\qquad\qquad \textbf{[2]}
```

Fitted in R (`R/05_run_regressions.R`; the formula is built with `reformulate()`, and covariates that are constant in the sample are left out):

```math
\begin{aligned}
&\mathtt{lmer(outcome }\sim\mathtt{ age\_c + age2 + C1 + C2 + ... + C20} \\
&\qquad \mathtt{+ sex + genotype\_cohort + gPC1 + gPC2 + gPC3 + gPC4 + }{\color{red}\mathtt{(1 | idc)}}\mathtt{, REML = FALSE)}
\end{aligned}
```

```math
\textbf{H1:}\;\; \theta_1 = \theta_2 = \dots = \theta_K = 0 \qquad (K\text{-df } F\text{-test, Satterthwaite})
```

### 3. Outcomes measured once (ordinary least squares) — equation [3]

For SES, diet and online use (one row per individual):

```math
\begin{aligned}
\underbrace{y_{i}}_{\text{outcome}}
={}& \underbrace{\beta_0}_{\text{intercept}}
 + \underbrace{\textstyle\sum_{k}\theta_k\,S_{ik}}_{\text{GPS (one, or a set)}}
 + \underbrace{\gamma_{\mathrm{sex}}\,\mathrm{sex}_i}_{\text{sex}}
 + \underbrace{\boldsymbol{\gamma}_{\mathrm{coh}}^{\top}\mathbf{c}_i}_{\text{genotyping cohort}}
 + \underbrace{\textstyle\sum_{p=1}^{4}\gamma_{p}\,\mathrm{gPC}_{p,i}}_{\text{ancestry: gPC 1 to 4}}\\[10pt]
&+ \underbrace{\delta\,\mathrm{age}_i}_{\text{age at measurement (if available)}}
 + \underbrace{\varepsilon_{i}}_{\text{residual}},
\qquad \varepsilon_i \sim \mathcal{N}(0, \sigma^2)
\end{aligned}
\qquad\qquad \textbf{[3]}
```

Fitted in R with ordinary least squares (`R/05_run_regressions.R`; the formula is built with `reformulate()`, and covariates that are constant in the sample are left out):

```math
\begin{aligned}
&\mathtt{lm(outcome }\sim\mathtt{ [age] + sex + genotype\_cohort + gPC1 + gPC2 + gPC3 + gPC4 + score)}\qquad\text{(no random effect)}
\end{aligned}
```

- The sum has either one score or a whole set.
- $`\mathrm{age}_i`$ (age at measurement) is included only for outcomes that have one, e.g. `AgeChild_GR1203` for online use.
- H1 ($`\boldsymbol{\theta} = \mathbf{0}`$) is the exact nested-model *F*-test. For a single score, the reported estimate is $`\theta_1`$.

### Hypotheses and multiple testing

Every hypothesis tested, with the equation used to test it:

- **H1, one GPS, repeated measures** (CBCL, MRI, DTI). Equation [1], model $`\mathcal{M}_1`$: $`\beta_3 = 0`$. Satterthwaite *F*-test with 1 df.
- **H2, one GPS, repeated measures** (CBCL, MRI, DTI): the GPS association changes with age. Equation [1], model $`\mathcal{M}_2`$: $`\beta_4 = \beta_5 = 0`$. Satterthwaite *F*-test with 2 df.
- **H1, one GPS, single-time outcomes** (SES, diet, online). Equation [3] with one score: $`\theta_1 = 0`$. Nested-model *F*-test.
- **H1, a GPS set jointly, repeated measures** (CBCL, MRI, DTI). Equation [2]: $`\theta_1 = \dots = \theta_K = 0`$. Satterthwaite *F*-test with $`K`$ df.
- **H1, a GPS set jointly, single-time outcomes** (SES, diet, online). Equation [3] with all $`K`$ scores of the set: $`\theta_1 = \dots = \theta_K = 0`$. Nested-model *F*-test.

**Multiple testing.** p-values are corrected with the Benjamini–Hochberg false
discovery rate (FDR), separately within each score set × hypothesis:

- For each single-score set (`PCS`, `UK_PCS`, `PGS`), H1 is one family. It
  covers every score in the set and every outcome, both repeated measures [1]
  and single-time outcomes [3].
- For each single-score set, H2 is another family. It covers every score in
  the set and every repeated measure [1].
- Each jointly tested set (`PCS_combined`, `UK_PCS_combined`) is its own H1
  family, covering every outcome, from [2] and [3].

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

## Figures (step 6)

`R/06_plot_results.R` reads `results.csv` from step 5 and writes, for each
score set tested one score at a time, to `<results_dir>/figures/`:

- `<set>_phewas.png`: −log10 p of every score × outcome test, grouped and
  coloured by domain.
  - Panel a shows H1. Triangles point up or down with the sign of the estimate.
  - Panel b shows H2.
  - Each outcome has one point per component (20 for C1–C20).
  - Component scores `C1`…`C20` are labelled `GPS1`…`GPS20` in the figures
    (`label_score()` in `R/06_plot_results.R`). The pipeline keeps the names
    `C1`…`C20` internally.
  - The red dashed line marks p = 0.05 for a single test. The purple dashed
    line marks the FDR 5% threshold across all scores and outcomes of the
    panel. Points passing FDR are labelled.
  - Domains appear in the order CBCL, SES, Diet, Online, T1w MRI, dMRI, on a white
    background separated by thin dashed grey vertical lines.
- `<set>_heatmap_H1.png`: score × outcome heatmap of the standardised estimate,
  on a diverging blue–grey–red scale.
- `<set>_heatmap_H2.png`: score × outcome heatmap of −log10 p for the GPS × age
  interaction.

In both heatmaps, `*` marks p < 0.05 and `**` marks FDR-corrected p < 0.05.

To try the figures on simulated results (20 components, all outcomes, no model
fitting):

```bash
Rscript toy/make_toy_results.R
```

```bash
Rscript R/06_plot_results.R config/toy_viz.R
```

Requires `ggplot2` and `ggrepel`.

## Verification

```bash
bash toy/run_toy.sh
```

It does three things:

1. It simulates toy data (`toy/make_toy_data.R`): 500 GenR individuals, 480
   genotyped individuals and 3,000 SNPs. The data include SNPs that QC must remove,
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

## Link to the genomica pipeline

The GenomICA components used here are produced by
[sourena-s/genomica](https://github.com/sourena-s/genomica). What that code
does, as far as it affects GenRomica (points to confirm are marked **(check)**):

- **Input GWAS.** Pan-UK Biobank summary statistics, European-ancestry (EUR)
  columns only (`panukb_dl.sh`). For each trait, Z = `beta_EUR / se_EUR`.
  Variants flagged `low_confidence_EUR` are set to missing.
- **Trait selection.** Traits are selected by genetic correlation with
  psychiatric traits (|r_g| > 0.5, `stage_1_filter_gwas.sh`), plus a manually
  curated list.
- **Variant set.** A mask in `clumps.sorted.mask.in.var.catalogue` (column 2 = 1)
  selects the variants that enter the ICA (`h5_mask_concat.py`). The current
  script uses a HapMap3 mask; earlier runs used LD-clumped masks.
- **Missing values.** Missing Z-scores, including low-confidence variants, are
  set to 0 before the ICA (`np.nan_to_num`).
- **ICA.** FSL MELODIC with 20 components, `pow3` non-linearity, variance
  normalisation (`--vn`) and no mixture modelling (`--no_mm`)
  (`stage_2_melodic.sh`). It runs on two trait sets: "manual" and "r2psy".
- **Export.** `stage_6_melodic_to_text.sh` writes each component of
  `melodic_oIC.nii.gz` as one value per masked variant, in catalogue order.

What this means for GenRomica:

1. **Two masks, two score sets.** Both GPS sets probably come from the
   "manual" trait set, decomposed twice with different variant masks:
   - **HapMap3:** common variants, about 1.3 million; used in `h5_mask_concat.py`.
   - **Common and rare variants:** about 5 million. `stage_6_melodic_to_text.sh`
     writes this list as `<catalogue>.5million`.

   Record which mask gives `C1–C20` and which gives `UK_PCS_1–20`. **(check)**
2. **One catalogue per mask.** Step 2 matches component rows to catalogue rows
   by position (`include == 1`). Each score set therefore needs the catalogue of
   its own mask. `config/project.R` has a single `cfg$variant_catalogue` and
   builds only the `C` set. To compute the second set, run steps 2–3 again with
   the other catalogue and its component files, write the output to separate
   files, and add those files to `cfg$score_weights`. The pipeline stops if the
   row counts do not match, so a wrong catalogue is caught.
3. **Rare variants in the 5-million mask.** For rare variants, GenR imputation
   quality is lower and fewer of them will match in `snp_match`. Allele
   frequencies also differ in a multi-ethnic cohort, and their weights come from
   European-ancestry GWAS. No MAF or imputation-quality filter is applied in
   GenRomica (see [Not applied](#not-applied)). Compare the number of matched
   variants per score in `*_snp_matching.csv`, and consider an INFO or MAF
   threshold for the 5-million set.
4. **Open issue: component weights are `melodic_oIC`, not Z-scores.**
   `stage_6_melodic_to_text.sh` exports `melodic_oIC`, the original component
   maps. MELODIC's Z-statistic maps are `melodic_IC`: each value divided by that
   variant's residual noise SD. GenRomica applies no z-transform. Step 2 uses
   the `melodic_oIC` values as the weights (`BETA`) and computes
   P = 2Φ(−|Z|) from them. Consequences:
   - The GPS uses the original weights. A single scaling factor per component
     would not change any result, because the GPS is z-scored before modelling.
     A per-variant conversion such as `melodic_IC` changes the relative SNP
     weights and therefore the GPS.
   - P is only valid on a Z scale. P is used only in QC filter 1 (P > 0), so on
     another scale that filter removes a different set of variants.
   - `component_z_summary.csv` (SD and λ per component) shows how far the
     values are from a Z scale.

   To decide: use `melodic_oIC` or `melodic_IC` as weights, and whether to
   z-transform before computing P. See concern B6.
5. **Component sign is arbitrary.** An ICA component and its negative are the
   same solution. The direction of a GPS association must be read together with
   the sign of each trait's loading in `melodic_mix` (the trait-by-component
   matrix). Use the components and loadings from the same MELODIC run.
6. **Component order.** MELODIC orders components by explained variance, and a
   rerun can change order and sign. Keep the exact run used for scoring.
7. **Effect allele (confirmed).** Pan-UKB variant IDs are `chr_pos_ref_alt`,
   and `beta` refers to the alt allele. The catalogue's A2, the second allele
   column, is the alt allele and sets the sign of the weights.
8. **Genome build.** Pan-UKB positions are GRCh37, so the GenR genotypes must
   also be GRCh37 for step 4 to match on position (concern A5).
9. **Ancestry.** The weights come from European-ancestry GWAS, which supports
   concern B5 (GenR is multi-ethnic).
10. **Missing set to zero.** Variants missing for some traits were given Z = 0
   in the ICA input, which shrinks their component weights towards zero.

## Concerns

These are open for discussion. None has been acted on yet, so that the
original results can be replicated first.

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
   GenR genotypes must use the same build. Pan-UKB, the source of the weights,
   is GRCh37 (see [Link to the genomica pipeline](#link-to-the-genomica-pipeline)).
   This is not checked.
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
9. **Warnings were not recorded.** The original code captured errors but not
   warnings, such as non-convergence. They are now logged in
   `model_warnings.csv` without changing any result. Singular fits are flagged
   in `results.csv` (`singular`), but their p-values are still reported.
10. **Silent skips.** A model with too little variation returns no row and is
    not logged in `failed_models.csv`. Compare the number of result rows with
    the number expected.

### B. Summary statistics and scores

1. **p = 0 is removed.** QC keeps `P > 0`. Because P is computed from Z, P is
   0 only when |Z| > ~37.5 (double-precision underflow). Those SNPs are the
   strongest signals and are dropped. Consider a filter on Z instead of P.
2. **Palindromic SNPs are removed by default.** This costs roughly 15% of SNPs.
   With well-defined strands (same reference panel), `snp_match(strand_flip = FALSE)`
   could keep them. The strand-flip switch in step 4 only helps if they are not
   already removed in step 3.
3. **Lowercase alleles are dropped** rather than upper-cased.
4. **Score construction.** The score is Σ dosage × Z over all SNPs, with no LD
   clumping or shrinkage. Correlated SNPs therefore count several times, and
   Z-scores scale with the discovery study's sample size. The score projects
   genotypes onto the ICA component; LD clumping or shrinkage (e.g. LDpred2)
   would be needed for an optimised PGS.
5. **Ancestry.** GenR is multi-ethnic. If the ICA study used European-ancestry
   data, score performance will differ by ancestry. Four gPCs may not fully
   control stratification; consider restricting to, or stratifying by,
   genetic ancestry.
6. **Open issue: component weights are not Z-transformed.** The component
   files are `melodic_oIC` (original maps), not MELODIC's Z-statistic maps
   (`melodic_IC`), and GenRomica does not z-transform them. They are used as
   weights as they are, and P is computed from them as if they were Z-scores.
   Details in [Link to the genomica pipeline](#link-to-the-genomica-pipeline),
   point 4.

### C. Statistical models

1. **Age centred at 17, the oldest wave.** See the caution box below.

> [!CAUTION]
> **Age centred at 17, the oldest wave.** The centring age changes the
> interpretation and precision of single coefficients in M2. The H1 and H2
> tests, fitted values and R² stay the same.
>
> - **Single coefficients.** In M2, $`\beta_3`$ is the GPS association at the
>   centring age, and $`\beta_4`$ (`age_interaction` in `h2_details.csv`) is the
>   linear GPS × age slope at that age. At 17, the edge of the age range,
>   `age_c` and `age_c²` are strongly correlated and these two estimates are
>   imprecise. $`\beta_5`$ is the same for any centring age.
> - **H1** ($`\beta_3`$ in M1). M1 has no GPS × age terms, so $`\beta_3`$ and its
>   p-value are the same for any centring age.
> - **H2** (joint test $`\beta_4 = \beta_5 = 0`$ in M2). Moving the centre *d*
>   years earlier ($`a = a' + d`$) turns $`(\beta_4, \beta_5)`$ into
>   $`(\beta_4 + 2d\,\beta_5,\ \beta_5)`$. Both pairs are zero together, so the
>   joint *F*-test gives the same result.
> - **Fitted model.** $`(a - c)^2 = a^2 - 2ca + c^2`$, so the age terms and their
>   GPS interactions describe the same model for any centring age *c*. Fitted
>   values, likelihood, variance components and R² are identical.
> - **Toy data check.** Centring at 11 instead of 17 changed H1 and H2
>   p-values, estimates and R² by less than 2×10⁻⁷. The SE of
>   `age_interaction` was 3.4 times larger at 17, and its p-value changed by
>   0.28 (median).
> - **Centring on the average age.** Uncomment `# centre <- mean(d$age)` in
>   `centre_age()` in `R/05_run_regressions.R` to centre on the average age of
>   each model's analysis sample. On the toy data this changed H1 and H2
>   p-values by less than 2×10⁻⁶.
> - **Other standardisation.** z-scoring the outcome, the GPS and the gPCs
>   changes the units of the estimates and leaves p-values unchanged.

2. **Different CBCL instruments.** Wave f05 used the CBCL/1½-5, and later waves
   used the CBCL/6-18. The scales differ in items and range, but are pooled and
   z-scored together. Standardising within wave, or modelling the instrument,
   would be more defensible. See [D. Comparability of scores across waves](#d-comparability-of-scores-across-waves).
3. **Outcome standardisation is pooled over waves.** It does not change the
   p-values, but estimates are in pooled-SD units.
4. **Skewed outcomes.** CBCL sums are counts with many zeros, and the models
   assume normal residuals. Consider a transform, or a count model, for robustness.
5. **Siblings and twins.** GenR includes related individuals. The random intercept
   is per individual only, so related individuals are treated as independent. Keep one
   individual per family or add a family random effect.
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
    income are parent characteristics. An association with the individual's score
    reflects shared parental genetics (genetic nurture, or confounding),
    not an effect on the individual.
13. **Ordinal outcomes treated as continuous** (`INCOME5`, `EDUCM`, `EDUCP`).
    Step 1 stops if they are stored as factors, so their coding has to be
    chosen explicitly.
14. **Timing.** `DietScore_8y` is measured at about age 8 and the online measure
    at about 17. Single-time outcomes are adjusted for age only where an age
    variable is given.

### D. Comparability of scores across waves

Scores whose items, scale version, informant or scanner differ between waves.
Points not yet confirmed are marked **(check)**.

#### Used in the current analysis

**CBCL at f05 vs f09, f13 and f17: not comparable.**
- In Generation R, the f05 wave (about age 6, range 5–7) used the
  preschool version, CBCL/1½–5 (99 items), confirmed in published papers
  ([stability of problem behaviour](https://pmc.ncbi.nlm.nih.gov/articles/PMC4729812),
  [age-6 mental health](https://pmc.ncbi.nlm.nih.gov/articles/PMC4628613)). The
  preschool form was chosen because most individuals were expected to be younger
  than 6 at that visit. The later waves used the school-age version, CBCL/6–18.
  The two versions have different syndrome scales built from different items.
- Item counts are from the CBCL manuals:

| Scale in the pipeline | CBCL/1½–5 (f05) | CBCL/6–18 (f09–f17) | Comparable? |
|---|---|---|---|
| `anx`: anxious/depressed | 8 items | 13 items | No: different items and range |
| `som`: somatic complaints | 7 items | 11 items | No |
| `wit`: withdrawn | "Withdrawn", 8 items | "Withdrawn/Depressed", 8 items | No: same item count, different content |
| `att`: attention problems | 5 items | 10 items | No |
| `agg`: aggressive behaviour | 19 items | 18 items | No: overlapping but not identical items |
| `sop`, `tho`, `rul`: social, thought, rule-breaking | **do not exist** | 11, 15 and 17 items | Not measured at f05 |

- The preschool version also has scales that the school-age version lacks:
  emotionally reactive and sleep problems.
- **The `sum_wid_f05` → `sum_wit_f05` rename is suspect.** The two names
  probably reflect two different scales: preschool "Withdrawn" and school-age
  "Withdrawn/Depressed". If so, the original code's rename merges different
  scales under one name, and the pipeline does the same. **(check** what `wid`
  stands for in the GenR data dictionary.)
- **f09, f13 and f17 are the same instrument**, so their raw sums are
  comparable, provided the informant is the same. **(check** that f17 is still
  parent-reported CBCL. If f17 is self-report (YSR) or an adult form, it is not
  comparable either.)
- **Consequence for the models (see C2).** The outcome is z-scored over all
  waves pooled, and age is modelled as a smooth trajectory. A change of
  instrument at f05 therefore looks like a developmental change, and it is
  confounded with the age terms and the GPS × age interaction (H2).

**MRI and DTI: possibly not comparable across waves (check).**
- The first imaging wave (2009, ages 6–9, ~1,070 individuals) used a different
  scanner from later waves. A dedicated 3T scanner was installed in 2013
  ([White et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5803295)).
- Scanner and acquisition differences shift grey-matter volume, cortical
  thickness and especially FA/MD.
- Check the scanner, sequence, FreeSurfer version and DTI pipeline for each
  wave. If they differ, harmonise the data (e.g. ComBat) or add wave/scanner as
  a covariate.

**Single-time variables: no cross-wave issue.** These are the diet score at 8,
online activity at 17, income, and parental education.

Cross-wave comparability of candidate variables not used yet is covered in
[Candidate GenR variables not used yet](#comparability-of-candidate-variables-across-waves).

#### Options for the CBCL (none applied yet, to keep replicating the original)

1. **TODO — same instrument only:** use only f09–f17 for the longitudinal CBCL
   models and drop f05 (analyse it separately). Uncomment line 106 in
   `config/project.R`:
   `cfg$longitudinal <- cfg$longitudinal[!(domain == "CBCL" & wave == "f05")]`.
2. **Standardise within each wave before pooling.** This removes level
   differences between instruments, but also removes real change with age.
3. **Common items only:** use only the items shared by both versions,
   re-scored at every wave.
4. **Age-normed T-scores:** comparable in meaning, but not on the raw scale.

## Candidate GenR variables not used yet

Each variable has been reported in at least one published Generation R paper,
cited in square brackets (see [References for candidate variables](#references-for-candidate-variables)).
Each is tagged with its data type: *continuous*, *count*, *ordinal*, *binary*
or *categorical*. The Generation R data dictionary is not public, so check exact
variable names, waves and sample sizes with Generation R data management before
adding a variable. A new variable is added as a row in `cfg$longitudinal`
(repeated measures) or `cfg$single_outcomes` (measured once) in
`config/project.R`.

### Behaviour and psychiatric traits

- **Autistic traits:** Social Responsiveness Scale (SRS), 18-item short form, at about 6 years and at about 13.5 years. The short form correlates 0.93–0.99 with the full 65-item SRS. *Continuous.* [1, 2, 3]
- **Executive function:** BRIEF-P (parent report) at age 4, plus ADHD diagnoses from the DISC-YC interview in individuals who screened positive. *Continuous* and *binary.* [4]
- **Bullying and peer victimisation:** teacher report at about 6.6 years; peer nomination and self-report at about 7.5 years; parent and teacher report at about 8 years. *Count* and *binary.* [5, 6]
- **Hallucinatory (psychotic-like) experiences:** self-report at about 10 years. *Ordinal.* [7]

### Cognition and education

- **Nonverbal IQ:** Snijders-Oomen nonverbal intelligence test (SON-R) at about 6 years. *Continuous.* [8]
- **Cognition at about 13 years:** four WISC-V subtests (vocabulary, matrix reasoning, digit span, coding) and an estimated full-scale IQ. *Continuous.* [9]
- **Cito end-of-primary-school test** (language and mathematics) at about 12 years. *Continuous;* left-skewed with a ceiling effect, so a transform may be needed. [9]

### Brain imaging

- **Subcortical and regional volumes:** e.g. hippocampus. *Continuous.* [10]
- **Total brain volume and other global measures** from the same scans. *Continuous.* [11]
- **Resting-state functional connectivity (fMRI):** at ages 6–10, and longitudinally at ages 10 and 14. *Continuous.* [12, 13, 14]

### Body composition and physical health

- **Body composition:** total and visceral fat at about 13 years, and body-fat indicators at about 6 years. *Continuous.* [15, 16]
- **Lung function:** spirometry at about 10 years. *Continuous.* [17]
- **Blood pressure and cardiometabolic markers:** blood pressure at ages 6 and 9–10, a cardiometabolic risk score at age 8, and lipids and insulin at age 6. *Continuous.* [18, 19, 20]
- **Oral and dental health:** collected at ages 3, 6, 9, 13 and 17. *Count* and *binary.* [21]

### Lifestyle

- **Screen time at about 6 years:** TV and computer-game use (parent report). *Continuous* (hours/day). [22, 16]
- **Screen time at about 13 years:** self-report. *Continuous* (hours/day). [15]
- **Physical activity:** sports participation, outdoor play and active travel to school at about 6 years (parent report); self-report at about 13 years. *Continuous* and *binary.* [16, 15]
- **Accelerometer-measured physical activity:** only found at age 2 (uniaxial ActiGraph, n = 347). *Continuous.* No accelerometer data were found for later ages. [23]
- **Diet components:** snacks and sugar-sweetened drinks at about 6 years, and the individual food groups of the age-8 food-frequency questionnaire. *Continuous* (portions/day). [22, 24]
- **Substance use:** alcohol use (self-report, adolescence) and age at first alcohol or tobacco use (at about 14 years). *Binary* and *continuous.* [25, 26]
- **Sleep:** sleep problems in early childhood; sleep and circadian rhythm in the age-9 protocol. *Continuous* and *ordinal.* [27, 28]

### Socioeconomic status and family

- **Parental psychopathology:** Brief Symptom Inventory (BSI) for mothers and fathers, in pregnancy and in childhood, including the hostility subscale. *Continuous.* [29]
- **Harsh parental discipline:** maternal and paternal report. *Continuous.* [30]
- **Family functioning** (parent report). *Continuous.* [29]
- **Neighbourhood socioeconomic status**, used in the bullying analyses. *Continuous.* [5]

### Prenatal exposures

- **Maternal smoking in pregnancy:** never, first trimester only, or continued. *Categorical.* [31, 32]
- **Maternal alcohol use in pregnancy**, coded in the same way. *Categorical.* [33]
- **Paternal smoking.** *Binary.* [32]

### Environmental exposures (urban exposome)

About 30 exposures, from pregnancy and childhood, at the home address. *Continuous* unless noted. [17, 10]

- air pollution;
- road-traffic noise over 24 hours, and night-time noise at home;
- green space: NDVI within 100 m, and major green space within 300 m (*binary*);
- blue space;
- walkability;
- land use.

### Comparability of candidate variables across waves

See also concern D for the variables already used.

- **Autistic traits (SRS):** the same 18-item short form was used at about 6 and at about 13.5 years, so scores are comparable in items. Check the informant at each wave. [2, 3]
- **Cognition:** not comparable (different tests).
  - At about 6y: SON-R nonverbal IQ. [8]
  - At about 13y: four WISC-V subtests (vocabulary, matrix reasoning, digit span, coding), combined into an estimated full-scale IQ. [9]
- **Executive function:** BRIEF-P (preschool version, parent report) at 4y only. No later BRIEF wave was found in published Generation R papers, so there is nothing to compare it with. [4]
- **Screen time and physical activity:** different informants.
  - At about 6y: parent questionnaire on TV viewing, computer-game use, sports participation, outdoor play and active travel to school. [16]
  - At about 13y: self-report questionnaires, n = 3,258. [15]
  - Whether the items match could not be checked from the published papers; check the questionnaires.
- **Bullying:** different informants at different ages; no adolescent wave was found.
  - About 6.6y: teacher report, n = 3,192.
  - About 7.5y: peer nomination and self-report, n ≈ 1,100–1,300. [5]
  - About 8y: parent and teacher report, n = 2,602. [6]

### References for candidate variables

1. Autistic traits and peer relationships in Generation R (dataset record). [https://pure.eur.nl/en/datasets/childrens-autistic-traits-and-peer-relationships-do-non-verbal-iq/](https://pure.eur.nl/en/datasets/childrens-autistic-traits-and-peer-relationships-do-non-verbal-iq/)
2. SRS short form at about 6 years in Generation R. *Molecular Autism*, 2022. [https://link.springer.com/10.1186/s13229-022-00507-4](https://link.springer.com/10.1186/s13229-022-00507-4)
3. SRS short form at about 13.5 years in Generation R. *Molecular Autism*, 2023. [https://link.springer.com/10.1186/s13229-023-00578-x](https://link.springer.com/10.1186/s13229-023-00578-x)
4. Preschool executive functioning (BRIEF-P) and ASD/ADHD traits. *Child and Adolescent Psychiatry and Mental Health*, 2019. [https://link.springer.com/10.1186/s13034-019-0299-7](https://link.springer.com/10.1186/s13034-019-0299-7)
5. Verlinden-Bondaruk M. Bullying in early elementary school (PhD thesis, Erasmus University Rotterdam), 2014. [https://repub.eur.nl/pub/77215/141202_Verlinden-Bondaruk-Maryna.pdf](https://repub.eur.nl/pub/77215/141202_Verlinden-Bondaruk-Maryna.pdf)
6. Bullying involvement and brain structure in Generation R. *Frontiers in Psychiatry*, 2019. [https://www.frontiersin.org/journals/psychiatry/articles/10.3389/fpsyt.2019.00696/pdf](https://www.frontiersin.org/journals/psychiatry/articles/10.3389/fpsyt.2019.00696/pdf)
7. Steenkamp LR, *et al.* Bullying, peer rejection and hallucinatory experiences, 2021. [https://repub.eur.nl/pub/134983](https://repub.eur.nl/pub/134983)
8. Mental health at 6 and IQ–achievement discrepancy (SON-R nonverbal IQ). [https://pmc.ncbi.nlm.nih.gov/articles/PMC9527951/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9527951/)
9. Neurocognition (WISC-V at 13) and Cito achievement among bereaved individuals in Generation R, 2024. [https://pmc.ncbi.nlm.nih.gov/articles/PMC11391039](https://pmc.ncbi.nlm.nih.gov/articles/PMC11391039)
10. Urban environment in early life and brain morphology in preadolescents (Generation R). [https://pure.eur.nl/en/publications/urban-environment-in-early-life-and-brain-morphology-in-preadoles/](https://pure.eur.nl/en/publications/urban-environment-in-early-life-and-brain-morphology-in-preadoles/)
11. White T, *et al.* Paediatric population neuroimaging and the Generation R Study: the second wave. *Eur J Epidemiol*, 2018. [https://pmc.ncbi.nlm.nih.gov/articles/PMC5803295](https://pmc.ncbi.nlm.nih.gov/articles/PMC5803295)
12. Developmental changes in dynamic functional connectivity from childhood into adolescence. *Frontiers in Systems Neuroscience*, 2021. [https://www.doi.org/10.3389/FNSYS.2021.724805](https://www.doi.org/10.3389/FNSYS.2021.724805)
13. Subcortical functional connectivity at ages 6–10 in Generation R. [https://core.ac.uk/download/160828747.pdf](https://core.ac.uk/download/160828747.pdf)
14. Functional connectivity and psychiatric symptoms in ABCD and Generation R. [https://pmc.ncbi.nlm.nih.gov/articles/PMC11332032](https://pmc.ncbi.nlm.nih.gov/articles/PMC11332032)
15. Physical activity, screen time and body composition in 13-year-old adolescents: the Generation R Study. *Pediatric Obesity*, 2023. [https://pure.eur.nl/en/publications/physical-activity-screen-time-and-body-composition-in-13-year-old/](https://pure.eur.nl/en/publications/physical-activity-screen-time-and-body-composition-in-13-year-old/)
16. Wijtzes AI. Lifestyle behaviours and body fat in childhood (PhD thesis, Erasmus University Rotterdam), 2015. [https://repub.eur.nl/pub/78112/150513_Wijtzes-Anne-Irene.pdf](https://repub.eur.nl/pub/78112/150513_Wijtzes-Anne-Irene.pdf)
17. Abellan A, *et al.* Urban environment during pregnancy and lung function, wheezing and asthma in school-age children. *Environmental Pollution*, 2024. [https://repositori.upf.edu/items/75d98624-ed04-4297-a6ca-e67dbb27bad5/full](https://repositori.upf.edu/items/75d98624-ed04-4297-a6ca-e67dbb27bad5/full)
18. Prediction of high blood pressure at 9–10 years from birth to 6 years (Generation R). [https://repub.eur.nl/pub/112552/Repub_112552_O-A.pdf](https://repub.eur.nl/pub/112552/Repub_112552_O-A.pdf)
19. Diet quality and cardiometabolic health in childhood: the Generation R Study. [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8854322/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8854322/)
20. Infant fat mass and cardiometabolic risk factors at 6 years (Generation R). [https://pmc.ncbi.nlm.nih.gov/articles/PMC5426532](https://pmc.ncbi.nlm.nih.gov/articles/PMC5426532)
21. Oral and craniofacial research in the Generation R Study. *Clinical Oral Investigations*, 2023. [https://link.springer.com/article/10.1007/s00784-023-05076-1](https://link.springer.com/article/10.1007/s00784-023-05076-1)
22. Lifestyle behaviour clusters at age 6 (screen time, activity, snacks, sugary drinks). [https://pmc.ncbi.nlm.nih.gov/articles/PMC7352876](https://pmc.ncbi.nlm.nih.gov/articles/PMC7352876)
23. Correlates of physical activity in 2-year-old toddlers: the Generation R Study. *J Pediatr*, 2013. [https://research.tilburguniversity.edu/en/publications/correlates-of-physical-activity-in-2-year-old-toddlers-the-genera/](https://research.tilburguniversity.edu/en/publications/correlates-of-physical-activity-in-2-year-old-toddlers-the-genera/)
24. Diet quality in childhood (age-8 food-frequency questionnaire), Generation R. [https://repub.eur.nl/pub/121454/121454_Chapter_2_4-Diet_quality_in_childhood.pdf](https://repub.eur.nl/pub/121454/121454_Chapter_2_4-Diet_quality_in_childhood.pdf)
25. Alcohol use and brain morphology in adolescence: a longitudinal study in three cohorts. [https://dspace.library.uu.nl/handle/1874/443370](https://dspace.library.uu.nl/handle/1874/443370)
26. Brain morphology in late childhood and initiation of alcohol and tobacco use (Generation R). [https://pure.eur.nl/en/publications/examining-associations-between-brain-morphology-in-late-childhood/](https://pure.eur.nl/en/publications/examining-associations-between-brain-morphology-in-late-childhood/)
27. Early sleep problems and autistic traits in Generation R. *Child and Adolescent Psychiatry and Mental Health*, 2021. [https://capmh.biomedcentral.com/track/pdf/10.1186/s13034-021-00421-2](https://capmh.biomedcentral.com/track/pdf/10.1186/s13034-021-00421-2)
28. Generation R Focus at 9 study protocol (onderzoekmetmensen.nl). [https://onderzoekmetmensen.nl/en/node/41521/pdf](https://onderzoekmetmensen.nl/en/node/41521/pdf)
29. Parental psychopathology, family functioning and child problems in Generation R. [https://pmc.ncbi.nlm.nih.gov/articles/PMC3135831](https://pmc.ncbi.nlm.nih.gov/articles/PMC3135831)
30. Mackenbach MD, *et al.* Harsh parental discipline and child problem behaviour. *PLoS ONE*, 2014. [https://doi.org/10.1371/JOURNAL.PONE.0104793](https://doi.org/10.1371/JOURNAL.PONE.0104793)
31. Zou R, *et al.* Maternal tobacco use during pregnancy and preadolescent brain morphology. [https://dspace.library.uu.nl/handle/1874/421990](https://dspace.library.uu.nl/handle/1874/421990)
32. Roza SJ, *et al.* Maternal smoking during pregnancy and child behaviour problems: the Generation R Study. *Int J Epidemiol*, 2009. [https://pure.eur.nl/en/publications/maternal-smoking-during-pregnancy-and-child-behaviour-problems-th/](https://pure.eur.nl/en/publications/maternal-smoking-during-pregnancy-and-child-behaviour-problems-th/)
33. Kooijman MN, *et al.* The Generation R Study: design and cohort update 2017. *Eur J Epidemiol*, 2016. [https://pmc.ncbi.nlm.nih.gov/articles/PMC5233749](https://pmc.ncbi.nlm.nih.gov/articles/PMC5233749)
