# microbiome-power-prototype

## Introduction

This repository contains two Shiny applets for planning studies of associations between metagenomic features—species-level genome bins (SGBs) or enzyme commission (EC) abundances—and outcomes of interest. One applet compares independent case and control groups; the other models associations with a continuous outcome, with options for covariate adjustment and repeated measurements.

The applets use Human Microbiome Project phase 1-II profiles to estimate feature variability at the body site of interest. After feature filtering, relative abundances are arcsine-square-root transformed, and feature-specific standard deviations (SDs) are calculated among positive observations. Power can be explored across the resulting SD distribution, with the 75th percentile providing a higher-variability planning scenario.

For a specified sample size and target power, the applets estimate the minimum detectable difference between group means or the minimum detectable regression slope. Bonferroni correction is based on the specified number of features and tests per feature. Covariates and repeated measurements enter the continuous-outcome power calculation separately.

To make transformed effect sizes easier to interpret, the applets provide illustrative relative-abundance changes at a reference abundance. These are back-transformed fitted values, not exact differences in arithmetic mean abundance. The continuous-outcome applet also allows users to explore assumptions about how much outcome variation is explained by covariates. All estimates are conditional on the selected reference variability and modeling assumptions.

## Development

This prototype was developed after I identified limitations in existing approaches used to plan microbiome studies, including outdated reference distributions, stool-only assumptions, inappropriate handling of unequal groups and longitudinal designs, and inflexible multiple-testing assumptions. I selected the reference populations and feature-specific filtering criteria, specified the statistical behavior and user requirements, and directed and validated the implementation through iterative AI-assisted development.

## How to set up

Install R 4.5 or later, then open a terminal in this repository's directory. Install the required R packages once:

```sh
Rscript scripts/install_dependencies.R
```

Start the applet for your study design:

| Applet | Command | Open in your browser |
|---|---|---|
| Independent cases and controls | `Rscript scripts/run_app.R case-control` | http://127.0.0.1:3850 |
| Continuous outcome | `Rscript scripts/run_app.R continuous` | http://127.0.0.1:3851 |

Keep the terminal running while using the applet; press Ctrl+C to stop it. To run both applets simultaneously, use a separate terminal for each. Precomputed training summaries are included, so you do not need to extract or rebuild the training data to use either applet.

In the applet:

1. Select the body site and feature type that best represent your planned study.
2. Enter the number of independent participants. For the repeated-measures option, enter the number of participants with both visits, not the total number of specimens.
3. Set the number of features you plan to test and the number of tests per feature in the correction family. The default feature count is the number retained in the selected reference dataset; you can replace it with the size of a prespecified hypothesis set.
4. Choose your target power (default 80%) and family-wise significance level (default 0.05).
5. For continuous outcomes, specify the adjustment covariates and explore the outcome-variance and residual-SD assumptions. For repeated measurements, also specify the assumed correlations between visits.

The results and copyable grant text update with your settings.

## How to interpret

### What effect is being estimated?

**Case-control:** the detectable difference between group means on the arcsine-square-root abundance scale. This calculation assumes independent groups and a common within-group SD. The case-control applet does not currently model covariate adjustment or repeated measurements.

**Continuous outcome:** the detectable regression slope, expressed as a change in transformed feature abundance per one-standard-deviation increase in the outcome. The model treats feature abundance as the response and the continuous outcome of interest as a predictor, alongside adjustment covariates.

**Repeated-measures continuous outcome:** an overall association across two visits, accounting for a participant-specific random intercept. This is not a test of within-person change alone or a treatment-by-time interaction. Visit correlations are user-specified assumptions, not estimates obtained from the reference data.

### Choosing a variability scenario

Results are shown at the 25th, 50th, and 75th percentiles of the reference feature-SD distribution. The 75th percentile is a higher-variability scenario, not a worst-case bound: some features have larger SDs. These reference distributions are useful for planning across a set of features, but they do not establish the variability of a particular feature in your study population.

SDs are calculated among positive observations. They are planning proxies, not validated zero-inclusive regression residual SDs. Participant counts must correspond to the observations actually included in the fitted model. For sparse features or hypotheses about acquisition or loss of carriage, a binary-outcome or simulation-based analysis may be more appropriate.

### Multiple testing

The per-test significance threshold is:

```text
family-wise alpha / number of features / tests per feature
```

For example, five prespecified features with one test each give a threshold of 0.05/5 = 0.01. The covariate count and number of visits do not automatically multiply this correction in these applets. If you plan separate tests at each visit and want to correct them together, include those tests explicitly in the tests-per-feature setting.

Reported power is per feature, not the probability of detecting every feature or at least one feature in the testing family.

### Covariates and repeated measurements

In the continuous-outcome applet, the covariate count is the number of adjustment coefficients; a categorical covariate with K levels usually contributes K−1 coefficients. The covariate R² setting describes the fraction of variation in the continuous outcome explained by adjustment covariates. Higher R² leaves less independent outcome variation with which to estimate its association with feature abundance. This is not a microbiome PERMANOVA R² or a measure of variance explained in the microbial feature.

The residual-SD multiplier separately describes unexplained feature variability relative to the training SD. Repeated-measures calculations additionally depend on the assumed feature intraclass correlation and correlation between adjusted outcome measurements. Explore plausible settings when these quantities are unknown.

### Relative-abundance examples and grant text

The abundance examples translate a transformed difference or slope into an illustrative change from a specified reference abundance. The same transformed effect corresponds to different absolute and fold changes at different reference abundances. These examples are not exact changes in arithmetic mean abundance on the original scale.

The generated grant paragraph records the selected training data, filtering criteria, study design, correction, and power assumptions. Review these against the planned analysis before using the text.

## Reference data

Reference metagenomes are from the expanded Human Microbiome Project: Lloyd-Price et al. (2017), *Strains, functions and dynamics in the expanded Human Microbiome Project*, [PMID: 28953883](https://pubmed.ncbi.nlm.nih.gov/28953883/).

Taxonomic profiles were generated with MetaPhlAn 4.0.6 and functional profiles with HUMAnN 4, as specified by the data provider. EC analyses use unstratified EC abundances.

| Distribution | Profiles | Detected features | Retained features |
|---|---:|---:|---:|
| Stool EC | 547 | 2,487 | 1,017 |
| Stool SGB | 548 | 1,943 | 390 |
| Oral SGB | 1,284 | 968 | 303 |
| Skin SGB | 320 | 951 | 59 |
| Vaginal SGB | 233 | 543 | 31 |

SGB features require prevalence ≥10% and mean relative abundance >0.01% among positive observations. Stool ECs require detection in ≥80% of samples, with no abundance threshold. Full profiles are normalized before filtering; the retained subset is not renormalized. Excluded features are omitted from the SD distribution and default testing family. These thresholds are modeling choices to reduce dimensionality, not universal requirements for linear models.

Profiles may include repeated visits and are not counts of independent participants. The supplied processed tables and their provenance are described in [data provenance](docs/PROVENANCE.md).

## Rebuilding and testing

To reproduce training summaries from the compressed inputs:

```sh
python3 scripts/unpack_training.py
Rscript scripts/rebuild_training.R
```

Python uses only its standard library. Rebuilding reads large extracted tables and may take several minutes. Original inputs remain zipped; SHA-256 checks are in `data/manifest.json`. Extracted files are ignored by Git.

To run the calculation and applet checks:

```sh
Rscript tests/run_tests.R
```

See [methods](docs/METHODS.md) and [validation](docs/VALIDATION.md) for calculation details and test coverage.

## Research-use disclaimer

**Research prototype—not a validated statistical product.** These tools are under active development and may contain errors. Power estimates depend on the selected reference data and modeling assumptions and may not generalize to a particular study. Users are responsible for independently verifying calculations, assessing whether the methods suit their study, and interpreting the results. The software is provided “as is,” without warranty; its use does not imply endorsement of any study design or conclusions.
