# microbiome-power-prototype

Two Shiny applets for feature-level microbiome power planning using only Human Microbiome Project phase 1-II reference profiles.

## Run

From this directory, with R 4.5 or later:

```sh
Rscript scripts/install_dependencies.R
Rscript scripts/run_app.R case-control
Rscript scripts/run_app.R continuous
```

Open http://127.0.0.1:3850 for case-control or http://127.0.0.1:3851 for continuous outcome. Run each app in its own terminal. Precomputed summaries are included; extraction is not needed to run an app.

- **Case-control:** independent two-group comparison, adjustable number of hypotheses, Bonferroni correction, SD quartiles, detectable transformed mean differences, and copyable grant text. This app does not model covariate adjustment.
- **Continuous outcome:** feature abundance regressed on a standardized continuous outcome plus covariates. Includes cross-sectional and explicitly assumption-based two-visit repeated-measures modes. The latter estimates an overall association, not a within-person change-only effect.

## Reference data

The app offers five distributions: stool EC, stool SGB, oral SGB, skin SGB, and vaginal SGB. All use the same feature filters: prevalence ≥10% and mean relative abundance >0.01% among positive observations. Profiles are normalized before filtering. SDs exclude zeros. The app and grant summaries report detected and retained feature counts; excluded features are omitted from the default testing family. These thresholds are modeling choices, not universal requirements for linear models.

The reference includes 548 stool taxonomy, 1,284 oral taxonomy, 320 skin taxonomy, 233 vaginal taxonomy, and 547 stool EC profiles. Profiles may include repeated visits and are not counts of independent participants.

## Rebuild and test

```sh
python3 scripts/unpack_training.py
Rscript scripts/rebuild_training.R
Rscript tests/run_tests.R
```

Python uses only its standard library. Rebuilding reads large extracted tables and may take several minutes. Original inputs remain zipped; SHA-256 checks are in `data/manifest.json`. Extraction is ignored by Git.

## Interpretation

SDs are calculated on arcsine-square-root relative abundance among positive observations. They are planning proxies, not validated zero-inclusive regression residual SDs. Sample sizes must correspond to observations actually analyzed. Sparse features may require a carriage or simulation-based model. Back-transformed abundance examples are fitted-value illustrations, not arithmetic mean differences. Power is per feature, not joint power across the testing family. Covariate R² in the continuous app concerns the outcome's association with covariates; it is not a microbiome PERMANOVA R².

See [methods](docs/METHODS.md) and [data provenance](docs/PROVENANCE.md).
