# power-prototype

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

| Profile | Samples | Common features (≥20% prevalence) |
|---|---:|---:|
| Stool taxonomy | 548 | 260 |
| Oral taxonomy | 1,284 | 231 |
| Skin taxonomy | 320 | 26 |
| Vaginal taxonomy | 233 | 9 |
| Stool EC functions | 547 | 1,513 |

Counts describe the supplied processed training tables. Profiles may include repeated visits and are not counts of independent participants. Common and more stringently filtered reference choices are available. Taxonomy filtering requires >10% prevalence and >0.1% mean relative abundance among positives; EC filtering requires >80% prevalence.

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
