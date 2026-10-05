# Validation

The test runner checks:

- Both applets across every HMP-only training option, including generated grant summaries.
- Case-control, cross-sectional regression, and repeated-measures detectable-effect inversion.
- Regression agreement with the two-sided noncentral-t coefficient test.
- Explicit rejection of zero or one control observation.
- HMP-only input selection and exclusion of rare features from the default test family.

Archive contents and SHA-256 checks are validated by the extraction script. Training reproducibility is checked by rebuilding the common and filtered summaries from the archived inputs and comparing them numerically with the exported summaries. Reference sample counts and feature filters are described in README.md.

These checks verify implementation consistency, not the empirical adequacy of the positive-only SD assumption or the normal-model approximation for sparse microbiome data.
