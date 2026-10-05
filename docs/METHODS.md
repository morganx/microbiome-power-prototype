# Methods

Training profiles are matched to HMP metadata by public sample IDs. Taxonomy uses terminal SGB rows by body area; stool ECs use unstratified rows. Nonempty full profiles are normalized before feature filtering. Positive-only arcsine-square-root SDs are summarized at the 25th, 50th, and 75th percentiles. The common reference uses prevalence ≥20%; rare features are excluded from the default hypothesis family. The filtered reference uses >10% prevalence plus nonzero mean >0.1% for SGBs, and >80% prevalence for ECs. Excluded-feature counts remain available.

Case-control power uses the two-sided noncentral-t two-sample calculation through pwr.t2n.test, with separate case and control counts and Bonferroni alpha = family alpha / features / tests per feature. Zero controls are rejected. No covariate adjustment is modeled in this app.

Cross-sectional regression power uses a noncentral-t coefficient test with n − k − 2 residual degrees of freedom and noncentrality beta × sqrt((n−1)(1−R²)) / residual SD. Outcome R² describes variance in the standardized continuous predictor explained by adjustment covariates. Residual feature SD is specified by a training-SD multiplier.

Repeated-measures mode uses a balanced two-visit random-intercept model and large-sample Wald approximation with separately assumed feature ICC and adjusted outcome correlation. Those correlations are not estimated from HMP by this app. It tests the overall outcome association, assumes a common within/between coefficient, and does not test a treatment-by-time interaction. See the app's calculation details and R/regression_power.R.

Grant text is conditional on these assumptions. The 75th percentile is a higher-variability scenario, not a worst-case bound. Inverse transformations use sin(asin(sqrt(p)) + effect)^2 only within the valid transformed range.
