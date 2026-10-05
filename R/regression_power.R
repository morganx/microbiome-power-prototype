# Fixed-design Gaussian linear regression, one two-sided coefficient test.
# x is the continuous outcome standardized to sample SD 1 before adjustment.
validate_regression_design <- function(n, covariates, alpha, outcome_r2, sd_multiplier) {
  stopifnot(length(n) == 1, is.finite(n), n == floor(n),
    length(covariates) == 1, is.finite(covariates), covariates >= 0,
    covariates == floor(covariates), n > covariates + 2,
    is.finite(alpha), alpha > 0, alpha < 1,
    is.finite(outcome_r2), outcome_r2 >= 0, outcome_r2 < 1,
    is.finite(sd_multiplier), sd_multiplier > 0)
}
regression_power <- function(beta, feature_sd, n, covariates = 0, alpha = 0.05,
                             outcome_r2 = 0, sd_multiplier = 1) {
  validate_regression_design(n, covariates, alpha, outcome_r2, sd_multiplier)
  stopifnot(all(is.finite(beta)), all(is.finite(feature_sd)), all(feature_sd > 0))
  df <- n - covariates - 2
  sxx <- (n - 1) * (1 - outcome_r2)
  sigma <- feature_sd * sd_multiplier
  ncp_squared <- (beta / sigma)^2 * sxx
  # F(1, df, lambda^2) is the square of the noncentral coefficient t statistic.
  critical <- qf(alpha, df1 = 1, df2 = df, lower.tail = FALSE)
  pf(critical, df1 = 1, df2 = df, ncp = ncp_squared, lower.tail = FALSE)
}
regression_detectable <- function(feature_sd, n, covariates = 0, alpha = 0.05,
                                  target = 0.8, outcome_r2 = 0, sd_multiplier = 1) {
  validate_regression_design(n, covariates, alpha, outcome_r2, sd_multiplier)
  stopifnot(is.finite(target), target > alpha, target < 1,
            all(is.finite(feature_sd)), all(feature_sd > 0))
  f <- function(standardized_beta) regression_power(standardized_beta, 1, n,
    covariates, alpha, outcome_r2, 1) - target
  upper <- 1
  while (f(upper) < 0) {
    upper <- upper * 2
    if (upper > 1e6) stop("The requested effect cannot be resolved for this design.")
  }
  effect <- uniroot(f, c(0, upper), tol = 1e-10)$root
  feature_sd * sd_multiplier * effect
}

# Balanced, two-visit random-intercept model. Large-sample GLS/Wald planning
# approximation: variance components are assumed, not estimated in the power formula.
# feature_sd*sd_multiplier is marginal per-visit SD after fixed effects (b_i + e_it).
# outcome_rho is correlation of residualized outcome between visits, after nuisance terms.
mixed_information <- function(n, covariates = 2, outcome_r2 = .05,
                              feature_rho = .5, outcome_rho = .7) {
  validate_regression_design(n, covariates, .05, outcome_r2, 1)
  stopifnot(n > covariates + 3, is.finite(feature_rho), feature_rho >= 0,
    feature_rho < 1, is.finite(outcome_rho), abs(outcome_rho) <= 1)
  # Information for the outcome coefficient after GLS projection off intercept,
  # time and covariates, under the documented balanced predictor design.
  ((2*n - 1)/2) * (1 - outcome_r2) *
    ((1 + outcome_rho)/(1 + feature_rho) + (1 - outcome_rho)/(1 - feature_rho))
}
mixed_power <- function(beta, feature_sd, n, covariates = 2, alpha = .05,
                        outcome_r2 = .05, sd_multiplier = 1,
                        feature_rho = .5, outcome_rho = .7) {
  validate_regression_design(n,covariates,alpha,outcome_r2,sd_multiplier)
  stopifnot(all(is.finite(beta)), all(is.finite(feature_sd)), all(feature_sd > 0))
  info <- mixed_information(n,covariates,outcome_r2,feature_rho,outcome_rho)
  lambda <- abs(beta) * sqrt(info) / (feature_sd * sd_multiplier)
  critical <- qnorm(alpha/2, lower.tail=FALSE)
  pnorm(-critical-lambda) + pnorm(critical-lambda, lower.tail=FALSE)
}
mixed_detectable <- function(feature_sd, n, covariates = 2, alpha = .05,
                             target = .8, outcome_r2 = .05, sd_multiplier = 1,
                             feature_rho = .5, outcome_rho = .7) {
  validate_regression_design(n,covariates,alpha,outcome_r2,sd_multiplier)
  stopifnot(target > alpha, target < 1, all(is.finite(feature_sd)), all(feature_sd > 0))
  info <- mixed_information(n,covariates,outcome_r2,feature_rho,outcome_rho)
  critical <- qnorm(alpha/2, lower.tail=FALSE)
  lambda <- uniroot(function(z) pnorm(-critical-z) + pnorm(critical-z,lower.tail=FALSE)-target,
    c(0, critical+qnorm(target)+10),tol=1e-10)$root
  lambda * feature_sd * sd_multiplier / sqrt(info)
}
# Shared entry points used by UI, summary, plots, sensitivity tables and downloads.
analysis_power <- function(beta, feature_sd, n, covariates=2, alpha=.05,
                           outcome_r2=.05, sd_multiplier=1, mode="cross",
                           feature_rho=.5, outcome_rho=.7) {
  stopifnot(mode %in% c("cross","mixed"))
  if (mode == "mixed") mixed_power(beta,feature_sd,n,covariates,alpha,outcome_r2,
    sd_multiplier,feature_rho,outcome_rho) else
    regression_power(beta,feature_sd,n,covariates,alpha,outcome_r2,sd_multiplier)
}
analysis_detectable <- function(feature_sd,n,covariates=2,alpha=.05,target=.8,
                                outcome_r2=.05,sd_multiplier=1,mode="cross",
                                feature_rho=.5,outcome_rho=.7) {
  stopifnot(mode %in% c("cross","mixed"))
  if(mode == "mixed") mixed_detectable(feature_sd,n,covariates,alpha,target,
    outcome_r2,sd_multiplier,feature_rho,outcome_rho) else
    regression_detectable(feature_sd,n,covariates,alpha,target,outcome_r2,sd_multiplier)
}
