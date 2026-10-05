make_regression_grant_text <- function(t, s, target, features, analyses, family_alpha,
                                       outcome_name, covariate_names, reference_percent,
                                       mode = "cross", rho = .5) {
  stopifnot(nrow(t) == 1, mode %in% c("cross","mixed"), s$mode == mode,
    is.finite(reference_percent), reference_percent > 0, reference_percent < 100)
  beta <- do.call(analysis_detectable, c(list(feature_sd = t$sd_75, target = target), s))
  type <- if (grepl("ECs",t$dataset)) "EC" else "taxonomic"
  intro <- training_grant_intro(t)
  design <- if (mode == "cross") sprintf("In a cross-sectional analysis of %s participants, adjusting for %s (%s covariate coefficients)",s$n,covariate_names,s$covariates) else
    sprintf("In a longitudinal analysis of %s participants measured at two visits (%s observations), using a linear mixed model adjusted for %s and time with a participant-specific random intercept to account for repeated measurements",s$n,2*s$n,covariate_names)
  correction <- sprintf("Bonferroni correction across %s features and %s test%s per feature (family-wise significance level %s)",
    features,analyses,if(analyses==1) "" else "s",format(family_alpha))
  effect <- sprintf("%s and applying %s, we estimate approximately %.0f%% power to detect a difference in average transformed %s feature abundance of %.5f units per one-standard-deviation higher %s.",
    design,correction,100*target,type,beta,outcome_name)
  angle <- asin(sqrt(reference_percent/100)) + beta
  if (angle <= pi/2) {
    p1 <- 100*sin(angle)^2
    effect <- sprintf("%s and applying %s, we estimate approximately %.0f%% power to detect an association between a one-standard-deviation higher %s and a difference in average transformed %s feature abundance equivalent to a fitted relative abundance of approximately %.3f%% versus a reference of %g%% (beta = %.5f on the transformed scale).",
      design,correction,100*target,outcome_name,type,p1,reference_percent,beta)
    if(mode == "cross") effect <- sprintf("%s and applying %s, we estimate %.0f%% power to detect an association equivalent to a %s feature increasing from %g%% to approximately %.3f%% relative abundance per one-standard-deviation increase in %s (beta = %.5f on the transformed scale).",
      design,correction,100*target,type,reference_percent,p1,outcome_name,beta)
    illustration <- "The abundance example is an illustrative back-transformed fitted value, not an arithmetic mean on the raw abundance scale."
  } else illustration <- "The positive shift at the reference abundance exceeds the valid transformed-abundance range; no abundance conversion is reported."
  assumption <- sprintf("This estimate uses the 75th-percentile training SD to specify an assumed %s SD of %.5f and assumes the fixed adjustment terms explain %.0f%% of variation in %s.",
    if(mode=="mixed") "unexplained per-visit" else "residual",t$sd_75*s$sd_multiplier,100*s$outcome_r2,outcome_name)
  if(mode=="mixed") assumption <- paste(assumption,sprintf("The large-sample power approximation assumes a feature residual intraclass correlation of %.2f and a correlation of %.2f between repeated %s measurements after fixed-effect adjustment, equal per-visit variances, and two complete visits per participant. The coefficient represents an overall repeated-measures association, not specifically within-person change.",s$feature_rho,s$outcome_rho,outcome_name))
  paste(intro,effect,assumption,illustration)
}
