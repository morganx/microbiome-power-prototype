# Interpret a transformed mean difference at a fixed 1% reference abundance.
case_control_abundance_example <- function(delta) {
  stopifnot(all(is.finite(delta)),all(delta>=0))
  angle<-asin(sqrt(.01))+delta
  ifelse(angle<=pi/2,100*sin(angle)^2,NA_real_)
}
make_case_control_grant_text <- function(t,s,target,features,tests_per_feature,family_alpha,quartile="75%") {
  stopifnot(nrow(t)==1,quartile %in% c("25%","50%","75%"))
  index<-match(quartile,c("25%","50%","75%"))
  qs<-as.numeric(t[1,c("sd_25","sd_50","sd_75")]);sd<-qs[index]
  delta<-do.call(case_control_detectable,c(list(feature_sd=sd,target=target),s))
  example<-case_control_abundance_example(delta)
  type<-if(grepl("ECs",t$dataset,fixed=TRUE))"enzyme-function (EC)" else "taxonomic"
  source<-sprintf("%d %s profiles from the expanded Human Microbiome Project (Lloyd-Price et al.; PMID: 28953883)",t$samples,t$dataset)
  intro<-sprintf("Using variability estimated from %s, we evaluated %s features. %s The 25th, 50th and 75th percentiles of retained feature SDs were %s on the arcsine-square-root relative-abundance scale, calculated among positive samples.",source,type,t$definition,paste(sprintf("%.5g",qs),collapse=", "))
  design<-sprintf("For a cross-sectional comparison of %d independent cases and %d independent controls, a two-sided two-sample t test with Bonferroni correction across %d features and %d test%s per feature (family-wise alpha %.3g; per-test alpha %.6g)",s$cases,s$controls,features,tests_per_feature,if(tests_per_feature==1)"" else "s",family_alpha,s$alpha)
  percentile<-c("25th","50th (median)","75th")[index]
  effect<-sprintf("%s has an estimated %.0f%% power to detect a difference between group means of %.5f on the transformed scale, assuming a common within-group SD of %.5f corresponding to the %s percentile of training feature SDs.",design,100*target,delta,sd,percentile)
  illustration<-if(is.na(example))"A positive shift of this size from a 1% reference exceeds the valid transformed-abundance range, so no relative-abundance example is reported." else sprintf("For an individual feature with a reference fitted relative abundance of 1%% in controls, this difference corresponds to approximately %.3f%% in cases. This is an illustrative back-transformed fitted value, not a difference between arithmetic means on the raw relative-abundance scale.",example)
  caveat<-"Power is per feature, not the probability of detecting all features. Positive-only training SDs are planning proxies; these estimates do not model carriage changes, covariate adjustment or repeated measurements."
  paste(intro,effect,illustration,caveat)
}
