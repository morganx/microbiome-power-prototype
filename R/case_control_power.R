validate_case_control <- function(cases, controls, alpha) {
  if (length(cases)!=1 || !is.finite(cases) || cases!=floor(cases) || cases<2)
    stop("Enter at least two real cases.")
  if (length(controls)!=1 || !is.finite(controls) || controls!=floor(controls) || controls<2)
    stop("Enter at least two real controls. No control group is imputed; use the continuous-outcome app for regression.")
  stopifnot(is.finite(alpha), alpha>0, alpha<1)
}
case_control_power <- function(delta, feature_sd, cases, controls, alpha=.05) {
  validate_case_control(cases,controls,alpha)
  stopifnot(all(is.finite(delta)),all(is.finite(feature_sd)),all(feature_sd>0))
  d <- abs(delta)/feature_sd
  vapply(d,function(effect) pwr::pwr.t2n.test(n1=cases,n2=controls,d=effect,
    sig.level=alpha,alternative="two.sided")$power,numeric(1))
}
case_control_detectable <- function(feature_sd,cases,controls,alpha=.05,target=.8) {
  validate_case_control(cases,controls,alpha)
  stopifnot(all(is.finite(feature_sd)),all(feature_sd>0),target>alpha,target<1)
  f <- function(d) case_control_power(d,1,cases,controls,alpha)-target
  upper <- 1
  while(f(upper)<0) {upper<-upper*2; if(upper>1e6) stop("Unable to resolve detectable effect.")}
  uniroot(f,c(0,upper),tol=1e-10)$root*feature_sd
}
