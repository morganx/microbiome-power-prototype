source('R/training_data.R');source('R/case_control_power.R');source('R/regression_power.R')
t<-load_applet_training();stopifnot(nrow(t)==10,all(grepl('^HMP1-2',t$dataset)),!any(t$feature_class=='Rare'))
for(n0 in c(0,1))stopifnot(inherits(try(case_control_detectable(.03,100,n0),silent=TRUE),'try-error'))
sds<-c(.01,.02,.03);d<-case_control_detectable(sds,100,120,.01,.8)
stopifnot(max(abs(case_control_power(d,sds,100,120,.01)-.8))<1e-7)
b<-regression_detectable(sds,200,5,.01,.8,.1)
stopifnot(max(abs(regression_power(b,sds,200,5,.01,.1)-.8))<1e-7)
crit<-qt(.995,193);nc<-b[1]*sqrt(199*.9)/sds[1]
stopifnot(abs(pt(-crit,193,ncp=nc)+pt(crit,193,ncp=nc,lower.tail=FALSE)-.8)<1e-7)
b<-mixed_detectable(sds,200,5,.01,.8,.1,1,.5,.7)
stopifnot(max(abs(mixed_power(b,sds,200,5,.01,.1,1,.5,.7)-.8))<1e-7)
e<-new.env();sys.source('apps/case-control/app.R',e)
shiny::testServer(e$server,{
 session$setInputs(dataset='Original: HMP1-2 taxonomy stool',cases=100,controls=100,target=.8,family_alpha=.05,features=260,tests_per_feature=1,grant_quartile='75%')
 for(key in unique(e$training$key)){
  session$setInputs(dataset=key)
  stopifnot(nrow(estimates())==3,all(is.finite(estimates()$Detectable_difference)),grepl('28953883',grant_text()),nchar(output$grantSummary$html)>0)
 }
})
e<-new.env();sys.source('apps/continuous/app.R',e)
shiny::testServer(e$server,{
 session$setInputs(mode='cross',dataset='Original: HMP1-2 taxonomy stool',n=200,covariates=2,target=.8,family_alpha=.05,features=260,analyses=1,outcome_r2='0.05',sd_multiplier=1,beta=.005,feature_rho=.5,outcome_rho=.7,outcome_name='outcome',covariate_names='age and sex')
 for(mode in c('cross','mixed'))for(key in unique(e$training$key)){
  session$setInputs(mode=mode,dataset=key)
  stopifnot(nrow(estimates())==3,all(is.finite(estimates()$Detectable_slope)),grepl('28953883',grant_text()),nchar(output$grantSummary$html)>0)
 }
})
cat('PASS: HMP-only training, numerical power inversions, independent t identity, both applets and every training option.\n')
