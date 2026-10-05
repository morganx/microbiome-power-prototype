source('R/training_data.R');source('R/case_control_power.R');source('R/regression_power.R')
t<-load_applet_training();stopifnot(nrow(t)==5,all(grepl('^HMP1-2',t$dataset)),!any(t$feature_class=='Rare'))
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
 session$setInputs(dataset='Stool SGB',cases=100,controls=100,target=.8,family_alpha=.05,features=260,tests_per_feature=1,grant_quartile='75%')
 for(key in unique(e$training$key)){
  session$setInputs(dataset=key)
  stopifnot(nrow(estimates())==3,all(is.finite(estimates()$Detectable_difference)),grepl('28953883',grant_text()),nchar(output$grantSummary$html)>0)
 }
})
e<-new.env();sys.source('apps/continuous/app.R',e)
shiny::testServer(e$server,{
 session$setInputs(mode='cross',dataset='Stool SGB',n=200,covariates=2,target=.8,family_alpha=.05,features=260,analyses=1,outcome_r2='0.05',sd_multiplier=1,beta=.005,feature_rho=.5,outcome_rho=.7,outcome_name='outcome',covariate_names='age and sex')
 for(mode in c('cross','mixed'))for(key in unique(e$training$key)){
  session$setInputs(mode=mode,dataset=key)
  stopifnot(nrow(estimates())==3,all(is.finite(estimates()$Detectable_slope)),grepl('28953883',grant_text()),nchar(output$grantSummary$html)>0)
 }
})


stopifnot(identical(t$key,c("Stool EC","Stool SGB","Oral SGB","Skin SGB","Vaginal SGB")),all(t$input_features==t$class_features+t$excluded),all(grepl("0.01%",t$definition[t$key!="Stool EC"],fixed=TRUE)),grepl("≥80%",t$definition[t$key=="Stool EC"],fixed=TRUE))
source("scripts/calculate_training.R")
x<-matrix(c(.0001,.000100001,.999799999),3,10)
rownames(x)<-c("at_threshold","above_threshold","other")
f<-filter_features(x,"SGB")
stopifnot(!f$details$retained[1],f$details$retained[2])
x<-matrix(0,2,10);rownames(x)<-c("at_prevalence_threshold","other");x[2,]<-1;x[1,1]<-.01
stopifnot(filter_features(x,"SGB")$details$retained[1])

cat('PASS: both applets, five reference choices, numerical power checks, filter thresholds, and count reconciliation.\n')

x<-matrix(0,3,10);rownames(x)<-c("at_80","below_80","other");x[3,]<-1;x[1,1:8]<-1e-8;x[2,1:7]<-1e-8
f<-filter_features(x,"EC")
stopifnot(f$details$retained[1],!f$details$retained[2],all(f$details$abundance_pass))
cat("PASS: EC 80% boundary and no abundance threshold.\n")
v<-t[t$key=='Vaginal SGB',];vtext<-training_grant_intro(v)
etext<-training_grant_intro(t[t$key=='Stool EC',])
stopifnot(grepl('233 vaginal metagenomes',vtext,fixed=TRUE),grepl('543 detected features',vtext,fixed=TRUE),grepl('retaining 31 features',vtext,fixed=TRUE),grepl('MetaPhlAn 4.0.6',vtext,fixed=TRUE),grepl('among positive samples',vtext,fixed=TRUE),grepl('547 stool metagenomes',etext,fixed=TRUE),grepl('2487 detected unstratified EC features',etext,fixed=TRUE),grepl('retaining 1017 features',etext,fixed=TRUE),grepl('HUMAnN 4',etext,fixed=TRUE),grepl('at least 80%',etext,fixed=TRUE),!grepl('0.01%',etext,fixed=TRUE),!grepl('MetaPhlAn',etext,fixed=TRUE))
cat('PASS: dataset-specific grant prose, profiling methods, and filters.\n')
