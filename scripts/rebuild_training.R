for(script in c("scripts/calculate_common_training.R","scripts/calculate_realistic_training.R")) {
 status<-system2(file.path(R.home("bin"),"Rscript"),script)
 if(status!=0)stop("Rebuild failed: ",script)
}
