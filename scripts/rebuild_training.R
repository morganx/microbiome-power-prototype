status <- system2(file.path(R.home("bin"), "Rscript"), "scripts/calculate_training.R")
if(status != 0) stop("Training rebuild failed")
