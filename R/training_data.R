load_applet_training <- function() {
  t <- read.csv("data/derived/filtered_training.csv", stringsAsFactors=FALSE)
  labels <- c(hmp_ec_stool="Stool EC", hmp_tax_stool="Stool SGB", hmp_tax_oral="Oral SGB", hmp_tax_skin="Skin SGB", hmp_tax_vaginal="Vaginal SGB")
  t <- t[match(names(labels),t$id),]
  data.frame(dataset=t$dataset,feature_class="Retained",samples=t$samples,
    total_features=t$retained_features,class_features=t$retained_features,
    input_features=t$observed_features,sd_25=t$sd_25,sd_50=t$sd_50,sd_75=t$sd_75,
    key=unname(labels),excluded=t$excluded_observed_features,
    definition=sprintf("Of %d detected features, %d were retained using %s. These filters reduce the feature set for linear modeling. Excluded features are omitted from the SD distribution and default testing family.", t$observed_features,t$retained_features, ifelse(t$feature_type=="EC", "detection in ≥80% of samples, with no abundance threshold", "prevalence ≥10% and mean relative abundance >0.01% among positive samples")))
}
