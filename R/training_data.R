load_regression_training <- function() {
  original <- read.csv("data/derived/app_training_data_sd_quartiles.csv", stringsAsFactors = FALSE)
  original$key <- paste0("Original: ", original$dataset)
  original$definition <- "Original classes: common >=20% prevalence; rare <20%. Both classes remain in the tested feature family."
  original$excluded <- 0L
  realistic <- read.csv("data/derived/power.analysis.realistic_training.csv", stringsAsFactors = FALSE)
  filtered <- data.frame(dataset = realistic$dataset, feature_class = "Retained",
    samples = realistic$samples, total_features = realistic$retained_features,
    class_features = realistic$retained_features, sd_25 = realistic$sd_25,
    sd_50 = realistic$sd_50, sd_75 = realistic$sd_75,
    key = paste0("Realistic: ", realistic$dataset),
    definition = ifelse(realistic$feature_type == "EC", "ECs: >80% prevalence.",
      "SGBs: >10% prevalence AND >0.1% mean relative abundance among positive samples."),
    excluded = realistic$excluded_total_features)
  rbind(original, filtered)
}

load_applet_training <- function() {
  t <- load_regression_training()
  t <- t[t$feature_class != "Rare", , drop = FALSE]
  original <- t$feature_class == "Common"
  t$excluded[original] <- t$total_features[original] - t$class_features[original]
  t$total_features <- t$class_features
  t$definition[original] <- "Common features: present in at least 20% of training samples. Rare features are excluded from analysis and correction."
  t
}
