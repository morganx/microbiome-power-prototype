source("R/training_inputs.R")
common_prevalence_threshold <- 0.20
normalize_samples <- function(matrix_data) {
  matrix_data <- matrix_data[rowSums(matrix_data > 0, na.rm = TRUE) > 0, , drop = FALSE]
  sample_sums <- colSums(matrix_data, na.rm = TRUE)
  keep_samples <- sample_sums > 0 & !is.na(sample_sums)
  sweep(matrix_data[, keep_samples, drop = FALSE], 2, sample_sums[keep_samples], "/")
}

sd_ignore_sparse <- function(values) {
  values <- values[values > 0 & is.finite(values)]
  if (length(values) <= 1) {
    return(NA_real_)
  }
  sd(asin(sqrt(values)))
}

summarize_matrix <- function(dataset, matrix_data) {
  normalized <- normalize_samples(matrix_data)
  prevalence <- rowSums(normalized > 0, na.rm = TRUE) / ncol(normalized)
  feature_class <- ifelse(prevalence >= common_prevalence_threshold, "Common", "Rare")
  feature_sd <- apply(normalized, 1, sd_ignore_sparse)

  do.call(rbind, lapply(c("Common", "Rare"), function(class_name) {
    qs <- quantile(feature_sd[feature_class == class_name], c(0.25, 0.5, 0.75), na.rm = TRUE)
    data.frame(
      dataset = dataset,
      feature_class = class_name,
      samples = ncol(normalized),
      total_features = nrow(normalized),
      class_features = sum(feature_class == class_name),
      sd_25 = as.numeric(qs[1]),
      sd_50 = as.numeric(qs[2]),
      sd_75 = as.numeric(qs[3])
    )
  }))
}

results <- rbind(
  summarize_matrix("HMP1-2 taxonomy stool", load_hmp_taxonomy_area("Gut")),
  summarize_matrix("HMP1-2 taxonomy oral", load_hmp_taxonomy_area("Oral")),
  summarize_matrix("HMP1-2 taxonomy skin", load_hmp_taxonomy_area("Skin")),
  summarize_matrix("HMP1-2 taxonomy vaginal", load_hmp_taxonomy_area("Vaginal")),
  summarize_matrix("HMP1-2 ECs stool", load_hmp_ec_stool())
)

if (!dir.exists("data/derived")) {
  dir.create("data/derived")
}
write.csv(results, "data/derived/app_training_data_sd_quartiles.csv", row.names = FALSE)
print(results)
