source("R/training_inputs.R")
# Treat floating-point roundoff at an exact threshold as equality.
strictly_above <- function(x, threshold) {
  x > threshold & (x - threshold) > 8 * .Machine$double.eps * pmax(abs(x), abs(threshold))
}

filter_features <- function(matrix_data, feature_type) {
  stopifnot(feature_type %in% c("SGB", "EC"),
            all(is.finite(matrix_data)), all(matrix_data >= 0))
  totals <- colSums(matrix_data)
  # Match the existing training workflow: omit empty sample profiles.
  matrix_data <- matrix_data[, totals > 0, drop = FALSE]
  totals <- totals[totals > 0]
  stopifnot(length(totals) > 0)
  # Normalize before feature filtering; do not renormalize the retained subset.
  normalized <- sweep(matrix_data, 2, totals, "/")
  positive_n <- rowSums(normalized > 0)
  prevalence <- positive_n / ncol(normalized)
  nonzero_mean <- rowSums(normalized) / positive_n
  nonzero_mean[positive_n == 0] <- NA_real_
  prevalence_pass <- prevalence >= if (feature_type == "EC") 0.80 else 0.10
  abundance_pass <- if (feature_type == "EC") rep(TRUE, nrow(normalized)) else
    !is.na(nonzero_mean) & strictly_above(nonzero_mean, 0.0001)
  retained <- prevalence_pass & abundance_pass
  list(normalized = normalized, details = data.frame(
    feature = rownames(matrix_data), positive_samples = positive_n,
    prevalence = prevalence, nonzero_mean_relative_abundance = nonzero_mean,
    prevalence_pass = prevalence_pass, abundance_pass = abundance_pass,
    retained = retained
  ))
}

summarize_filtered <- function(id, label, feature_type, matrix_data) {
  filtered <- filter_features(matrix_data, feature_type)
  details <- filtered$details
  # Compute variability only for retained features, never for excluded features.
  retained_sd <- apply(filtered$normalized[details$retained, , drop = FALSE], 1,
    function(x) sd(asin(sqrt(x[x > 0]))))
  stopifnot(length(retained_sd) > 0, all(is.finite(retained_sd)), all(retained_sd > 0))
  qs <- quantile(retained_sd, c(0.25, 0.50, 0.75))
  details$transformed_sd <- NA_real_
  details$transformed_sd[details$retained] <- retained_sd
  details$dataset <- label
  write.csv(details, paste0("data/derived/filtered_", id, "_features.csv"), row.names = FALSE)
  data.frame(id = id, dataset = label, feature_type = feature_type,
    samples = ncol(filtered$normalized), empty_samples_excluded = ncol(matrix_data) - ncol(filtered$normalized),
    input_features = nrow(matrix_data),
    observed_features = sum(details$positive_samples > 0),
    retained_features = sum(details$retained),
    excluded_observed_features = sum(!details$retained & details$positive_samples > 0),
    absent_features = sum(details$positive_samples == 0),
    excluded_total_features = sum(!details$retained),
    sd_25 = unname(qs[1]), sd_50 = unname(qs[2]), sd_75 = unname(qs[3]))
}

if (sys.nframe() == 0) {
  dir.create("data/derived", showWarnings = FALSE)
  results <- rbind(
    summarize_filtered("hmp_tax_stool", "HMP1-2 taxonomy stool", "SGB", load_hmp_taxonomy_area("Gut")),
    summarize_filtered("hmp_tax_oral", "HMP1-2 taxonomy oral", "SGB", load_hmp_taxonomy_area("Oral")),
    summarize_filtered("hmp_tax_skin", "HMP1-2 taxonomy skin", "SGB", load_hmp_taxonomy_area("Skin")),
    summarize_filtered("hmp_tax_vaginal", "HMP1-2 taxonomy vaginal", "SGB", load_hmp_taxonomy_area("Vaginal")),
    summarize_filtered("hmp_ec_stool", "HMP1-2 ECs stool", "EC", load_hmp_ec_stool())
  )
  write.csv(results, "data/derived/filtered_training.csv", row.names = FALSE)
  print(results, row.names = FALSE)
}
