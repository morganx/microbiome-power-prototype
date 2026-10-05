# Source after running scripts/unpack_training.py.
metadata <- read.delim(
  "data/extracted/hmp12_data/hmp1-II_public_metadata.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("#N/A", "NA", "")
)

clean_sample_id <- function(cols) {
  ids <- sub("_taxonomic$", "", cols)
  ids <- sub("_Abundance-RPKs$", "", ids)
  ids <- sub("_Abundance$", "", ids)
  sub("_relab$", "", ids)
}

load_hmp_taxonomy_area <- function(area) {
  input <- read.delim(
    "data/extracted/hmp12_data/all_tsv.tsv",
    skip = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(input)[1] <- "feature"
  sample_ids <- clean_sample_id(names(input)[-1])
  keep_srs <- metadata$SRS[metadata$STArea == area & !is.na(metadata$SRS)]
  keep_sample <- sample_ids %in% keep_srs

  terminal_taxa <- input[grepl("t__", input$feature, fixed = TRUE), c(TRUE, keep_sample), drop = FALSE]
  matrix_data <- as.matrix(terminal_taxa[, -1, drop = FALSE])
  mode(matrix_data) <- "numeric"
  rownames(matrix_data) <- terminal_taxa$feature
  matrix_data
}

load_hmp_ec_stool <- function() {
  input <- read.delim(
    "data/extracted/hmp12_data/ecs_relab.tsv",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(input)[1] <- "feature"
  input <- input[!grepl("\\|", input$feature), , drop = FALSE]

  sample_ids <- clean_sample_id(names(input)[-1])
  keep_srs <- metadata$SRS[metadata$STSite == "Stool" & !is.na(metadata$SRS)]
  keep_sample <- sample_ids %in% keep_srs
  input <- input[, c(TRUE, keep_sample), drop = FALSE]

  matrix_data <- as.matrix(input[, -1, drop = FALSE])
  mode(matrix_data) <- "numeric"
  rownames(matrix_data) <- input$feature
  matrix_data
}

read_feature_matrix <- function(file) {
  input <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  features <- input[[1]]
  matrix_data <- as.matrix(input[, -1, drop = FALSE])
  mode(matrix_data) <- "numeric"
  rownames(matrix_data) <- features
  rowsum(matrix_data, group = rownames(matrix_data), reorder = FALSE)
}

combine_feature_matrices <- function(files) {
  matrices <- lapply(files, read_feature_matrix)
  all_features <- unique(unlist(lapply(matrices, rownames)))
  padded <- lapply(matrices, function(matrix_data) {
    output <- matrix(
      0,
      nrow = length(all_features),
      ncol = ncol(matrix_data),
      dimnames = list(all_features, colnames(matrix_data))
    )
    output[rownames(matrix_data), ] <- matrix_data
    output
  })
  do.call(cbind, padded)
}
