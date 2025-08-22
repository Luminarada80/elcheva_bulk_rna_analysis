args <- commandArgs(TRUE)
input_dir <- args[1]
output_path_and_name <- args[2]

merge_feature_counts <- function(samples, read_count_files) {
  message("Merging ", length(read_count_files), " files")
  # Read first file to set reference order
  fc0 <- readRDS(read_count_files[1])
  if (is.null(fc0$counts)) stop("First RDS has no $counts")
  ref_ids <- rownames(fc0$counts)
  if (is.null(ref_ids)) stop("Counts object has no rownames (gene IDs)")

  out_list <- vector("list", length(read_count_files))
  names(out_list) <- samples

  for (i in seq_along(read_count_files)) {
    f <- read_count_files[i]
    s <- samples[i]
    message("Reading: ", f)
    x <- readRDS(f)
    if (is.null(x$counts)) stop("File ", f, " has no $counts")
    mat <- x$counts
    # Reorder rows to match the reference gene order
    idx <- match(ref_ids, rownames(mat))
    if (anyNA(idx)) {
      missing <- sum(is.na(idx))
      stop("File ", f, " is missing ", missing, " genes from the reference set")
    }
    vec <- mat[idx, , drop = FALSE]
    out_list[[s]] <- vec
  }

  # Column-bind in the same sample order
  merged <- do.call(cbind, out_list)
  return(merged)
}

# ---- Collect files safely ----
files <- list.files(path = input_dir,
                    pattern = "\\.[Rr]ds$",
                    full.names = TRUE)
if (length(files) == 0) stop("No .rds files found in: ", input_dir)
# Optional: sort deterministically
files <- sort(files)

# Derive sample names from basenames
samples <- tools::file_path_sans_ext(basename(files))

# Extra sanity check: all files exist
missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop("These files do not exist:\n", paste(missing, collapse = "\n"))
}

merged_feature_counts <- merge_feature_counts(samples, files)
print(merged_feature_counts[1:5, 1:min(5, ncol(merged_feature_counts))])

# Save outputs
output_rds <- paste0(output_path_and_name, ".rds")
saveRDS(merged_feature_counts, output_rds)
cat(output_rds, "\n")

merged_with_id <- data.frame(id = rownames(merged_feature_counts),
                             merged_feature_counts,
                             check.names = FALSE)
print(merged_with_id[1:5, 1:min(5, ncol(merged_with_id))])

output_txt <- paste0(output_path_and_name, ".txt")
write.table(merged_with_id, file = output_txt,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cat(output_txt, "\n")
