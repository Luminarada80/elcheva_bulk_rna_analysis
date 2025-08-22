# merge_read_counts.R
args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
input_dir <- args[1]
output_path_and_name <- args[2]

suppressWarnings(suppressMessages({
  library(tools)
}))

# ---- Helpers ----
read_featurecounts_txt <- function(f) {
  # featureCounts writes commented header lines starting with '#'
  # and a table with columns: Geneid Chr Start End Strand Length <count columns...>
  con <- file(f, open = "r")
  on.exit(close(con))
  first_lines <- readLines(con, n = 100L, warn = FALSE)
  # rewind
  close(con); con <- file(f, open = "r"); on.exit(close(con))
  skip <- which.max(grepl("^Geneid\\t", first_lines)) - 1L
  if (skip < 0) skip <- 0L
  df <- tryCatch(
    read.delim(f, header = TRUE, sep = "\t", comment.char = "", skip = skip, check.names = FALSE),
    error = function(e) stop("Failed to read featureCounts file: ", f, "\n", conditionMessage(e))
  )
  if (!"Geneid" %in% names(df)) stop("No 'Geneid' column in: ", f)
  # Take the LAST column as the counts (works for per-sample outputs)
  cnt_col <- ncol(df)
  data.frame(Geneid = df$Geneid, counts = df[[cnt_col]], stringsAsFactors = FALSE)
}

merge_featurecounts_txt_dir <- function(files) {
  # derive sample names from filenames (strip .txt)
  samples <- file_path_sans_ext(basename(files))
  # read first to get gene order
  x0 <- read_featurecounts_txt(files[1])
  ref_ids <- x0$Geneid
  mat <- matrix(NA_integer_, nrow = length(ref_ids), ncol = length(files),
                dimnames = list(ref_ids, samples))
  mat[, 1] <- x0$counts
  for (i in seq_along(files)[-1]) {
    xi <- read_featurecounts_txt(files[i])
    m <- match(ref_ids, xi$Geneid)
    if (anyNA(m)) {
      stop("File ", files[i], " missing ", sum(is.na(m)), " genes found in the first file.")
    }
    mat[, i] <- xi$counts[m]
  }
  mat
}

merge_rds_dir <- function(files) {
  # Your original RDS merge (expects $counts and matching rownames)
  samples <- file_path_sans_ext(basename(files))
  fc0 <- readRDS(files[1]); stopifnot(!is.null(fc0$counts))
  ref_ids <- rownames(fc0$counts); stopifnot(!is.null(ref_ids))
  out_list <- vector("list", length(files)); names(out_list) <- samples
  out_list[[1]] <- fc0$counts[ref_ids, , drop = FALSE]
  for (i in seq_along(files)[-1]) {
    x <- readRDS(files[i]); stopifnot(!is.null(x$counts))
    mat <- x$counts
    idx <- match(ref_ids, rownames(mat))
    if (anyNA(idx)) stop("File ", files[i], " missing ", sum(is.na(idx)), " genes from reference")
    out_list[[i]] <- mat[idx, , drop = FALSE]
  }
  do.call(cbind, out_list)
}

# ---- Detect inputs ----
rds_files <- sort(list.files(input_dir, pattern = "\\.[Rr]ds$", full.names = TRUE))
txt_files <- sort(list.files(input_dir, pattern = "\\.txt$",   full.names = TRUE))

if (length(rds_files) > 0) {
  message("Merging RDS files (n=", length(rds_files), ") in: ", input_dir)
  merged <- merge_rds_dir(rds_files)
} else if (length(txt_files) > 0) {
  # filter to typical featureCounts outputs
  fc_txt <- txt_files[grepl("featureCounts", basename(txt_files))]
  if (length(fc_txt) == 0) fc_txt <- txt_files
  message("Merging featureCounts TXT files (n=", length(fc_txt), ") in: ", input_dir)
  merged <- merge_featurecounts_txt_dir(fc_txt)
} else {
  stop("No .rds or .txt files found in: ", input_dir)
}

# ---- Save outputs ----
dir.create(dirname(output_path_and_name), showWarnings = FALSE, recursive = TRUE)

saveRDS(merged, paste0(output_path_and_name, ".rds"))
write.table(
  data.frame(id = rownames(merged), merged, check.names = FALSE),
  file = paste0(output_path_and_name, ".txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Wrote: ", paste0(output_path_and_name, ".rds"))
message("Wrote: ", paste0(output_path_and_name, ".txt"))
