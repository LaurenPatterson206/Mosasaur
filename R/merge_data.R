#' Merge bedMethyl files into a single unified table
#'
#' Reads multiple bedMethyl (`*.mods.bed`) files, extracts coverage and
#' methylation-related columns, renames them with sample-specific prefixes,
#' and merges all samples by genomic coordinate.
#'
#' @param path Character string. Directory containing bedMethyl files.
#'    Defaults to the working directory (`"."`).
#' @param pattern Regular expression pattern used to identify bedMethyl files.
#'    Default matches files ending in `"\\mods.bed$"`.
#'
#' @returns A merged `data.frame` with one row per CpG site and columns for
#' coverage, methylation, unmethylation, and hydroxymethylation for each sample.
#'
#' @details
#' This function:
#' \itemize{
#' \item Identifies all files in `path` matching `pattern`
#' \item Reads each file using `data.table::fread()`
#' \item Keeps only rows where `V4 == "m"` (methylated sites)
#' \item Constructs a `Coord` variable in the format `"chr:start-end"`
#' \item Extracts and renames relevant measurement columns
#' \item Merges all samples using a full join on `Coord`
#' }
#'
#' @examples
#' \dontrun{
#' merged_df <- merge_data(path = "data/", pattern = "\\mods.bed$")
#' }
#'
#'@importFrom purrr reduce
#'
#'@export
#'
merge_data <- function(
    path = ".",
    pattern = "\\mods.bed$"
) {

  bed.files <- list.files(path = path, pattern = pattern, full.names = TRUE)

  if (length(bed.files) == 0) {
    stop ("No bedmethyl files found. Check path and/or pattern")

  }

  sample_list <- list()

  for (file in bed.files) {
    raw <- as.data.frame(data.table::fread(file))
    sample_name <- strsplit(basename(file), "[.]")[[1]][1]

    raw$Coord <- paste0(raw$V1, ":", raw$V2, "-", raw$V3)
    raw <- raw[raw$V4 == "m", ]

    colnames(raw)[10] <- paste0("Coverage_", sample_name)
    colnames(raw)[12] <- paste0("Methyl_", sample_name)
    colnames(raw)[13] <- paste0("Unmethyl_", sample_name)
    colnames(raw)[14] <- paste0("Hydroxymethyl_", sample_name)

    raw <- raw[, c("Coord",
                   paste0("Coverage_", sample_name),
                   paste0("Methyl_", sample_name),
                   paste0("Unmethyl_", sample_name),
                   paste0("Hydroxymethyl_", sample_name))]

    sample_list[[sample_name]] <- raw

  }

  merged <- purrr::reduce(sample_list, dplyr::full_join, by = "Coord")

  return(merged)

}
utils::globalVariables(c("Coord", "Mean_Coverage", "Missing", "numIDV", "pos"))
