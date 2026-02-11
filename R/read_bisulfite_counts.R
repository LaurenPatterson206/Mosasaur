#' Read bisulfite methylation count files and merge into Mosasaur format
#'
#' Reads bisulfite sequencing count tables containing methylated and unmethylated
#' read counts, standardizes them, and merges across samples by genomic coordinate.
#' This function is intended to provide an entry point for bisulfite-based data
#' analogous to \code{\link{merge_data}} for nanopore bedMethyl.
#'
#' Supported inputs:
#' \itemize{
#'   \item \strong{Per-sample files (recommended)} with columns like:
#'         \code{chr, start, end, meth, unmeth} (column names configurable)
#'   \item \strong{A single wide file} containing multiple samples as columns
#'         (use \code{wide = TRUE})
#' }
#'
#' @param path Directory containing per-sample files. Ignored if \code{files} provided.
#' @param pattern Regex used to match per-sample files (e.g., "\\\\.txt$").
#' @param files Optional character vector of file paths to read.
#' @param sample_names Optional sample names (same length as \code{files}).
#'   If NULL, inferred from filename prefix before first dot.
#' @param chr_col,start_col,end_col Column names for genomic coordinates
#'   in per-sample files.
#' @param meth_col,unmeth_col Column names for methylated/unmethylated counts
#'   in per-sample files.
#' @param coord_style Either \code{"chr:start-end"} or \code{"chr:start"}.
#' @param wide Logical; if TRUE, reads a single wide table with sample columns
#'   already present (see Details).
#' @param wide_file Path to a single wide-format table (required if \code{wide=TRUE}).
#' @param return_coord_cols Logical; if TRUE, also return \code{chr,start,end}.
#'
#' @details
#' The returned merged table follows the same repeating 4-column block layout used
#' by Mosasaur's nanopore merge step:
#' \preformatted{
#' Coord | Coverage_sample | Methyl_sample | Unmethyl_sample | Hydroxymethyl_sample
#' }
#' For bisulfite data, Hydroxymethyl columns are filled with zeros, and Coverage is
#' computed as \code{meth + unmeth}.
#'
#' @return A merged \code{data.frame} with one row per CpG coordinate and
#'   sample-specific columns for coverage/methyl/unmethyl/hydroxy (hydroxy=0).
#'
#' @examples
#' \dontrun{
#' merged_bs <- read_bisulfite_counts(
#'   path = "bisulfite_counts/",
#'   pattern = "\\\\.counts\\.txt$",
#'   chr_col = "chr",
#'   start_col = "pos",
#'   end_col = "pos_end",
#'   meth_col = "mC",
#'   unmeth_col = "uC"
#' )
#' }
#'
#' @export
read_bisulfite_counts <- function(
    path = ".",
    pattern = "\\.txt$",
    files = NULL,
    sample_names = NULL,
    chr_col = "chr",
    start_col = "start",
    end_col = "end",
    meth_col = "meth",
    unmeth_col = "unmeth",
    coord_style = c("chr:start-end", "chr:start"),
    wide = FALSE,
    wide_file = NULL,
    return_coord_cols = FALSE
) {
  coord_style <- match.arg(coord_style)

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }

  #wide format mode
  if (isTRUE(wide)) {
    if (is.null(wide_file)) stop("If wide=TRUE, you must provide wide_file.")
    df <- as.data.frame(data.table::fread(wide_file))

    #require coords
    needed <- c(chr_col, start_col)
    if (!all(needed %in% names(df))) {
      stop("Wide file must include at least: ", paste(needed, collapse = ", "))
    }
    if (!end_col %in% names(df)) {
      df[[end_col]] <- df[[start_col]] + 1
    }

    if (coord_style == "chr:start-end") {
      df$Coord <- paste0(df[[chr_col]], ":", df[[start_col]], "-", df[[end_col]])
    } else {
      df$Coord <- paste0(df[[chr_col]], ":", df[[start_col]])
    }

    #expect user already has Methyl_ and Unmethyl_ columns in wide file.
    #if not, they can rename before calling
    if (!any(grepl("^Methyl_", names(df))) || !any(grepl("^Unmethyl_", names(df)))) {
      stop("Wide file must contain columns starting with 'Methyl_' and 'Unmethyl_'.")
    }

    #build coverage and hydroxy=0 for each sample
    methyl_cols <- grep("^Methyl_", names(df), value = TRUE)
    unmeth_cols <- grep("^Unmethyl_", names(df), value = TRUE)

    #infer matching samples
    samps_m <- sub("^Methyl_", "", methyl_cols)
    samps_u <- sub("^Unmethyl_", "", unmeth_cols)
    if (!setequal(samps_m, samps_u)) stop("Methyl_ and Unmethyl_ sample sets do not match.")

    for (s in samps_m) {
      df[[paste0("Coverage_", s)]] <- df[[paste0("Methyl_", s)]] + df[[paste0("Unmethyl_", s)]]
      df[[paste0("Hydroxymethyl_", s)]] <- 0
    }

    keep <- c("Coord", unlist(lapply(samps_m, function(s) {
      c(paste0("Coverage_", s), paste0("Methyl_", s), paste0("Unmethyl_", s), paste0("Hydroxymethyl_", s))
    })))

    out <- df[, keep, drop = FALSE]
    if (return_coord_cols) out <- cbind(df[, c(chr_col, start_col, end_col), drop = FALSE], out)
    return(out)
  }

  #per-sample files mode
  if (is.null(files)) {
    files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  }
  if (length(files) == 0) stop("No bisulfite count files found. Check path/pattern/files.")

  if (is.null(sample_names)) {
    sample_names <- vapply(files, function(f) strsplit(basename(f), "[.]")[[1]][1], character(1))
  }
  if (length(sample_names) != length(files)) {
    stop("sample_names must be same length as files.")
  }

  sample_list <- vector("list", length(files))
  names(sample_list) <- sample_names

  for (i in seq_along(files)) {
    f <- files[i]
    s <- sample_names[i]

    raw <- as.data.frame(data.table::fread(f))
    needed <- c(chr_col, start_col, meth_col, unmeth_col)
    if (!all(needed %in% names(raw))) {
      stop("File ", f, " is missing required columns: ",
           paste(setdiff(needed, names(raw)), collapse = ", "))
    }

    if (!end_col %in% names(raw)) {
      raw[[end_col]] <- raw[[start_col]] + 1
    }

    if (coord_style == "chr:start-end") {
      raw$Coord <- paste0(raw[[chr_col]], ":", raw[[start_col]], "-", raw[[end_col]])
    } else {
      raw$Coord <- paste0(raw[[chr_col]], ":", raw[[start_col]])
    }

    meth <- raw[[meth_col]]
    unm  <- raw[[unmeth_col]]

    df <- data.frame(
      Coord = raw$Coord,
      stringsAsFactors = FALSE
    )
    df[[paste0("Coverage_", s)]] <- meth + unm
    df[[paste0("Methyl_", s)]] <- meth
    df[[paste0("Unmethyl_", s)]] <- unm
    df[[paste0("Hydroxymethyl_", s)]] <- 0

    sample_list[[s]] <- df
  }

  merged <- purrr::reduce(sample_list, dplyr::full_join, by = "Coord")
  merged
}
