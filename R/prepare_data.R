#' Prepare count matrices for PQLseq2 from merged bedMethyl data
#'
#' Builds success (MCounts) and trial (TCounts) matrices for PQLseq2.
#' By default, returns raw counts without pseudocounts (recommended for modeling).
#'
#' @param filtered A data.frame containing a \code{Coord} column ("chr:start-end")
#'   and count columns with prefixes: \code{Methyl_}, \code{Hydroxymethyl_}, \code{Unmethyl_}.
#' @param compare One of \code{"modified_vs_unmodified"}, \code{"methyl_vs_unmodified"},
#'   \code{"hydroxy_vs_unmodified"}.
#' @param add_pseudocounts Logical; if TRUE adds pseudocounts (not recommended for PQLseq2).
#' @param success_pseudocount Added to successes if \code{add_pseudocounts=TRUE}.
#' @param trial_pseudocount Added to trials if \code{add_pseudocounts=TRUE}.
#' @param strip_sample_prefix_regex Optional regex applied to sample names (e.g. \code{"^SS"}).
#' @param write_outputs Logical; if TRUE writes \code{<prefix>_MCounts.txt} and \code{<prefix>_TCounts.txt}.
#' @param prefix Output prefix.
#'
#' @return A list with \code{MCounts}, \code{TCounts}, and metadata.
#' @export
prepare_data <- function(
    filtered,
    compare = c("modified_vs_unmodified", "methyl_vs_unmodified", "hydroxy_vs_unmodified"),
    add_pseudocounts = FALSE,
    success_pseudocount = 1,
    trial_pseudocount = 0,
    strip_sample_prefix_regex = NULL,
    write_outputs = TRUE,
    prefix = "PQLseq2"
) {
  compare <- match.arg(compare)

  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required for prepare_data().")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required for prepare_data().")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required for prepare_data().")

  if (!("Coord" %in% names(filtered))) {
    stop("filtered must contain a 'Coord' column formatted as 'chr:start-end'.")
  }

  # Pull component matrices
  meth_df   <- dplyr::select(filtered, dplyr::starts_with("Methyl_"))
  hydro_df  <- dplyr::select(filtered, dplyr::starts_with("Hydroxymethyl_"))
  unmeth_df <- dplyr::select(filtered, dplyr::starts_with("Unmethyl_"))

  if (ncol(meth_df) == 0) stop("No columns with prefix 'Methyl_'.")
  if (ncol(unmeth_df) == 0) stop("No columns with prefix 'Unmethyl_'.")
  if (compare %in% c("modified_vs_unmodified", "hydroxy_vs_unmodified") && ncol(hydro_df) == 0) {
    stop("compare='", compare, "' requires 'Hydroxymethyl_' columns, but none were found.")
  }

  strip_component <- function(nm) sub("^(Methyl_|Hydroxymethyl_|Unmethyl_)", "", nm)

  samp_m <- strip_component(colnames(meth_df))
  samp_u <- strip_component(colnames(unmeth_df))
  if (!setequal(samp_m, samp_u)) stop("Sample IDs in Methyl_ and Unmethyl_ columns do not match.")

  if (compare %in% c("modified_vs_unmodified", "hydroxy_vs_unmodified")) {
    samp_h <- strip_component(colnames(hydro_df))
    if (!setequal(samp_m, samp_h)) stop("Sample IDs in Hydroxymethyl_ do not match Methyl_.")
  }

  sample_ids <- sort(unique(samp_m))

  reorder_component <- function(df, prefix) {
    df[, paste0(prefix, sample_ids), drop = FALSE]
  }

  meth_df   <- reorder_component(meth_df, "Methyl_")
  unmeth_df <- reorder_component(unmeth_df, "Unmethyl_")
  if (compare %in% c("modified_vs_unmodified", "hydroxy_vs_unmodified")) {
    hydro_df <- reorder_component(hydro_df, "Hydroxymethyl_")
  }

  to_num <- function(df) {
    df[is.na(df)] <- 0
    m <- as.matrix(df)
    storage.mode(m) <- "numeric"
    m
  }

  M_m <- to_num(meth_df)
  M_u <- to_num(unmeth_df)
  M_h <- if (compare %in% c("modified_vs_unmodified", "hydroxy_vs_unmodified")) to_num(hydro_df) else NULL

  if (compare == "modified_vs_unmodified") {
    Y <- M_m + M_h
    T <- M_m + M_h + M_u
  } else if (compare == "methyl_vs_unmodified") {
    Y <- M_m
    T <- M_m + M_u
  } else {
    Y <- M_h
    T <- M_h + M_u
  }

  if (add_pseudocounts) {
    Y <- Y + success_pseudocount
    T <- T + trial_pseudocount
  }

  if (any(T < Y, na.rm = TRUE)) stop("Found TCounts < MCounts. Check inputs.")

  # Coord parsing
  coord_df <- data.frame(Coord = filtered$Coord, stringsAsFactors = FALSE)
  coord_df <- tidyr::separate(coord_df, Coord, into = c("chr", "pos"), sep = ":", remove = TRUE)
  coord_df <- tidyr::separate(coord_df, pos, into = c("start", "end"), sep = "-", remove = TRUE)
  coord_df$start <- suppressWarnings(as.integer(coord_df$start))

  MCounts <- data.frame(chr = coord_df$chr, start = coord_df$start, Y, check.names = FALSE)
  TCounts <- data.frame(chr = coord_df$chr, start = coord_df$start, T, check.names = FALSE)

  colnames(MCounts)[-(1:2)] <- sample_ids
  colnames(TCounts)[-(1:2)] <- sample_ids

  if (!is.null(strip_sample_prefix_regex)) {
    new_ids <- sub(strip_sample_prefix_regex, "", sample_ids)
    colnames(MCounts)[-(1:2)] <- new_ids
    colnames(TCounts)[-(1:2)] <- new_ids
    sample_ids <- new_ids
  }

  mcounts_file <- paste0(prefix, "_MCounts.txt")
  tcounts_file <- paste0(prefix, "_TCounts.txt")

  if (isTRUE(write_outputs)) {
    data.table::fwrite(MCounts, mcounts_file, sep = "\t")
    data.table::fwrite(TCounts, tcounts_file, sep = "\t")
  }

  list(
    MCounts = MCounts,
    TCounts = TCounts,
    compare = compare,
    add_pseudocounts = add_pseudocounts,
    success_pseudocount = success_pseudocount,
    trial_pseudocount = trial_pseudocount,
    sample_ids = sample_ids,
    mcounts_file = mcounts_file,
    tcounts_file = tcounts_file
  )
}
