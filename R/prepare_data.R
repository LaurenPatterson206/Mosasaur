#' Prepare count matrices for PQLseq2 based on a chosen comparison
#'
#' Constructs successes (MCounts) and trials (TCounts) matrices for binomial mixed
#' modeling using PQLseq2. The constructed matrices depend on the selected
#' comparison mode.
#'
#' Supported comparison modes:
#' \itemize{
#'   \item \code{"modified_vs_unmodified"}: (methyl + hydroxy) vs unmethyl
#'   \item \code{"methyl_vs_unmodified"}: methyl vs unmethyl
#'   \item \code{"hydroxy_vs_unmodified"}: hydroxy vs unmethyl
#' }
#'
#' @param filtered A data.frame typically output from \code{\link{filter_data}}.
#'   Must contain \code{Coord} formatted as "chr:start-end" and count columns
#'   beginning with Methyl_, Hydroxymethyl_, and Unmethyl_.
#' @param compare Which comparison mode to prepare.
#' @param write_outputs Logical; if TRUE, writes `<prefix>_MCounts.txt` and
#'   `<prefix>_TCounts.txt`.
#' @param prefix Output filename prefix.
#' @param strip_sample_prefix_regex Optional regex to strip from sample column
#'   names (e.g. "^SSS").
#' @param warn_on_mismatch Logical; if TRUE, warns when upstream filter criteria
#'   suggest a different comparison than \code{compare}.
#'
#' @return A list containing MCounts, TCounts, compare mode, and output filenames.
#'
#' @examples
#' \dontrun{
#' prep <- prepare_data(filtered, compare = "methyl_vs_unmodified",
#'                      prefix = "EtOH6m_methyl_vs_unmod",
#'                      strip_sample_prefix_regex = "^SSS")
#' }
#'
#' @export
prepare_data <- function(
    filtered,
    compare = c("modified_vs_unmodified", "methyl_vs_unmodified", "hydroxy_vs_unmodified"),
    write_outputs = TRUE,
    prefix = "PQLseq2",
    strip_sample_prefix_regex = NULL,
    warn_on_mismatch = TRUE
) {
  compare <- match.arg(compare)

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for prepare_data().")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for prepare_data().")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for prepare_data().")
  }

  if (!("Coord" %in% names(filtered))) {
    stop("filtered must contain a 'Coord' column in the format 'chr:start-end'.")
  }

  if (warn_on_mismatch) {
    crit <- attr(filtered, "filter_criteria")
    if (!is.null(crit) && length(crit) > 0) {
      cn <- names(crit)
      hydroxy_related <- any(grepl("Hydroxy", cn, fixed = TRUE))
      methyl_related <- any(grepl("Methyl", cn, fixed = TRUE))

      if (compare == "methyl_vs_unmodified" && hydroxy_related && !methyl_related) {
        warning("Filter criteria include hydroxy-related metrics but compare is methyl_vs_unmodified.")
      }
      if (compare == "hydroxy_vs_unmodified" && methyl_related && !hydroxy_related) {
        warning("Filter criteria include methyl-related metrics but compare is hydroxy_vs_unmodified.")
      }
    }
  }

  meth_mat   <- dplyr::select(filtered, dplyr::starts_with("Methyl_"))
  hydro_mat  <- dplyr::select(filtered, dplyr::starts_with("Hydroxymethyl_"))
  unmeth_mat <- dplyr::select(filtered, dplyr::starts_with("Unmethyl_"))

  if (compare == "modified_vs_unmodified") {
    mod <- meth_mat + hydro_mat + 1
    cov <- meth_mat + hydro_mat + unmeth_mat + 2
  } else if (compare == "methyl_vs_unmodified") {
    mod <- meth_mat + 1
    cov <- meth_mat + unmeth_mat + 2
  } else {
    mod <- hydro_mat + 1
    cov <- hydro_mat + unmeth_mat + 2
  }

  rownames(mod) <- filtered$Coord
  rownames(cov) <- filtered$Coord
  mod$Coord <- rownames(mod)
  cov$Coord <- rownames(cov)

  mod <- tidyr::separate(mod, Coord, into = c("chr", "pos"), sep = ":")
  mod <- tidyr::separate(mod, pos, into = c("drop", "start"), sep = "-")
  mod <- mod[, c("chr", "start", setdiff(names(mod), c("chr", "start", "drop")))]

  cov <- tidyr::separate(cov, Coord, into = c("chr", "pos"), sep = ":")
  cov <- tidyr::separate(cov, pos, into = c("drop", "start"), sep = "-")
  cov <- cov[, c("chr", "start", setdiff(names(cov), c("chr", "start", "drop")))]

  colnames(mod) <- gsub("^(Methyl_|Hydroxymethyl_|Unmethyl_)", "", colnames(mod))
  colnames(cov) <- gsub("^(Methyl_|Hydroxymethyl_|Unmethyl_)", "", colnames(cov))

  if (!is.null(strip_sample_prefix_regex)) {
    colnames(mod) <- sub(strip_sample_prefix_regex, "", colnames(mod))
    colnames(cov) <- sub(strip_sample_prefix_regex, "", colnames(cov))
  }

  mcounts_file <- paste0(prefix, "_MCounts.txt")
  tcounts_file <- paste0(prefix, "_TCounts.txt")

  if (write_outputs) {
    data.table::fwrite(mod, mcounts_file, sep = "\t")
    data.table::fwrite(cov, tcounts_file, sep = "\t")
  }

  list(
    MCounts = mod,
    TCounts = cov,
    compare = compare,
    mcounts_file = mcounts_file,
    tcounts_file = tcounts_file
  )
}

