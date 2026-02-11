#' Filter merged bedMethyl data using flexible, user-defined criteria
#'
#' Computes per-site QC and modification metrics from a merged bedMethyl table
#' (output of \code{\link{merge_data}}) and filters rows using a user-provided
#' set of criteria. Criteria are named ranges over derived metrics, enabling
#' filtering for coverage/missingness and multiple comparisons such as methyl vs
#' unmodified, hydroxy vs unmodified, and methyl vs hydroxy.
#'
#' @param merged A data.frame produced by \code{\link{merge_data}} containing
#'   \code{Coord} and repeating blocks of columns per sample:
#'   Coverage_*, Methyl_*, Unmethyl_*, Hydroxymethyl_*.
#' @param criteria Optional named list of filtering rules. Each entry name must
#'   be one of the supported metric names (see Details). Each entry may be either:
#'   \itemize{
#'     \item numeric \code{c(min, max)}
#'     \item list \code{list(min = <min>, max = <max>)}
#'   }
#'   If \code{criteria = NULL} (default), no filtering is performed.
#' @param add_metrics Logical; if TRUE (default), append computed metric columns
#'   to the output.
#'
#' @return A data.frame of filtered rows (or all rows if \code{criteria = NULL}).
#'   If \code{add_metrics = TRUE}, includes computed metric columns. Also stores:
#'   \itemize{
#'     \item \code{attr(x, "filter_criteria")}
#'     \item \code{attr(x, "filter_metrics_available")}
#'   }
#'
#' @details
#' Supported metric names for \code{criteria}:
#' \describe{
#'   \item{Mean_Coverage}{Mean coverage across samples.}
#'   \item{Missing}{Fraction of samples with 0 coverage (0 to 1).}
#'   \item{Percent_Methyl}{Mean methyl / mean coverage (0 to 1).}
#'   \item{Percent_Hydroxy}{Mean hydroxy / mean coverage (0 to 1).}
#'   \item{Percent_Modified}{(Mean methyl + mean hydroxy) / mean coverage (0 to 1).}
#'   \item{Ratio_Methyl_Hydroxy}{Mean methyl / mean hydroxy (0 to Inf).}
#'   \item{Ratio_Methyl_Unmodified}{Mean methyl / mean unmethyl (0 to Inf).}
#'   \item{Ratio_Hydroxy_Unmodified}{Mean hydroxy / mean unmethyl (0 to Inf).}
#'   \item{Fraction_Hydroxy_of_Modified}{Mean hydroxy / (mean methyl + mean hydroxy) (0 to 1).}
#' }
#'
#' Ratios use a small epsilon in denominators to avoid division by zero.
#'
#' @examples
#' \dontrun{
#' # QC-only filtering
#' filtered <- filter_data(
#'   merged,
#'   criteria = list(
#'     Mean_Coverage = c(20, Inf),
#'     Missing = c(0, 0.30)
#'   )
#' )
#'
#' # Methyl vs unmodified (ratio), plus QC
#' filtered_m_u <- filter_data(
#'   merged,
#'   criteria = list(
#'     Mean_Coverage = c(20, Inf),
#'     Missing = c(0, 0.30),
#'     Ratio_Methyl_Unmodified = c(2, Inf)
#'   )
#' )
#' }
#'
#' @export
filter_data <- function(
    merged,
    criteria = NULL,
    add_metrics = TRUE
) {
  if (!("Coord" %in% names(merged))) {
    stop("merged must contain a 'Coord' column (output of merge_data()).")
  }

  coverage <- merged[, seq(2, ncol(merged), 4), drop = FALSE]
  methyl   <- merged[, seq(3, ncol(merged), 4), drop = FALSE]
  unmethyl <- merged[, seq(4, ncol(merged), 4), drop = FALSE]
  hydroxy  <- merged[, seq(5, ncol(merged), 4), drop = FALSE]

  coverage[is.na(coverage)] <- 0
  methyl[is.na(methyl)] <- 0
  unmethyl[is.na(unmethyl)] <- 0
  hydroxy[is.na(hydroxy)] <- 0

  mean_cov <- rowMeans(coverage)
  missing_amt <- rowMeans(coverage == 0)

  mean_m <- rowMeans(methyl)
  mean_u <- rowMeans(unmethyl)
  mean_h <- rowMeans(hydroxy)

  eps <- 1e-9
  denom_cov <- mean_cov + eps

  percent_m <- mean_m / denom_cov
  percent_h <- mean_h / denom_cov
  percent_mod <- (mean_m + mean_h) / denom_cov

  ratio_m_h <- mean_m / (mean_h + eps)
  ratio_m_u <- mean_m / (mean_u + eps)
  ratio_h_u <- mean_h / (mean_u + eps)

  frac_h_of_mod <- mean_h / ((mean_m + mean_h) + eps)

  metrics <- data.frame(
    Mean_Coverage = mean_cov,
    Missing = missing_amt,
    Percent_Methyl = percent_m,
    Percent_Hydroxy = percent_h,
    Percent_Modified = percent_mod,
    Ratio_Methyl_Hydroxy = ratio_m_h,
    Ratio_Methyl_Unmodified = ratio_m_u,
    Ratio_Hydroxy_Unmodified = ratio_h_u,
    Fraction_Hydroxy_of_Modified = frac_h_of_mod,
    stringsAsFactors = FALSE
  )

  out <- merged
  if (add_metrics) out <- cbind(out, metrics)

  attr(out, "filter_metrics_available") <- names(metrics)

  if (is.null(criteria)) {
    attr(out, "filter_criteria") <- NULL
    return(out)
  }

  bad <- setdiff(names(criteria), names(metrics))
  if (length(bad) > 0) {
    stop(
      "Unknown metric(s) in criteria: ",
      paste(bad, collapse = ", "),
      "\nValid metrics are:\n  ",
      paste(names(metrics), collapse = ", ")
    )
  }

  keep <- rep(TRUE, nrow(out))

  for (nm in names(criteria)) {
    crit <- criteria[[nm]]

    if (is.numeric(crit) && length(crit) == 2) {
      minv <- crit[1]
      maxv <- crit[2]
    } else if (is.list(crit)) {
      minv <- if (!is.null(crit$min)) crit$min else 0
      maxv <- if (!is.null(crit$max)) crit$max else Inf
    } else {
      stop("criteria[['", nm, "']] must be c(min,max) or list(min=, max=).")
    }

    keep <- keep & (metrics[[nm]] >= minv) & (metrics[[nm]] <= maxv)
  }

  res <- out[keep, , drop = FALSE]
  attr(res, "filter_criteria") <- criteria
  attr(res, "filter_metrics_available") <- names(metrics)
  res
}
