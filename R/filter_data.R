#' Filter merged bedMethyl data using flexible, user-defined criteria
#'
#' Computes per-site QC and modification metrics from a merged bedMethyl table
#' (output of \code{\link{merge_data}}) and filters rows using a user-provided
#' set of criteria. Percent/ratio metrics are computed only among samples with a coverage
#' greater than zero.
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
#'   \item{Mean_Coverage}{Mean coverage across all samples (zeros included).}
#'   \item{Missing}{Fraction of samples with 0 coverage (0 to 1).}
#'   \item{Percent_Methyl}{Mean methyl / mean coverage among covered samples (0 to 1+).}
#'   \item{Percent_Hydroxy}{Mean hydroxy / mean coverage among covered samples (0 to 1+).}
#'   \item{Percent_Modified}{(Mean methyl + mean hydroxy) / mean coverage among covered samples (0 to 1+).}
#'   \item{Ratio_Methyl_Hydroxy}{Mean methyl / mean hydroxy among covered samples (0 to Inf).}
#'   \item{Ratio_Methyl_Unmodified}{Mean methyl / mean unmethyl among covered samples (0 to Inf).}
#'   \item{Ratio_Hydroxy_Unmodified}{Mean hydroxy / mean unmethyl among covered samples (0 to Inf).}
#'   \item{Fraction_Hydroxy_of_Modified}{Mean hydroxy / (mean methyl + mean hydroxy) among covered samples (0 to 1).}
#' }
#'
#' Ratios use a small epsilon in denominators to avoid division by zero.
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

  # ---- Column structure checks ----
  if ((ncol(merged) - 1) %% 4 != 0) {
    stop(
      "Expected merged to have 1 Coord column plus repeating blocks of 4 columns per sample.\n",
      "But (ncol(merged) - 1) is not divisible by 4."
    )
  }

  cn <- names(merged)
  cn_data <- cn[-1]

  # Check expected order within each 4-col block
  ok_cov <- all(startsWith(cn_data[seq(1, length(cn_data), 4)], "Coverage_"))
  ok_met <- all(startsWith(cn_data[seq(2, length(cn_data), 4)], "Methyl_"))
  ok_unm <- all(startsWith(cn_data[seq(3, length(cn_data), 4)], "Unmethyl_"))
  ok_hyd <- all(startsWith(cn_data[seq(4, length(cn_data), 4)], "Hydroxymethyl_"))

  if (!(ok_cov && ok_met && ok_unm && ok_hyd)) {
    stop(
      "Column order does not match expected repeating blocks after 'Coord':\n",
      "Coverage_*, Methyl_*, Unmethyl_*, Hydroxymethyl_*\n\n",
      "Tip: ensure merge_data() returns columns in this exact order and you haven't reordered/added columns."
    )
  }

  # ---- Extract matrices (assumes Coord then repeating blocks of 4 columns) ----
  coverage <- merged[, seq(2, ncol(merged), 4), drop = FALSE]
  methyl   <- merged[, seq(3, ncol(merged), 4), drop = FALSE]
  unmethyl <- merged[, seq(4, ncol(merged), 4), drop = FALSE]
  hydroxy  <- merged[, seq(5, ncol(merged), 4), drop = FALSE]

  # Coerce to numeric defensively (fread should already be numeric)
  coverage[] <- lapply(coverage, as.numeric)
  methyl[]   <- lapply(methyl, as.numeric)
  unmethyl[] <- lapply(unmethyl, as.numeric)
  hydroxy[]  <- lapply(hydroxy, as.numeric)

  # Replace NA with 0 for coverage determination and sums
  coverage[is.na(coverage)] <- 0
  methyl[is.na(methyl)] <- 0
  unmethyl[is.na(unmethyl)] <- 0
  hydroxy[is.na(hydroxy)] <- 0

  eps <- 1e-9

  # ---- QC metrics (zeros included) ----
  mean_cov <- rowMeans(coverage)
  missing_amt <- rowMeans(coverage == 0)

  # ---- Compute modification metrics ONLY where there is coverage ----
  covered <- coverage > 0
  n_cov <- rowSums(covered)

  # Mean coverage among covered samples (used as denominator for percents)
  mean_cov_covered <- rowSums(coverage * covered) / (n_cov + eps)

  mean_m_cov <- rowSums(methyl * covered) / (n_cov + eps)
  mean_u_cov <- rowSums(unmethyl * covered) / (n_cov + eps)
  mean_h_cov <- rowSums(hydroxy * covered) / (n_cov + eps)

  denom_cov_cov <- mean_cov_covered + eps

  percent_m <- mean_m_cov / denom_cov_cov
  percent_h <- mean_h_cov / denom_cov_cov
  percent_mod <- (mean_m_cov + mean_h_cov) / denom_cov_cov

  ratio_m_h <- mean_m_cov / (mean_h_cov + eps)
  ratio_m_u <- mean_m_cov / (mean_u_cov + eps)
  ratio_h_u <- mean_h_cov / (mean_u_cov + eps)

  frac_h_of_mod <- mean_h_cov / ((mean_m_cov + mean_h_cov) + eps)

  # If a site has zero covered samples, set coverage-based metrics to NA
  zero_cov_sites <- n_cov == 0
  percent_m[zero_cov_sites] <- NA_real_
  percent_h[zero_cov_sites] <- NA_real_
  percent_mod[zero_cov_sites] <- NA_real_
  ratio_m_h[zero_cov_sites] <- NA_real_
  ratio_m_u[zero_cov_sites] <- NA_real_
  ratio_h_u[zero_cov_sites] <- NA_real_
  frac_h_of_mod[zero_cov_sites] <- NA_real_

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
  if (isTRUE(add_metrics)) out <- cbind(out, metrics)

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

  keep <- rep(TRUE, nrow(metrics))

  for (nm in names(criteria)) {
    crit <- criteria[[nm]]

    if (is.numeric(crit) && length(crit) == 2) {
      minv <- crit[1]
      maxv <- crit[2]
    } else if (is.list(crit)) {
      minv <- if (!is.null(crit$min)) crit$min else -Inf
      maxv <- if (!is.null(crit$max)) crit$max else Inf
    } else {
      stop("criteria[['", nm, "']] must be c(min,max) or list(min=, max=).")
    }

    x <- metrics[[nm]]
    keep <- keep & !is.na(x) & (x >= minv) & (x <= maxv)
  }

  res <- out[keep, , drop = FALSE]
  attr(res, "filter_criteria") <- criteria
  attr(res, "filter_metrics_available") <- names(metrics)
  res
}
