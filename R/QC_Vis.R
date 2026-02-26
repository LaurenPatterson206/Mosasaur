#' Visualize QC metrics for merged modification count tables
#'
#' Generates histograms and/or density plots for key QC summaries derived from a
#' merged Mosasaur table (e.g., output of \code{\link{merge_data}} or
#' \code{\link{read_bisulfite_counts}}). This visualization is intended to guide
#' selection of filtering thresholds (e.g., minimum coverage and maximum missingness).
#'
#' @param merged A data.frame with a \code{Coord} column followed by repeating
#'   blocks of sample columns in the order:
#'   \code{Coverage_}, \code{Methyl_}, \code{Unmethyl_}, \code{Hydroxymethyl_}.
#' @param plot_type Character. One of \code{"hist"}, \code{"density"}, or \code{"both"}.
#' @param bins Integer. Number of bins for histograms (ignored for density-only).
#' @param save Logical; if TRUE, save all plots to a single PDF.
#' @param prefix Character; filename prefix when \code{save=TRUE}.
#'
#' @return A data.frame with per-site QC summaries:
#'   \itemize{
#'     \item \code{Mean_Coverage}
#'     \item \code{Missing} (proportion of samples with zero coverage)
#'     \item \code{Percent_Methyl}
#'     \item \code{Percent_Hydroxy}
#'   }
#'
#' @examples
#' \dontrun{
#' merged <- merge_data(path = "bedmethyl/")
#' qc <- QC_Vis(merged, plot_type = "both", save = TRUE, prefix = "MyStudy")
#' head(qc)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density theme_classic labs
#' @importFrom grDevices pdf dev.off
#' @importFrom rlang .data
#'
#' @export
QC_Vis <- function(
    merged,
    plot_type = c("hist", "density", "both"),
    bins = 50,
    save = FALSE,
    prefix = "QC"
) {
  plot_type <- match.arg(plot_type)

  # Extract matrices (assumes Coord then repeating blocks of 4 columns)
  coverage <- merged[, seq(2, ncol(merged), 4), drop = FALSE]
  methyl   <- merged[, seq(3, ncol(merged), 4), drop = FALSE]
  unmethyl <- merged[, seq(4, ncol(merged), 4), drop = FALSE]
  hydroxy  <- merged[, seq(5, ncol(merged), 4), drop = FALSE]

  coverage[is.na(coverage)] <- 0
  methyl[is.na(methyl)] <- 0
  unmethyl[is.na(unmethyl)] <- 0
  hydroxy[is.na(hydroxy)] <- 0

  mean_coverage <- rowMeans(coverage)
  missing_prop  <- apply(coverage, 1, function(x) mean(x == 0))

  percent_methyl  <- rowMeans(methyl)  / (mean_coverage + 1e-9)
  percent_hydroxy <- rowMeans(hydroxy) / (mean_coverage + 1e-9)

  qc_df <- data.frame(
    Mean_Coverage = mean_coverage,
    Missing = missing_prop,
    Percent_Methyl = percent_methyl,
    Percent_Hydroxy = percent_hydroxy
  )

  make_plot <- function(df, var, xlab) {
    var_sym <- rlang::sym(var)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = !!var_sym))

    if (plot_type %in% c("hist", "both")) {
      p <- p + ggplot2::geom_histogram(
        bins = bins,
        fill = "steelblue",
        color = "black",
        alpha = ifelse(plot_type == "both", 0.6, 1)
      )
    }

    if (plot_type %in% c("density", "both")) {
      p <- p + ggplot2::geom_density(color = "firebrick", linewidth = 1)
    }

    p + ggplot2::theme_classic() +
      ggplot2::labs(x = xlab, y = "Frequency / Density")
  }

  p1 <- make_plot(qc_df, "Mean_Coverage", "Mean Coverage")
  p2 <- make_plot(qc_df, "Missing", "Proportion Missing (zero coverage)")
  p3 <- make_plot(qc_df, "Percent_Methyl", "Percent Methylated")
  p4 <- make_plot(qc_df, "Percent_Hydroxy", "Percent Hydroxymethylated")

  if (isTRUE(save)) {
    grDevices::pdf(paste0(prefix, "_QC_Plots.pdf"), width = 8, height = 6)
    print(p1); print(p2); print(p3); print(p4)
    grDevices::dev.off()
  } else {
    print(p1); print(p2); print(p3); print(p4)
  }

  qc_df
}
