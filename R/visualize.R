#' Ternary visualization of composition for annotated DMCs and DMRs
#'
#' Creates ternary plots showing the relative composition of methylated,
#' unmodified, and hydroxymethylated signal for:
#' \itemize{
#'   \item Annotated DMCs (site-level features)
#'   \item Annotated DMRs (region-level features)
#' }
#'
#' Composition is computed from the merged bedMethyl count table produced by
#' \code{\link{merge_data}} by averaging counts across samples within each group
#' defined in a phenotype table. DMR composition is computed by aggregating
#' (mean) CpG site compositions within each region.
#'
#' @param merged A merged bedMethyl data.frame produced by \code{\link{merge_data}}.
#'   Must contain columns beginning with \code{"Methyl_"}, \code{"Unmethyl_"},
#'   and \code{"Hydroxymethyl_"}, as well as a \code{Coord} column.
#' @param dmcs Optional data.frame of annotated DMC results. Must include \code{chr}
#'   and \code{start}. If provided, ternary data are computed for those sites.
#' @param dmrs Optional data.frame of annotated DMR results. Must include \code{chr},
#'   \code{start}, and \code{end}. If provided, ternary data are computed by
#'   aggregating CpG sites that overlap each region.
#' @param pheno_file Path to a CSV phenotype file.
#' @param subject_id_var Column name in phenotype file containing sample IDs
#'   that match the merged count column suffixes (after optional stripping).
#' @param group_var Column name in phenotype file defining the grouping variable.
#'   Typically "group" or "Group".
#' @param group_levels Character vector of length 2 giving the two groups to plot,
#'   in the desired order (e.g., c("control","case") or c("C","HVHD")).
#' @param strip_sample_prefix_regex Optional regex to strip from merged sample
#'   column suffixes before matching to phenotype IDs (e.g. "^SSS").
#' @param output_pdf Logical; if TRUE writes plots to a PDF.
#' @param pdf_file Output PDF filename (used if output_pdf = TRUE).
#' @param point_alpha Numeric in (0,1]. Alpha for points.
#' @param point_size Numeric. Point size.
#' @param color_by Either "group" (default) or a column name present in dmcs/dmrs
#'   (e.g., "context") to color by annotation category. If set to a column name,
#'   group is used for facetting.
#'
#' @return A list with:
#' \describe{
#'   \item{dmc_data}{data.frame of ternary composition for DMCs (or NULL)}
#'   \item{dmr_data}{data.frame of ternary composition for DMRs (or NULL)}
#'   \item{dmc_plot}{ggtern plot for DMCs (or NULL)}
#'   \item{dmr_plot}{ggtern plot for DMRs (or NULL)}
#' }
#'
#' @details
#' This function requires the \code{ggtern} package. Install with:
#' \preformatted{
#' install.packages("ggtern")
#' }
#'
#' @examples
#' \dontrun{
#' out <- ternary_vis(
#'   merged = merged,
#'   dmcs = dmcs_anno,
#'   dmrs = dmrs_anno,
#'   pheno_file = "pheno.filtered.csv",
#'   subject_id_var = "MATRR_ID",
#'   group_var = "group",
#'   group_levels = c("control","case"),
#'   strip_sample_prefix_regex = "^SSS",
#'   output_pdf = TRUE,
#'   pdf_file = "ternary_plots.pdf"
#' )
#' }
#'
#' @importFrom rlang .data
#'
#' @export
ternary_vis <- function(
    merged,
    dmcs = NULL,
    dmrs = NULL,
    pheno_file,
    subject_id_var = "subject_id",
    group_var = "group",
    group_levels = c("control", "case"),
    strip_sample_prefix_regex = NULL,
    output_pdf = TRUE,
    pdf_file = "ternary_plots.pdf",
    point_alpha = 0.5,
    point_size = 1.2,
    color_by = "group"
) {
  #dependencies
  if (!requireNamespace("ggtern", quietly = TRUE)) {
    stop("Package 'ggtern' is required for ternary_vis(). Install with install.packages('ggtern').")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for ternary_vis().")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required for ternary_vis().")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("Package 'IRanges' is required for ternary_vis().")
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required for ternary_vis().")
  }

  if (!("Coord" %in% names(merged))) stop("merged must contain a Coord column.")
  if (length(group_levels) != 2) stop("group_levels must be length 2 (two groups).")

  pheno <- utils::read.csv(pheno_file, header = TRUE)
  if (!(subject_id_var %in% names(pheno))) stop("subject_id_var not found in phenotype table.")
  if (!(group_var %in% names(pheno))) stop("group_var not found in phenotype table.")

  pheno[[subject_id_var]] <- as.character(pheno[[subject_id_var]])
  pheno[[group_var]] <- as.character(pheno[[group_var]])

  pheno <- pheno[pheno[[group_var]] %in% group_levels, , drop = FALSE]
  pheno[[group_var]] <- factor(pheno[[group_var]], levels = group_levels)

  coord_split1 <- strsplit(as.character(merged$Coord), ":", fixed = TRUE)
  chr <- vapply(coord_split1, function(x) x[1], character(1))
  pos <- vapply(coord_split1, function(x) x[2], character(1))
  coord_split2 <- strsplit(pos, "-", fixed = TRUE)
  start <- as.integer(vapply(coord_split2, function(x) x[1], character(1)))
  end <- as.integer(vapply(coord_split2, function(x) x[2], character(1)))

  meth_cols <- grep("^Methyl_", names(merged), value = TRUE)
  unmod_cols <- grep("^Unmethyl_", names(merged), value = TRUE)
  hydro_cols <- grep("^Hydroxymethyl_", names(merged), value = TRUE)

  if (length(meth_cols) == 0 || length(unmod_cols) == 0 || length(hydro_cols) == 0) {
    stop("merged must contain Methyl_*, Unmethyl_*, and Hydroxymethyl_* columns.")
  }

  suffix <- sub("^Methyl_", "", meth_cols)
  if (!is.null(strip_sample_prefix_regex)) {
    suffix <- sub(strip_sample_prefix_regex, "", suffix)
  }

  idx <- match(pheno[[subject_id_var]], suffix)
  if (any(is.na(idx))) {
    missing_ids <- pheno[[subject_id_var]][is.na(idx)]
    stop("These phenotype sample IDs were not found in merged columns after stripping: ",
         paste(missing_ids, collapse = ", "))
  }

  meth_mat <- as.matrix(merged[, meth_cols, drop = FALSE])
  unmod_mat <- as.matrix(merged[, unmod_cols, drop = FALSE])
  hydro_mat <- as.matrix(merged[, hydro_cols, drop = FALSE])

  meth_mat[is.na(meth_mat)] <- 0
  unmod_mat[is.na(unmod_mat)] <- 0
  hydro_mat[is.na(hydro_mat)] <- 0

  eps <- 1e-9

  group_prop <- function(group_name) {
    sids <- pheno[[subject_id_var]][pheno[[group_var]] == group_name]
    j <- match(sids, suffix)

    m_mean <- rowMeans(meth_mat[, j, drop = FALSE])
    u_mean <- rowMeans(unmod_mat[, j, drop = FALSE])
    h_mean <- rowMeans(hydro_mat[, j, drop = FALSE])

    tot <- m_mean + u_mean + h_mean + eps

    data.frame(
      chr = chr,
      start = start,
      end = end,
      methyl = m_mean / tot,
      unmodified = u_mean / tot,
      hydroxy = h_mean / tot,
      group = group_name,
      stringsAsFactors = FALSE
    )
  }

  site_g1 <- group_prop(as.character(group_levels[1]))
  site_g2 <- group_prop(as.character(group_levels[2]))
  site_props <- rbind(site_g1, site_g2)

  site_gr <- GenomicRanges::GRanges(
    seqnames = site_props$chr,
    ranges = IRanges::IRanges(start = site_props$start, end = site_props$end)
  )

  dmc_data <- NULL
  dmc_plot <- NULL

  if (!is.null(dmcs)) {
    dmcs_df <- as.data.frame(dmcs)
    if (!all(c("chr", "start") %in% names(dmcs_df))) {
      stop("dmcs must contain columns chr and start.")
    }
    dmcs_df$chr <- as.character(dmcs_df$chr)
    dmcs_df$start <- as.integer(dmcs_df$start)


    key <- paste(site_props$chr, site_props$start, site_props$group, sep = "|")
    site_props$key <- key

    dkey <- paste(dmcs_df$chr, dmcs_df$start, sep = "|")

    dmc_expand <- rbind(
      transform(dmcs_df, group = as.character(group_levels[1]), key = paste(chr, start, as.character(group_levels[1]), sep = "|")),
      transform(dmcs_df, group = as.character(group_levels[2]), key = paste(chr, start, as.character(group_levels[2]), sep = "|"))
    )

    dmc_data <- merge(
      dmc_expand,
      site_props[, c("key", "methyl", "unmodified", "hydroxy")],
      by = "key",
      all.x = TRUE
    )


    dmc_data <- dmc_data[stats::complete.cases(dmc_data[, c("methyl", "unmodified", "hydroxy")]), , drop = FALSE]

    #plotting
    aes_color <- if (color_by %in% names(dmc_data)) color_by else "group"
    p <- ggtern::ggtern(
      data = dmc_data,
      ggtern::aes(x = .data$methyl, y = .data$unmodified, z = .data$hydroxy)
    ) +
      ggplot2::geom_point(ggplot2::aes(color = .data[[aes_color]]), alpha = point_alpha, size = point_size) +
      ggplot2::labs(
        title = "DMC composition (ternary)",
        x = "Methylated",
        y = "Unmodified",
        z = "Hydroxymethylated"
      ) +
      ggplot2::theme_bw()

    if (aes_color != "group") {
      p <- p + ggplot2::facet_wrap(~group)
    }

    dmc_plot <- p
  }

  #DMR ternary
  dmr_data <- NULL
  dmr_plot <- NULL

  if (!is.null(dmrs)) {
    dmrs_df <- as.data.frame(dmrs)
    if (!all(c("chr", "start", "end") %in% names(dmrs_df))) {
      stop("dmrs must contain columns chr, start, end.")
    }
    dmrs_df$chr <- as.character(dmrs_df$chr)
    dmrs_df$start <- as.integer(dmrs_df$start)
    dmrs_df$end <- as.integer(dmrs_df$end)

    dmr_gr <- GenomicRanges::GRanges(
      seqnames = dmrs_df$chr,
      ranges = IRanges::IRanges(start = dmrs_df$start, end = dmrs_df$end)
    )

    dmr_rows <- list()
    for (g in as.character(group_levels)) {
      spg <- site_props[site_props$group == g, , drop = FALSE]
      sgr <- GenomicRanges::GRanges(
        seqnames = spg$chr,
        ranges = IRanges::IRanges(start = spg$start, end = spg$end)
      )
      hits <- GenomicRanges::findOverlaps(dmr_gr, sgr)

      qh <- S4Vectors::queryHits(hits)
      sh <- S4Vectors::subjectHits(hits)

      if (length(qh) == 0) {
        next
      }

      tmp <- data.frame(
        region_idx = qh,
        methyl = spg$methyl[sh],
        unmodified = spg$unmodified[sh],
        hydroxy = spg$hydroxy[sh],
        stringsAsFactors = FALSE
      )

      agg <- stats::aggregate(
        tmp[, c("methyl", "unmodified", "hydroxy")],
        by = list(region_idx = tmp$region_idx),
        FUN = mean
      )

      rr <- dmrs_df[agg$region_idx, , drop = FALSE]
      rr$group <- g
      rr$methyl <- agg$methyl
      rr$unmodified <- agg$unmodified
      rr$hydroxy <- agg$hydroxy

      dmr_rows[[g]] <- rr
    }

    if (length(dmr_rows) > 0) {
      dmr_data <- do.call(rbind, dmr_rows)

      aes_color <- if (color_by %in% names(dmr_data)) color_by else "group"
      p <- ggtern::ggtern(
        data = dmr_data,
        ggtern::aes(x = .data$methyl, y = .data$unmodified, z = .data$dhydroxy)
      ) +
        ggplot2::geom_point(ggplot2::aes(color = .data[[aes_color]]), alpha = point_alpha, size = point_size) +
        ggplot2::labs(
          title = "DMR composition (ternary)",
          x = "Methylated",
          y = "Unmodified",
          z = "Hydroxymethylated"
        ) +
        ggplot2::theme_bw()

      if (aes_color != "group") {
        p <- p + ggplot2::facet_wrap(~group)
      }

      dmr_plot <- p
    } else {
      dmr_data <- data.frame()
      dmr_plot <- NULL
    }
  }

  #write PDF if requested by user
  if (output_pdf) {
    grDevices::pdf(pdf_file)
    if (!is.null(dmc_plot)) print(dmc_plot)
    if (!is.null(dmr_plot)) print(dmr_plot)
    grDevices::dev.off()
  }

  list(
    dmc_data = dmc_data,
    dmr_data = dmr_data,
    dmc_plot = dmc_plot,
    dmr_plot = dmr_plot
  )
}
