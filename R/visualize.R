#' Compute composition for DMCs and/or DMRs from merged bedMethyl counts
#'
#' Computes group-wise composition (methyl, unmodified, hydroxy) at CpG sites
#' from a merged bedMethyl table, then optionally aggregates site composition
#' into DMRs by averaging CpGs that overlap each region.
#'
#' @param merged A data.frame from \code{\link{merge_data}} containing \code{Coord}
#'   and per-sample columns \code{Methyl_*}, \code{Unmethyl_*}, \code{Hydroxymethyl_*}.
#' @param pheno_file Path to a CSV phenotype file.
#' @param subject_id_var Phenotype column giving sample IDs (must match merged suffixes).
#' @param group_var Phenotype column defining groups (e.g., Cohort, Group).
#' @param group_levels Character vector of length 2 specifying the two groups to include,
#'   in plotting order.
#' @param dmcs Optional data.frame with \code{chr} and \code{start} to subset site results.
#' @param dmrs Optional data.frame with \code{chr}, \code{start}, \code{end} defining regions.
#' @param strip_sample_prefix_regex Optional regex to strip from merged sample suffixes
#'   (e.g., "^SS" to turn "SS10208" into "10208").
#'
#' @return A list with:
#' \describe{
#'   \item{site_props}{All site-level compositions for both groups (data.frame).}
#'   \item{dmc_data}{Site compositions restricted to dmcs (or NULL).}
#'   \item{dmr_data}{Region compositions restricted to dmrs (or NULL).}
#' }
#'
#' @export
compute_region_composition <- function(
    merged,
    pheno_file,
    subject_id_var = "subject_id",
    group_var = "group",
    group_levels = c("control", "case"),
    dmcs = NULL,
    dmrs = NULL,
    strip_sample_prefix_regex = NULL
) {

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) stop("Need GenomicRanges.")
  if (!requireNamespace("IRanges", quietly = TRUE)) stop("Need IRanges.")
  if (!requireNamespace("S4Vectors", quietly = TRUE)) stop("Need S4Vectors.")

  if (!("Coord" %in% names(merged))) stop("merged must contain a Coord column.")
  if (length(group_levels) != 2) stop("group_levels must be length 2.")

  pheno <- utils::read.csv(pheno_file, stringsAsFactors = FALSE)
  if (!(subject_id_var %in% names(pheno))) stop("subject_id_var not found in phenotype.")
  if (!(group_var %in% names(pheno))) stop("group_var not found in phenotype.")

  pheno[[subject_id_var]] <- as.character(pheno[[subject_id_var]])
  pheno[[group_var]] <- as.character(pheno[[group_var]])
  pheno <- pheno[pheno[[group_var]] %in% group_levels, , drop = FALSE]
  pheno[[group_var]] <- factor(pheno[[group_var]], levels = group_levels)


  coord_split1 <- strsplit(as.character(merged$Coord), ":", fixed = TRUE)
  chr <- vapply(coord_split1, `[[`, character(1), 1)
  pos <- vapply(coord_split1, `[[`, character(1), 2)
  coord_split2 <- strsplit(pos, "-", fixed = TRUE)
  start <- as.integer(vapply(coord_split2, `[[`, character(1), 1))
  end   <- as.integer(vapply(coord_split2, `[[`, character(1), 2))

  meth_cols  <- grep("^Methyl_", names(merged), value = TRUE)
  unmod_cols <- grep("^Unmethyl_", names(merged), value = TRUE)
  hydro_cols <- grep("^Hydroxymethyl_", names(merged), value = TRUE)

  if (length(meth_cols) == 0 || length(unmod_cols) == 0 || length(hydro_cols) == 0) {
    stop("merged must contain Methyl_*, Unmethyl_*, Hydroxymethyl_* columns.")
  }


  suffix <- sub("^Methyl_", "", meth_cols)
  if (!is.null(strip_sample_prefix_regex)) {
    suffix <- sub(strip_sample_prefix_regex, "", suffix)
  }


  idx <- match(pheno[[subject_id_var]], suffix)
  if (any(is.na(idx))) {
    missing_ids <- pheno[[subject_id_var]][is.na(idx)]
    stop("These phenotype sample IDs were not found in merged columns: ",
         paste(missing_ids, collapse = ", "))
  }


  meth_mat  <- as.matrix(merged[, meth_cols, drop = FALSE])
  unmod_mat <- as.matrix(merged[, unmod_cols, drop = FALSE])
  hydro_mat <- as.matrix(merged[, hydro_cols, drop = FALSE])

  meth_mat[is.na(meth_mat)] <- 0
  unmod_mat[is.na(unmod_mat)] <- 0
  hydro_mat[is.na(hydro_mat)] <- 0

  eps <- 1e-9

  group_prop <- function(g) {
    sids <- pheno[[subject_id_var]][pheno[[group_var]] == g]
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
      group = g,
      stringsAsFactors = FALSE
    )
  }

  site_props <- rbind(group_prop(as.character(group_levels[1])),
                      group_prop(as.character(group_levels[2])))


  dmc_data <- NULL
  if (!is.null(dmcs)) {
    dmcs_df <- as.data.frame(dmcs)
    if (!all(c("chr", "start") %in% names(dmcs_df))) stop("dmcs must have chr and start.")
    dmcs_df$chr <- as.character(dmcs_df$chr)
    dmcs_df$start <- as.integer(dmcs_df$start)

    key_site <- paste(site_props$chr, site_props$start, site_props$group, sep = "|")

    dmc_expand <- rbind(
      transform(dmcs_df, group = as.character(group_levels[1])),
      transform(dmcs_df, group = as.character(group_levels[2]))
    )
    dmc_expand$key <- paste(dmc_expand$chr, dmc_expand$start, dmc_expand$group, sep = "|")

    dmc_data <- merge(
      dmc_expand,
      data.frame(key = key_site,
                 methyl = site_props$methyl,
                 unmodified = site_props$unmodified,
                 hydroxy = site_props$hydroxy,
                 stringsAsFactors = FALSE),
      by = "key",
      all.x = TRUE
    )
    dmc_data <- dmc_data[stats::complete.cases(dmc_data[, c("methyl","unmodified","hydroxy")]), , drop = FALSE]
  }


  dmr_data <- NULL
  if (!is.null(dmrs)) {
    dmrs_df <- as.data.frame(dmrs)
    if (!all(c("chr", "start", "end") %in% names(dmrs_df))) stop("dmrs must have chr, start, end.")
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
      if (length(hits) == 0) next

      qh <- S4Vectors::queryHits(hits)
      sh <- S4Vectors::subjectHits(hits)

      tmp <- data.frame(
        region_idx = qh,
        methyl = spg$methyl[sh],
        unmodified = spg$unmodified[sh],
        hydroxy = spg$hydroxy[sh],
        stringsAsFactors = FALSE
      )

      agg <- stats::aggregate(tmp[, c("methyl","unmodified","hydroxy")],
                              by = list(region_idx = tmp$region_idx),
                              FUN = mean)

      rr <- dmrs_df[agg$region_idx, , drop = FALSE]
      rr$group <- g
      rr$methyl <- agg$methyl
      rr$unmodified <- agg$unmodified
      rr$hydroxy <- agg$hydroxy

      dmr_rows[[g]] <- rr
    }

    dmr_data <- if (length(dmr_rows)) do.call(rbind, dmr_rows) else data.frame()
  }

  list(site_props = site_props, dmc_data = dmc_data, dmr_data = dmr_data)
}


#' Visualize differential methylation results with selectable plot types
#'
#' Creates one or more plots for site-level (DMC) or region-level (DMR) results.
#' Supports ternary composition (if methyl/unmodified/hydroxy are present),
#' volcano plots (effect vs significance), and Manhattan plots (genome-wide view).
#'
#' @param x A data.frame of results. Must include \code{chr} and \code{start}
#'   (and \code{end} for region-length plots).
#' @param plots Character vector of plot names to make. Options:
#'   \code{"ternary"}, \code{"len_vs_p"}, \code{"p_along_genome"},
#'   \code{"volcano"}, \code{"manhattan"}.
#' @param chr_col,start_col,end_col Column names for coordinates.
#' @param p_col Column name for p-values. If NULL, tries \code{"pvalue"} then \code{"p"}.
#' @param effect_col Column name for effect size (volcano). If NULL, tries \code{"beta"}.
#' @param color_by Column name used for color mapping (e.g., \code{"group"}, \code{"context"}).
#'   If NULL, uses \code{"group"} if present.
#' @param facet_by Optional column name to facet by (often \code{"group"}).
#' @param point_alpha Alpha for points.
#' @param point_size Point size.
#' @param highlight_p Optional threshold; if provided, adds a horizontal line at -log10(threshold).
#' @param output_pdf Logical; if TRUE writes plots to PDF.
#' @param pdf_file Output PDF filename (if output_pdf=TRUE).
#'
#' @return A named list of ggplot objects.
#'
#' @export
visualize_diffmeth <- function(
    x,
    plots = c("manhattan", "volcano"),
    chr_col = "chr",
    start_col = "start",
    end_col = "end",
    p_col = NULL,
    effect_col = NULL,
    color_by = NULL,
    facet_by = NULL,
    point_alpha = 0.6,
    point_size = 1.6,
    highlight_p = NULL,
    output_pdf = FALSE,
    pdf_file = "diffmeth_plots.pdf"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Need ggplot2.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Need dplyr.")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Need rlang.")

  df <- as.data.frame(x)


  if (!(chr_col %in% names(df))) stop("Missing chromosome column: ", chr_col)
  if (!(start_col %in% names(df))) stop("Missing start/position column: ", start_col)


  if (is.null(p_col)) {
    if ("pvalue" %in% names(df)) p_col <- "pvalue"
    else if ("p" %in% names(df)) p_col <- "p"
  }
  if (!is.null(p_col) && !(p_col %in% names(df))) stop("p_col not found: ", p_col)


  if (is.null(effect_col)) {
    if ("beta" %in% names(df)) effect_col <- "beta"
  }
  if ("volcano" %in% plots && (is.null(effect_col) || !(effect_col %in% names(df)))) {
    stop("Volcano plot requested but no effect column found. Provide effect_col= (e.g., 'beta').")
  }


  df[[chr_col]] <- as.character(df[[chr_col]])
  df[[start_col]] <- suppressWarnings(as.numeric(df[[start_col]]))

  if (end_col %in% names(df)) {
    df[[end_col]] <- suppressWarnings(as.numeric(df[[end_col]]))
    df$length_bp <- pmax(1, df[[end_col]] - df[[start_col]])
  }


  if (!is.null(p_col)) {
    df[[p_col]] <- suppressWarnings(as.numeric(df[[p_col]]))
    df$neglog10p <- -log10(pmax(df[[p_col]], 1e-300))
  }


  if (is.null(color_by)) {
    color_by <- if ("group" %in% names(df)) "group" else NULL
  }
  color_aes <- if (!is.null(color_by) && color_by %in% names(df)) ggplot2::aes(color = .data[[color_by]]) else NULL

  out <- list()


  if ("ternary" %in% plots) {
    if (!requireNamespace("ggtern", quietly = TRUE)) stop("Need ggtern for ternary.")
    need_comp <- c("methyl","unmodified","hydroxy")
    miss <- setdiff(need_comp, names(df))
    if (length(miss)) stop("Ternary requires: ", paste(miss, collapse = ", "))

    p_tern <- ggtern::ggtern(df, ggtern::aes(x = .data$methyl, y = .data$unmodified, z = .data$hydroxy)) +
      ggplot2::geom_point(mapping = color_aes, alpha = point_alpha, size = point_size) +
      ggplot2::labs(title = "Composition (ternary)", x = "Methyl", y = "Unmodified", z = "Hydroxy") +
      ggplot2::theme_bw()

    if (!is.null(facet_by) && facet_by %in% names(df)) {
      p_tern <- p_tern + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
    }
    out$ternary <- p_tern
  }


  if ("len_vs_p" %in% plots) {
    if (!("length_bp" %in% names(df))) stop("len_vs_p needs end_col present to compute length.")
    if (!("neglog10p" %in% names(df))) stop("len_vs_p needs p-values (p_col).")

    p_len <- ggplot2::ggplot(df, ggplot2::aes(x = .data$length_bp, y = .data$neglog10p)) +
      ggplot2::geom_point(mapping = color_aes, alpha = point_alpha, size = point_size) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Length vs -log10(p)", x = "Length (bp, log10)", y = "-log10(p)") +
      ggplot2::theme_bw()

    if (!is.null(facet_by) && facet_by %in% names(df)) {
      p_len <- p_len + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
    }
    if (!is.null(highlight_p)) {
      p_len <- p_len + ggplot2::geom_hline(yintercept = -log10(highlight_p), linetype = 2)
    }
    out$len_vs_p <- p_len
  }


  if ("volcano" %in% plots) {
    if (!("neglog10p" %in% names(df))) stop("Volcano needs p-values (p_col).")

    df[[effect_col]] <- suppressWarnings(as.numeric(df[[effect_col]]))

    p_vol <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[effect_col]], y = .data$neglog10p)) +
      ggplot2::geom_point(mapping = color_aes, alpha = point_alpha, size = point_size) +
      ggplot2::labs(
        title = "Volcano plot",
        x = effect_col,
        y = "-log10(p)"
      ) +
      ggplot2::theme_bw()

    if (!is.null(facet_by) && facet_by %in% names(df)) {
      p_vol <- p_vol + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
    }
    if (!is.null(highlight_p)) {
      p_vol <- p_vol + ggplot2::geom_hline(yintercept = -log10(highlight_p), linetype = 2)
    }
    out$volcano <- p_vol
  }


  if ("manhattan" %in% plots) {
    if (!("neglog10p" %in% names(df))) stop("Manhattan needs p-values (p_col).")


    chr_raw <- df[[chr_col]]
    chr_num <- suppressWarnings(as.numeric(gsub("^chr", "", chr_raw, ignore.case = TRUE)))
    chr_is_num <- !is.na(chr_num)

    chr_levels <- unique(chr_raw)
    if (any(chr_is_num)) {

      chr_levels <- c(
        unique(chr_raw[order(chr_num, na.last = NA)]),
        setdiff(unique(chr_raw), unique(chr_raw[order(chr_num, na.last = NA)]))
      )
    }
    df[[chr_col]] <- factor(df[[chr_col]], levels = chr_levels)


    df <- df[order(df[[chr_col]], df[[start_col]]), , drop = FALSE]
    chr_max <- tapply(df[[start_col]], df[[chr_col]], max, na.rm = TRUE)
    chr_offsets <- c(0, cumsum(as.numeric(chr_max))[-length(chr_max)])
    names(chr_offsets) <- names(chr_max)

    df$cum_pos <- df[[start_col]] + chr_offsets[as.character(df[[chr_col]])]


    chr_min <- tapply(df$cum_pos, df[[chr_col]], min, na.rm = TRUE)
    chr_max2 <- tapply(df$cum_pos, df[[chr_col]], max, na.rm = TRUE)
    chr_mid <- (chr_min + chr_max2) / 2
    axis_df <- data.frame(chr = names(chr_mid), mid = as.numeric(chr_mid))

    p_man <- ggplot2::ggplot(df, ggplot2::aes(x = .data$cum_pos, y = .data$neglog10p)) +
      ggplot2::geom_point(mapping = color_aes, alpha = point_alpha, size = point_size) +
      ggplot2::scale_x_continuous(breaks = axis_df$mid, labels = axis_df$chr) +
      ggplot2::labs(title = "Manhattan plot", x = "Chromosome", y = "-log10(p)") +
      ggplot2::theme_bw()

    if (!is.null(highlight_p)) {
      p_man <- p_man + ggplot2::geom_hline(yintercept = -log10(highlight_p), linetype = 2)
    }
    out$manhattan <- p_man
  }


  if (isTRUE(output_pdf) && length(out)) {
    grDevices::pdf(pdf_file, width = 10, height = 6)
    for (nm in names(out)) print(out[[nm]])
    grDevices::dev.off()
  }

  out
}
