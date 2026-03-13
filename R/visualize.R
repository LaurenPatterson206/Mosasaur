#' Compute modification composition for DMCs and/or DMRs from merged bedMethyl counts
#'
#' Computes group-wise composition (methyl, unmodified, hydroxy) at CpG sites
#' from a merged bedMethyl table, then optionally aggregates site composition
#' into DMRs by averaging CpGs that overlap each region.
#'
#' This function is only needed for composition-based visualizations such as
#' ternary plots. It is not required for volcano or Manhattan plots from
#' \code{\link{run_pqlseq}} output.
#'
#' @param merged A data.frame from \code{\link{merge_data}} containing \code{Coord}
#'   and per-sample columns beginning with \code{Methyl_}, \code{Unmethyl_},
#'   and \code{Hydroxymethyl_}.
#' @param pheno_file Path to a CSV phenotype file.
#' @param subject_id_var Phenotype column containing sample IDs that match the
#'   merged count column suffixes.
#' @param group_var Phenotype column defining the grouping variable.
#' @param group_levels Character vector of length 2 specifying the two groups to include.
#' @param dmcs Optional data.frame with \code{chr} and \code{start} to subset
#'   site compositions for DMC visualization.
#' @param dmrs Optional data.frame with \code{chr}, \code{start}, and \code{end}
#'   defining regions for DMR composition.
#' @param strip_sample_prefix_regex Optional regex to strip from merged sample
#'   suffixes before matching to phenotype IDs (e.g. \code{"^SS"}).
#'
#' @return A list with:
#' \describe{
#'   \item{site_props}{All site-level compositions for both groups.}
#'   \item{dmc_data}{Site compositions restricted to \code{dmcs}, or NULL.}
#'   \item{dmr_data}{Region compositions restricted to \code{dmrs}, or NULL.}
#' }
#'
#' @examples
#' \dontrun{
#' comp <- compute_modification_composition(
#'   merged = merged,
#'   dmrs = dmrs_anno,
#'   pheno_file = "NHP_Blood_Samples_Modified.csv",
#'   subject_id_var = "MATRR_ID",
#'   group_var = "Cohort",
#'   group_levels = c("10", "14"),
#'   strip_sample_prefix_regex = "^SS"
#' )
#' }
#'
#' @export
compute_modification_composition <- function(
    merged,
    pheno_file,
    subject_id_var = "subject_id",
    group_var = "group",
    group_levels = c("control", "case"),
    dmcs = NULL,
    dmrs = NULL,
    strip_sample_prefix_regex = NULL
) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required for compute_modification_composition().")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("Package 'IRanges' is required for compute_modification_composition().")
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required for compute_modification_composition().")
  }

  if (!("Coord" %in% names(merged))) {
    stop("merged must contain a Coord column.")
  }
  if (length(group_levels) != 2) {
    stop("group_levels must be length 2.")
  }

  pheno <- utils::read.csv(pheno_file, stringsAsFactors = FALSE)
  if (!(subject_id_var %in% names(pheno))) {
    stop("subject_id_var not found in phenotype table.")
  }
  if (!(group_var %in% names(pheno))) {
    stop("group_var not found in phenotype table.")
  }

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
    stop(
      "These phenotype sample IDs were not found in merged columns after stripping: ",
      paste(missing_ids, collapse = ", ")
    )
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

  site_props <- rbind(
    group_prop(as.character(group_levels[1])),
    group_prop(as.character(group_levels[2]))
  )

  dmc_data <- NULL
  if (!is.null(dmcs)) {
    dmcs_df <- as.data.frame(dmcs)
    if (!all(c("chr", "start") %in% names(dmcs_df))) {
      stop("dmcs must contain columns chr and start.")
    }
    dmcs_df$chr <- as.character(dmcs_df$chr)
    dmcs_df$start <- as.integer(dmcs_df$start)

    site_props$key <- paste(site_props$chr, site_props$start, site_props$group, sep = "|")

    dmc_expand <- rbind(
      transform(dmcs_df, group = as.character(group_levels[1])),
      transform(dmcs_df, group = as.character(group_levels[2]))
    )
    dmc_expand$key <- paste(dmc_expand$chr, dmc_expand$start, dmc_expand$group, sep = "|")

    dmc_data <- merge(
      dmc_expand,
      site_props[, c("key", "methyl", "unmodified", "hydroxy")],
      by = "key",
      all.x = TRUE
    )

    dmc_data <- dmc_data[
      stats::complete.cases(dmc_data[, c("methyl", "unmodified", "hydroxy")]),
      ,
      drop = FALSE
    ]
  }

  dmr_data <- NULL
  if (!is.null(dmrs)) {
    dmrs_df <- as.data.frame(dmrs)
    if (!all(c("chr", "start", "end") %in% names(dmrs_df))) {
      stop("dmrs must contain columns chr, start, and end.")
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
      if (length(hits) == 0) {
        next
      }

      qh <- S4Vectors::queryHits(hits)
      sh <- S4Vectors::subjectHits(hits)

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

    dmr_data <- if (length(dmr_rows) > 0) {
      do.call(rbind, dmr_rows)
    } else {
      data.frame()
    }
  }

  list(
    site_props = site_props,
    dmc_data = dmc_data,
    dmr_data = dmr_data
  )
}


#' Visualize differential methylation results with selectable plot types
#'
#' Creates one or more plots for site-level (DMC) or region-level (DMR) results.
#'
#' Supported plots:
#' \itemize{
#'   \item \code{"volcano"}: effect size vs significance
#'   \item \code{"manhattan"}: genome-wide significance across chromosomes
#'   \item \code{"density"}: distribution of beta values
#'   \item \code{"ternary"}: methyl / unmodified / hydroxy composition
#'   \item \code{"len_vs_p"}: region length vs significance
#'   \item \code{"p_along_genome"}: per-chromosome significance vs position
#' }
#'
#' For PQLseq2 site-level results, common usage is:
#' \preformatted{
#' plots = c("density", "volcano", "manhattan")
#' }
#'
#' For comb-p DMR results, common usage is:
#' \preformatted{
#' plots = c("manhattan", "len_vs_p")
#' }
#'
#' For ternary plots, input must already contain \code{methyl}, \code{unmodified},
#' and \code{hydroxy} columns (e.g. from \code{\link{compute_modification_composition}}).
#'
#' @param x A data.frame of differential methylation results.
#' @param plots Character vector of plot names to make.
#' @param chr_col Column name for chromosome.
#' @param start_col Column name for genomic start/position.
#' @param end_col Column name for genomic end.
#' @param p_col Column name for p-values. If NULL, tries \code{"pvalue"} then \code{"p"}.
#' @param effect_col Column name for effect size. If NULL, tries \code{"beta"}.
#' @param color_by Optional column name used for color.
#' @param facet_by Optional column name used for facetting.
#' @param point_alpha Point transparency.
#' @param point_size Point size.
#' @param highlight_p Optional p-value threshold for a horizontal line.
#' @param output_pdf Logical; if TRUE, save all requested plots to a PDF.
#' @param pdf_file Output PDF filename.
#'
#' @return A named list of ggplot objects.
#'
#' @examples
#' \dontrun{
#' plots <- visualize_diffmeth(
#'   x = res,
#'   plots = c("density", "volcano", "manhattan"),
#'   p_col = "pvalue",
#'   effect_col = "beta",
#'   highlight_p = 0.05
#' )
#'
#' comp <- compute_modification_composition(
#'   merged = merged,
#'   dmrs = dmrs_anno,
#'   pheno_file = "NHP_Blood_Samples_Modified.csv",
#'   subject_id_var = "MATRR_ID",
#'   group_var = "Cohort",
#'   group_levels = c("10", "14"),
#'   strip_sample_prefix_regex = "^SS"
#' )
#'
#' tern <- visualize_diffmeth(
#'   x = comp$dmr_data,
#'   plots = "ternary",
#'   color_by = "group",
#'   facet_by = "group"
#' )
#' }
#'
#' @importFrom rlang .data
#' @export
visualize_diffmeth <- function(
    x,
    plots = c("manhattan", "volcano", "density"),
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
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for visualize_diffmeth().")
  }

  df <- as.data.frame(x)

  if (!(chr_col %in% names(df))) {
    stop("Missing chromosome column: ", chr_col)
  }
  if (!(start_col %in% names(df))) {
    stop("Missing start/position column: ", start_col)
  }

  if (is.null(p_col)) {
    if ("pvalue" %in% names(df)) {
      p_col <- "pvalue"
    } else if ("p" %in% names(df)) {
      p_col <- "p"
    }
  }
  if (!is.null(p_col) && !(p_col %in% names(df))) {
    stop("p_col not found in input data: ", p_col)
  }

  if (is.null(effect_col)) {
    if ("beta" %in% names(df)) {
      effect_col <- "beta"
    }
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

  if (!is.null(effect_col) && effect_col %in% names(df)) {
    df[[effect_col]] <- suppressWarnings(as.numeric(df[[effect_col]]))
  }

  if (is.null(color_by)) {
    if ("group" %in% names(df)) {
      color_by <- "group"
    }
  }

  color_mapping <- !is.null(color_by) && color_by %in% names(df)
  facet_mapping <- !is.null(facet_by) && facet_by %in% names(df)

  out <- list()


  if ("volcano" %in% plots) {
    if (is.null(effect_col) || !(effect_col %in% names(df))) {
      stop("Volcano plot requires an effect column (e.g., effect_col = 'beta').")
    }
    if (!("neglog10p" %in% names(df))) {
      stop("Volcano plot requires p-values (set p_col).")
    }

    if (color_mapping) {
      p_vol <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data[[effect_col]], y = .data$neglog10p, color = .data[[color_by]])
      )
    } else {
      p_vol <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data[[effect_col]], y = .data$neglog10p)
      )
    }

    p_vol <- p_vol +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::labs(
        title = "Volcano plot",
        x = effect_col,
        y = "-log10(p-value)"
      ) +
      ggplot2::theme_bw()

    if (!is.null(highlight_p)) {
      p_vol <- p_vol +
        ggplot2::geom_hline(yintercept = -log10(highlight_p), linetype = 2)
    }

    if (facet_mapping) {
      p_vol <- p_vol +
        ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
    }

    out$volcano <- p_vol
  }


  if ("manhattan" %in% plots) {
    if (!("neglog10p" %in% names(df))) {
      stop("Manhattan plot requires p-values (set p_col).")
    }

    df_man <- df

    chr_num <- suppressWarnings(as.numeric(gsub("^chr", "", df_man[[chr_col]], ignore.case = TRUE)))
    chr_is_num <- !is.na(chr_num)

    chr_levels <- unique(df_man[[chr_col]])
    if (any(chr_is_num)) {
      chr_levels <- c(
        unique(df_man[[chr_col]][order(chr_num, na.last = NA)]),
        setdiff(unique(df_man[[chr_col]]), unique(df_man[[chr_col]][order(chr_num, na.last = NA)]))
      )
    }

    df_man[[chr_col]] <- factor(df_man[[chr_col]], levels = chr_levels)
    df_man <- df_man[order(df_man[[chr_col]], df_man[[start_col]]), , drop = FALSE]

    chr_max <- tapply(df_man[[start_col]], df_man[[chr_col]], max, na.rm = TRUE)
    chr_offsets <- c(0, cumsum(as.numeric(chr_max))[-length(chr_max)])
    names(chr_offsets) <- names(chr_max)

    df_man$cum_pos <- df_man[[start_col]] + chr_offsets[as.character(df_man[[chr_col]])]

    chr_min <- tapply(df_man$cum_pos, df_man[[chr_col]], min, na.rm = TRUE)
    chr_max2 <- tapply(df_man$cum_pos, df_man[[chr_col]], max, na.rm = TRUE)
    chr_mid <- (chr_min + chr_max2) / 2
    axis_df <- data.frame(chr = names(chr_mid), mid = as.numeric(chr_mid), stringsAsFactors = FALSE)

    if (color_mapping) {
      p_man <- ggplot2::ggplot(
        df_man,
        ggplot2::aes(x = .data$cum_pos, y = .data$neglog10p, color = .data[[color_by]])
      )
    } else {
      p_man <- ggplot2::ggplot(
        df_man,
        ggplot2::aes(x = .data$cum_pos, y = .data$neglog10p)
      )
    }

    p_man <- p_man +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::scale_x_continuous(breaks = axis_df$mid, labels = axis_df$chr) +
      ggplot2::labs(
        title = "Manhattan plot",
        x = "Chromosome",
        y = "-log10(p-value)"
      ) +
      ggplot2::theme_bw()

    if (!is.null(highlight_p)) {
      p_man <- p_man +
        ggplot2::geom_hline(yintercept = -log10(highlight_p), linetype = 2)
    }

    out$manhattan <- p_man
  }


  if ("density" %in% plots) {
    if (is.null(effect_col) || !(effect_col %in% names(df))) {
      stop("Density plot requires an effect column (e.g., effect_col = 'beta').")
    }

    if (color_mapping) {
      p_density <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data[[effect_col]], fill = .data[[color_by]])
      ) +
        ggplot2::geom_density(alpha = 0.4, color = "black")
    } else {
      p_density <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data[[effect_col]])
      ) +
        ggplot2::geom_density(alpha = 0.4, fill = "steelblue", color = "black")
    }

    p_density <- p_density +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::labs(
        title = "Distribution of differential methylation effect sizes",
        x = "Effect size (beta)",
        y = "Density"
      ) +
      ggplot2::theme_bw()

    if (facet_mapping) {
      p_density <- p_density +
        ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
    }

    out$density <- p_density
  }


  if ("ternary" %in% plots) {
    if (!requireNamespace("ggtern", quietly = TRUE)) {
      stop("Plot 'ternary' requested but package 'ggtern' is not installed.")
    }

    need_comp <- c("methyl", "unmodified", "hydroxy")
    miss <- setdiff(need_comp, names(df))
    if (length(miss) > 0) {
      stop("Ternary plot requires columns: ", paste(miss, collapse = ", "))
    }

    if (color_mapping) {
      p_tern <- ggtern::ggtern(
        df,
        ggtern::aes(
          x = .data$methyl,
          y = .data$unmodified,
          z = .data$hydroxy,
          color = .data[[color_by]]
        )
      )
    } else {
      p_tern <- ggtern::ggtern(
        df,
        ggtern::aes(
          x = .data$methyl,
          y = .data$unmodified,
          z = .data$hydroxy
        )
      )
    }

    p_tern <- p_tern +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::labs(
        title = "Composition (ternary)",
        x = "Methyl",
        y = "Unmodified",
        z = "Hydroxy"
      ) +
      ggplot2::theme_bw()

    if (facet_mapping) {
      p_tern <- p_tern +
        ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
    }

    out$ternary <- p_tern
  }


  if ("len_vs_p" %in% plots) {
    if (!("length_bp" %in% names(df))) {
      stop("len_vs_p requires an end column so region length can be computed.")
    }
    if (!("neglog10p" %in% names(df))) {
      stop("len_vs_p requires p-values (set p_col).")
    }

    if (color_mapping) {
      p_len <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$length_bp, y = .data$neglog10p, color = .data[[color_by]])
      )
    } else {
      p_len <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$length_bp, y = .data$neglog10p)
      )
    }

    p_len <- p_len +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(
        title = "Length vs -log10(p)",
        x = "Length (bp, log10 scale)",
        y = "-log10(p)"
      ) +
      ggplot2::theme_bw()

    if (!is.null(highlight_p)) {
      p_len <- p_len +
        ggplot2::geom_hline(yintercept = -log10(highlight_p), linetype = 2)
    }

    if (facet_mapping) {
      p_len <- p_len +
        ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
    }

    out$len_vs_p <- p_len
  }


  if ("p_along_genome" %in% plots) {
    if (!("neglog10p" %in% names(df))) {
      stop("p_along_genome requires p-values (set p_col).")
    }

    df_pg <- df
    df_pg$mid <- if (end_col %in% names(df_pg)) {
      (df_pg[[start_col]] + df_pg[[end_col]]) / 2
    } else {
      df_pg[[start_col]]
    }

    if (color_mapping) {
      p_pg <- ggplot2::ggplot(
        df_pg,
        ggplot2::aes(x = .data$mid, y = .data$neglog10p, color = .data[[color_by]])
      )
    } else {
      p_pg <- ggplot2::ggplot(
        df_pg,
        ggplot2::aes(x = .data$mid, y = .data$neglog10p)
      )
    }

    p_pg <- p_pg +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", chr_col)), scales = "free_x") +
      ggplot2::labs(
        title = "Signal along genome",
        x = "Position (bp)",
        y = "-log10(p)"
      ) +
      ggplot2::theme_bw()

    if (!is.null(highlight_p)) {
      p_pg <- p_pg +
        ggplot2::geom_hline(yintercept = -log10(highlight_p), linetype = 2)
    }

    out$p_along_genome <- p_pg
  }

  if (isTRUE(output_pdf) && length(out) > 0) {
    grDevices::pdf(pdf_file, width = 10, height = 6)
    for (nm in names(out)) {
      print(out[[nm]])
    }
    grDevices::dev.off()
  }

  out
}
