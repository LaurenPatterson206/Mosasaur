#' Read and convert EPIC methylation array beta values to Mosasaur count matrices
#'
#' Reads EPIC (Infinium) methylation array beta values and converts them into
#' pseudo-count matrices compatible with Mosasaur's downstream modeling
#' (e.g., PQLseq2). Probe coordinates are obtained from a manifest/annotation file
#' and can optionally be lifted over to a target genome using a UCSC chain file.
#'
#' Conversion uses a constant pseudo-coverage per probe per sample
#' (\code{coverage_per_probe}):
#' \itemize{
#'   \item \code{TCounts} = \code{coverage_per_probe}
#'   \item \code{MCounts} = \code{round(beta * coverage_per_probe)}
#' }
#'
#' @param beta_file Path to a tab-delimited table of beta values with one row per
#'   probe and one column per sample, plus a probe ID column.
#' @param sample_info_file Path to a CSV of sample metadata containing a column
#'   with sample IDs matching beta table column names.
#' @param annotation_file Path to the EPIC manifest/annotation CSV.
#' @param perfect_matches_file Optional CSV of probe IDs to retain (e.g., perfect
#'   rhesus matches). If provided, only probes in this file are kept.
#' @param beta_probe_id_col Column name in \code{beta_file} containing probe IDs.
#' @param anno_probe_id_col Column name in \code{annotation_file} for probe IDs.
#' @param anno_chr_col Column name in \code{annotation_file} for chromosome.
#' @param anno_pos_col Column name in \code{annotation_file} for probe position.
#' @param anno_gene_col Optional column name in \code{annotation_file} for gene IDs.
#' @param sample_id_col Column name in \code{sample_info_file} listing sample IDs.
#' @param liftover_chain Optional path to a UCSC chain file used for liftover.
#' @param liftover_seqstyle Sequence naming style (e.g., \code{"UCSC"}).
#' @param probe_width Integer probe window size (bp) used for liftover. Defaults to 40.
#' @param coverage_per_probe Integer pseudo-coverage used to convert beta to counts.
#' @param write_outputs Logical; if TRUE, writes \code{<prefix>_MCounts.txt} and
#'   \code{<prefix>_TCounts.txt}.
#' @param prefix Output prefix used when writing files.
#' @param return_annotation Logical; if TRUE, return probe annotation used.
#'
#' @return A list with:
#' \describe{
#'   \item{MCounts}{data.frame with \code{chr}, \code{start}, and sample columns}
#'   \item{TCounts}{data.frame with \code{chr}, \code{start}, and sample columns}
#'   \item{annotation}{(optional) data.frame of probe coordinates and gene IDs}
#' }
#'
#' @examples
#' \dontrun{
#' out <- read_epic_array(
#'   beta_file = "EPIC_beta.txt",
#'   sample_info_file = "EPIC_samples.csv",
#'   annotation_file = "EPIC_manifest.csv",
#'   perfect_matches_file = "Perfect_MM_Matches.csv",
#'   liftover_chain = "hg19ToRheMac10.over.chain",
#'   prefix = "EPIC"
#' )
#' }
#'
#' @importFrom utils read.csv
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom methods as
#'
#' @export
read_epic_array <- function(
    beta_file,
    sample_info_file,
    annotation_file,
    perfect_matches_file = NULL,
    beta_probe_id_col = "TargetID",
    anno_probe_id_col = "IlmnID",
    anno_chr_col = "CHR",
    anno_pos_col = "MAPINFO",
    anno_gene_col = "UCSC_RefGene_Name",
    sample_id_col = "Sample_ID",
    liftover_chain = NULL,
    liftover_seqstyle = "UCSC",
    probe_width = 40L,
    coverage_per_probe = 500L,
    write_outputs = TRUE,
    prefix = "EPIC",
    return_annotation = TRUE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) stop("GenomicRanges is required.")
  if (!requireNamespace("IRanges", quietly = TRUE)) stop("IRanges is required.")
  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) stop("GenomeInfoDb is required.")
  if (!is.null(liftover_chain) && !requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("rtracklayer is required when liftover_chain is provided.")
  }

  beta <- as.data.frame(data.table::fread(beta_file))
  anno <- as.data.frame(data.table::fread(annotation_file))
  samp <- utils::read.csv(sample_info_file, header = TRUE)

  if (!(beta_probe_id_col %in% names(beta))) stop("beta_probe_id_col not found in beta table.")
  if (!(sample_id_col %in% names(samp))) stop("sample_id_col not found in sample info.")
  if (!all(c(anno_probe_id_col, anno_chr_col, anno_pos_col) %in% names(anno))) {
    stop("annotation_file missing required columns.")
  }

  #optional probe whitelist
  keep_probes <- NULL
  if (!is.null(perfect_matches_file)) {
    pm <- utils::read.csv(perfect_matches_file, header = FALSE)
    if (ncol(pm) < 1) stop("perfect_matches_file does not appear to contain probe IDs.")
    keep_probes <- unique(pm[[1]])
  }

  anno$ProbeID <- anno[[anno_probe_id_col]]
  gene_vec <- if (anno_gene_col %in% names(anno)) anno[[anno_gene_col]] else NA_character_

  pos <- suppressWarnings(as.numeric(anno[[anno_pos_col]]))
  chr <- as.character(anno[[anno_chr_col]])

  half <- as.integer(probe_width / 2)
  anno2 <- data.frame(
    seqnames = chr,
    start = pos - half,
    end = pos + half - 1,
    ProbeID = anno$ProbeID,
    GeneID = gene_vec,
    stringsAsFactors = FALSE
  )
  anno2 <- stats::na.omit(anno2)

  if (!is.null(keep_probes)) {
    anno2 <- anno2[anno2$ProbeID %in% keep_probes, , drop = FALSE]
  }

  #optional liftOver
  if (!is.null(liftover_chain)) {
    gr <- GenomicRanges::makeGRangesFromDataFrame(anno2, keep.extra.columns = TRUE)
    GenomeInfoDb::seqlevelsStyle(gr) <- liftover_seqstyle

    chain <- rtracklayer::import.chain(liftover_chain)
    lifted_list <- rtracklayer::liftOver(gr, chain)


    lifted_gr <- as(unlist(lifted_list), "GenomicRanges::GRanges")
    lifted_df <- as.data.frame(lifted_gr)


    lifted_df <- lifted_df[lifted_df$width == probe_width, , drop = FALSE]

    lifted_df$chr <- sub("^chr", "", as.character(lifted_df$seqnames))

    lifted_df$start <- lifted_df$start + half

    anno_final <- lifted_df[, c("ProbeID", "chr", "start", "GeneID"), drop = FALSE]
  } else {
    anno2$chr <- sub("^chr", "", as.character(anno2$seqnames))
    anno2$start <- anno2$start + half
    anno_final <- anno2[, c("ProbeID", "chr", "start", "GeneID"), drop = FALSE]
  }

  beta_merged <- merge(
    anno_final,
    beta,
    by.x = "ProbeID",
    by.y = beta_probe_id_col
  )

  samp_ids <- as.character(samp[[sample_id_col]])
  missing_cols <- setdiff(samp_ids, colnames(beta_merged))
  if (length(missing_cols) > 0) {
    warning("Dropping sample IDs not found in beta table: ", paste(missing_cols, collapse = ", "))
    samp_ids <- intersect(samp_ids, colnames(beta_merged))
  }
  if (length(samp_ids) == 0) stop("No sample columns matched between sample info and beta table.")

  chr_info <- beta_merged[, c("chr", "start"), drop = FALSE]
  beta_mat <- beta_merged[, samp_ids, drop = FALSE]

  TCounts <- beta_mat
  TCounts[,] <- as.integer(coverage_per_probe)
  TCounts <- cbind(chr_info, TCounts)

  MCounts <- round(beta_mat * as.integer(coverage_per_probe), 0)
  MCounts <- cbind(chr_info, MCounts)

  if (isTRUE(write_outputs)) {
    data.table::fwrite(TCounts, paste0(prefix, "_TCounts.txt"), sep = "\t")
    data.table::fwrite(MCounts, paste0(prefix, "_MCounts.txt"), sep = "\t")
  }

  out <- list(MCounts = as.data.frame(MCounts), TCounts = as.data.frame(TCounts))
  if (isTRUE(return_annotation)) out$annotation <- anno_final
  out
}
