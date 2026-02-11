#' Annotate differentially methylated cytosines (DMCs) and regions (DMRs)
#'
#' Adds basic genomic context annotations to either site-level (DMC) or
#' region-level (DMR) results using overlaps with features derived from a
#' transcript database (TxDb).
#'
#' Two annotation levels are supported:
#' \itemize{
#'   \item \code{level = "DMC"}: annotates single cytosines/CpGs. If \code{end}
#'     is missing, it is created as \code{start + 1}.
#'   \item \code{level = "DMR"}: annotates regions. Requires \code{chr}, \code{start},
#'     and \code{end}.
#' }
#'
#' Added annotations include overlap flags (gene/exon/intron/promoter), a context
#' label, nearest gene identifier, and distance to nearest gene.
#'
#' @param x A data.frame containing genomic coordinates. For DMCs must include
#'   \code{chr} and \code{start}. For DMRs must include \code{chr}, \code{start},
#'   and \code{end}.
#' @param level One of \code{"DMC"} or \code{"DMR"}.
#' @param txdb Optional TxDb object. If provided, \code{gtf_file} is ignored.
#' @param gtf_file Optional path to a GTF/GFF file used to build a TxDb via
#'   \code{GenomicFeatures::makeTxDbFromGFF()}. Required if \code{txdb} not provided.
#' @param chr_col Column name for chromosome.
#' @param start_col Column name for start coordinate.
#' @param end_col Column name for end coordinate.
#' @param promoter_upstream Bases upstream used to define promoter regions.
#' @param promoter_downstream Bases downstream used to define promoter regions.
#'
#' @return A data.frame with original columns plus:
#' \describe{
#'   \item{context}{Context string: gene/intergenic plus feature tags.}
#'   \item{overlaps_gene}{Logical overlap with gene bodies.}
#'   \item{overlaps_exon}{Logical overlap with exons.}
#'   \item{overlaps_intron}{Logical overlap with introns.}
#'   \item{overlaps_promoter}{Logical overlap with promoters.}
#'   \item{nearest_gene_id}{Nearest gene ID from TxDb (may be NA).}
#'   \item{distance_to_nearest_gene}{Distance to nearest gene (bp).}
#' }
#'
#' @examples
#' \dontrun{
#' dmcs <- read.table("PQLseq2_EtOH6m_methyl_vs_unmod.txt", header = TRUE)
#' dmcs_anno <- annotate(dmcs, level = "DMC", gtf_file = "my.gtf")
#'
#' combp_res <- run_combp("PQLseq2_EtOH6m_methyl_vs_unmod.txt", nCores = 1)
#' dmrs <- as.data.frame(combp_res$regions)
#' dmrs_anno <- annotate(dmrs, level = "DMR", gtf_file = "my.gtf")
#' }
#'
#' @export
annotate <- function(
    x,
    level = c("DMC", "DMR"),
    txdb = NULL,
    gtf_file = NULL,
    chr_col = "chr",
    start_col = "start",
    end_col = "end",
    promoter_upstream = 2000,
    promoter_downstream = 0
) {
  level <- match.arg(level)

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required for annotate().")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("Package 'IRanges' is required for annotate().")
  }
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop("Package 'GenomicFeatures' is required for annotate().")
  }

  df <- as.data.frame(x)

  if (!(chr_col %in% names(df)) || !(start_col %in% names(df))) {
    stop("Input must contain columns: ", chr_col, " and ", start_col)
  }

  if (level == "DMC") {
    if (!(end_col %in% names(df))) {
      df[[end_col]] <- as.numeric(df[[start_col]]) + 1
    }
  } else {
    if (!(end_col %in% names(df))) {
      stop("DMR annotation requires an end column (default name 'end').")
    }
  }

  if (is.null(txdb)) {
    if (is.null(gtf_file)) {
      stop("Provide either txdb= or gtf_file= to build a TxDb.")
    }
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file)
  }

  gr <- GenomicRanges::GRanges(
    seqnames = as.character(df[[chr_col]]),
    ranges = IRanges::IRanges(
      start = as.numeric(df[[start_col]]),
      end   = as.numeric(df[[end_col]])
    )
  )

  gene_gr <- GenomicFeatures::genes(txdb)
  exon_gr <- GenomicFeatures::exons(txdb)
  intron_gr <- unlist(GenomicFeatures::intronsByTranscript(txdb), use.names = FALSE)
  prom_gr <- GenomicFeatures::promoters(txdb, upstream = promoter_upstream, downstream = promoter_downstream)

  in_gene <- GenomicRanges::countOverlaps(gr, gene_gr) > 0
  in_exon <- GenomicRanges::countOverlaps(gr, exon_gr) > 0
  in_intron <- GenomicRanges::countOverlaps(gr, intron_gr) > 0
  in_prom <- GenomicRanges::countOverlaps(gr, prom_gr) > 0

  context <- ifelse(in_gene, "gene", "intergenic")
  context[in_prom] <- paste0(context[in_prom], ";promoter")
  context[in_exon] <- paste0(context[in_exon], ";exon")
  context[in_intron] <- paste0(context[in_intron], ";intron")

  nearest <- GenomicRanges::distanceToNearest(gr, gene_gr)
  nearest_gene_idx <- S4Vectors::subjectHits(nearest)
  nearest_dist <- S4Vectors::mcols(nearest)$distance

  nearest_gene_idx <- S4Vectors::subjectHits(nearest)
  nearest_gene_id[S4Vectors::queryHits(nearest)] <- names(gene_gr)[nearest_gene_idx]

  df$context <- context
  df$overlaps_gene <- in_gene
  df$overlaps_exon <- in_exon
  df$overlaps_intron <- in_intron
  df$overlaps_promoter <- in_prom
  df$nearest_gene_id <- nearest_gene_id
  df$distance_to_nearest_gene <- nearest_dist

  df
}
