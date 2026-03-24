#' Annotate DMCs or DMRs using TxDb overlaps
#'
#' Adds overlap-based annotations (gene/exon/intron/promoter) and nearest-gene
#' information to site-level (DMC) or region-level (DMR) results.
#'
#' @param x data.frame with coordinates. DMC: chr/start (end optional). DMR: chr/start/end.
#' @param level "DMC" or "DMR".
#' @param txdb Optional TxDb object. If provided, gtf_file is ignored.
#' @param gtf_file Path to GTF/GFF to build a TxDb (uses txdbmaker::makeTxDbFromGFF).
#' @param chr_col,start_col,end_col Column names for coordinates.
#' @param promoter_upstream,promoter_downstream Promoter definition.
#' @param map_gene_names Logical; if TRUE, attempts to map nearest gene IDs to gene symbols.
#' @param biomart_dataset biomaRt dataset name. Example: "mmulatta_gene_ensembl".
#' @param biomart_id_col ID type expected in TxDb gene IDs for biomaRt mapping.
#'   Common values: "entrezgene_id", "ensembl_gene_id".
#'
#' @return data.frame with added annotation columns.
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
    promoter_downstream = 0,
    map_gene_names = TRUE,
    biomart_dataset = "mmulatta_gene_ensembl",
    biomart_id_col = "entrezgene_id"
) {
  level <- match.arg(level)

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) stop("Package 'GenomicRanges' is required.")
  if (!requireNamespace("IRanges", quietly = TRUE)) stop("Package 'IRanges' is required.")
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) stop("Package 'GenomicFeatures' is required.")
  if (!requireNamespace("S4Vectors", quietly = TRUE)) stop("Package 'S4Vectors' is required.")
  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) stop("Package 'GenomeInfoDb' is required.")
  if (is.null(txdb) && is.null(gtf_file)) stop("Provide either txdb= or gtf_file=.")

  df <- as.data.frame(x)

  if (!(chr_col %in% names(df)) || !(start_col %in% names(df))) {
    stop("Input must contain columns: ", chr_col, " and ", start_col)
  }

  if (level == "DMC") {
    if (!(end_col %in% names(df))) df[[end_col]] <- suppressWarnings(as.integer(df[[start_col]]) + 1L)
  } else {
    if (!(end_col %in% names(df))) stop("DMR annotation requires an end column (default 'end').")
  }

  # Build TxDb
  if (is.null(txdb)) {
    if (!requireNamespace("txdbmaker", quietly = TRUE)) {
      stop("Package 'txdbmaker' is required to build TxDb from GTF/GFF.\n",
           "Install with BiocManager::install('txdbmaker').")
    }
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_file)
  }

  #Input GRanges
  gr <- GenomicRanges::GRanges(
    seqnames = as.character(df[[chr_col]]),
    ranges = IRanges::IRanges(
      start = suppressWarnings(as.integer(df[[start_col]])),
      end   = suppressWarnings(as.integer(df[[end_col]]))
    )
  )

  tx_seqs <- GenomeInfoDb::seqlevels(txdb)
  if (length(tx_seqs) > 0) {
    tx_has_chr <- any(grepl("^chr", tx_seqs, ignore.case = TRUE))
    gr_has_chr <- any(grepl("^chr", as.character(GenomicRanges::seqnames(gr)), ignore.case = TRUE))

    if (tx_has_chr && !gr_has_chr) {
      GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    } else if (!tx_has_chr && gr_has_chr) {
      GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
    }

    gr <- GenomeInfoDb::keepSeqlevels(
      gr,
      intersect(GenomeInfoDb::seqlevels(gr), tx_seqs),
      pruning.mode = "coarse"
    )
  }

  if (length(gr) == 0) stop("No ranges remain after chromosome harmonization. Check chr naming vs TxDb.")

  gene_gr   <- GenomicFeatures::genes(txdb)
  exon_gr   <- GenomicFeatures::exons(txdb)
  intron_gr <- unlist(GenomicFeatures::intronsByTranscript(txdb), use.names = FALSE)
  prom_gr   <- GenomicFeatures::promoters(txdb, upstream = promoter_upstream, downstream = promoter_downstream)

  in_gene   <- GenomicRanges::countOverlaps(gr, gene_gr) > 0
  in_exon   <- GenomicRanges::countOverlaps(gr, exon_gr) > 0
  in_intron <- GenomicRanges::countOverlaps(gr, intron_gr) > 0
  in_prom   <- GenomicRanges::countOverlaps(gr, prom_gr) > 0

  context <- ifelse(in_gene, "gene", "intergenic")
  context[in_prom]   <- paste0(context[in_prom],   ";promoter")
  context[in_exon]   <- paste0(context[in_exon],   ";exon")
  context[in_intron] <- paste0(context[in_intron], ";intron")

  #Nearest gene mapped back to all rows
  nearest_hits <- GenomicRanges::distanceToNearest(gr, gene_gr)
  qh <- S4Vectors::queryHits(nearest_hits)
  sh <- S4Vectors::subjectHits(nearest_hits)
  dist_vec <- S4Vectors::mcols(nearest_hits)$distance

  nearest_gene_id <- rep(NA_character_, length(gr))
  distance_to_nearest_gene <- rep(NA_real_, length(gr))
  nearest_gene_name <- rep(NA_character_, length(gr))

  gene_ids <- names(gene_gr)
  if (is.null(gene_ids) || all(is.na(gene_ids)) || all(gene_ids == "")) {
    gene_ids <- as.character(seq_along(gene_gr))
  }

  nearest_gene_id[qh] <- gene_ids[sh]
  distance_to_nearest_gene[qh] <- dist_vec

  #Map gene IDs to gene symbols with biomaRt
  if (map_gene_names) {
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      warning("Package 'biomaRt' not installed; nearest_gene_name will be NA.")
    } else {
      unique_ids <- unique(stats::na.omit(nearest_gene_id))

      if (length(unique_ids) > 0) {
        mart_df <- tryCatch({
          mart <- biomaRt::useEnsembl(
            biomart = "genes",
            dataset = biomart_dataset
          )

          biomaRt::getBM(
            attributes = c(biomart_id_col, "external_gene_name"),
            filters = biomart_id_col,
            values = unique_ids,
            mart = mart
          )
        }, error = function(e) {
          warning("biomaRt mapping failed: ", conditionMessage(e))
          NULL
        })

        if (!is.null(mart_df) && nrow(mart_df) > 0) {
          mart_df <- mart_df[mart_df$external_gene_name != "" & !is.na(mart_df$external_gene_name), , drop = FALSE]
          id_to_name <- stats::setNames(mart_df$external_gene_name, as.character(mart_df[[biomart_id_col]]))
          nearest_gene_name <- unname(id_to_name[nearest_gene_id])
        }
      }
    }
  }

  #If mapping failed, fall back to gene ID
  nearest_gene_name[is.na(nearest_gene_name) | nearest_gene_name == ""] <-
    nearest_gene_id[is.na(nearest_gene_name) | nearest_gene_name == ""]

  #Attach annotations
  df2 <- df
  key_df <- paste(df2[[chr_col]], df2[[start_col]], df2[[end_col]], sep = "|")
  key_gr <- paste(
    as.character(GenomicRanges::seqnames(gr)),
    IRanges::start(gr),
    IRanges::end(gr),
    sep = "|"
  )
  df2 <- df2[match(key_gr, key_df), , drop = FALSE]

  df2$context <- context
  df2$overlaps_gene <- in_gene
  df2$overlaps_exon <- in_exon
  df2$overlaps_intron <- in_intron
  df2$overlaps_promoter <- in_prom
  df2$nearest_gene_id <- nearest_gene_id
  df2$nearest_gene_name <- nearest_gene_name
  df2$distance_to_nearest_gene <- distance_to_nearest_gene

  df2
}
