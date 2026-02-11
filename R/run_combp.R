#' Run comb-p on converged PQLseq2 site-level results
#'
#' Formats site-level results (DMCs) into the input required by \code{ENmix::combp()}
#' and identifies spatially correlated differentially methylated regions (DMRs).
#'
#' On Windows, \code{nCores > 1} is not supported by \code{mclapply}; this function
#' automatically sets \code{nCores = 1} and warns.
#'
#' If \code{write_output = TRUE}, writes:
#' \itemize{
#'   \item \code{<prefix>.combp.bed} (region-level results)
#'   \item \code{<prefix>.combp.sites.bed} (site-level results returned by combp)
#'   \item \code{<prefix>.combp.input.bed} (formatted input to combp)
#' }
#'
#' @param pqlseq_results Either a filepath (tab-delimited) or a data.frame containing
#'   at least \code{chr}, \code{start}, and \code{pvalue}.
#' @param dist_cutoff Maximum distance (bp) allowed between adjacent sites when forming regions.
#' @param bin_size Bin size (bp) used for the Stouffer-Liptak-Kechris correction.
#' @param seed Seed p-value threshold for candidate regions.
#' @param nCores Number of cores. Windows: forced to 1.
#' @param write_output Logical; if TRUE, write output files.
#' @param outfile_prefix Prefix for output files.
#' @param require_converged Logical; if TRUE and \code{converged} exists, keep only TRUE.
#'   If \code{sigma2} exists, keep only \code{sigma2 > 0}.
#'
#' @return The object returned by \code{ENmix::combp()}, typically a list with
#'   \code{regions} and \code{data}.
#'
#' @examples
#' \dontrun{
#' combp_res <- run_combp(
#'   pqlseq_results = "PQLseq2_EtOH6m_methyl_vs_unmod.txt",
#'   dist_cutoff = 1000,
#'   bin_size = 50,
#'   seed = 0.05,
#'   nCores = 1,
#'   outfile_prefix = "EtOH6m_combP"
#' )
#' dmrs <- as.data.frame(combp_res$regions)
#' }
#'
#' @export
run_combp <- function(
    pqlseq_results,
    dist_cutoff = 500,
    bin_size = 50,
    seed = 0.05,
    nCores = parallel::detectCores(),
    write_output = TRUE,
    outfile_prefix = "combp_results",
    require_converged = TRUE
) {
  if (!requireNamespace("ENmix", quietly = TRUE)) {
    stop("Package 'ENmix' is required for run_combp(). Install with BiocManager::install('ENmix').")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for run_combp().")
  }

  if (.Platform$OS.type == "windows" && nCores > 1) {
    warning("Windows does not support mclapply; setting nCores = 1 for combp.")
    nCores <- 1
  }

  if (is.character(pqlseq_results) && length(pqlseq_results) == 1) {
    message("Reading PQLseq results...")
    df0 <- as.data.frame(data.table::fread(pqlseq_results))
  } else {
    df0 <- as.data.frame(pqlseq_results)
  }

  req <- c("chr", "start", "pvalue")
  miss <- setdiff(req, names(df0))
  if (length(miss) > 0) {
    stop("PQLseq results missing required columns: ", paste(miss, collapse = ", "))
  }

  df <- df0
  if (require_converged && "converged" %in% names(df)) {
    df <- df[df$converged == TRUE, , drop = FALSE]
  }
  if ("sigma2" %in% names(df)) {
    df <- df[df$sigma2 > 0, , drop = FALSE]
  }

  df$chr <- as.character(df$chr)
  df$start <- as.numeric(as.character(df$start))
  df$pvalue <- as.numeric(as.character(df$pvalue))

  df$end <- df$start + 1
  bed <- df[, c("chr", "start", "end", "pvalue")]
  bed <- bed[order(bed$chr, bed$start), , drop = FALSE]
  bed$Probe <- paste0("chr", bed$chr, ":", bed$start, "-", bed$end)

  message("Running combp...")
  combp_res <- ENmix::combp(
    bed,
    dist.cutoff = dist_cutoff,
    bin.size = bin_size,
    seed = seed,
    region_plot = FALSE,
    mht_plot = FALSE,
    nCores = nCores,
    verbose = TRUE
  )

  if (write_output) {
    message("Writing combp output...")

    regions_df <- tryCatch(as.data.frame(combp_res$regions), error = function(e) data.frame())
    sites_df <- tryCatch(as.data.frame(combp_res$data), error = function(e) data.frame())

    data.table::fwrite(regions_df, paste0(outfile_prefix, ".combp.bed"), sep = "\t", quote = FALSE)
    data.table::fwrite(sites_df, paste0(outfile_prefix, ".combp.sites.bed"), sep = "\t", quote = FALSE)
    data.table::fwrite(bed, paste0(outfile_prefix, ".combp.input.bed"), sep = "\t", quote = FALSE)
  }

  message("combp finished.")
  combp_res
}

