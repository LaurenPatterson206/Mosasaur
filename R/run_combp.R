#' Identify DMRs with ENmix::combp() from PQLseq2 site-level results
#'
#' ENmix::combp() prints results and writes region-level DMRs to a file named
#' "resu_combp.csv" in the current working directory, but it does not return an
#' object (it returns NULL). This wrapper runs combp in a controlled output
#' directory, then reads "resu_combp.csv" back into R and returns it.
#'
#' @param pqlseq_results Either:
#'   (1) a data.frame containing columns chr, start, pvalue (and optionally converged), or
#'   (2) a path to a tab-delimited file readable by data.table::fread().
#' @param prefix Filename prefix for outputs (default "combp").
#' @param dist_cutoff Maximum distance (bp) to merge nearby significant sites into regions.
#' @param bin_size Bin size (bp) used for autocorrelation calculation.
#' @param seed FDR threshold used by combp to seed candidate regions (typical 0.05 or 0.01).
#' @param keep_converged_only If TRUE (default), and a converged column exists, keep only converged == TRUE.
#' @param n_cores Number of cores (Windows is forced to 1 because combp uses mclapply).
#'
#' @return A list with:
#'   \itemize{
#'     \item regions: data.frame of DMRs read from resu_combp.csv (may be empty)
#'     \item input_bed: the combp input bed used (chr, start, end, p, probe)
#'     \item out_dir: the directory containing outputs
#'     \item regions_file: the written regions csv path
#'   }
#'
#' @examples
#' \dontrun{
#' combp_out <- run_combp(res, prefix="EtOH6m_ModVsUnmod", dist_cutoff=1000, bin_size=50, seed=0.05)
#' combp_out$regions
#' }
#'
#' @export
run_combp <- function(
    pqlseq_results,
    prefix = "combp",
    dist_cutoff = 1000,
    bin_size = 50,
    seed = 0.05,
    keep_converged_only = TRUE,
    n_cores = 1
) {
  if (!requireNamespace("ENmix", quietly = TRUE)) {
    stop("Package 'ENmix' is required. Install with BiocManager::install('ENmix').")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }

  # Windows safety: ENmix combp uses mclapply internally
  if (.Platform$OS.type == "windows" && n_cores > 1) {
    warning("Windows detected: setting n_cores = 1 (combp uses mclapply).")
    n_cores <- 1
  }

  # Read input
  df0 <- NULL
  if (is.character(pqlseq_results) && length(pqlseq_results) == 1) {
    df0 <- as.data.frame(data.table::fread(pqlseq_results))
  } else {
    df0 <- as.data.frame(pqlseq_results)
  }

  req <- c("chr", "start", "pvalue")
  miss <- setdiff(req, names(df0))
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "),
         "\nNeed at least: chr, start, pvalue.")
  }

  df <- df0


  if (isTRUE(keep_converged_only) && "converged" %in% names(df)) {
    df <- df[df$converged == TRUE, , drop = FALSE]
  }


  df$chr <- as.character(df$chr)
  df$start <- suppressWarnings(as.numeric(df$start))
  df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
  df <- df[!is.na(df$start) & !is.na(df$pvalue), , drop = FALSE]

  if (nrow(df) < 2) {
    stop("Too few sites after filtering to run combp (need >= 2).")
  }

  # ENmix wants columns: chr, start, end, p, probe
  bed <- data.frame(
    chr   = df$chr,
    start = df$start,
    end   = df$start,
    p     = df$pvalue,
    probe = paste0(df$chr, ":", df$start),
    stringsAsFactors = FALSE
  )
  bed <- bed[order(bed$chr, bed$start), , drop = FALSE]

  message("Running combp on ", nrow(bed), " sites...")


  out_dir <- file.path(getwd(), paste0(prefix, "_combp"))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)


  if (file.exists("resu_combp.csv")) file.remove("resu_combp.csv")


  ENmix::combp(
    data = bed,
    dist.cutoff = dist_cutoff,
    bin.size = bin_size,
    seed = seed,
    region_plot = FALSE,
    mht_plot = FALSE,
    nCores = n_cores,
    verbose = TRUE
  )

  regions_file <- file.path(out_dir, "resu_combp.csv")
  if (!file.exists(regions_file)) {
    stop("combp finished but did not create resu_combp.csv in: ", out_dir,
         "\nThis usually indicates an ENmix/working-directory issue.")
  }

  regions <- tryCatch(
    utils::read.csv(regions_file, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )


  input_file <- file.path(out_dir, paste0(prefix, ".combp.input.bed"))
  data.table::fwrite(bed, input_file, sep = "\t")


  friendly_regions_file <- file.path(out_dir, paste0(prefix, ".combp.regions.csv"))
  file.copy(regions_file, friendly_regions_file, overwrite = TRUE)

  message("combp done. Regions file: ", friendly_regions_file)

  list(
    regions = regions,
    input_bed = bed,
    out_dir = out_dir,
    regions_file = friendly_regions_file
  )
}
