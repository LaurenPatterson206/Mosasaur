#' Run PQLseq2
#'
#' Runs \code{PQLseq2::pqlseq2()} for a binomial mixed model using count matrices
#' from \code{\link{prepare_data}} and a phenotype CSV.
#'
#' Designed to be stable for small cohorts by automatically filtering
#' problematic sites and preventing known PQLseq2 crashes.
#'
#' @param pheno_file Path to a CSV phenotype file.
#' @param prep List output from \code{\link{prepare_data}} containing MCounts and TCounts.
#' @param predictor_var Phenotype column name for predictor.
#' @param covariates Character vector of covariate column names.
#' @param subject_id_var Phenotype column with sample IDs matching count columns.
#' @param output_file Output path for results.
#' @param keep_only_converged Logical; if TRUE (default) return only converged sites.
#' @param kinship_file Optional path to a tab-delimited square kinship matrix.
#'   If NULL (default), an identity matrix is used.
#'
#' @return A data.frame of results (also written to disk).
#' @export
run_pqlseq <- function(
    pheno_file,
    prep,
    predictor_var,
    covariates,
    subject_id_var,
    output_file = "PQLseq2_Results.txt",
    keep_only_converged = TRUE,
    kinship_file = NULL
) {

  if (!requireNamespace("PQLseq2", quietly = TRUE))
    stop("Package 'PQLseq2' is required.")

  pheno <- utils::read.csv(pheno_file, stringsAsFactors = FALSE)

  if (!subject_id_var %in% names(pheno))
    stop("subject_id_var not found in phenotype.")

  if (!predictor_var %in% names(pheno))
    stop("predictor_var not found in phenotype.")

  missing_cov <- setdiff(covariates, names(pheno))
  if (length(missing_cov) > 0)
    stop("Missing covariates: ", paste(missing_cov, collapse = ", "))

  pheno[[subject_id_var]] <- as.character(pheno[[subject_id_var]])

  for (cv in covariates) {
    if (is.numeric(pheno[[cv]]) || is.integer(pheno[[cv]])) {
      u <- unique(pheno[[cv]][!is.na(pheno[[cv]])])
      if (length(u) > 1 && length(u) <= 10) {
        pheno[[cv]] <- as.factor(pheno[[cv]])
        message("Treating numeric covariate as factor: ", cv,
                " (", length(u), " levels: ", paste(sort(u), collapse = ", "), ")")
      }
    }
    if (is.character(pheno[[cv]]))
      pheno[[cv]] <- as.factor(pheno[[cv]])
  }

  if (is.null(prep$MCounts) || is.null(prep$TCounts))
    stop("prep must contain $MCounts and $TCounts.")

  MC_all <- as.data.frame(prep$MCounts, check.names = FALSE)
  TC_all <- as.data.frame(prep$TCounts, check.names = FALSE)

  if (!all(c("chr", "start") %in% names(MC_all)))
    stop("prep$MCounts must include 'chr' and 'start'.")

  if (!all(c("chr", "start") %in% names(TC_all)))
    stop("prep$TCounts must include 'chr' and 'start'.")

  coords <- MC_all[, c("chr", "start"), drop = FALSE]
  MC <- MC_all[, setdiff(names(MC_all), c("chr", "start")), drop = FALSE]
  TC <- TC_all[, setdiff(names(TC_all), c("chr", "start")), drop = FALSE]

  if (!identical(colnames(MC), colnames(TC)))
    stop("MCounts and TCounts sample columns do not match.")

  keep <- pheno[[subject_id_var]] %in% colnames(MC)

  if (!all(keep)) {
    warning("Dropping phenotype samples without counts: ",
            paste(pheno[[subject_id_var]][!keep], collapse = ", "))
  }

  pheno <- pheno[keep, , drop = FALSE]
  sample_ids <- pheno[[subject_id_var]]

  MC <- MC[, sample_ids, drop = FALSE]
  TC <- TC[, sample_ids, drop = FALSE]

  x <- as.numeric(pheno[[predictor_var]])
  if (anyNA(x))
    stop("Predictor contains NA after coercion.")
  x <- matrix(x, ncol = 1)

  W <- stats::model.matrix(
    stats::as.formula(paste("~", paste(covariates, collapse = "+"))),
    data = pheno
  )

  n <- nrow(pheno)

  if (is.null(kinship_file)) {
    K <- diag(n)
  } else {
    kin_df <- utils::read.table(kinship_file, header = TRUE, check.names = FALSE, sep = "\t")
    kin_mat <- as.matrix(kin_df)

    if (!is.numeric(kin_mat[1,1]) && ncol(kin_mat) > 1) {
      rn <- kin_mat[,1]
      kin_mat <- kin_mat[,-1, drop = FALSE]
      rownames(kin_mat) <- rn
    }

    storage.mode(kin_mat) <- "numeric"

    if (!is.null(rownames(kin_mat)) && !is.null(colnames(kin_mat)) &&
        all(sample_ids %in% rownames(kin_mat)) && all(sample_ids %in% colnames(kin_mat))) {
      kin_mat <- kin_mat[sample_ids, sample_ids, drop = FALSE]
    }

    if (!all(dim(kin_mat) == c(n, n)))
      stop("Kinship matrix must be ", n, " x ", n,
           " after alignment. Got: ", paste(dim(kin_mat), collapse = " x "))

    K <- kin_mat
  }

  Y <- as.matrix(MC)
  storage.mode(Y) <- "numeric"

  Tmat <- as.matrix(TC)
  storage.mode(Tmat) <- "numeric"

  if (any(Tmat < Y, na.rm = TRUE))
    stop("Some trials < successes.")

  nonboundary <- (Tmat > 0) & (Y > 0) & (Y < Tmat)
  keep_sites <- rowSums(nonboundary) >= max(6, ceiling(ncol(Y) / 2))

  Y <- Y[keep_sites, , drop = FALSE]
  Tmat <- Tmat[keep_sites, , drop = FALSE]
  coords <- coords[keep_sites, , drop = FALSE]

  message("Stability filter kept ", nrow(Y), " sites.")
  if (nrow(Y) == 0)
    stop("No sites remain after stability filtering.")


  pqlseq_warnings <- character(0)
  shown_warnings <- 0L

  fit <- withCallingHandlers(
    PQLseq2::pqlseq2(
      Y = Y,
      x = x,
      K = K,
      W = W,
      lib_size = Tmat,
      model = "BMM",
      ncores = 1,
      filter = FALSE,
      verbose = FALSE
    ),
    warning = function(w) {
      msg <- conditionMessage(w)
      pqlseq_warnings <<- c(pqlseq_warnings, msg)

      if (shown_warnings < 50L) {
        shown_warnings <<- shown_warnings + 1L
        message("Warning [", shown_warnings, "]: ", msg)
      } else if (shown_warnings == 50L) {
        shown_warnings <<- shown_warnings + 1L
        message("Additional warnings suppressed. Type pqlseq_error to view all warnings.")
      }

      invokeRestart("muffleWarning")
    }
  )

  assign("pqlseq_error", pqlseq_warnings, envir = parent.frame())

  if (length(pqlseq_warnings) > 0) {
    message("Captured ", length(pqlseq_warnings),
            " total warnings. Type pqlseq_error to view all warnings.")
  }


  res <- cbind(coords, as.data.frame(fit, stringsAsFactors = FALSE))
  rownames(res) <- NULL

  if (keep_only_converged && "converged" %in% names(res)) {
    res$converged <- as.logical(res$converged)
    res <- res[res$converged == TRUE, , drop = FALSE]
    rownames(res) <- NULL
  }

  if ("pvalue" %in% names(res)) {
    res$pvalue <- suppressWarnings(as.numeric(res$pvalue))
    res <- res[order(res$pvalue), , drop = FALSE]
    rownames(res) <- NULL
  }

  utils::write.table(res, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  res
}
