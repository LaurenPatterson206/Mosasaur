#' Run PQLseq2 on prepared count matrices (beta-binomial mixed model)
#'
#' Fits a beta-binomial mixed model at each CpG site using \code{PQLseq2::pqlseq2()},
#' supporting arbitrary predictors and covariates with an optional kinship matrix.
#' If no kinship matrix is provided, an identity matrix is used. Results are
#' filtered to keep only \code{converged == TRUE} (if present).
#'
#' This function can take either:
#' \itemize{
#'   \item \code{prep}: the list output of \code{\link{prepare_data}} (preferred)
#'   \item or file paths \code{mcounts_file} and \code{tcounts_file}
#' }
#'
#' @param pheno_file Path to a CSV phenotype file.
#' @param prep Optional list output from \code{\link{prepare_data}} containing
#'   \code{MCounts} and \code{TCounts}. If provided, \code{mcounts_file} and
#'   \code{tcounts_file} are ignored.
#' @param kinship_file Optional path to a kinship matrix file (tab-delimited, square).
#'   If \code{NULL}, an identity matrix is used.
#' @param mcounts_file Path to an MCounts file (ignored if \code{prep} provided).
#' @param tcounts_file Path to a TCounts file (ignored if \code{prep} provided).
#' @param output_file Path to write tab-delimited PQLseq2 results.
#' @param predictor_var Column name in phenotype file used as predictor (W).
#'   Must be numeric or coercible to numeric without introducing NAs.
#' @param covariates Character vector of covariate column names to include in the
#'   design matrix (x). The design matrix includes an intercept by default.
#' @param subject_id_var Column name in phenotype file containing sample IDs that
#'   match count matrix column names.
#' @param n_cores Number of CPU cores passed to \code{PQLseq2::pqlseq2()} (\code{ncores}).
#' @param keep_only_converged Logical; if TRUE and a \code{converged} column exists,
#'   keep only \code{converged == TRUE}.
#' @param filter_in_pqlseq2 Logical; passed to \code{PQLseq2::pqlseq2(filter=...)}.
#' @param check_K Logical; passed to \code{PQLseq2::pqlseq2(check_K=...)}.
#' @param verbose Logical; passed to \code{PQLseq2::pqlseq2(verbose=...)}.
#'
#' @return A data.frame containing \code{chr}, \code{start}, and PQLseq2 output columns
#'   (optionally converged-only), sorted by ascending p-value when present.
#'
#' @examples
#' \dontrun{
#' res <- run_pqlseq(
#'   pheno_file = "pheno.filtered.csv",
#'   prep = prep,
#'   kinship_file = NULL,
#'   predictor_var = "EtOH_6m",
#'   covariates = c("Cohort","Age"),
#'   subject_id_var = "MATRR_ID",
#'   n_cores = 8,
#'   output_file = "PQLseq2_results.txt"
#' )
#' }
#'
#' @export
run_pqlseq <- function(
    pheno_file,
    prep = NULL,
    kinship_file = NULL,
    mcounts_file = "MCounts.txt",
    tcounts_file = "TCounts.txt",
    output_file = "PQLseq2_Results.txt",
    predictor_var,
    covariates = c("Cohort", "Age"),
    subject_id_var = "subject_id",
    n_cores = 4,
    keep_only_converged = TRUE,
    filter_in_pqlseq2 = TRUE,
    check_K = FALSE,
    verbose = FALSE
) {
  #dependencies
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for run_pqlseq().")
  }
  if (!requireNamespace("PQLseq2", quietly = TRUE)) {
    stop("Package 'PQLseq2' is required for run_pqlseq().")
  }

  #read phenotype
  pheno <- utils::read.csv(pheno_file, header = TRUE)
  if (!(subject_id_var %in% names(pheno))) {
    stop("subject_id_var not found in phenotype file: ", subject_id_var)
  }
  if (!(predictor_var %in% names(pheno))) {
    stop("predictor_var not found in phenotype file: ", predictor_var)
  }
  pheno[[subject_id_var]] <- as.character(pheno[[subject_id_var]])

  #load counts from prep or files
  if (!is.null(prep)) {
    MC_all <- as.data.frame(prep$MCounts)
    TC_all <- as.data.frame(prep$TCounts)
  } else {
    MC_all <- as.data.frame(data.table::fread(mcounts_file))
    TC_all <- as.data.frame(data.table::fread(tcounts_file))
  }

  #expect chr/start columns in MCounts
  if (!all(c("chr", "start") %in% names(MC_all))) {
    stop("MCounts must contain columns 'chr' and 'start'.")
  }
  coords <- MC_all[, c("chr", "start")]
  MC <- MC_all[, !(names(MC_all) %in% c("chr", "start")), drop = FALSE]
  TC <- TC_all[, !(names(TC_all) %in% c("chr", "start")), drop = FALSE]

  #align phenotype samples to count columns
  sample_order <- pheno[[subject_id_var]]
  valid <- sample_order %in% colnames(MC)
  if (!all(valid)) {
    warning("Dropping phenotype samples without counts: ", paste(sample_order[!valid], collapse = ", "))
  }
  pheno <- pheno[valid, , drop = FALSE]
  sample_order <- pheno[[subject_id_var]]

  MC <- MC[, sample_order, drop = FALSE]
  TC <- TC[, sample_order, drop = FALSE]


  W_vec <- pheno[[predictor_var]]
  if (is.factor(W_vec)) W_vec <- as.character(W_vec)
  W_vec <- as.numeric(W_vec)
  if (anyNA(W_vec)) {
    stop("Predictor could not be coerced to numeric without NA. Check predictor_var: ", predictor_var)
  }
  W <- matrix(W_vec, ncol = 1)

  #make sure covariates exist
  missing_cov <- setdiff(covariates, names(pheno))
  if (length(missing_cov) > 0) {
    stop("These covariate columns are missing in phenotype: ", paste(missing_cov, collapse = ", "))
  }

  for (cv in covariates) {
    if (is.character(pheno[[cv]])) pheno[[cv]] <- as.factor(pheno[[cv]])
  }

  #build design matrix (x) with intercept
  x_formula <- stats::as.formula(paste("~", paste(covariates, collapse = " + ")))
  x <- stats::model.matrix(x_formula, data = pheno)

  #kinship matrix (K)
  if (is.null(kinship_file)) {
    kin <- diag(nrow(pheno))
  } else {
    kin <- utils::read.table(kinship_file, header = TRUE, check.names = FALSE)
    kin <- as.matrix(kin)
  }

  #convert counts to numeric matrices
  Y <- as.matrix(MC)
  lib_size <- as.matrix(TC)

  if (nrow(x) != nrow(W)) stop("Design matrix x and W have mismatched rows.")
  if (nrow(kin) != nrow(W) || ncol(kin) != nrow(W)) stop("Kinship matrix K dimensions do not match samples.")
  if (ncol(Y) != nrow(W)) stop("Y columns do not match number of samples.")
  if (ncol(lib_size) != nrow(W)) stop("lib_size columns do not match number of samples.")

  tmp_out <- tempfile(pattern = "pqlseq2_", fileext = ".txt")

  #fir model
  fit <- PQLseq2::pqlseq2(
    Y = Y,
    x = x,
    K = kin,
    W = W,
    lib_size = lib_size,
    model = "BMM",
    ncores = n_cores,
    filter = filter_in_pqlseq2,
    check_K = check_K,
    outfile = tmp_out,
    verbose = verbose
  )

  if (file.exists(tmp_out)) unlink(tmp_out)

  fit <- as.data.frame(fit)

  #attach coords
  df <- cbind(coords, fit)

  #keep converged only
  if (keep_only_converged && "converged" %in% names(df)) {
    df <- df[df$converged == TRUE, , drop = FALSE]
  }

  #sort by pvalue if present
  if ("pvalue" %in% names(df)) {
    df$pvalue <- as.numeric(as.character(df$pvalue))
    df <- df[order(df$pvalue), , drop = FALSE]
  }

  #write final output
  utils::write.table(df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  df
}

