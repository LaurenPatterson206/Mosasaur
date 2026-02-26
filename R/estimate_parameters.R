#' Estimate Dirichlet alpha parameters from raw bedMethyl counts (fast streaming)
#'
#' Estimates Dirichlet parameters for three categories (methylated,
#' unmodified, hydroxymethylated) from raw bedMethyl count data using
#' a method-of-moments approach based on empirical means and SDs of proportions
#'
#' @param path Directory containing bedMethyl files. Defaults to `"."`.
#' @param pattern Regex to match bedMethyl files. Defaults to `"\\\\.bed$"`.
#' @param keep_mod_code Value in column V4 to keep. Defaults to `"m"`.
#' @param pseudocount Numeric added to each count before forming proportions.
#'   Helps stabilize estimation when many zeros exist. Default `0`.
#' @param subsample Integer or NULL. If provided, randomly samples up to this
#'   many observations per file (after filtering `V4`). Default `NULL`.
#' @param seed Optional RNG seed used with `subsample`.
#' @param categories Names for the three returned alpha parameters.
#'
#' @return A list with:
#' \describe{
#'   \item{alpha}{Named numeric vector of length 3.}
#'   \item{alpha1,alpha2,alpha3}{Convenience scalars.}
#'   \item{p_mean}{Mean proportions for each category.}
#'   \item{p_sd}{SD of proportions for each category.}
#'   \item{n_obs}{Number of observations used.}
#'   \item{alpha0}{Estimated concentration (sum of alphas).}
#' }
#'
#' @export
Estimate_Parameters <- function(
    path = ".",
    pattern = "\\\\.bed$",
    keep_mod_code = "m",
    pseudocount = 0,
    subsample = NULL,
    seed = NULL,
    categories = c("methylated", "unmodified", "hydroxymethylated")
) {
  stopifnot(length(categories) == 3)

  bed.files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  if (length(bed.files) == 0) {
    stop("No bedMethyl files found. Check 'path' and/or 'pattern'.")
  }


  n <- 0L
  mean_vec <- c(0, 0, 0)
  M2_vec <- c(0, 0, 0)

  if (!is.null(seed)) set.seed(seed)

  for (f in bed.files) {

    dt <- data.table::fread(
      f,
      select = c(4, 12, 13, 14),
      showProgress = FALSE
    )
    if (ncol(dt) < 4) next

    dt <- dt[dt[[1]] == keep_mod_code, ]
    if (nrow(dt) == 0) next


    if (!is.null(subsample) && is.finite(subsample) && subsample > 0 &&
        nrow(dt) > subsample) {
      dt <- dt[sample.int(nrow(dt), subsample), ]
    }

    m <- as.numeric(dt[[2]])
    u <- as.numeric(dt[[3]])
    h <- as.numeric(dt[[4]])


    m[is.na(m)] <- 0
    u[is.na(u)] <- 0
    h[is.na(h)] <- 0

    if (pseudocount != 0) {
      m <- m + pseudocount
      u <- u + pseudocount
      h <- h + pseudocount
    }

    denom <- m + u + h
    keep <- denom > 0
    if (!any(keep)) next


    p1 <- m[keep] / denom[keep]
    p2 <- u[keep] / denom[keep]
    p3 <- h[keep] / denom[keep]


    for (i in seq_along(p1)) {
      n <- n + 1L
      x <- c(p1[i], p2[i], p3[i])
      delta <- x - mean_vec
      mean_vec <- mean_vec + delta / n
      delta2 <- x - mean_vec
      M2_vec <- M2_vec + delta * delta2
    }
  }

  if (n < 2) {
    stop("Not enough observations to estimate parameters (n_obs < 2).")
  }

  var_vec <- M2_vec / (n - 1L)
  sd_vec <- sqrt(var_vec)


  mu <- mean_vec
  alpha0_candidates <- (mu * (1 - mu) / var_vec) - 1
  alpha0_candidates <- alpha0_candidates[
    is.finite(alpha0_candidates) & alpha0_candidates > 0
  ]

  if (length(alpha0_candidates) == 0) {
    stop(
      "Cannot estimate Dirichlet parameters (variance too small/invalid). ",
      "Try setting pseudocount to 0.5 or 1."
    )
  }

  alpha0 <- stats::median(alpha0_candidates)

  alpha_hat <- alpha0 * mu
  if (any(!is.finite(alpha_hat)) || any(alpha_hat <= 0)) {
    stop("Estimated alpha parameters are invalid. Try pseudocount > 0.")
  }

  names(alpha_hat) <- categories
  names(mu) <- categories
  names(sd_vec) <- categories

  list(
    alpha = alpha_hat,
    alpha1 = unname(alpha_hat[1]),
    alpha2 = unname(alpha_hat[2]),
    alpha3 = unname(alpha_hat[3]),
    p_mean = mu,
    p_sd = sd_vec,
    n_obs = n,
    alpha0 = alpha0
  )
}
