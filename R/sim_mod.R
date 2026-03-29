#' Dirichlet-multinomial simulation of modification counts
#'
#' Simulates per-site, per-sample counts across modification categories using a
#' Dirichlet-multinomial model. An effect can be applied to a subset of the CpG
#' sites in cases by shifting the Dirichlet parameters.
#'
#' @param n_controls Number of control samples.
#' @param n_cases Number of case samples.
#' @param n_cpg Number of CpG sites.
#' @param m Total count (or coverage) per CpG site per sample.
#' @param alpha Numeric vector of length \code{d} (Dirichlet parameters) for controls.
#' @param delta Numeric vector of length \code{d} (effect shift added to \code{alpha} in cases).
#' @param d Number of categories (defaults to 3: methylated, unmodified, hydroxymethylated).
#' @param effect_prop Proportion of CpG sites that are truly affected in cases.
#'   Defaults to 1 (all sites affected).
#' @param seed Optional integer for reproducibility.
#'
#' @return A data.frame with columns \code{methylated}, \code{unmodified},
#'   \code{hydroxymethylated}, \code{sample_id}, \code{group}, \code{cpg_site},
#'   and \code{is_effect}.
#'
#' @importFrom MGLM rdirmn
#' @exportcheck()
sim_mod <- function(
    n_controls,
    n_cases,
    n_cpg,
    m,
    alpha,
    delta,
    d = 3,
    effect_prop = 1,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  if (length(alpha) != d || length(delta) != d) stop("alpha and delta must have length d.")
  if (effect_prop <= 0 || effect_prop > 1) stop("effect_prop must be in (0, 1].")

  n_total <- n_controls + n_cases
  sample_ids <- paste0("S", seq_len(n_total))
  group_labels <- factor(c(rep("control", n_controls), rep("case", n_cases)))

  n_eff <- max(1, round(n_cpg * effect_prop))
  effect_sites <- paste0("CpG", sort(sample(seq_len(n_cpg), n_eff)))
  cpg_ids <- paste0("CpG", seq_len(n_cpg))

  all_data <- vector("list", n_total)

  for (i in seq_len(n_total)) {
    group <- group_labels[i]
    a_use <- alpha

    out_i <- vector("list", n_cpg)
    for (j in seq_len(n_cpg)) {
      cpg <- cpg_ids[j]
      if (group == "case" && (cpg %in% effect_sites)) {
        a_use <- alpha + delta
      } else {
        a_use <- alpha
      }

      counts <- MGLM::rdirmn(1, m, a_use)
      df <- as.data.frame(counts)
      colnames(df) <- c("methylated", "unmodified", "hydroxymethylated")
      df$cpg_site <- cpg
      df$is_effect <- cpg %in% effect_sites
      out_i[[j]] <- df
    }

    out_i <- do.call(rbind, out_i)
    out_i$sample_id <- sample_ids[i]
    out_i$group <- group
    all_data[[i]] <- out_i
  }

  result <- do.call(rbind, all_data)
  rownames(result) <- NULL
  result
}
