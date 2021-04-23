#' Differential Expression Analysis on Mzroll List
#'
#' @inheritParams test_mzroll_list
#' @inheritParams plot_heatmap
#' @param test_model a RHS formula for regression
#' @param null_model if provided a null RHS formula to compare to
#'   \code{test_model} using the likelihood-ratio test.
#'
#' @return a tibble of significance tests for each feature
#'
#' @examples
#' library(dplyr)
#'
#' mzroll_list <- readRDS(
#'   authutils::get_clamr_assets("X0083-mzroll-list.Rds")
#' ) %>%
#'   floor_peaks(12)
#'
#' diffex_mzroll(mzroll_list, test_model = "mutation + sex")
#' @export
diffex_mzroll <- function(mzroll_list,
                          value.var = "centered_log2_abundance",
                          test_model,
                          null_model = NULL) {
  test_mzroll_list(mzroll_list)

  quant_vars <- colnames(mzroll_list$peaks)[
    purrr::map_chr(mzroll_list$peaks, class) == "numeric"
  ]
  valid_quant_vars <- setdiff(quant_vars, c("groupId", "sampleId"))
  if (!(value.var %in% valid_quant_vars)) {
    stop(
      value.var,
      " not present in \"mzroll_list$peaks\", valid numeric variables are:",
      paste(valid_quant_vars, collapse = ", ")
    )
  }

  # setup data

  nested_peaks <- mzroll_list$peaks %>%
    dplyr::left_join(mzroll_list$samples, "sampleId") %>%
    tidyr::nest(one_peak_data = -groupId)

  # setup and validate formulas
  viable_sample_fields <- setdiff(
    colnames(mzroll_list$samples),
    c("sampleId", "tube label", "MS ID string", "MS ID string alternative")
  )

  test_model_formula <- validate_formulas(
    test_model,
    value.var,
    viable_sample_fields
  )

  # flag peakgroups which cannot be fit due to too few samples

  n_conditions <- length(unique(mzroll_list$samples$`condition #`))

  underspecified_groups <- nested_peaks %>%
    dplyr::filter(purrr::map_dbl(one_peak_data, nrow) < n_conditions + 1)

  if (nrow(underspecified_groups) != 0) {
    warning(
      nrow(underspecified_groups),
      " peakgroups were discarded
      consider using floor_peaks() or fitting a simpler model"
    )
    nested_peaks <- nested_peaks %>%
      dplyr::anti_join(underspecified_groups, by = "groupId")
  }

  # carry-out feature-wise significance testing

  if ("NULL" %in% class(null_model)) {
    peakgroup_signif <- nested_peaks %>%
      dplyr::mutate(lm_fits = purrr::map(
        one_peak_data,
        diffex_lm,
        test_model_formula
      )) %>%
      dplyr::select(-one_peak_data) %>%
      tidyr::unnest(lm_fits)
  } else {
    null_model_formula <- validate_formulas(
      null_model,
      value.var,
      viable_sample_fields
    )

    peakgroup_signif <- nested_peaks %>%
      dplyr::mutate(anova_signif = purrr::map(
        one_peak_data,
        diffex_anova,
        test_model_formula,
        null_model_formula
      )) %>%
      dplyr::select(-one_peak_data) %>%
      tidyr::unnest(anova_signif)
  }

  # FDR control

  peakgroup_signif %>%
    tidyr::nest(term_data = -term) %>%
    dplyr::mutate(fdr_summary = purrr::map(term_data, diffex_fdr)) %>%
    dplyr::select(-term_data) %>%
    tidyr::unnest(fdr_summary)
}


validate_formulas <- function(model, value.var, viable_sample_fields) {
  if (stringr::str_detect(model, "~")) {
    stop("your regression model: ", model, " should not include an \"~\"")
  }

  # ensure that all covariates are in viable_sample_fields
  RHS_formula_vars <- stats::as.formula(paste("~", model)) %>%
    all.vars()

  inviable_RHS_vars <- setdiff(RHS_formula_vars, viable_sample_fields)
  if (length(inviable_RHS_vars) != 0) {
    stop(
      "some regression fields in your model (~ ",
      model,
      ") are not valid fields from the samples table: ",
      paste(inviable_RHS_vars, collapse = ", ")
    )
  }

  full_formula <- stats::as.formula(paste(value.var, "~", model))

  return(full_formula)
}

diffex_lm <- function(one_peak_data, test_model_formula) {
  stats::lm(
    formula = test_model_formula,
    data = one_peak_data
  ) %>%
    broom::tidy()
}

diffex_anova <- function(one_peak_data,
                         test_model_formula,
                         null_model_formula) {
  test_lm <- stats::lm(test_model_formula, one_peak_data)
  null_lm <- stats::lm(null_model_formula, one_peak_data)

  stats::anova(null_lm, test_lm) %>%
    broom::tidy() %>%
    dplyr::slice(-1) %>%
    dplyr::mutate(term = "LRT")
}

diffex_fdr <- function(term_data) {
  p_values <- term_data$p.value
  q_values <- try(qvalue::qvalue(p_values)$qvalues, silent = TRUE)

  if ("try-error" %in% class(q_values)) {
    # if qvalue fails this is probably because there are no p-values greater
    # than 0.95 (the highest lambda value)
    # if so add a single p-value of 1 to try to combat the problem
    q_values <- qvalue::qvalue(c(p_values, 1))$qvalues
    q_values <- q_values[-length(q_values)]

    warning(
      "q-value calculation initially failed but calicomics
        was able to recover results"
    )
  }

  term_data %>%
    dplyr::mutate(qvalue = q_values)
}

#' Volcano plot
#'
#' @param regression_significance output of diffex_mzroll
#' @param max_p_trans maximum value of -log10 pvalues to plot
#' @param FDR_cutoff FDR cutoff to label for significance
#'
#' @export
volcano_plot <- function(regression_significance,
                         max_p_trans = 10,
                         FDR_cutoff = 0.1) {
  checkmate::assertDataFrame(regression_significance)
  stopifnot("term" %in% colnames(regression_significance))

  effect_var <- dplyr::case_when(
    "estimate" %in% colnames(regression_significance) ~ "estimate",
    "sumsq" %in% colnames(regression_significance) ~ "sumsq",
    TRUE ~ NA_character_
  )

  if (is.na(effect_var)) {
    stop("volcano plot cannot be generated due to unknown test")
  }

  regression_significance %>%
    dplyr::filter(!is.na(p.value)) %>%
    dplyr::mutate(
      p.value.trans = trans_pvalues(p.value, max_p_trans = max_p_trans),
      is_discovery = qvalue < FDR_cutoff
    ) %>%
    ggplot(aes_string(x = effect_var)) +
    geom_point(aes(y = p.value.trans, color = is_discovery)) +
    facet_wrap(~term, scales = "free_x") +
    scale_y_continuous(expression(-log[10] ~ "pvalue")) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "RED")) +
    theme_bw() +
    ggtitle("Volcano plot")
}

trans_pvalues <- function(p, max_p_trans = 10) {
  pmin(-1 * log10(p), max_p_trans)
}

#' P-value Histogram
#'
#' @inheritParams volcano_plot
#'
#' @export
pvalue_histograms <- function(regression_significance) {
  stopifnot("data.frame" %in% class(regression_significance))
  stopifnot("term" %in% colnames(regression_significance))

  regression_significance %>%
    dplyr::filter(!is.na(p.value)) %>%
    ggplot(aes(x = p.value)) +
    geom_histogram(bins = 100) +
    facet_wrap(~term, scales = "free_y") +
    expand_limits(x = c(0, 1)) +
    theme_bw() +
    ggtitle("P-value histograms")
}
