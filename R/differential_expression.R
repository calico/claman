#' Differential Expression Analysis on Mzroll List
#'
#' @inheritParams test_mzroll_list
#' @param value_var measurement variable
#' @param test_model a RHS formula for regression
#' @param null_model if provided a null RHS formula to compare to
#'   \code{test_model} using the likelihood-ratio test.
#' @param additional_grouping_vars sample-, or measurement-level variables
#'   to groupby when performing regression in addition to groupId.
#' @param fdr_by_groupings if TRUE then calculate FDR separately across
#'  \code{additional_grouping_vars} (if they exist).
#'
#' @return a tibble of significance tests for each feature
#'
#' @examples
#' 
#' # standard feature-wise regression
#' diffex_mzroll(
#'   nplug_mzroll_normalized,
#'   "normalized_log2_abundance",
#'   "limitation + limitation:DR + 0"
#'   )
#' 
#' # separate regression for each feature and limitatioon
#' diffex_mzroll(
#'   nplug_mzroll_normalized,
#'   "normalized_log2_abundance",
#'   "DR + 0",
#'   additional_grouping_vars = "limitation"
#'   )
#'   
#' # feature-wise ANOVA
#' diffex_mzroll(
#'   nplug_mzroll_normalized,
#'   "normalized_log2_abundance",
#'   "limitation + 0",
#'   "0"
#'   )
#' @export
diffex_mzroll <- function(
  mzroll_list,
  value_var,
  test_model,
  null_model = NULL,
  additional_grouping_vars = NULL,
  fdr_by_groupings = FALSE
  ) {
  
  test_mzroll_list(mzroll_list)
  design_tbl <- romic::get_design_tbl(mzroll_list)

  quant_vars <- colnames(mzroll_list$measurements)[
    purrr::map_chr(mzroll_list$measurements, class) == "numeric"
  ]
  valid_quant_vars <- setdiff(quant_vars, c("groupId", "sampleId"))
  if (!(value_var %in% valid_quant_vars)) {
    stop(
      value_var,
      " is not a valid quant variable in \"mzroll_list$measurements\",
      - valid numeric variables are:
      - ",
      paste(valid_quant_vars, collapse = ", ")
    )
  }
  
  grouping_vars <- format_grouping_vars(
    mzroll_list,
    additional_grouping_vars
  )
  
  # setup and validate formulas
  viable_sample_fields <- setdiff(colnames(mzroll_list$samples), "sampleId")

  test_model_formula <- validate_formulas(
    test_model,
    value_var,
    viable_sample_fields
  )

  checkmate::assertLogical(fdr_by_groupings, len = 1)
  
  # setup data so each row in a table corresponds to a single regression
  
  nested_vars <- setdiff(unique(design_tbl$variable), grouping_vars)
  
  # warning is suppressed because I could only get this to work with the
  #   deprecated key argument
  nested_peaks <- suppressWarnings(
    mzroll_list %>%
      romic::triple_to_tidy() %>%
      {.$data} %>%
      tidyr::nest(!!!rlang::syms(nested_vars), .key = "one_peak_data")
  )
  
  # determine how many distinct conditions there are from a feature
  # with the most samples
  n_conditions <- ncol(stats::model.matrix(
    test_model_formula,
    nested_peaks %>%
      dplyr::mutate(n = purrr::map_int(one_peak_data, nrow)) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::slice(1) %>%
      tidyr::unnest(one_peak_data)
    ))
  
  # flag peakgroups which cannot be fit due to too few samples

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
      value_var,
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

  # nest by term and by additional_grouping_vars
  
  if (fdr_by_groupings) {
    nested_vars <- setdiff(
      colnames(peakgroup_signif),
      c("term", additional_grouping_vars)
    )
  } else {
    nested_vars <- setdiff(colnames(peakgroup_signif), "term")
  }
  
  output <- suppressWarnings(
    peakgroup_signif %>%
      tidyr::nest(!!!rlang::syms(nested_vars), .key = "term_data")
    ) %>%
    dplyr::mutate(fdr_summary = purrr::map(term_data, diffex_fdr)) %>%
    dplyr::select(-term_data) %>%
    tidyr::unnest(fdr_summary) %>%
    dplyr::mutate(
      signif_qual = dplyr::case_when(
        qvalue < 0.001 ~ " ***",
        qvalue < 0.01 ~ " **",
        qvalue < 0.1 ~ " *",
        TRUE ~ ""
      ),
      diffex_label = glue::glue("{round(statistic, 3)}{signif_qual}")
    ) %>%
    dplyr::select(-signif_qual)
  
  return(output)
}

format_grouping_vars <- function(
  mzroll_list,
  additional_grouping_vars = NULL
  ) {
  
  design_tbl <- romic::get_design_tbl(mzroll_list)
  
  if (is.null(additional_grouping_vars)) {
    grouping_vars <- mzroll_list$design$feature_pk
  } else {
    sample_measurement_vars <- design_tbl %>%
      dplyr::filter(table %in% c("samples", "measurements")) %>%
      dplyr::filter(type != "sample_primary_key")
    
    purrr::walk(
      additional_grouping_vars,
      checkmate::assertChoice,
      sample_measurement_vars$variable
    )
    
    grouping_vars <- c(mzroll_list$design$feature_pk, additional_grouping_vars)
  }
  
  return(grouping_vars)
}

validate_formulas <- function(model, value_var, viable_sample_fields) {
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

  full_formula <- stats::as.formula(paste(value_var, "~", model))

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
      "q-value calculation initially failed due to too many small p-values
        but claman was able to recover results"
    )
  }

  term_data %>%
    dplyr::mutate(qvalue = q_values)
}

#' Volcano plot
#'
#' @param regression_significance returned by \code{\link{diffex_mzroll}};
#'   a tibble of tests performed.
#' @param max_p_trans maximum value of -log10 pvalues to plot
#' @param FDR_cutoff FDR cutoff to label for significance
#'
#' @returns a grob
#'
#' @examples 
#' regression_significance <- diffex_mzroll(
#'   nplug_mzroll_normalized,
#'   "normalized_log2_abundance",
#'   "limitation + limitation:DR + 0"
#'   )
#'
#' plot_volcano(regression_significance, 10, 0.1)
#'
#' @export
plot_volcano <- function(
  regression_significance,
  max_p_trans = 10,
  FDR_cutoff = 0.1
  ) {
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
#' @inheritParams plot_volcano
#'
#' @examples 
#' regression_significance <- diffex_mzroll(
#'   nplug_mzroll_normalized,
#'   "normalized_log2_abundance",
#'   "limitation + limitation:DR + 0"
#'   )
#'
#' plot_pvalues(regression_significance)
#'
#' @export
plot_pvalues <- function(regression_significance) {
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
