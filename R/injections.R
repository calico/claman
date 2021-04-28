#' Collapse Injections
#'
#' @inheritParams test_mzroll_list
#' @param grouping_vars character vector of sample variables to use for
#'   grouping
#' @param peak_quant_vars character vector of quantification variables to
#'   average
#' @param collapse_fxn function name to use for collapse
#'
#' @return a \code{\link{process_mzroll}} collapsed across grouping_vars
#'
#' @examples
#' collapse_injections(
#'   nplug_mzroll_augmented,
#'   grouping_vars = "condition",
#'   peak_quant_vars = "log2_abundance"
#'   )
#' @export
collapse_injections <- function(
  mzroll_list,
  grouping_vars,
  peak_quant_vars,
  collapse_fxn = "mean"
) {
  
  checkmate::assertClass(mzroll_list, "tomic")
  checkmate::assertClass(mzroll_list, "mzroll")
  checkmate::assertCharacter(grouping_vars, min.len = 1)
  purrr::walk(
    grouping_vars,
    checkmate::assertChoice,
    choices = mzroll_list$design$samples$variable
  )
  checkmate::assertCharacter(peak_quant_vars, min.len = 1)
  purrr::walk(
    peak_quant_vars,
    checkmate::assertChoice,
    choices = mzroll_list$design$measurements$variable
  )
  checkmate::assertString(collapse_fxn)
  stopifnot(length(find(collapse_fxn, mode="function")) >= 1)
  
  # reduce samples to fields which are defined unambiguously with grouping vars
  
  found_injections <- find_injections(mzroll_list, grouping_vars)
  new_samples <- found_injections$new_samples
  collapse_dict <- found_injections$collapse_dict
  
  new_measurements <- mzroll_list$measurements %>%
    dplyr::rename(old_sampleId = sampleId) %>%
    dplyr::left_join(collapse_dict, by = "old_sampleId") %>%
    dplyr::select(-old_sampleId) %>%
    dplyr::group_by(groupId, sampleId) %>%
    dplyr::summarize_each(funs = collapse_fxn, peak_quant_vars) %>%
    dplyr::ungroup()
  
  # update sample and measurements and the schema
  mzroll_list <- romic::update_tomic(mzroll_list, new_measurements)
  mzroll_list <- romic::update_tomic(mzroll_list, new_samples)
  
  return(mzroll_list)
}

n_unique <- function(x) {
  length(unique(x))
}

#' Compare Injections
#' 
#' Create a scatterplot of injection correlations.
#' 
#' @inheritParams collapse_injections
#' @param peak_quant_var variable to plot in comparison
#' 
#' @examples
#' plot_compare_injection(
#'   nplug_mzroll_augmented,
#'   grouping_vars = "condition",
#'   peak_quant_var = "log2_abundance"
#'   )
#'   
#' @export
plot_compare_injection <- function(
  mzroll_list,
  grouping_vars = "condition",
  peak_quant_var = "log2_abundance"
  ) {
  
  checkmate::assertClass(mzroll_list, "tomic")
  checkmate::assertClass(mzroll_list, "mzroll")
  checkmate::assertCharacter(grouping_vars, min.len = 1)
  purrr::walk(
    grouping_vars,
    checkmate::assertChoice,
    choices = mzroll_list$design$samples$variable
  )
  checkmate::assertChoice(
    peak_quant_var,
    mzroll_list$design$measurements$variable
    )
  
  found_injections <- find_injections(mzroll_list, grouping_vars)
  collapse_dict <- found_injections$collapse_dict
  
  tidy_mzroll_list_data <- romic::triple_to_tidy(mzroll_list)$data %>%
    dplyr::select(
      groupId,
      old_sampleId = sampleId,
      !!rlang::sym(peak_quant_var)
      ) %>%
    dplyr::left_join(collapse_dict, by = "old_sampleId") %>%
    dplyr::mutate(.entry = 1:dplyr::n())
  
  quant_comparison <- tidy_mzroll_list_data %>%
    # generate all x all comparisons of sampe peakgroup and same unique sample
    dplyr::rename(quant_1 = !!rlang::sym(peak_quant_var)) %>%
    dplyr::inner_join(
      tidy_mzroll_list_data %>%
        dplyr::rename(quant_2 = !!rlang::sym(peak_quant_var)),
      by = c("groupId", "sampleId")) %>%
    # remove self comparisons, and only retain one pair for each {A-B, B-A}
    dplyr::filter(.entry.x < .entry.y)
  
  R_squared <- round(cor(
    quant_comparison$quant_1,
    quant_comparison$quant_2
    )^2, 3)
  
  label_dat <- tibble(
    quant_1 = min(quant_comparison$quant_1),
    quant_2 = max(quant_comparison$quant_2),
    )

  bivariate_hex <- ggplot(quant_comparison, aes(x = quant_1, y = quant_2)) +
    geom_hex(bins = 50) +
    geom_abline(intercept = 0, slope = 1, size = 1, color = "gray50") +
    geom_label(hjust = 0, vjust = 1, parse = TRUE,
      data = label_dat, 
      label = glue::glue("R^2 ~ \"=\" ~ {R_squared}")
    ) +
    scale_fill_viridis_c(trans = "log2", option = "plasma") +
    coord_equal() +
    scale_x_continuous(paste0(stringr::str_replace_all(
      peak_quant_var, "[-_]", " "
      ), " A"), expand = c(0,0)) +
    scale_y_continuous(paste0(stringr::str_replace_all(
      peak_quant_var, "[-_]", " "
    ), " B"), expand = c(0,0)) +
    theme_minimal() +
    theme(axis.line = element_line(color = "black"))
 
  return(bivariate_hex)
}

find_injections <- function(mzroll_list, grouping_vars) {
  
  checkmate::assertClass(mzroll_list, "tomic")
  checkmate::assertClass(mzroll_list, "mzroll")
  checkmate::assertCharacter(grouping_vars, min.len = 1)
  purrr::walk(
    grouping_vars,
    checkmate::assertChoice,
    choices = mzroll_list$design$samples$variable
  )
  
  not_grouping_vars <- setdiff(colnames(mzroll_list$samples), grouping_vars)
  
  consistent_reduced_fields <- mzroll_list$samples %>%
    dplyr::group_by(!!!syms(grouping_vars)) %>%
    dplyr::summarize_all(n_unique) %>%
    tidyr::gather(field, n_unique_records, !!!(not_grouping_vars)) %>%
    dplyr::group_by(field) %>%
    dplyr::filter(all(n_unique_records == 1)) %>%
    dplyr::distinct(field) %>%
    dplyr::ungroup()
  
  # names aren't consistent but should still be maintained
  reduced_names <- mzroll_list$samples %>%
    dplyr::select(name, !!!syms(grouping_vars)) %>%
    dplyr::group_by(!!!syms(grouping_vars)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  new_samples <- mzroll_list$samples %>%
    dplyr::select(c(!!!rlang::syms(grouping_vars), consistent_reduced_fields$field)) %>%
    dplyr::group_by(!!!syms(grouping_vars)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      sampleId = as.character(1:dplyr::n()),
      sampleId = factor(sampleId, levels = sampleId)
      ) %>%
    dplyr::left_join(reduced_names, by = grouping_vars)
  
  dropped_fields <- setdiff(colnames(mzroll_list$samples), colnames(new_samples))
  if (length(dropped_fields) > 0) {
    message(glue::glue(
      "{length(dropped_fields)} sample variables will be dropped since they
        - vary for the same grouping_vars:
        - {paste(dropped_fields, collapse = ', ')}"
    ))
  }
  
  collapse_dict <- mzroll_list$samples %>%
    dplyr::select(old_sampleId = sampleId, !!!syms(grouping_vars)) %>%
    dplyr::left_join(new_samples %>%
                       dplyr::select(sampleId, grouping_vars),
                     by = grouping_vars
    ) %>%
    dplyr::select(old_sampleId, sampleId)
  
  return(list(
    new_samples = new_samples,
    collapse_dict = collapse_dict
    ))
}
