#' Filter GroupIds
#'
#' Remove all data besides those from provided groupIds
#'
#' @inheritParams test_mzroll_list
#' @param groupIds groupIds to retain
#' @param invert if TRUE then remove provided groupIds; if FALSE then retain
#'   provided groupIds
#'
#' @examples
#' filter_groupIds(nplug_mzroll_normalized, factor(1:10))
#' @export
filter_groupIds <- function(mzroll_list, groupIds, invert = FALSE) {
  checkmate::assertFactor(groupIds)
  checkmate::assertLogical(invert, len = 1)

  if (invert) {
    filter_groupIds <- setdiff(mzroll_list$measurements$groupId, groupIds)
  } else {
    filter_groupIds <- groupIds
  }

  mzroll_list <- romic::filter_tomic(
    mzroll_list,
    filter_type = "quo",
    filter_table = "features",
    filter_value = rlang::quo(groupId %in% filter_groupIds)
  )

  return(mzroll_list)
}





#' Check groupIds for infinite or NA values
#'
#' Extract groupIds and optionally print compoundNames of peaks containing
#' non-finite values
#'
#' @param mzroll_list output of \link{process_mzroll} or
#' \link{process_mzroll_multi}
#' @param infinite_type one of "infinite" or "NA" to indicate whether to check
#' for NA values or non-NA infinite (Inf, -Inf, and NaN) values
#' @param quant_var variable in measurements to use for abundance
#' @param threshold threshold above which groupIds will be returned. A value of
#' 0 is recommended for \code{infinite_type = "infinite"} to return _any_
#' groupId containing any infinite value. A higher value or a value of 1 is
#' commended for \code{infinite_type = "NA"} to return groupIds with a high
#' fraction or complete fraction of NA values.
#' @inheritParams extract_ids_from_metadata
#' @param verbose when \code{TRUE} prints the compoundNames of groupIds that
#' exceed \code{threshold}
#'
#' @returns a character vector of groupIds recommended for checking or filtering
#'
#' @export
check_infinite_values <- function(mzroll_list,
                                  infinite_type = "infinite",
                                  quant_var = "centered_log2_abundance",
                                  threshold = 0,
                                  filter_var = NULL,
                                  filter_ids = NULL,
                                  verbose = TRUE) {
  # run checks
  romic:::check_triple_omic(mzroll_list)
  checkmate::assertChoice(quant_var, colnames(mzroll_list$measurements))
  checkmate::assertNumber(threshold, lower = 0, upper = 1)

  # if metadata filter desired
  if (!is.null(filter_ids) && !is.null(filter_var)) {
    is_filter <- mzroll_list %>%
      claman::extract_ids_from_metadata(
        filter_var = filter_var,
        filter_ids = filter_ids
      )
    mzroll_list <- mzroll_list %>%
      romic::filter_tomic(
        filter_type = "category",
        filter_table = is_filter$filter_var,
        filter_value = is_filter$filter_ids
      )
  }

  # run type
  check_methods <- tibble::tribble(
    ~method_name, ~function_name,
    "infinite", "check_peaks_infinite",
    "NA", "check_peaks_NA"
  )

  checkmate::assertChoice(infinite_type, check_methods$method_name)

  check_method_call <- check_methods %>%
    dplyr::filter(method_name == infinite_type) %>%
    dplyr::pull(function_name)

  check_method_output <- do.call(
    check_method_call,
    list(
      mzroll_list = mzroll_list,
      quant_var = quant_var,
      threshold = threshold
    )
  )

  if (verbose) {
    if (length(check_method_output) > 0) {
      problem_compounds <- mzroll_list$features %>%
        dplyr::filter(groupId %in% check_method_output)
      cat("\nNon-finite measurements associated with compoundNames: ")
      cat(paste(unique(problem_compounds$compoundName), collapse = ", "))
    } else {
      cat("\nAll peaks below non-finite value threshold")
    }
  }
  return(check_method_output)
}



#' Extract peaks containing non-NA, non-finite values
#'
#' @inheritParams check_infinite_values
#'
#' @rdname check_infinite_values
check_peaks_infinite <- function(mzroll_list, quant_var, threshold) {
  # get sample n
  unique_sampleId_n <- mzroll_list$samples %>%
    dplyr::n_distinct(.$sampleId)

  non_finite_groupIds <- mzroll_list %>%
    romic::filter_tomic(
      filter_type = "quo",
      filter_table = "measurements",
      filter_value = rlang::quo(is.infinite(!!rlang::sym(quant_var)) | is.nan(!!rlang::sym(quant_var)))
    ) %>%
    purrr::pluck("measurements") %>%
    dplyr::group_by(groupId) %>%
    dplyr::summarise(percent_inf = dplyr::n() / unique_sampleId_n) %>%
    dplyr::filter(percent_inf >= threshold) %>%
    dplyr::mutate(groupId = as.character(groupId)) %>%
    dplyr::pull(groupId)
  return(non_finite_groupIds)
}


#' Extract peaks containing non-NA, non-finite values
#'
#' @inheritParams check_infinite_values
#'
#' @rdname check_infinite_values
check_peaks_NA <- function(mzroll_list, quant_var, threshold) {
  # get sample n
  unique_sampleId_n <- mzroll_list$samples %>%
    dplyr::n_distinct(.$sampleId)

  high_na_groupIds <- mzroll_list %>%
    claman::expand_peaks() %>%
    romic::filter_tomic(
      filter_type = "quo",
      filter_table = "measurements",
      filter_value = rlang::quo(is.na(!!rlang::sym(quant_var)))
    ) %>%
    romic::filter_tomic(
      filter_type = "quo",
      filter_table = "measurements",
      filter_value = rlang::quo(!is.infinite(!!rlang::sym(quant_var)))
    ) %>%
    romic::filter_tomic(
      filter_type = "quo",
      filter_table = "measurements",
      filter_value = rlang::quo(!is.nan(!!rlang::sym(quant_var)))
    ) %>%
    purrr::pluck("measurements") %>%
    dplyr::group_by(groupId) %>%
    dplyr::summarise(percent_na = dplyr::n() / unique_sampleId_n) %>%
    dplyr::filter(percent_na >= threshold) %>%
    dplyr::mutate(groupId = as.character(groupId)) %>%
    dplyr::pull(groupId)
  return(high_na_groupIds)
}





#' Extract groupIds or sampleIds per metadata
#'
#' @description
#' For a given filter_var and filter_ids, extract the associated groupIds and
#' sampleIds for easy filtering
#'
#' @param mzroll_list output of \link{process_mzroll} or
#' \link{process_mzroll_multi}
#' @param filter_ids sample or peak type ids on which to filter data from column
#' provided by \code{filter_var}
#' @param filter_var column name on which to filter \code{filter_ids}; must be
#' a column name present in the \code{samples} or \code{features} tables
#'
#' @returns a list of 2 objects, where the first object is one of "sampleId" or
#' "groupId" and the second is a character vector of associated Ids
#'
#' @export
extract_ids_from_metadata <- function(mzroll_list,
                                      filter_var,
                                      filter_ids) {
  # run checks
  romic:::check_triple_omic(mzroll_list)
  checkmate::assertString(filter_var)
  checkmate::assertString(filter_ids)
  checkmate::assertChoice(
    filter_var,
    union(
      colnames(mzroll_list$samples),
      colnames(mzroll_list$features)
    )
  )

  if (filter_var %in% colnames(mzroll_list$samples)) {
    if (any(filter_ids %in% unique(mzroll_list$samples[[filter_var]]))) {
      filter_var_use <- "sampleId"
      filter_ids_use <- mzroll_list$samples %>%
        dplyr::filter(!!rlang::sym(filter_var) %in% filter_ids) %>%
        dplyr::distinct(sampleId) %>%
        dplyr::pull()
    } else {
      stop("\nfilter_var identified as a column in samples table, but no filter_ids matched values in the column")
    }
  }

  if (filter_var %in% colnames(mzroll_list$features)) {
    if (any(filter_ids %in% unique(mzroll_list$features[[filter_var]]))) {
      filter_var_use <- "groupId"
      filter_ids_use <- mzroll_list$features %>%
        dplyr::filter(!!rlang::sym(filter_var) %in% filter_ids) %>%
        dplyr::distinct(groupId) %>%
        dplyr::pull()
    } else {
      stop("\nfilter_var identified as a column in features table, but no filter_ids matched values in the column")
    }
  }

  return(list(
    filter_var = filter_var_use,
    filter_ids = filter_ids_use
  ))
}
