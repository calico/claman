#' Expand all samples and features
#'
#' This function takes an mzroll_list (mzroll_tomic) object and expands
#' the features and samples in the measurements dataframe. All original peak
#' area values are retained. If \code{log2_floor_value} is not provided, all
#' missing values are enumerated with NA values. If a \code{log2_floor_value} is
#' provided as a numeric value, all NA values will be replaced with that value.
#' The \code{log2_floor_value} argument can also accept strings "min" or
#' "halfmin" which will replace NA values with groupId-specific log2 minimum or
#' half-minimum values, respectively.
#'
#' This function differs from \code{claman::floor_peaks} as values below the
#' \code{log2_floor_value} are _not_ replaced with the floor values.
#'
#' @param mzroll_list an \code{mzroll_list} tomic object
#' @param quant_var peak measurement variable on which to determine missing
#' values
#' @param log2_floor_value minimum value to set for missing peaks only. Will
#' also accept "min" or "halfmin" to impute NA values with groupId-specific
#' minimum or half-minimum values, respectively
#'
#' @returns an mzroll_list object with groupId and sampleId values expanded
#'
#' @export
expand_peaks <- function(mzroll_list,
                         quant_var = "log2_abundance",
                         log2_floor_value = NULL) {
  claman:::test_mzroll_list(mzroll_list)
  missing_peaks <- tidyr::expand_grid(
    groupId = mzroll_list$features$groupId,
    sampleId = mzroll_list$samples$sampleId
  ) %>%
    dplyr::anti_join(mzroll_list$measurements,
      by = c("groupId", "sampleId")
    ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(`:=`(!!rlang::sym(quant_var), NA))

  completed_peaks <- dplyr::bind_rows(
    mzroll_list$measurements,
    missing_peaks
  )

  if (!is.null(log2_floor_value)) {
    # replace NAs with floor value
    if (is.numeric(log2_floor_value) && length(log2_floor_value) == 1) {
      completed_peaks <- completed_peaks %>%
        dplyr::mutate(!!rlang::sym(quant_var) :=
          ifelse(is.na(!!rlang::sym(quant_var)),
            log2_floor_value,
            !!rlang::sym(quant_var)
          ))

      # replace NAs with groupId minimum
    } else if (log2_floor_value == "min") {
      completed_peaks <- completed_peaks %>%
        dplyr::group_by(groupId) %>%
        dplyr::mutate(!!rlang::sym(quant_var) :=
          ifelse(is.na(!!rlang::sym(quant_var)),
            min(!!rlang::sym(quant_var), na.rm = TRUE),
            !!rlang::sym(quant_var)
          )) %>%
        dplyr::ungroup()

      # replace NAs with groupId half minimum
    } else if (log2_floor_value == "halfmin") {
      completed_peaks <- completed_peaks %>%
        dplyr::group_by(groupId) %>%
        dplyr::mutate(!!rlang::sym(quant_var) :=
          ifelse(is.na(!!rlang::sym(quant_var)),
            min(!!rlang::sym(quant_var), na.rm = TRUE) - 1,
            !!rlang::sym(quant_var)
          )) %>%
        dplyr::ungroup()
    } else {
      warning("log2_floor_value must be a single number, 'min', or 'halfmin'. Returning expanded mzroll_list with NA values.")
    }
  }

  mzroll_list$measurements <- completed_peaks
  return(mzroll_list)
}




#' Floor Peaks
#'
#' Set a minimum peak abundance of floor_value for low abundance and
#'   undetected peaks.
#'
#' @inheritParams test_mzroll_list
#' @param log2_floor_value minimum value to set for low abundance or
#'   missing peaks
#' @param floor_var measurement variable to floor to \code{log2_floor_value}.
#'
#' @return \code{\link{process_mzroll}}
#'
#' @examples
#' floored_peaks <- floor_peaks(nplug_mzroll_augmented, 12)
#' @export
floor_peaks <- function(mzroll_list,
                        log2_floor_value = 12,
                        floor_var = "log2_abundance") {
  test_mzroll_list(mzroll_list)

  valid_floor_var <- setdiff(
    mzroll_list$design$measurements$variable,
    c(mzroll_list$design$feature_pk, mzroll_list$design$sample_pk)
  )

  checkmate::assertChoice(floor_var, valid_floor_var)
  checkmate::assertNumeric(mzroll_list$measurements[[floor_var]])
  checkmate::assertNumber(log2_floor_value)

  # summarize peaks associated with each peakgroup

  missing_peaks <- tidyr::expand_grid(
    groupId = mzroll_list$features$groupId,
    sampleId = mzroll_list$samples$sampleId
  ) %>%
    dplyr::anti_join(
      mzroll_list$measurements,
      by = c("groupId", "sampleId")
    ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(!!rlang::sym(floor_var) := log2_floor_value)

  # combine detected peaks with peaks that were missing for some samples
  completed_peaks <- dplyr::bind_rows(
    mzroll_list$measurements %>%
      dplyr::mutate(!!rlang::sym(floor_var) := dplyr::case_when(
        is.na(!!rlang::sym(floor_var)) ~ log2_floor_value,
        !!rlang::sym(floor_var) < log2_floor_value ~ log2_floor_value,
        !!rlang::sym(floor_var) >= log2_floor_value ~ !!rlang::sym(floor_var)
      )),
    missing_peaks
  )

  mzroll_list$measurements <- completed_peaks

  return(mzroll_list)
}

#' Impute missing peaks with provided feature-level imputation values
#'
#' @param mzroll_list data in triple omic structure
#' @param lod_values a tibble that maps groupId to log2 feature-level imputation values
#' @param quant_var column to use for peak values
#' @param imputation_sd standard deviation of Gaussian distribution to use for missing peak imputation
#'
#' @return triple omic data with imputed missing peaks
#'
#' @examples
#' library(dplyr)
#'
#' lod_values <- nplug_mzroll_augmented[["measurements"]] %>%
#'   dplyr::select(groupId, log2_abundance) %>%
#'   dplyr::distinct(groupId, .keep_all = TRUE)
#'
#' mzroll_list <- nplug_mzroll_augmented
#' mzroll_list$measurements <- mzroll_list$measurements %>% filter(groupId != 2 & sampleId != 1)
#'
#' mzroll_list_imputed <- impute_missing_peaks(mzroll_list, lod_values, "log2_abundance", 0.15)
#'
#' @export
impute_missing_peaks <- function(mzroll_list,
                                 lod_values,
                                 quant_var = "log2_abundance",
                                 imputation_sd = 0.15) {
  test_mzroll_list(mzroll_list)

  stopifnot(colnames(lod_values) %in% c("groupId", rlang::sym(quant_var)))

  if (nrow(lod_values) > nrow(lod_values %>% dplyr::distinct(groupId, .keep_all = TRUE))) {
    stop("only one value per feature must be provided to impute missing peaks")
  }

  features <- mzroll_list$features %>% dplyr::select(groupId)

  if (!all(features$groupId %in% lod_values$groupId)) {
    stop("groupId values of lod_value table and feature table of triple omic data must match")
  }

  valid_quant_var <- setdiff(
    mzroll_list$design$measurements$variable,
    c(mzroll_list$design$feature_pk, mzroll_list$design$sample_pk)
  )

  checkmate::assertChoice(quant_var, valid_quant_var)
  checkmate::assertNumeric(mzroll_list$measurements[[quant_var]])

  ## find missing peaks
  missing_peaks <- tidyr::expand_grid(
    groupId = mzroll_list$features$groupId,
    sampleId = mzroll_list$samples$sampleId
  ) %>%
    dplyr::anti_join(
      mzroll_list$measurements,
      by = c("groupId", "sampleId")
    )

  ## impute missing peaks
  missing_peaks_imputed <- dplyr::left_join(
    missing_peaks,
    lod_values,
    by = c("groupId")
  ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(!!rlang::sym(quant_var) := stats::rnorm(1, mean = !!rlang::sym(quant_var) + 1, sd = imputation_sd))

  ## merge measured peaks with imputed peaks
  completed_peaks <- dplyr::bind_rows(
    mzroll_list$measurements,
    missing_peaks_imputed
  )

  mzroll_list$measurements <- completed_peaks
  return(mzroll_list)
}

#' Fill in missing peaks
#'
#' @description
#' If \code{fill_values} is a data frame, this function calls \code{impute_missing_peaks()}.
#' If it is a numeric vector, 'this function calls \code{floor_peaks()}. Other types are not currently supported
#'
#' @param mzroll_list data in triple omic structure
#' @param fill_values either a numeric constant or a tibble that maps groupId to log2 feature-level imputation values
#' @param quant_var column to use for peak values
#' @param imputation_sd standard deviation of Gaussian distribution to use for missing peak imputation
#'
#' @return triple omic data with imputed missing peaks
#'
#' @examples
#' library(dplyr)
#'
#' lod_values <- nplug_mzroll_augmented[["measurements"]] %>%
#'   dplyr::select(groupId, log2_abundance) %>%
#'   dplyr::distinct(groupId, .keep_all = TRUE)
#'
#' mzroll_list <- nplug_mzroll_augmented
#' mzroll_list$measurements <- mzroll_list$measurements %>% filter(groupId != 2 & sampleId != 1)
#'
#' mzroll_list_imputed <- fill_in_missing_peaks(mzroll_list, lod_values, "log2_abundance", 0.15)
#' mzroll_list_imputed <- fill_in_missing_peaks(mzroll_list, 12, "log2_abundance")
#'
#' @export
fill_in_missing_peaks <- function(mzroll_list,
                                  fill_values,
                                  quant_var = "log2_abundance",
                                  imputation_sd = 0.15) {
  if (is.data.frame(fill_values)) {
    output <- impute_missing_peaks(mzroll_list, fill_values, quant_var, imputation_sd)
  } else if (is.numeric(fill_values)) {
    output <- floor_peaks(mzroll_list, fill_values, quant_var)
  } else {
    stop("fill_values should either be a single numeric value or a tibble that maps groupId to feature-level imputation values")
  }
  return(output)
}

#' Normalize Peaks
#'
#' @inheritParams test_mzroll_list
#' @param normalization_method Normalization method to apply
#' \itemize{
#'   \item{\code{median polish}: column normalization based on average signal
#'   (adds \code{median_polish_scaling_factor} as feature variable in addition to
#'     \code{norm_peak_varname})}
#'   \item{\code{loading value}: column normalization using a sample-level
#'     value
#'     }
#'   \item{\code{center batches}: batch centering}
#'   \item{\code{center}: zero-center each compound}
#'   \item{\code{reference sample}: compare to reference sample}
#'   \item{\code{reference condition}: compare each sample to its specified
#'     reference sample}
#'   \item{\code{loess}: weighted smoothing of IC over time (adds
#'     \code{.loess_fit} as a peaks variable in addition to
#'     \code{norm_peak_varname})}
#'   \item{\code{lm}: linear regression smoothing of IC over time (adds
#'     \code{lm_estimate} as feature variable in addition to
#'     \code{norm_peak_varname})}
#'   }
#' @param quant_peak_varname variable in measurements to use for abundance
#' @param norm_peak_varname variable in measurements to add for normalized
#'   abundance
#' @param ... additional arguments to pass to normalization methods
#'
#' @return a \code{mzroll_list} as generated by \code{\link{process_mzroll}}
#'   a \code{norm_peak_varname} variable added to measuremnets
#'
#' @examples
#'
#' normalize_peaks(
#'   nplug_mzroll_augmented,
#'   "median polish",
#'   quant_peak_varname = "log2_abundance",
#'   norm_peak_varname = "normalized_log2_abundance"
#' )
#'
#' # this examples doesn't make biological sense but its syntactically correct
#' normalize_peaks(
#'   nplug_mzroll_augmented,
#'   "loading value",
#'   quant_peak_varname = "log2_abundance",
#'   norm_peak_varname = "normalized_log2_abundance",
#'   loading_varname = "DR"
#' )
#'
#' normalize_peaks(
#'   nplug_mzroll_augmented,
#'   "center batches",
#'   quant_peak_varname = "log2_abundance",
#'   norm_peak_varname = "normalized_log2_abundance",
#'   batch_varnames = "month",
#'   centering_fxn = median
#' )
#'
#' normalize_peaks(
#'   nplug_mzroll_augmented,
#'   normalization_method = "reference condition",
#'   quant_peak_varname = "log2_abundance",
#'   norm_peak_varname = "normalized_log2_abundance",
#'   condition_varname = "condition",
#'   reference_varname = "reference"
#' )
#'
#' normalize_peaks(
#'   nplug_mzroll_augmented,
#'   normalization_method = "reference sample",
#'   quant_peak_varname = "log2_abundance",
#'   norm_peak_varname = "normalized_log2_abundance",
#'   batch_varnames = c("month", "extraction"),
#'   reference_varname = "exp_ref",
#'   reference_values = "ref"
#' )
#' @export
normalize_peaks <- function(mzroll_list,
                            normalization_method,
                            quant_peak_varname = "log2_abundance",
                            norm_peak_varname = "normalized_log2_IC",
                            ...) {
  dots <- list(...)
  test_mzroll_list(mzroll_list)

  checkmate::assertString(quant_peak_varname)
  if (!(quant_peak_varname %in% colnames(mzroll_list$measurements))) {
    stop(
      "\"quant_peak_varname\":",
      quant_peak_varname,
      ", not present in measurements table"
    )
  }
  checkmate::assertString(norm_peak_varname)

  normalization_methods <- tibble::tribble(
    ~method_name, ~function_name,
    "median polish", "normalize_peaks_median_polish",
    "loading value", "normalize_peaks_loading_value",
    "center batches", "normalize_peaks_batch_center",
    "reference sample", "normalize_peaks_reference_sample",
    "reference condition", "normalize_peaks_reference_condition",
    "loess", "normalize_peaks_loess",
    "lm", "normalize_peaks_lm",
    "center", "normalize_peaks_center"
  )

  checkmate::assertChoice(
    normalization_method,
    normalization_methods$method_name
  )

  normalization_method_call <- normalization_methods$function_name[
    normalization_methods$method_name == normalization_method
  ]

  method_args <- dots[intersect(
    names(dots),
    names(formals(normalization_method_call))
  )]

  unused_args <- setdiff(names(dots), names(formals(normalization_method_call)))
  if (length(unused_args) > 0) {
    warning(glue::glue(
      "{length(unused_args)} passed arguments could not be used by
       {normalization_method_call}: {paste(unused_args, collapse = ', ')} "
    ))
  }

  normalization_output <- do.call(
    normalization_method_call,
    append(
      list(
        mzroll_list = mzroll_list,
        quant_peak_varname = quant_peak_varname,
        norm_peak_varname = norm_peak_varname
      ),
      method_args
    )
  )

  normalization_output
}

#' Normalize Peaks - Median Polish
#'
#' Using a robust metabolomics normalization rescale ion counts based on the
#'   consensus signal of each sample.
#'
#' @details The robust median polish was reported in Kamphorst et al. 2015, see
#'   \url{https://github.com/shackett/Pancreatic_tumor_metabolomics}.
#'
#' @inheritParams normalize_peaks
#' @inheritParams floor_peaks
#' @param filter_ids sample type ids on which to filter data for calculation of
#' median polish quotient from column provided by \code{filter_var};
#' defaults to NULL
#' @param filter_var column name on which to filter \code{filter_ids}; should be
#' a column name present in the \code{samples} or \code{features} tables; defaults to NULL
#'
#' @rdname normalize_peaks
normalize_peaks_median_polish <- function(mzroll_list,
                                          quant_peak_varname,
                                          norm_peak_varname,
                                          filter_ids = NULL,
                                          filter_var = NULL,
                                          log2_floor_value = NA) {
  stopifnot(length(log2_floor_value) == 1)
  if (!is.na(log2_floor_value)) {
    stopifnot(class(log2_floor_value) == "numeric")

    normalization_peaks <- mzroll_list$measurements %>%
      dplyr::filter(
        !!rlang::sym(quant_peak_varname) > log2_floor_value + 0.001
      )
  } else {
    normalization_peaks <- mzroll_list$measurements
  }

  if (any(is.na(mzroll_list$measurements[[quant_peak_varname]]))) {
    stop(
      "NAs are not allowed in ",
      quant_peak_varname,
      " try calling floor_peaks() first"
    )
  }

  ## Check for any filter matching
  filter_var_use <- NULL
  if (!is.null(filter_ids) && !is.null(filter_var)) {
    if (filter_var %in% colnames(mzroll_list$samples)) {
      if (any(filter_ids %in% unique(mzroll_list$samples[[filter_var]]))) {
        filter_var_use <- "sampleId"
        filter_ids_use <- mzroll_list$samples %>%
          dplyr::filter(!!rlang::sym(filter_var) %in% filter_ids) %>%
          dplyr::distinct(sampleId) %>%
          dplyr::pull()
      }
    }
    if (filter_var %in% colnames(mzroll_list$features)) {
      if (any(filter_ids %in% unique(mzroll_list$features[[filter_var]]))) {
        filter_var_use <- "groupId"
        filter_ids_use <- mzroll_list$features %>%
          dplyr::filter(!!rlang::sym(filter_var) %in% filter_ids) %>%
          dplyr::distinct(groupId) %>%
          dplyr::pull()
      }
    }
  }

  ## No filter option (most common)
  if (is.null(filter_var_use)) {
    sample_scaling_factors <- normalization_peaks %>%
      dplyr::group_by(groupId) %>%
      dplyr::mutate(
        median_abund = stats::median(!!rlang::sym(quant_peak_varname))
      ) %>%
      dplyr::ungroup()

    ## Filter on groupIds
    ## This is a version of PQN that will only calculate median differences on a
    ## subset of groupIds (peaks) -- useful if someone wants to subset on "good"
    ## or robust peaks only (or even just a few known, low variability peaks)
  } else if (filter_var_use == "groupId") {
    cat(paste0(
      "PQN median calculations filtered to peaks with ",
      filter_var, " values matching to ",
      paste(
        intersect(
          filter_ids,
          unique(mzroll_list$features[[filter_var]])
        ),
        collapse = ", "
      ), "\n"
    ))

    sample_scaling_factors <- normalization_peaks %>%
      dplyr::filter(groupId %in% filter_ids_use) %>%
      dplyr::group_by(groupId) %>%
      dplyr::mutate(
        median_abund = stats::median(!!rlang::sym(quant_peak_varname))
      ) %>%
      dplyr::ungroup()

    ## Filter on sampleIds
    ## This is a version of PQN that will only calculate median values of each
    ## compound based on a subset of sampleIds. This is useful if it is desired
    ## for the median reference value of each compound to only be calculated on
    ## QC or biological samples. Note that due to selecting a median value on
    ## default, filtering in this manner does not tend to change results
    ## significantly, and may only be useful for very small datasets with a
    ## large proportion of blank or standard samples
  } else if (filter_var_use == "sampleId") {
    cat(paste0(
      "PQN median calculations filtered to samples with ",
      filter_var, " values matching to ",
      paste(
        intersect(
          filter_ids,
          unique(mzroll_list$samples[[filter_var]])
        ),
        collapse = ", "
      ), "\n"
    ))

    sample_scaling_factors_temp <- normalization_peaks %>%
      dplyr::filter(sampleId %in% filter_ids_use) %>%
      dplyr::group_by(groupId) %>%
      dplyr::mutate(
        median_abund = stats::median(!!rlang::sym(quant_peak_varname))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(groupId, median_abund)

    sample_scaling_factors <- normalization_peaks %>%
      dplyr::left_join(., sample_scaling_factors_temp, by = "groupId")
  } else {
    stop("Issue filtering measurements dataframe for median polish")
  }

  ## Calculate median difference of each groupId
  sample_scaling_factors <- sample_scaling_factors %>%
    dplyr::mutate(
      diff_to_median = !!rlang::sym(quant_peak_varname) - median_abund
    ) %>%
    dplyr::group_by(sampleId) %>%
    dplyr::summarize(
      median_polish_scaling_factor = stats::median(diff_to_median),
      .groups = "drop"
    )

  missing_sample_scaling_factors <- mzroll_list$samples %>%
    dplyr::anti_join(sample_scaling_factors, by = "sampleId")

  if (nrow(missing_sample_scaling_factors) != 0) {
    stop(
      "sample scaling factors could not be calculated for ",
      nrow(missing_sample_scaling_factors),
      " samples"
    )
  }

  if (!is.na(log2_floor_value)) {
    updated_measurements <- mzroll_list$measurements %>%
      dplyr::filter(
        !!rlang::sym(quant_peak_varname) > log2_floor_value + 0.001
      ) %>%
      dplyr::left_join(sample_scaling_factors, by = "sampleId") %>%
      dplyr::mutate(
        !!rlang::sym(norm_peak_varname) :=
          !!rlang::sym(quant_peak_varname) - median_polish_scaling_factor
      ) %>%
      # measurements starting at limit of detection are reset to
      #   log2_floor_value
      dplyr::bind_rows(
        mzroll_list$measurements %>%
          dplyr::filter(
            !!rlang::sym(quant_peak_varname) <= log2_floor_value + 0.001
          ) %>%
          dplyr::mutate(!!rlang::sym(norm_peak_varname) := log2_floor_value)
      ) %>%
      dplyr::select(-median_polish_scaling_factor) %>%
      # measurements pushed below limit of detection are reset to
      #   log2_floor_value
      dplyr::mutate(
        !!rlang::sym(norm_peak_varname) :=
          pmax(!!rlang::sym(norm_peak_varname), log2_floor_value)
      )
  } else {
    updated_measurements <- mzroll_list$measurements %>%
      dplyr::left_join(sample_scaling_factors, by = "sampleId") %>%
      dplyr::mutate(
        !!rlang::sym(norm_peak_varname) :=
          !!rlang::sym(quant_peak_varname) - median_polish_scaling_factor
      ) %>%
      dplyr::select(-median_polish_scaling_factor)
  }

  mzroll_list <- romic::update_tomic(mzroll_list, updated_measurements)

  updated_samples <- mzroll_list$samples %>%
    dplyr::left_join(., sample_scaling_factors, by = "sampleId")

  mzroll_list <- romic::update_tomic(mzroll_list, updated_samples)

  return(mzroll_list)
}


#' Predict Dilutions from Median Polish Scaling Factor
#'
#' @description
#' Using `median_polish_scaling_factor` output from `normalize_peaks_median_polish`,
#' predict sample-wise dilutions
#'
#' @details This function performs an inverse-log transformation on the scaling factors,
#' and scales the samples based on the observed maximum scaling factor. Assumes that
#' median polish is perform on log2 transformed data
#'
#' @param mzroll_list data in triple omic structure
#' @param scaling_factor scaling factor, defaults to `median_polish_scaling_factor`
#' and must be a `samples` column in `mzroll_list`
#' @param norm_scale_varname variable in samples to add for dilution predictions
#' @param group_var optional grouping variable on which to calculate maximum dilution,
#' must be a `samples` column in `mzroll_list`
#'
#' @return a \code{mzroll_list} with \code{norm_scale_varname} variable added to
#' samples
#'
#' @export
median_polish_predict_dilutions <- function(mzroll_list,
                                            scaling_factor = "median_polish_scaling_factor",
                                            norm_scale_varname = "median_polish_predicted_dilutions",
                                            group_var = NULL) {
  test_mzroll_list(mzroll_list)

  checkmate::assertString(scaling_factor)
  if (!(scaling_factor %in% colnames(mzroll_list$samples))) {
    stop(
      "\"scaling_factor\":",
      scaling_factor,
      ", not present in samples table"
    )
  }
  checkmate::assertString(norm_scale_varname)

  updated_samples <- mzroll_list$samples %>%
    dplyr::mutate(temp_scaling_factor = !!rlang::sym(scaling_factor)) %>%
    dplyr::mutate(inverse_log_scaling_factor = 2^temp_scaling_factor)

  if (!is.null(group_var) && any(group_var %in% colnames(mzroll_list$samples))) {
    updated_samples <- updated_samples %>%
      dplyr::group_by_at(group_var) %>%
      dplyr::mutate(m = max(inverse_log_scaling_factor))
  } else {
    max_temp <- max(updated_samples$inverse_log_scaling_factor, na.rm = T)
    updated_samples <- updated_samples %>%
      dplyr::mutate(m = .env$max_temp)
  }

  updated_samples <- updated_samples %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(`:=`(
      !!rlang::sym(norm_scale_varname),
      inverse_log_scaling_factor / m
    )) %>%
    dplyr::select(-temp_scaling_factor, -inverse_log_scaling_factor, -m)

  mzroll_list <- romic::update_tomic(mzroll_list, updated_samples)
  return(mzroll_list)
}

#' Normalize Peaks - Loading Value
#'
#' Using a sample-level summary, such as number of cells, adjust all values.
#'
#' @details Note, log2 intensities will be shifted down by this value, so
#'   make sure the provided values are appropriately transformed.
#'
#' @inheritParams normalize_peaks
#' @inheritParams floor_peaks
#' @param loading_varname sample variable used for adjustment
#'
#' @rdname normalize_peaks
normalize_peaks_loading_value <- function(mzroll_list,
                                          quant_peak_varname,
                                          norm_peak_varname,
                                          loading_varname,
                                          log2_floor_value = NA) {
  stopifnot(length(log2_floor_value) == 1)
  if (!is.na(log2_floor_value)) {
    stopifnot(class(log2_floor_value) == "numeric")

    normalization_peaks <- mzroll_list$measurements %>%
      dplyr::filter(
        !!rlang::sym(quant_peak_varname) > log2_floor_value + 0.001
      )
  } else {
    normalization_peaks <- mzroll_list$measurements
  }

  if (any(is.na(mzroll_list$measurements[[quant_peak_varname]]))) {
    stop(
      "NAs are not allowed in ",
      quant_peak_varname,
      " try calling floor_peaks() first"
    )
  }

  # validate loading_varname

  checkmate::assertChoice(loading_varname, mzroll_list$design$samples$variable)

  loading_values <- mzroll_list$samples %>%
    dplyr::select("sampleId", scaling_factor = !!rlang::sym(loading_varname))

  checkmate::assertNumeric(
    loading_values$scaling_factor,
    any.missing = FALSE,
    finite = TRUE
  )

  if (!is.na(log2_floor_value)) {
    updated_measurements <- mzroll_list$measurements %>%
      dplyr::filter(
        !!rlang::sym(quant_peak_varname) > log2_floor_value + 0.001
      ) %>%
      dplyr::left_join(loading_values, by = "sampleId") %>%
      dplyr::mutate(
        !!rlang::sym(norm_peak_varname) :=
          !!rlang::sym(quant_peak_varname) - scaling_factor
      ) %>%
      # measurements starting at limit of detection are reset to
      #   log2_floor_value
      dplyr::bind_rows(
        mzroll_list$measurements %>%
          dplyr::filter(
            !!rlang::sym(quant_peak_varname) <= log2_floor_value + 0.001
          ) %>%
          dplyr::mutate(!!rlang::sym(norm_peak_varname) := log2_floor_value)
      ) %>%
      dplyr::select(-scaling_factor) %>%
      # measurements pushed below limit of detection are reset to
      #   log2_floor_value
      dplyr::mutate(
        !!rlang::sym(norm_peak_varname) :=
          pmax(!!rlang::sym(norm_peak_varname), log2_floor_value)
      )
  } else {
    updated_measurements <- mzroll_list$measurements %>%
      dplyr::left_join(loading_values, by = "sampleId") %>%
      dplyr::mutate(
        !!rlang::sym(norm_peak_varname) :=
          !!rlang::sym(quant_peak_varname) - scaling_factor
      ) %>%
      dplyr::select(-scaling_factor)
  }

  mzroll_list <- romic::update_tomic(mzroll_list, updated_measurements)

  return(mzroll_list)
}

#' Normalize Peaks - Batch Center
#'
#' @inheritParams normalize_peaks
#' @inheritParams floor_peaks
#' @param batch_varnames variables to use for grouping when removing batch
#'   effects
#' @param centering_fxn function to use when centering (mean, median, ...)
#'
#' @rdname normalize_peaks
normalize_peaks_batch_center <- function(mzroll_list,
                                         quant_peak_varname,
                                         norm_peak_varname,
                                         batch_varnames,
                                         log2_floor_value = NA,
                                         centering_fxn = mean) {
  checkmate::assertChoice(batch_varnames, colnames(mzroll_list$samples))
  stopifnot(length(log2_floor_value) == 1)
  if (!is.na(log2_floor_value)) {
    stopifnot(class(log2_floor_value) == "numeric")
  }
  checkmate::assertFunction(centering_fxn)

  centered_peaks <- mzroll_list$measurements %>%
    dplyr::left_join(
      mzroll_list$samples %>%
        dplyr::select(sampleId, !!!rlang::syms(batch_varnames)),
      by = "sampleId"
    ) %>%
    dplyr::group_by(groupId, !!!rlang::syms(batch_varnames)) %>%
    dplyr::mutate(
      .center = centering_fxn(!!rlang::sym(quant_peak_varname))
    ) %>%
    dplyr::group_by(groupId) %>%
    dplyr::mutate(
      .group_center = centering_fxn(!!rlang::sym(quant_peak_varname))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      !!rlang::sym(norm_peak_varname) :=
        !!rlang::sym(quant_peak_varname) - (.center - .group_center)
    ) %>%
    dplyr::select(
      !!!rlang::syms(c(colnames(mzroll_list$measurements), norm_peak_varname))
    )

  centered_refloor <- normalization_refloor(
    centered_peaks,
    log2_floor_value,
    norm_peak_varname,
    quant_peak_varname
  )

  mzroll_list <- romic::update_tomic(mzroll_list, centered_refloor)

  return(mzroll_list)
}

#' Normalize Peaks - Reference Sample
#'
#' Compare each sample to the reference samples within the same batch
#'
#' @inheritParams normalize_peaks_batch_center
#' @param reference_varname variable specifying which samples are references
#' @param reference_values values of \code{reference_varname} indicating
#'   reference samples
#'
#' @rdname normalize_peaks
normalize_peaks_reference_sample <- function(mzroll_list,
                                             quant_peak_varname,
                                             norm_peak_varname,
                                             batch_varnames,
                                             reference_varname = "sample",
                                             reference_values = "posctl",
                                             log2_floor_value = NA) {
  purrr::walk(
    batch_varnames,
    checkmate::assertChoice,
    colnames(mzroll_list$samples)
  )
  checkmate::assertChoice(reference_varname, colnames(mzroll_list$samples))
  checkmate::assertCharacter(reference_values)

  stopifnot(length(log2_floor_value) == 1)
  if (!is.na(log2_floor_value)) {
    stopifnot(class(log2_floor_value) == "numeric")
  }

  annotated_peaks <- mzroll_list$measurements %>%
    dplyr::left_join(
      mzroll_list$samples %>%
        dplyr::select(
          sampleId,
          !!!rlang::syms(c(batch_varnames, reference_varname))
        ),
      by = "sampleId"
    )

  reference_peaks <- annotated_peaks %>%
    dplyr::filter(
      !!rlang::quo(!!rlang::sym(reference_varname) %in% reference_values)
    ) %>%
    dplyr::group_by(groupId, !!!rlang::syms(batch_varnames)) %>%
    dplyr::summarize(
      .reference = median(!!rlang::sym(quant_peak_varname)),
      .groups = "drop"
    )

  if (is.na(log2_floor_value)) {
    # use the relative scale
    reference_peaks <- reference_peaks %>%
      dplyr::mutate(.reference_diff = .reference)
  } else {
    # retain the original scale
    reference_peaks <- reference_peaks %>%
      dplyr::group_by(groupId) %>%
      dplyr::mutate(.reference_diff = .reference - mean(.reference)) %>%
      dplyr::ungroup()
  }

  normalized_peaks <- annotated_peaks %>%
    dplyr::left_join(reference_peaks, by = c("groupId", batch_varnames)) %>%
    dplyr::mutate(
      !!rlang::sym(norm_peak_varname) :=
        !!rlang::sym(quant_peak_varname) - .reference_diff
    ) %>%
    dplyr::select(
      !!!rlang::syms(c(colnames(mzroll_list$measurements), norm_peak_varname))
    )

  # ensure that normalized values still respect the limit of detection
  reference_refloor <- normalization_refloor(
    normalized_peaks,
    log2_floor_value,
    norm_peak_varname,
    quant_peak_varname
  )

  mzroll_list <- romic::update_tomic(mzroll_list, reference_refloor)

  return(mzroll_list)
}

#' Normalize Peaks - Reference Condition
#'
#' Compare each sample to the reference samples within the same batch
#'
#' @inheritParams normalize_peaks_batch_center
#' @param condition_varname variable specifying each sample's condition
#' @inheritParams normalize_peaks_reference_sample
#'
#' @details
#' reference_varname specifies which condition to contrast a given sample to
#'   (using the condition values from condition_varname)
#'
#' @rdname normalize_peaks
normalize_peaks_reference_condition <- function(mzroll_list,
                                                quant_peak_varname,
                                                norm_peak_varname,
                                                condition_varname = "condition #",
                                                reference_varname = "reference condition #",
                                                log2_floor_value = NA) {
  # check for valid inputs

  checkmate::assertString(condition_varname)
  if (!(condition_varname %in% colnames(mzroll_list$samples))) {
    stop(
      "\"reference_varname\":",
      condition_varname,
      ", not present in samples table"
    )
  }

  checkmate::assertString(reference_varname)
  if (!(reference_varname %in% colnames(mzroll_list$samples))) {
    stop(
      "\"reference_varname\":",
      reference_varname,
      ", not present in samples table"
    )
  }

  # check for reference mis-specification

  if (!(
    all(mzroll_list$samples[[reference_varname]] != "") &&
      length(mzroll_list$samples[[reference_varname]]) ==
        nrow(mzroll_list$samples)
  )) {
    stop("not all samples possessed a reference condition")
  }

  missing_references <- setdiff(
    unique(mzroll_list$samples[[reference_varname]]),
    mzroll_list$samples[[condition_varname]]
  )
  if (length(missing_references) != 0) {
    stop(
      length(missing_references),
      " reference conditions were not present as conditions: ",
      paste(missing_references, collapse = ", ")
    )
  }

  stopifnot(length(log2_floor_value) == 1)
  if (!is.na(log2_floor_value)) {
    stopifnot(class(log2_floor_value) == "numeric")
  }

  # summarize reference conditions

  annotated_peaks <- mzroll_list$measurements %>%
    dplyr::left_join(
      mzroll_list$samples %>%
        dplyr::select(
          sampleId,
          !!!rlang::syms(c(condition_varname, reference_varname))
        ),
      by = "sampleId"
    )

  # estimate log abundance for each reference condition per group
  reference_peaks <- annotated_peaks %>%
    dplyr::filter(
      !!rlang::sym(condition_varname) %in%
        unique(mzroll_list$samples[[reference_varname]])
    ) %>%
    dplyr::group_by(groupId, !!rlang::sym(condition_varname)) %>%
    dplyr::summarize(
      .reference = median(!!rlang::sym(quant_peak_varname)),
      .groups = "drop"
    ) %>%
    dplyr::ungroup()

  if (is.na(log2_floor_value)) {
    # use the relative scale
    reference_peaks <- reference_peaks %>%
      dplyr::mutate(.reference_diff = .reference)
  } else {
    # retain the original scale
    reference_peaks <- reference_peaks %>%
      dplyr::group_by(groupId) %>%
      dplyr::mutate(.reference_diff = .reference - mean(.reference)) %>%
      dplyr::ungroup()
  }

  join_by <- rlang::set_names(
    rlang::quo_name(condition_varname),
    rlang::quo_name(reference_varname)
  )

  normalized_peaks <- annotated_peaks %>%
    # add reference to each peak
    dplyr::left_join(reference_peaks, by = c("groupId", join_by)) %>%
    dplyr::mutate(
      !!rlang::sym(norm_peak_varname) :=
        !!rlang::sym(quant_peak_varname) - .reference_diff
    ) %>%
    dplyr::select(
      !!!rlang::syms(c(colnames(mzroll_list$measurements), norm_peak_varname))
    )

  # ensure that normalized values still respect the limit of detection
  reference_refloor <- normalization_refloor(
    normalized_peaks,
    log2_floor_value,
    norm_peak_varname,
    quant_peak_varname
  )

  mzroll_list <- romic::update_tomic(mzroll_list, reference_refloor)

  return(mzroll_list)
}

#' @inheritParams normalize_peaks_batch_center
#' @param weights_tribble a table containing weights and sample variables to
#'   match them to.
#'
#' @rdname normalize_peaks
normalize_peaks_loess <- function(mzroll_list,
                                  quant_peak_varname,
                                  norm_peak_varname,
                                  weights_tribble = NULL,
                                  log2_floor_value = 12) {
  stopifnot("datetime" %in% colnames(mzroll_list$samples))

  stopifnot(length(log2_floor_value) == 1)
  if (!is.na(log2_floor_value)) {
    stopifnot(class(log2_floor_value) == "numeric")
  }

  if (!is.null(weights_tribble)) {
    # check weights_tribble validity
    stopifnot("data.frame" %in% class(weights_tribble))
    stopifnot("weights" %in% colnames(weights_tribble))

    weights_vars <- setdiff(colnames(weights_tribble), "weights")
    stopifnot(weights_vars %in% colnames(mzroll_list$samples))
    stopifnot(length(weights_vars) >= 1)

    sample_weights <- mzroll_list$samples %>%
      dplyr::inner_join(weights_tribble, by = weights_vars) %>%
      dplyr::select(sampleId, datetime, weights)

    if (nrow(sample_weights) == 0) {
      stop(
        "No samples with weights exist
      check configuration of \"weights_tribble\""
      )
    }
  } else {
    sample_weights <- mzroll_list$samples %>%
      dplyr::select(sampleId, datetime) %>%
      dplyr::mutate(weights = 1)
  }

  loess_fits <- mzroll_list$measurements %>%
    dplyr::left_join(sample_weights, by = "sampleId") %>%
    tidyr::nest(groupData = -groupId) %>%
    dplyr::mutate(loess_fits = purrr::map(
      groupData,
      fit_loess,
      quant_peak_varname = quant_peak_varname,
      norm_peak_varname = norm_peak_varname
    )) %>%
    dplyr::select(-groupData) %>%
    tidyr::unnest(loess_fits) %>%
    dplyr::select(!!!rlang::syms(c(
      colnames(mzroll_list$measurements),
      ".loess_fit", norm_peak_varname
    )))

  loess_floor <- normalization_refloor(
    loess_fits,
    log2_floor_value,
    norm_peak_varname,
    quant_peak_varname
  )

  mzroll_list <- romic::update_tomic(mzroll_list, loess_floor)

  return(mzroll_list)
}

fit_loess <- function(groupData, quant_peak_varname, norm_peak_varname) {
  loess_data <- groupData %>%
    dplyr::arrange(datetime) %>%
    dplyr::mutate(
      hours_elapsed = lubridate::interval(datetime[1], datetime) /
        lubridate::hours(1)
    )

  # weights should be inversely proportional to variance

  loess_model <- stats::loess(
    formula = stats::as.formula(
      glue::glue("{quant_peak_varname} ~ hours_elapsed")
    ),
    data = loess_data,
    weights = weights
  )

  loess_apply <- loess_data %>%
    dplyr::mutate(
      .loess_fit = stats::predict(loess_model),
      .loess_shift = .loess_fit - mean(.loess_fit),
      !!rlang::sym(norm_peak_varname) :=
        !!rlang::sym(quant_peak_varname) - .loess_shift
    ) %>%
    dplyr::select(!!!rlang::syms(c(
      colnames(groupData),
      ".loess_fit",
      norm_peak_varname
    )))

  return(loess_apply)
}

normalization_refloor <- function(normalized_peaks,
                                  log2_floor_value,
                                  norm_peak_varname,
                                  quant_peak_varname) {
  # If a floor value is supplied floor a peak to this value if it reached the
  # floor either pre- or post-normalization

  stopifnot(length(log2_floor_value) == 1)
  if (!is.na(log2_floor_value)) {
    # peaks starting or pushed below limit of detection are reset to
    # log2_floor_value
    stopifnot(class(log2_floor_value) == "numeric")

    normalized_peaks <- normalized_peaks %>%
      dplyr::mutate(!!rlang::sym(norm_peak_varname) := dplyr::case_when(
        !!rlang::sym(quant_peak_varname) <=
          log2_floor_value + 0.001 ~ log2_floor_value,
        !!rlang::sym(norm_peak_varname) <=
          log2_floor_value + 0.001 ~ log2_floor_value,
        TRUE ~ !!rlang::sym(norm_peak_varname)
      ))
  }

  return(normalized_peaks)
}

#' @inheritParams normalize_peaks_batch_center
#' @param time_col_varname variable in samples table to use for linear
#' correction
#' @param lm_drift_order polynomial degree for drift, as integer; defaults to
#' 1 (linear drift correction)
#' @param filter_ids sample type ids on which to filter data from column
#' provided by \code{filter_var}; defaults to NULL
#' @param filter_var column name on which to filter samples indicated by
#' \code{filter_ids}; defaults to NULL
#' @param remove_outliers logical that defines whether or not to remove outliers
#' from dataset before quantifying linear drift; defaults to TRUE
#' @param outlier_sd the standard deviation value above or below which samples
#' would be removes as outliers if \code{remove_outliers} is set to TRUE
#'
#' @rdname normalize_peaks
normalize_peaks_lm <- function(mzroll_list,
                               quant_peak_varname,
                               norm_peak_varname,
                               time_col_varname,
                               lm_drift_order = 1,
                               filter_ids = NULL,
                               filter_var = NULL,
                               remove_outliers = TRUE,
                               outlier_sd = 2) {
  stopifnot(time_col_varname %in% colnames(mzroll_list$samples))

  time_col_varname_add <- mzroll_list$samples %>%
    dplyr::select(c("sampleId", time_col_varname))

  lm_fit <- mzroll_list$measurements %>%
    dplyr::left_join(., time_col_varname_add, by = "sampleId") %>%
    tidyr::nest(groupData = -groupId) %>%
    dplyr::mutate(lm_fits = purrr::map(groupData,
      fit_lm,
      quant_peak_varname = quant_peak_varname,
      norm_peak_varname = norm_peak_varname,
      time_col_varname = time_col_varname,
      lm_drift_order = lm_drift_order,
      filter_ids = filter_ids,
      filter_var = filter_var,
      remove_outliers = remove_outliers,
      outlier_sd = outlier_sd
    )) %>%
    dplyr::select(-groupData) %>%
    tidyr::unnest(lm_fits)

  lm_fit_measurements <- lm_fit %>%
    dplyr::select(!!!rlang::syms(c(colnames(mzroll_list$measurements), norm_peak_varname)))

  lm_fit_features <- lm_fit %>%
    dplyr::distinct(groupId, lm_estimate) %>%
    dplyr::left_join(mzroll_list$features, ., by = "groupId")

  mzroll_list <- romic::update_tomic(mzroll_list, lm_fit_measurements)
  mzroll_list <- romic::update_tomic(mzroll_list, lm_fit_features)

  return(mzroll_list)
}

fit_lm <- function(groupData,
                   time_col_varname,
                   quant_peak_varname,
                   norm_peak_varname,
                   lm_drift_order,
                   filter_ids,
                   filter_var,
                   remove_outliers,
                   outlier_sd) {
  # Make data on which to apply drift correction
  lm_data <- groupData %>%
    dplyr::mutate(
      val_var = !!rlang::sym(quant_peak_varname),
      dri_var = !!rlang::sym(time_col_varname)
    )

  # Make data on which to calculate drift slope
  lm_predict <- lm_data %>%
    tidyr::drop_na(val_var)

  # If provided, filter fit dataset to only capture a specific data type
  if (!is.null(filter_ids) && filter_var %in% colnames(lm_predict)) {
    lm_predict <- lm_predict %>%
      dplyr::filter(!!rlang::sym(filter_var) %in% filter_ids)
  }

  # If provided, filter outliers out of the fit values
  if (remove_outliers) {
    outliers <- metstats::find_outliers(
      lm_predict,
      outlier_col_names = "val_var",
      outlier_sd = outlier_sd,
      sample_name_col = "sampleId"
    )

    lm_predict <- lm_predict %>%
      dplyr::filter(!sampleId %in% outliers)
  }

  # Run model
  lm_model <- stats::lm(val_var ~ poly(dri_var, degree = lm_drift_order), data = lm_predict)

  # Compute corrected values and assign new normalized variable norm_peak_varname
  # Predict values from linear fit are subtracted from quant_peak_varname,
  # which centers the data around zero
  # Median value from quant_peak_varname is added back to maintain abundance value
  lm_apply <- lm_data %>%
    dplyr::mutate(median_value = median(.data$val_var, na.rm = T)) %>%
    dplyr::mutate(`:=`(
      !!rlang::sym(norm_peak_varname),
      .data$val_var - .env$predict(lm_model, newdata = lm_data) + .data$median_value
    )) %>%
    dplyr::mutate(lm_estimate = summary(.env$lm_model)$coefficient[2, 1]) %>%
    dplyr::select(-c(val_var, dri_var, median_value))

  return(lm_apply)
}


#' @inheritParams normalize_peaks_batch_center
#'
#' @rdname normalize_peaks
normalize_peaks_center <- function(mzroll_list,
                                   quant_peak_varname,
                                   norm_peak_varname) {
  updated_measurements <- mzroll_list$measurements %>%
    dplyr::group_by(groupId) %>%
    dplyr::mutate(!!rlang::sym(norm_peak_varname) :=
      scale(!!rlang::sym(quant_peak_varname), scale = F, center = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(!!rlang::sym(norm_peak_varname) :=
      as.numeric(!!rlang::sym(norm_peak_varname)))

  mzroll_list <- romic::update_tomic(mzroll_list, updated_measurements)

  return(mzroll_list)
}
