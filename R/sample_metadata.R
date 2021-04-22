#' Add Samples Table
#'
#' Add a table of sample metadata to an existing mzroll_list
#'
#' @inheritParams test_mzroll_list
#' @param samples_tbl Table of sample metadata
#' @param id_strings one or more variables which will be used to match
#'   sample names. 
#' @param exact if true, an exact match between mzroll names and id_strings
#'   will be found; if false, then substring matches will be used.   
#'   
#' @examples 
#' mzroll_list <- process_mzroll(nplug_mzroll())
#' add_samples_tbl(mzroll_list, nplug_samples, "sample_name")
#' 
#' @export
add_samples_tbl <- function(
  mzroll_list,
  samples_tbl,
  id_strings,
  exact = TRUE
  ) {
  
  checkmate::assertClass(mzroll_list, "tomic")
  checkmate::assertClass(mzroll_list, "mzroll")
  checkmate::assertDataFrame(samples_tbl)
  checkmate::assertCharacter(id_strings)
  purrr::walk(id_strings, checkmate::assertChoice, colnames(samples_tbl))
  checkmate::assertLogical(exact)
  
  # match samples_tbl to mzroll samples
  
  samples_tbl <- samples_tbl %>%
    dplyr::mutate(.entry = 1:dplyr::n())
  
  # read all id_strings - substrings used to match sample names in mzroll
  
  ms_id_strings <- samples_tbl %>%
    dplyr::select(!!!rlang::syms(c(".entry", id_strings))) %>%
    dplyr::mutate_at(dplyr::vars(-.entry), as.character) %>%
    tidyr::gather(id_string_var, id_string_value, -.entry) %>%
    dplyr::filter(!is.na(id_string_value))
  
  sample_ms_id_string_matches <- mzroll_list$samples %>%
    tidyr::crossing(ms_id_strings)
  
  if (exact) {
    sample_ms_id_string_matches <- sample_ms_id_string_matches %>%
      dplyr::filter(name == id_string_value)
  } else {
    sample_ms_id_string_matches <- sample_ms_id_string_matches %>%
      dplyr::filter(stringr::str_detect(name, id_string_value))
  }
  
  # check whether 1 sample matches 2+ IDs
  
  sample_multimatch <- sample_ms_id_string_matches %>%
    dplyr::count(name) %>%
    dplyr::filter(n > 1)
  
  if (nrow(sample_multimatch) != 0) {
    multimatch_details <- sample_ms_id_string_matches %>%
      dplyr::semi_join(sample_multimatch, by = "name") %>%
      dplyr::arrange(name) %>%
      dplyr::mutate(out_string = glue::glue(
        "name: {name} matching row: {.entry}, id string column: {id_string_var}, id string entry: {id_string_value}"
      )) %>%
      {
        paste(.$out_string, collapse = "\n")
      }
    
    stop(
      nrow(sample_multimatch),
      " experimental samples matched multiple ID strings: ",
      paste(sample_multimatch$name, collapse = ", "),
      "\nDetails:\n",
      multimatch_details
    )
  }
  
  # check whether 1 ID matches 2+ samples
  
  condition_multimatch <- sample_ms_id_string_matches %>%
    dplyr::count(.entry) %>%
    dplyr::filter(n > 1)
  
  if (nrow(condition_multimatch) != 0) {
    multimatch_details <- sample_ms_id_string_matches %>%
      dplyr::semi_join(condition_multimatch, by = ".entry") %>%
      dplyr::arrange(.entry) %>%
      dplyr::mutate(out_string = glue::glue(
        "row: {.entry} matching name_match: {name} with id string column: {id_string_var}, id string entry: {id_string_value}"
      )) %>%
      {
        paste(.$out_string, collapse = "\n")
      }
    
    stop(
      nrow(condition_multimatch),
      " rows in sample_tbl matched multiple experimental samples: ",
      paste(condition_multimatch$.entry, collapse = ", "),
      "\nDetails:\n",
      multimatch_details
    )
  }
  
  matched_augmented_samples <- mzroll_list$samples %>%
    dplyr::inner_join(
      sample_ms_id_string_matches %>%
        dplyr::select(sampleId, .entry) %>%
        dplyr::left_join(samples_tbl, by = ".entry"),
      by = "sampleId"
    )
  
  unmatched_samples <- mzroll_list$samples %>%
    dplyr::anti_join(sample_ms_id_string_matches, by = "sampleId")
  
  if (nrow(unmatched_samples) != 0) {
    warning(
      nrow(unmatched_samples),
      " experimental samples were not matched to ID strings. Their measurements will be discarded.:\n  ",
      paste(unmatched_samples$name, collapse = "\n  ")
    )
  }
  
  mzroll_list <- romic::update_tomic(
    mzroll_list,
    matched_augmented_samples
    )
  
  return(mzroll_list)
}

#' Remove Constant Name
#'
#' Taking a set of filenames remove the constant leading &/or lagging portion
#'   of the names.
#'
#' @param name_set a character vector of names with some common leading or
#'   lagging substrings.
#'
#' @return name_set with constant substring removed.
#'
#' @examples
#' remove_constant_name(c("xxxxsample1yyyy", "xxxxcontrol1yyyy"))
#' @export
remove_constant_name <- function(name_set) {
  stopifnot(all(class(name_set) %in% c("character", "factor")))

  if (length(name_set) < 2) {
    stop(
      "only ",
      length(name_set),
      " names provided; this function can't determine the constant and
        variable porition of filenames"
    )
  }

  distinct_names <- tibble::tibble(name = name_set) %>%
    dplyr::distinct(name)

  name_set_tibble <- distinct_names %>%
    dplyr::mutate(name_tibble = purrr::map(name, tibble_name)) %>%
    tidyr::unnest(name_tibble)

  leading_constant_end <- name_set_tibble %>%
    dplyr::count(char, position) %>%
    dplyr::filter(n != nrow(distinct_names)) %>%
    dplyr::arrange(position) %>%
    {
      .$position[1] - 1
    }

  lagging_constant_r_end <- name_set_tibble %>%
    dplyr::count(char, r_position) %>%
    dplyr::filter(n != nrow(distinct_names)) %>%
    dplyr::arrange(r_position) %>%
    {
      .$r_position[1] - 1
    }

  name_reductions <- name_set_tibble %>%
    dplyr::filter(
      position > leading_constant_end,
      r_position > lagging_constant_r_end
    ) %>%
    dplyr::group_by(name) %>%
    dplyr::arrange(position) %>%
    dplyr::summarize(reduced_name = paste(char, collapse = "")) %>%
    dplyr::ungroup()

  tibble::tibble(name = name_set) %>%
    dplyr::left_join(name_reductions, by = "name") %>%
    {
      .$reduced_name
    }
}

tibble_name <- function(name) {
  tibble::tibble(char = strsplit(name, split = "")[[1]]) %>%
    dplyr::mutate(
      position = 1:dplyr::n(),
      r_position = dplyr::n():1
    )
}
