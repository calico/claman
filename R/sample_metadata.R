#' Merge Samples Table
#'
#' Merge a table of sample metadata with an existing mzroll_list
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
#' merge_samples_tbl(mzroll_list, nplug_samples, "sample_name", exact = TRUE)
#' @export
merge_samples_tbl <- function(mzroll_list,
                              samples_tbl,
                              id_strings,
                              exact = FALSE) {
  checkmate::assertClass(mzroll_list, "tomic")
  checkmate::assertClass(mzroll_list, "mzroll")
  checkmate::assertDataFrame(samples_tbl)

  if ("samples_tbl_row" %in% colnames(samples_tbl)) {
    stop(
      "\"samples_tbl_row\" is a reserved variable in samples_tbl,
        do not include it"
    )
  }

  checkmate::assertCharacter(id_strings)
  purrr::walk(id_strings, checkmate::assertChoice, colnames(samples_tbl))
  checkmate::assertLogical(exact)

  # match samples_tbl to mzroll samples
  samples_tbl <- samples_tbl %>%
    dplyr::mutate(samples_tbl_row = 1:dplyr::n())

  # read all id_strings - substrings used to match sample names in mzroll
  ms_id_strings <- samples_tbl %>%
    dplyr::select(!!!rlang::syms(c("samples_tbl_row", id_strings))) %>%
    dplyr::mutate_at(dplyr::vars(-samples_tbl_row), as.character) %>%
    tidyr::gather(id_string_var, id_string_value, -samples_tbl_row) %>%
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

  ## claman Issue-26: fix mis-match issues
  # first, check if any string ids are present in another string id
  id_internal_match <- expand.grid(
    ms_id_strings$id_string_value,
    ms_id_strings$id_string_value
  ) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), as.character)) %>%
    dplyr::filter(
      stringr::str_detect(Var1, Var2),
      Var1 != Var2
    )

  # if there are rows in id_internal_match, we need to filter out rows where
  # id_string_value == Var2 matches the same name that is being matched by
  # id_string_value == Var1
  if (nrow(id_internal_match) > 0) {
    # find problematic rows
    problematic_matches <- sample_ms_id_string_matches %>%
      dplyr::inner_join(id_internal_match,
        by = c("id_string_value" = "Var2")
      ) %>%
      dplyr::filter(stringr::str_detect(name, Var1))

    # remove problematic rows from main matching dataframe
    sample_ms_id_string_matches <- sample_ms_id_string_matches %>%
      dplyr::anti_join(problematic_matches,
        by = c("samples_tbl_row", "name", "id_string_value")
      )
  }

  # now, proceed with multi-match checking as before

  # check whether 1 sample matches 2+ IDs
  sample_multimatch <- sample_ms_id_string_matches %>%
    dplyr::count(name) %>%
    dplyr::filter(n > 1)

  if (nrow(sample_multimatch) != 0) {
    multimatch_details <- sample_ms_id_string_matches %>%
      dplyr::semi_join(sample_multimatch, by = "name") %>%
      dplyr::arrange(name) %>%
      dplyr::mutate(out_string = glue::glue(
        "name: {name} matching row: {samples_tbl_row}, id string column: {id_string_var}, id string entry: {id_string_value}"
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
    dplyr::count(samples_tbl_row) %>%
    dplyr::filter(n > 1)

  if (nrow(condition_multimatch) != 0) {
    multimatch_details <- sample_ms_id_string_matches %>%
      dplyr::semi_join(condition_multimatch, by = "samples_tbl_row") %>%
      dplyr::arrange(samples_tbl_row) %>%
      dplyr::mutate(out_string = glue::glue(
        "row: {samples_tbl_row} matching name_match: {name} with id string column: {id_string_var}, id string entry: {id_string_value}"
      )) %>%
      {
        paste(.$out_string, collapse = "\n")
      }

    stop(
      nrow(condition_multimatch),
      " rows in sample_tbl matched multiple experimental samples: ",
      paste(condition_multimatch$samples_tbl_row, collapse = ", "),
      "\nDetails:\n",
      multimatch_details
    )
  }

  matched_augmented_samples <- mzroll_list$samples %>%
    dplyr::inner_join(
      sample_ms_id_string_matches %>%
        dplyr::select(sampleId, samples_tbl_row) %>%
        dplyr::left_join(samples_tbl, by = "samples_tbl_row"),
      by = "sampleId"
    )

  unmatched_samples <- mzroll_list$samples %>%
    dplyr::anti_join(sample_ms_id_string_matches, by = "sampleId")

  if (nrow(unmatched_samples) != 0) {
    warning(
      nrow(unmatched_samples),
      " experimental samples were not matched to ID strings. Their measurements will be discarded.:\n  ",
      paste(unmatched_samples$name, collapse = "\n  "),
      "\n"
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

#' Generate relevant lipid components from compound name
#'
#' Produce data frame with lipid components from name
#'
#' @param compound_name vector of compound names
#'
#' @return compound_name vector expanded to large table, with structural
#'   components described as their own columns.
#'
#' @export
lipid_components <- function(compound_name) {
  compound_components <- tibble::tibble("compoundName" = compound_name) %>%
    # general compound information
    dplyr::mutate(
      lipidClass = stringr::str_extract(compoundName, "^.*(?=\\()")
    ) %>%
    dplyr::mutate(
      compound_name_adduct = stringr::str_extract(compoundName, "(?<=\\) ).*$")
    ) %>%
    dplyr::mutate(
      plasmalogen_type = stringr::str_extract(compoundName, "[op]-(?=[0-9]+:)")
    ) %>%
    dplyr::mutate(lipidClass_o_p = ifelse(
      is.na(plasmalogen_type),
      lipidClass,
      paste0(plasmalogen_type, lipidClass)
    )) %>%
    dplyr::mutate(lipidClass_plas = ifelse(
      is.na(plasmalogen_type),
      lipidClass,
      paste0("plas-", lipidClass)
    )) %>%
    # Parse all chains.  Also retrieves number of hydroxyl groups on each
    #   FA (m=1, d=2, t=3), if available
    dplyr::mutate(sn_chains = stringr::str_extract_all(
      compoundName,
      "(?<=[\\(/_])[A-Za-z]?-?[0-9]+:[0-9]+;?O?[0-9]?(?=[\\)/_])"
    )) %>%
    dplyr::mutate(num_sn_chains = sapply(sn_chains, function(x) {
      length(x)
    })) %>%
    dplyr::mutate(num_unique_sn_chains = sapply(sn_chains, function(x) {
      x %>%
        unique() %>%
        length()
    })) %>%
    # sn1
    dplyr::mutate(sn1 = sapply(sn_chains, function(x) {
      x[1]
    })) %>%
    dplyr::mutate(FA1 = stringr::str_extract(sn1, "[0-9]+:[0-9]+")) %>%
    dplyr::mutate(single_1 = as.numeric(
      stringr::str_extract(FA1, ".*(?=:)")
    )) %>%
    dplyr::mutate(single_1 = ifelse(is.na(single_1), 0, single_1)) %>%
    dplyr::mutate(double_1 = as.numeric(
      stringr::str_extract(FA1, "(?<=:).*")
    )) %>%
    dplyr::mutate(double_1 = ifelse(is.na(double_1), 0, double_1)) %>%
    dplyr::mutate(double_1 = ifelse(
      (!is.na(plasmalogen_type) & plasmalogen_type == "p-"),
      double_1 + 1,
      double_1
    )) %>%
    dplyr::mutate(OH_1 = dplyr::case_when(
      grepl("^m", sn1) ~ 1,
      grepl("^d", sn1) ~ 2,
      grepl("^t", sn1) ~ 3,
      grepl(";O4", sn1) ~ 4,
      grepl(";O3", sn1) ~ 3,
      grepl(";O2", sn1) ~ 2,
      grepl(";O1", sn1) ~ 1,
      grepl(";O", sn1) ~ 1,
      TRUE ~ 0
    )) %>%
    # sn2
    dplyr::mutate(sn2 = sapply(sn_chains, function(x) {
      x[2]
    })) %>%
    dplyr::mutate(FA2 = stringr::str_extract(sn2, "[0-9]+:[0-9]+")) %>%
    dplyr::mutate(single_2 = as.numeric(
      stringr::str_extract(FA2, ".*(?=:)")
    )) %>%
    dplyr::mutate(single_2 = ifelse(is.na(single_2), 0, single_2)) %>%
    dplyr::mutate(double_2 = as.numeric(
      stringr::str_extract(FA2, "(?<=:).*")
    )) %>%
    dplyr::mutate(double_2 = ifelse(is.na(double_2), 0, double_2)) %>%
    dplyr::mutate(OH_2 = dplyr::case_when(
      grepl("^m", sn2) ~ 1,
      grepl("^d", sn2) ~ 2,
      grepl("^t", sn2) ~ 3,
      grepl(";O4", sn2) ~ 4,
      grepl(";O3", sn2) ~ 3,
      grepl(";O2", sn2) ~ 2,
      grepl(";O1", sn2) ~ 1,
      grepl(";O", sn2) ~ 1,
      TRUE ~ 0
    )) %>%
    # sn3 (for TGs)
    dplyr::mutate(sn3 = sapply(sn_chains, function(x) {
      x[3]
    })) %>%
    dplyr::mutate(FA3 = stringr::str_extract(sn3, "[0-9]+:[0-9]+")) %>%
    dplyr::mutate(single_3 = as.numeric(
      stringr::str_extract(FA3, ".*(?=:)")
    )) %>%
    dplyr::mutate(single_3 = ifelse(is.na(single_3), 0, single_3)) %>%
    dplyr::mutate(double_3 = as.numeric(
      stringr::str_extract(FA3, "(?<=:).*")
    )) %>%
    dplyr::mutate(double_3 = ifelse(is.na(double_3), 0, double_3)) %>%
    dplyr::mutate(OH_3 = dplyr::case_when(
      grepl("^m", sn3) ~ 1,
      grepl("^d", sn3) ~ 2,
      grepl("^t", sn3) ~ 3,
      grepl(";O4", sn3) ~ 4,
      grepl(";O3", sn3) ~ 3,
      grepl(";O2", sn3) ~ 2,
      grepl(";O1", sn3) ~ 1,
      grepl(";O", sn3) ~ 1,
      TRUE ~ 0
    )) %>%
    # sn4 (for Cardiolipins CLs)
    dplyr::mutate(sn4 = sapply(sn_chains, function(x) {
      x[4]
    })) %>%
    dplyr::mutate(FA4 = stringr::str_extract(sn4, "[0-9]+:[0-9]+")) %>%
    dplyr::mutate(single_4 = as.numeric(
      stringr::str_extract(FA4, ".*(?=:)")
    )) %>%
    dplyr::mutate(single_4 = ifelse(is.na(single_4), 0, single_4)) %>%
    dplyr::mutate(double_4 = as.numeric(
      stringr::str_extract(FA4, "(?<=:).*")
    )) %>%
    dplyr::mutate(double_4 = ifelse(is.na(double_4), 0, double_4)) %>%
    dplyr::mutate(OH_4 = dplyr::case_when(
      grepl("^m", sn4) ~ 1,
      grepl("^d", sn4) ~ 2,
      grepl("^t", sn4) ~ 3,
      grepl(";O4", sn4) ~ 4,
      grepl(";O3", sn4) ~ 3,
      grepl(";O2", sn4) ~ 2,
      grepl(";O1", sn4) ~ 1,
      grepl(";O", sn4) ~ 1,
      TRUE ~ 0
    )) %>%
    # finally, put together all info from chains to describe summed
    #   composition. sometimes, the compound is already reported as a summed
    #   composition. In that case, return the original name.
    dplyr::mutate(total_single = ifelse(
      num_sn_chains == 0,
      as.numeric(stringr::str_extract(
        compoundName,
        "(?<=\\()[A-Za-z]?-?[0-9]+(?=:[0-9]+,)"
      )),
      (single_1 + single_2 + single_3 + single_4)
    )) %>%
    dplyr::mutate(total_single = ifelse(
      is.na(total_single),
      0,
      total_single
    )) %>%
    dplyr::mutate(total_double = ifelse(
      num_sn_chains == 0,
      as.numeric(stringr::str_extract(compoundName, "(?<=:)[0-9]+(?=,)")),
      (double_1 + double_2 + double_3 + double_4)
    )) %>%
    dplyr::mutate(total_double = ifelse(
      is.na(total_double),
      0,
      total_double
    )) %>%
    dplyr::mutate(total_OH = ifelse(
      num_sn_chains == 0,
      as.numeric(stringr::str_extract(compoundName, "(?<=,)[0-9]+(?=-)")),
      OH_1 + OH_2 + OH_3 + OH_4
    )) %>%
    dplyr::mutate(total_OH = ifelse(is.na(total_OH), 0, total_OH)) %>%
    dplyr::mutate(sumComposition = dplyr::case_when(
      num_sn_chains > 0 & total_OH == 0 ~
        glue::glue("{class}({plas}{single}:{double})",
          class = lipidClass,
          plas = ifelse(is.na(plasmalogen_type), "", plasmalogen_type),
          single = total_single,
          double = ifelse(
            (!is.na(plasmalogen_type) & plasmalogen_type == "p-"),
            (total_double - 1),
            total_double
          )
        ),
      num_sn_chains > 0 & total_OH == 1 ~
        glue::glue(
          "{class}({plas}{single}:{double};O)",
          class = lipidClass,
          plas = ifelse(is.na(plasmalogen_type), "", plasmalogen_type),
          single = total_single,
          double = ifelse(
            (!is.na(plasmalogen_type) & plasmalogen_type == "p-"),
            (total_double - 1),
            total_double
          )
        ),
      num_sn_chains > 0 & total_OH > 1 ~
        glue::glue(
          "{class}({plas}{single}:{double};O{num_OH})",
          class = lipidClass,
          plas = ifelse(is.na(plasmalogen_type), "", plasmalogen_type),
          single = total_single,
          double = ifelse(
            (!is.na(plasmalogen_type) & plasmalogen_type == "p-"),
            (total_double - 1),
            total_double
          ),
          num_OH = total_OH
        ),
      num_sn_chains == 0 ~ compoundName
    )) %>%
    dplyr::mutate(etherPlasmalogenCompoundName = ifelse(
      is.na(plasmalogen_type),
      compoundName,
      glue::glue(
        "{class}({single}:{double}e/{sn2_chain})",
        class = lipidClass,
        single = single_1,
        double = double_1,
        sn2_chain = FA2
      )
    )) %>%
    dplyr::mutate(etherPlasmalogenSumComposition = ifelse(
      is.na(plasmalogen_type),
      sumComposition,
      glue::glue(
        "{class}({single}:{double}e)",
        class = lipidClass,
        single = total_single,
        double = total_double
      )
    ))

  compound_components
}
