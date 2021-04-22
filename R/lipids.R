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