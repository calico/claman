#' Add Samples Table
#'
#' Add a table of sample metadata to an existing mzroll_list
#'
#' @param samples_tbl Table of sample metadata
add_samples_tbl <- function(mzroll_list, samples_tbl) {
  
  
  
  
  samples_tbl

  sample_list
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
