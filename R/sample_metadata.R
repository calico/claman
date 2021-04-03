#' Find Tracking Sheet
#'
#' Locate the meta data tracking sheet for the project
#'
#' @inheritParams stringr::str_detect
#'
#' @return tracking_sheet_id: if a unique tracking sheet is found, then
#'   return the googlesheets id
#'
#' @examples
#' \dontrun{
#' find_tracking_sheet(pattern = "X0083")
#' }
#'
#' @export
find_tracking_sheet <- function(pattern) {

  # search in Metabolomics Core > Experiments > Project ID
  experiments <- googledrive::drive_ls(
    googledrive::as_id("0B6szsSC6az3YN0dhZDZ4bzFLM00")
  )

  experiments_folder <- experiments %>%
    dplyr::filter(stringr::str_detect(name, pattern))

  if (nrow(experiments_folder) > 1) {
    stop(
      nrow(reduced_tracking_sheets),
      ' experiments folders matched "pattern": ',
      paste(reduced_tracking_sheets$name, collapse = ", ")
    )
  } else if (nrow(experiments_folder) == 1) {
    # check experiment folder for .gsheet
    experiment_tracking_sheet <- googledrive::as_id(experiments_folder$id) %>%
      googledrive::drive_ls() %>%
      # filter to just googlesheets with relevant worksheets
      dplyr::mutate(is_tracking_sheet = purrr::map_lgl(
        id,
        identify_tracking_sheet
      )) %>%
      dplyr::filter(is_tracking_sheet)

    if (nrow(experiment_tracking_sheet) > 1) {
      # too many sheets found
      stop(
        nrow(experiment_tracking_sheet),
        ' tracking sheets found in folder matching "pattern": ',
        paste(reduced_tracking_sheets$name, collapse = ", ")
      )
    } else if (nrow(experiment_tracking_sheet) == 1) {
      # return ID uniquely found in folder
      return(experiment_tracking_sheet$id)
    } else {
      # continue looking for sheets
      message(
        "No experiments folder matching pattern in Experiments folder
          checking in sample sheets"
      )
    }
  } else {
    # continue looking for sheets
    message(
      "No experiments folder matching pattern in Experiments folder,
        checking in sample sheets"
    )
  }

  # search in Metabolomics Core > Experimental Tracking Spreadsheets
  tracking_sheets <- googledrive::drive_ls(
    googledrive::as_id("1PSk5g-sYkwT5jMTOQRO9BUozoruAJJGI")
  )

  reduced_tracking_sheets <- tracking_sheets %>%
    dplyr::filter(stringr::str_detect(name, pattern))

  if (nrow(reduced_tracking_sheets) != 1) {
    stop(
      nrow(reduced_tracking_sheets),
      ' tracking sheets matched "pattern :"',
      paste(reduced_tracking_sheets$name, collapse = ", ")
    )
  } else {
    reduced_tracking_sheets$id
  }
}

identify_tracking_sheet <- function(id) {
  id_sheets <- try(googlesheets4::sheets_sheets(id), silent = TRUE)

  if (class(id_sheets) == "try-error") {
    return(FALSE)
  }

  if (all(
    c("USER SUBMISSION", "USER SAMPLE LIST", "Metabolomics User Sample List")
    %in% id_sheets
  )) {
    TRUE
  } else {
    FALSE
  }
}

#' Read Sample List
#'
#' @param tracking_sheet_id output of \link{find_tracking_sheet}
read_sample_list <- function(tracking_sheet_id) {
  sample_list <- googlesheets4::read_sheet(
    tracking_sheet_id,
    sheet = "Metabolomics User Sample List"
  )

  # require a few standard fields in the Metabolomics User Sample List
  required_sample_list_vars <- c(
    "tube label",
    "sample description",
    "condition #",
    "reference condition #",
    "MS ID string",
    "MS ID string alternative"
  )
  missing_required_vars <- setdiff(
    required_sample_list_vars,
    colnames(sample_list)
  )
  if (length(missing_required_vars) != 0) {
    stop(
      "\"Metabolomics User Sample List\" is missing required variables: ",
      paste(missing_required_vars, collapse = ", ")
    )
  }

  # check for tube label uniqueness

  duplicated_tube_labels <- sample_list %>%
    dplyr::count(`tube label`) %>%
    dplyr::filter(n > 1)

  if (nrow(duplicated_tube_labels) != 0) {
    stop(
      nrow(duplicated_tube_labels),
      " tube labels were not unique, duplicated labels: ",
      paste(duplicated_tube_labels$`tube label`, collapse = ", ")
    )
  }

  # check for uniqueness of every MS ID string and alternatives

  duplicated_id_strings <- sample_list %>%
    dplyr::select(`MS ID string`, `MS ID string alternative`) %>%
    tidyr::gather("ID string column", "ID string entry") %>%
    dplyr::filter(!is.na(`ID string entry`)) %>%
    dplyr::count(`ID string entry`) %>%
    dplyr::filter(n > 1)

  if (nrow(duplicated_id_strings) != 0) {
    stop(
      nrow(duplicated_id_strings),
      " MS ID strings were not unique: ",
      paste(duplicated_id_strings$`ID string entry`, collapse = ", ")
    )
  }

  # check that all reference conditions exist

  if (
    !all(sample_list$`reference condition #` %in% sample_list$`condition #`)
  ) {
    stop("some reference condition #s not found as condition #s")
  }

  sample_list
}

#' Import Sample List
#'
#' Find, read and process a sample meta sheet from googlesheets
#'
#' @inheritParams find_tracking_sheet
#'
#' @examples
#' \dontrun{
#' import_sample_sheet(pattern = "X0083")
#' }
#'
#' @export
import_sample_sheet <- function(pattern) {
  tracking_sheet_id <- find_tracking_sheet(pattern = pattern)

  sample_sheet <- read_sample_list(tracking_sheet_id)

  list(
    tracking_sheet_id = tracking_sheet_id,
    sample_sheet = sample_sheet
  )
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
