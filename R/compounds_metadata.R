#' Merge Compounds Table
#'
#' Merge a table of compound information with an existing mzroll_list
#'
#' @inheritParams test_mzroll_list
#' @param compounds_tbl Table of compound metadata
#' @param join_by variable shared by mzroll peakgroups and compounds_tbl
#'   to merge results by.
#'   
#' @examples 
#' mzroll_list <- process_mzroll(nplug_mzroll())
#' compounds_tbl <- nplug_compounds
#' merge_compounds_tbl(mzroll_list, compounds_tbl)
#' 
#' @export
merge_compounds_tbl <- function(
  mzroll_list,
  compounds_tbl,
  join_by = "compoundName"
) {
  
  checkmate::assertClass(mzroll_list, "tomic")
  checkmate::assertClass(mzroll_list, "mzroll")
  checkmate::assertDataFrame(compounds_tbl)
  checkmate::assertString(join_by)
  
  # ensure that the join_by variable is unique in compounds_tbl
  # (it doesn't have to be in the features table of the mzroll_list)
  duplicated_compounds <- compounds_tbl %>%
    dplyr::group_by(!!rlang::sym(join_by)) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::distinct(!!rlang::sym(join_by))
  
  if (nrow(duplicated_compounds) > 0) {
    stop(glue::glue(
      "{nrow(duplicated_compounds)} {join_by} entries were not unique:
        {paste(duplicated_compounds[[join_by]], collapse = ', ')}")
      )
    }
  
  unmatched_features <- mzroll_list$features %>%
    dplyr::anti_join(compounds_tbl, by = join_by)
  
  if (nrow(unmatched_features) == nrow(mzroll_list$features)) {
    stop(glue::glue(
      "No compounds could be matched to mzroll features using {join_by}"
      ))
  }
  if (nrow(unmatched_features) > 0) {
    warning(glue::glue(
      "{nrow(unmatched_features)} features did not match entries
        in the compound_tbl"
      ))
  }
  
  augmented_mzroll_list <- romic::update_tomic(
    mzroll_list,
    mzroll_list$features %>%
      dplyr::left_join(compounds_tbl, by = join_by)
    )
  
  return(augmented_mzroll_list)
}