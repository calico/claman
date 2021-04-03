#' Filter GroupIds
#'
#' Remove all data besides those from provided groupIds
#'
#' @inheritParams test_mzroll_list
#' @param groupIds groupIds to retain
#' @param invert if TRUE then remove provided groupIds; if FALSE then retain
#'   provided groupIds
#'
#' @export
filter_groupIds <- function(mzroll_list, groupIds, invert = FALSE) {
  checkmate::assertNumeric(groupIds)
  checkmate::assertLogical(invert, len = 1)

  if (invert) {
    filter_groupIds <- setdiff(mzroll_list$peakgroups$groupId, groupIds)
  } else {
    filter_groupIds <- groupIds
  }

  mzroll_list <- filter_peakgroups_quo(
    mzroll_list,
    rlang::quo(groupId %in% filter_groupIds)
  )

  return(mzroll_list)
}

#' Filter Calicomics Samples
#'
#' Remove samples based on sample table variables
#'
#' @inheritParams test_mzroll_list
#' @param filter_quosure quosure to use to filter samples
#'
#' @return an \code{mzroll_list} with samples (and corresponding peaks) removed
#'
#' @examples
#' mzroll_list <- readRDS(authutils::get_clamr_assets("X0083-mzroll-list.Rds"))
#' filter_quosure <- rlang::quo(sex == "Male")
#' filter_samples_quo(mzroll_list, filter_quosure)
#' @export
filter_samples_quo <- function(mzroll_list, filter_quosure) {
  test_mzroll_list(mzroll_list)
  stopifnot("formula" %in% class(filter_quosure))

  mzroll_list$samples <- mzroll_list$samples %>%
    dplyr::filter(!!filter_quosure)

  reconciled_list <- reconcile_mzroll_list(mzroll_list)

  reconciled_list
}

#' Filter Calicomics Peakgroups
#'
#' Remove peakgroups based on peakgroups table variables
#'
#' @inheritParams test_mzroll_list
#' @param filter_quosure quosure to use to filter peakgroups
#'
#' @return an \code{mzroll_list} with peakgroups (and corresponding peaks)
#'   removed
#'
#' @examples
#' mzroll_list <- readRDS(authutils::get_clamr_assets("X0083-mzroll-list.Rds"))
#' filter_quosure <- rlang::quo(focus_pathway %in% "Pyrimidine metabolism")
#' filter_peakgroups_quo(mzroll_list, filter_quosure)
#' @export
filter_peakgroups_quo <- function(mzroll_list, filter_quosure) {
  test_mzroll_list(mzroll_list)
  stopifnot("formula" %in% class(filter_quosure))

  mzroll_list$peakgroups <- mzroll_list$peakgroups %>%
    dplyr::filter(!!filter_quosure)

  reconciled_list <- reconcile_mzroll_list(mzroll_list)

  reconciled_list
}
#' Reconcile Mzroll List
#'
#' If some samples, peakgroups or peaks have been dropped; peaks and peakgroups
#'   which may now be irrelevant
#'
#' @inheritParams test_mzroll_list
reconcile_mzroll_list <- function(mzroll_list) {
  mzroll_list$peaks <- mzroll_list$peaks %>%
    dplyr::semi_join(mzroll_list$samples, by = "sampleId") %>%
    dplyr::semi_join(mzroll_list$peakgroups, by = "groupId")

  mzroll_list$peakgroups <- mzroll_list$peakgroups %>%
    dplyr::semi_join(mzroll_list$peaks, by = "groupId")

  mzroll_list
}
