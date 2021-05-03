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
#' filter_groupIds(nplug_mzroll_normalized, 1:10)
#'
#' @export
filter_groupIds <- function(mzroll_list, groupIds, invert = FALSE) {
  checkmate::assertNumeric(groupIds)
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
