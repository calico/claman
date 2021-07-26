#' Test MzRoll List
#'
#' @param mzroll_list output of \link{process_mzroll} or
#'   \link{process_mzroll_multi}
#'
#' \itemize{
#'   \item{features: one row per unique analyte (defined by a
#'     unique groupId)},
#'   \item{samples: one row per unique sample (defined by a unique sampleId)},
#'   \item{measurements: one row per peak (samples x peakgroups)}
#'   }
#'
#' @inheritParams romic:::check_triple_omic
#'
test_mzroll_list <- function(mzroll_list, fast_check = TRUE) {
  checkmate::assertClass(mzroll_list, "tomic")
  checkmate::assertClass(mzroll_list, "mzroll")

  # check that mzroll_list is a valid tomic

  romic::check_tomic(mzroll_list, fast_check)

  # check for claman-specific conventions

  if (mzroll_list$design$feature_pk != "groupId") {
    stop(glue::glue(
      "The mzroll feature primary key was {mzroll_list$design$feature_pk}
        this value must be \"groupId\""
    ))
  }

  if (mzroll_list$design$sample_pk != "sampleId") {
    stop(glue::glue(
      "The mzroll feature primary key was {mzroll_list$design$sample_pk}
        this value must be \"sampleId\""
    ))
  }

  # check for required variables

  check_required_variables(
    mzroll_list,
    "features",
    c(
      "groupId",
      "compoundName",
      "smiles",
      "tagString",
      "mz",
      "rt",
      "compoundDB",
      "searchTableName",
      "label"
    )
  )

  check_required_variables(
    mzroll_list,
    "measurements",
    c("groupId", "sampleId", "log2_abundance")
  )

  # the only required field is sampleId, other field will likely
  #   be discarded as samples are merged during normalization
  check_required_variables(mzroll_list, "samples", c("sampleId", "name"))

  # check for invalid variables

  checkmate::assertClass(mzroll_list[["features"]]$groupId, "factor")
  checkmate::assertClass(mzroll_list[["samples"]]$sampleId, "factor")

  unnamed_samples <- mzroll_list$samples %>% dplyr::filter(is.na(name))
  if (nrow(unnamed_samples) > 0) {
    stop(glue::glue(
      "{nrow(unnamed_samples)} samples were unnamed. All samples must be named"
    ))
  }

  duplicated_names <- mzroll_list$samples %>%
    dplyr::group_by(name) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::distinct(name)

  if (nrow(duplicated_names) > 0) {
    stop(glue::glue("{nrow(duplicated_names)} sample names were duplicated"))
  }

  return(invisible(0))
}

check_required_variables <- function(mzroll_list, table, required_variables) {
  checkmate::assertClass(mzroll_list, "tomic")
  checkmate::assertClass(mzroll_list, "mzroll")
  checkmate::assertChoice(table, c("features", "samples", "measurements"))
  checkmate::assertCharacter(required_variables)

  missing_measurements <- setdiff(
    required_variables,
    mzroll_list$design[[table]]$variable
  )

  if (length(missing_measurements) != 0) {
    stop(glue::glue(
      "required variable(s) {paste(missing_measurements, collapse = ', ')}
        missing from {table}
      "
    ))
  }

  return(invisible(0))
}


#' Util - Pretty Knitr Head
#'
#' @param tbl a data.frame or tibble
#' @param nrows the max number of rows to show
#' @inheritParams knitr::kable
#'
#' @return an html knitr table
#'
#' @export
#'
#' @examples
#' util_pretty_khead(mtcars, nrows = 5, caption = "cars!")
util_pretty_khead <- function(tbl, nrows = 10, caption = NULL) {
  checkmate::assertDataFrame(tbl)
  checkmate::assertNumber(nrows, lower = 1)

  tbl %>%
    dplyr::slice(1:nrows) %>%
    knitr::kable(caption = caption) %>%
    kableExtra::kable_styling(
      position = "left",
      bootstrap_options = "striped"
    )
}
