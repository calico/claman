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
  
  romic:::check_triple_omic(mzroll_list, fast_check)

  # check for claman-specific conventions
  
  if (mzroll_list$design$feature_pk != "groupId") {
    stop(glue::glue(
      "The mzroll feature primary key was {mzroll_list$design$feature_pk}
        this value must be \"groupId\""))
  }
  
  if (mzroll_list$design$sample_pk != "sampleId") {
    stop(glue::glue(
      "The mzroll feature primary key was {mzroll_list$design$sample_pk}
        this value must be \"sampleId\""))
  }
  
  # check for required variables
  
  check_required_variables(
    mzroll_list,
    "features", 
    c(
      "groupId",
      "compoundName",
      "smiles",
      "adductName",
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
    c("peakId", "groupId", "sampleId", "log2_abundance")
  )
  
  check_required_variables(
    mzroll_list,
    "samples", 
    c("sampleId", "name", "filename")
  )
  
  # check for invalid variables
  
  unnamed_samples <- mzroll_list$samples %>% dplyr::filter(is.na(name))
  if (nrow(unnamed_samples) > 0){
    stop(glue::glue(
      "{nrow(unnamed_samples)} samples were unnamed. All samples must be named"
      ))
  }
  
  duplicated_names <- mzroll_list$samples %>%
    dplyr::group_by(name) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::distinct(name)
  
  if (nrow(duplicated_names) > 0){
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
      "))
  }
  
  return (invisible(0))
}


#' Paths for X0106 examples
#'
#' methods and paths to X0106 example mzroll datasets.
#'
#' @examples
#' \dontrun{
#' X0106_paths <- tibble::tribble(
#'   ~method_tag, ~mzroll_db_path,
#'   "metabolites (neg)", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M001-peakdetector.mzrollDB",
#'   "metabolites (pos) - 1", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M002-peakdetector.mzrollDB",
#'   "metabolites (pos) - 2", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M002-peakdetector2.mzrollDB",
#'   "lipids (neg) - 1", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M004A-1.mzrollDB",
#'   "lipids (neg) - 2", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M004A-2.mzrollDB",
#'   "lipids (neg) - 3", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M004A-2b-AA-BDP-CL.mzrollDB",
#'   "lipids (neg) - 4", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M004A-2b-AA-BDP-CL+Stds.mzrollDB",
#'   "lipids (pos) - 1", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M005A-peakdetector.mzrollDB",
#'   "lipids (pos) - 2", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M005A-peakdetectorMGDGTGStds.mzrollDB"
#' )
#'
#' usethis::use_data(X0106_paths, overwrite = TRUE)
#' }
"X0106_paths"

#' Flatten mzroll list to single table
#'
#' Flatten mzroll list to single table using groupId and sampleId
#'
#' @param mzroll_list output of process_mzrollDB functions
#' \itemize{
#'   \item{peakgroups: one row per unique analyte (defined by a unique
#'     groupId)},
#'   \item{samples: one row per unique sample (defined by a unique sampleId)},
#'   \item{peaks: one row per peak (samples x peakgroups)}
#'   }
#'
#' @return mzroll_list (process_mzrollDB functions output) reformatted as a
#'   single table
#'
#' @export
mzroll_table <- function(mzroll_list) {

  # guards
  if (class(mzroll_list) != "list") {
    stop("input \"mzroll_list\" must be a list.")
  }
  if (length(mzroll_list) != 3) {
    stop("input \"mzroll_list\" must contain exactly 3 entries.")
  }
  if (names(mzroll_list[1]) != "peakgroups") {
    stop("input \"mzroll_list\" first element must be named \"peakgroups\"")
  }
  if (names(mzroll_list[2]) != "samples") {
    stop("input \"mzroll_list\" second element must be named \"samples\"")
  }
  if (names(mzroll_list[3]) != "peaks") {
    stop("input \"mzroll_list\" third element must be named \"peaks\"")
  }
  if (!"groupId" %in% colnames(mzroll_list[[1]])) {
    stop("peakgroups table missing required column \"groupId\"")
  }
  if (!"groupId" %in% colnames(mzroll_list[[3]])) {
    stop("peaks table missing required column \"groupId\"")
  }
  if (!"sampleId" %in% colnames(mzroll_list[[2]])) {
    stop("samples table missing required column \"sampleId\"")
  }
  if (!"sampleId" %in% colnames(mzroll_list[[3]])) {
    stop("peaks table missing required column \"sampleId\"")
  }

  mzroll_table <- dplyr::inner_join(
    mzroll_list$peakgroups,
    mzroll_list$peaks,
    by = c("groupId")
  ) %>%
    dplyr::inner_join(mzroll_list$samples, by = c("sampleId"))

  mzroll_table
}

#' Util Pretty Knitr Head
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
    slice(1:nrows) %>%
    knitr::kable(caption = caption) %>%
    kableExtra::kable_styling(
      position = "left",
      bootstrap_options = "striped"
    )
}
