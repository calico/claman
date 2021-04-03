#' Query Systematic Compounds
#'
#' @param systematicCompoundIds NULL to query all compounds or a numeric
#'   vector of systematic compound IDs
#' @param systematic_compounds_con A connection to an SQL database which
#'   stores the structure and systematic IDs of well-studied compounds.
#'
#' @return a list containing:
#' \itemize{
#'   \item{distinct_compounds: unique database entries}
#'   \item{database_ids: systematic IDs of compounds in KEGG, HMDB, ...}
#'   \item{compound_aliases: alternative names for a compound}
#'   \item{metabolic_pathways: KEGG pathways associated with a metabolite}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' systematic_compounds_con <- configure_db_access()$systematic_compounds_con
#' systematic_compounds <- query_systematic_compounds(
#'   systematic_compounds_con = systematic_compounds_con
#' )
#' }
query_systematic_compounds <- function(systematicCompoundIds = NULL,
                                       systematic_compounds_con) {
  stopifnot(all(
    class(systematicCompoundIds) %in% c("NULL", "numeric", "integer")
  ))

  distinct_compounds <- dplyr::tbl(
    systematic_compounds_con,
    "distinct_compounds"
  ) %>%
    dplyr::collect()
  database_ids <- dplyr::tbl(
    systematic_compounds_con,
    "database_ids"
  ) %>%
    dplyr::collect()
  compound_aliases <- dplyr::tbl(
    systematic_compounds_con,
    "compound_aliases"
  ) %>%
    dplyr::collect()
  metabolic_pathways <- dplyr::tbl(
    systematic_compounds_con,
    "metabolic_pathways"
  ) %>%
    dplyr::collect()

  if (is.null(systematicCompoundIds[1])) {
    systematicCompoundIds <- distinct_compounds$compoundId
  }

  list(
    distinct_compounds = distinct_compounds,
    database_ids = database_ids,
    compound_aliases = compound_aliases,
    metabolic_pathways = metabolic_pathways
  ) %>%
    purrr::map(
      .f = function(a_tbl, ids) {
        a_tbl %>%
          dplyr::rename(systematicCompoundId = compoundId) %>%
          dplyr::filter(systematicCompoundId %in% ids)
      },
      ids = systematicCompoundIds
    )
}

#' Summarize compound pathways
#'
#' Reduce compound pathway annotations to a set of core "focus" pathways and
#'   to well-represented (>= \code{min_pw_size} examples) general pathways.
#'
#' @param systematic_compounds The output of
#'   \code{\link{query_systematic_compounds}}.
#' @param min_pw_size Minimum number of examples of metabolites from a pathway
#'   for the pathway not to be filtered.
#' @param focus_pathways A character vector of names of KEGG pathways to be
#'   included as \code{focus_pathway} annotations.
#'
#' @return a tibble of systematicCompoundIds and pathway annotations
#'
#' @export
summarize_compound_pathways <- function(systematic_compounds,
                                        min_pw_size = 5L,
                                        focus_pathways = c(
                                          "Glycolysis / Gluconeogenesis",
                                          "Citrate cycle (TCA cycle)",
                                          "Pentose phosphate pathway",
                                          "Biosynthesis of amino acids",
                                          "Purine metabolism",
                                          "Pyrimidine metabolism"
                                        )) {
  checkmate::assertNumber(min_pw_size, lower = 0)
  checkmate::assertCharacter(focus_pathways)

  represented_pathways <- systematic_compounds$metabolic_pathways %>%
    dplyr::count(pathwayName) %>%
    dplyr::mutate(is_focus_pathway = ifelse(
      pathwayName %in% focus_pathways,
      TRUE,
      FALSE
    )) %>%
    dplyr::filter(n >= min_pw_size | is_focus_pathway)

  # label compounds in focus pathway with most specific annotation

  focus_annotations <- systematic_compounds$metabolic_pathways %>%
    dplyr::inner_join(represented_pathways %>%
      dplyr::filter(is_focus_pathway),
    by = "pathwayName"
    ) %>%
    dplyr::group_by(systematicCompoundId) %>%
    dplyr::arrange(n) %>%
    dplyr::slice(1) %>%
    dplyr::select(systematicCompoundId, focus_pathway = pathwayName) %>%
    dplyr::ungroup()

  general_annotations <- systematic_compounds$metabolic_pathways %>%
    dplyr::inner_join(represented_pathways, by = "pathwayName") %>%
    dplyr::group_by(systematicCompoundId) %>%
    dplyr::arrange(n) %>%
    dplyr::slice(1) %>%
    dplyr::select(systematicCompoundId, pathway = pathwayName) %>%
    dplyr::ungroup()

  general_annotations %>%
    dplyr::left_join(focus_annotations, by = "systematicCompoundId") %>%
    dplyr::mutate(focus_pathway = ifelse(
      is.na(focus_pathway),
      "Other",
      focus_pathway
    ))
}
