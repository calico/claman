#' Find Pathway Enrichments
#'
#' @inheritParams test_mzroll_list
#' @inheritParams plot_volcano
#' @param pathway_list a named list where names are pathways and entries
#'   are groupIds belonging to each pathway
#' @param ranking_measure variable in \code{regression_significance} to use
#'   when calculating enrichment.
#' @param test_absolute_effects If TRUE then only consider the magnitude of
#'   test-statistics when calculating rankings for enrichment; if FALSE then
#'   consider sign as well.
#'
#' @return a list containing:
#' \itemize{
#'   \item{enrichment_table: a tibble containing each term's enrichment for
#'     every pathway and an enrichment_plot}
#'   \item{enrichment_plots: pre-generated plots containing the most
#'     significant enrichments for each term}
#'   }
#'
#' @examples 
#' 
#' regression_significance <- diffex_mzroll(
#'   nplug_mzroll_normalized,
#'   "normalized_log2_abundance",
#'   "limitation + limitation:DR + 0"
#'   )
#' pathway_nest <- nplug_mzroll_normalized$features %>%
#'  dplyr::select(groupId, pathway) %>%
#'  tidyr::nest(pathway_members = groupId)
#'
#' pathway_list <- purrr::map(
#'  pathway_nest$pathway_members,
#'  function(x) {
#'    as.character(x$groupId)
#'    }
#'  )
#' names(pathway_list) <- pathway_nest$pathway
#'
#' find_pathway_enrichments(
#'   mzroll_list = nplug_mzroll_normalized,
#'   regression_significance,
#'   pathway_list,
#'   test_absolute_effects = FALSE
#'   )
#'
#' @export
find_pathway_enrichments <- function(
  mzroll_list,
  regression_significance,
  pathway_list,
  ranking_measure = "statistic",
  test_absolute_effects = TRUE
  ) {
  
  test_mzroll_list(mzroll_list)
  
  checkmate::assertDataFrame(regression_significance)
  checkmate::assertList(pathway_list, names = "unique")
  checkmate::assertChoice(ranking_measure, colnames(regression_significance))
  checkmate::assertLogical(test_absolute_effects, len = 1)

  regression_significance <- regression_significance %>%
    dplyr::mutate(ranking_measure = !!rlang::sym(ranking_measure))
  
  if (test_absolute_effects) {
    regression_significance <- regression_significance %>%
      dplyr::mutate(ranking_measure = abs(ranking_measure))
  }

  enrichment_by_term <- mzroll_list$features %>%
    dplyr::inner_join(regression_significance, by = "groupId") %>%
    tidyr::nest(term_data = -term) %>%
    dplyr::mutate(gsea_results = purrr::map(
      term_data,
      calculate_pathway_enrichment,
      pathway_list
    )) %>%
    dplyr::select(-term_data) %>%
    tidyr::unnest_wider(gsea_results)

  enrichment_table_grobs <- enrichment_by_term %>%
    dplyr::select(-fgsea_results) %>%
    tidyr::unnest(gsea_table_grob)

  term_plots <- enrichment_table_grobs$gsea_table_grob
  names(term_plots) <- enrichment_table_grobs$term

  enrichment_table <- enrichment_by_term %>%
    dplyr::select(-gsea_table_grob) %>%
    tidyr::unnest(fgsea_results) %>%
    dplyr::select(-leadingEdge)

  pathway_nest <- purrr::map2(
    names(pathway_list),
    unname(pathway_list),
    function(x,y){
      tibble::tibble(pathway = x) %>%
        dplyr::mutate(pathway_members = list(tibble::tibble(groupId = y)))
  }) %>%
    dplyr::bind_rows()
  
  # Issue 569: add back peak group and compound data, removed by fgsea::fgsea()
  enrichment_table_full <- dplyr::inner_join(
    enrichment_table,
    pathway_nest,
    by = "pathway"
    ) %>%
    dplyr::select(
      term,
      pathway,
      pathway_members,
      pval,
      padj,
      ES,
      NES,
      size,
      enrichment_plot
    )

  enrichment_table_full_expanded <- enrichment_table_full %>%
    tidyr::unnest(cols = c(pathway_members))

  pathway_enrichments_output <- list(
    enrichment_table = enrichment_table_full,
    enrichment_table_expanded = enrichment_table_full_expanded,
    enrichment_plots = term_plots
  )

  pathway_enrichments_output
}

#' Calculate Pathway Enrichments
#'
#' @param term_data test statistics for a single term
#' @param pathway_list a named list formatted for fgsea with a name
#'   corresponding to a pathway and entries matching the analytes included
#'   in the set.
#' @param padj_cutoff minimum adjusted p-value to use for creating a the
#'   pathway summary plot.
#'
#' @return a list containing a summary plot for the term and tibble
#'   containing data and plots for all comparisons performed.
calculate_pathway_enrichment <- function(
  term_data,
  pathway_list,
  padj_cutoff = 0.2
  ) {
  
  checkmate::assertDataFrame(term_data)
  checkmate::assertList(pathway_list, names = "unique")
  checkmate::assertNumber(padj_cutoff, lower = 0, upper = 1)
  
  ranked_coefs <- term_data$ranking_measure
  names(ranked_coefs) <- as.character(term_data$groupId)
  sorted_ranked_coefs <- sort(ranked_coefs, decreasing = TRUE)

  # choose an appropriate score type based on whether coefs can be negative
  score_type <- ifelse(all(sorted_ranked_coefs >= 0), "pos", "std")
  
  fgsea_results <- fgsea::fgsea(
    pathway_list,
    sorted_ranked_coefs,
    scoreType = score_type
  )
  
  # generate an enrichment plot for each pathway
  enrichment_plots <- tibble::tibble(
    pathway = names(pathway_list),
    members = unname(pathway_list)
  ) %>%
    dplyr::mutate(enrichment_plot = purrr::map(
      members,
      fgsea::plotEnrichment,
      ranked_coefs
  )) %>%
    dplyr::select(-members)
  
  fgsea_results <- fgsea_results %>%
    dplyr::left_join(enrichment_plots, by = "pathway")

  # generate a plot summarizing the largest pathway changes
  reported_pathways <- fgsea_results %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::arrange(padj)
    
  if (nrow(reported_pathways) == 0) {
    gsea_table_grob <- NULL
  } else {
    gsea_table_grob <- fgsea::plotGseaTable(
      pathway_list[reported_pathways$pathway],
      ranked_coefs,
      fgsea_results,
      gseaParam = 0.5,
      render = FALSE
    )
  }
  
  list(
    fgsea_results = fgsea_results,
    gsea_table_grob = list(plot = gsea_table_grob)
  )
}
