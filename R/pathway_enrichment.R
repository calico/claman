#' Find Pathway Enrichments
#'
#' @inheritParams test_mzroll_list
#' @inheritParams plot_volcano
#' @param pathway_list a named list where names are pathways and entries
#'   are groupIds belonging to each pathway.
#' @param enrichment_method method used to calculate pathway enrichments.
#' \describe{
#'   \item{gsea}{gene set enrichment analysis using \code{fgsea}}
#'   \item{fisher}{Fisher exact test to look for enrichment of differentially
#'     abundant metabolite among pathways.}
#' }
#' @param ranking_measure variable in \code{regression_significance} to use
#'   when calculating enrichment.
#' @param fdr_cutoff minimum adjusted p-value to use for creating a the
#'   pathway summary plot.
#' @param test_absolute_effects If TRUE then only consider the magnitude of
#'   test-statistics when calculating rankings for enrichment; if FALSE then
#'   consider sign as well.
#' @inheritParams diffex_mzroll
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
#' library(dplyr)
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
#' find_pathway_enrichments(
#'   mzroll_list = nplug_mzroll_normalized,
#'   regression_significance,
#'   pathway_list,
#'   test_absolute_effects = FALSE,
#'   enrichment_method = "fisher"
#'   )
#'
#' @export
find_pathway_enrichments <- function(
  mzroll_list,
  regression_significance,
  pathway_list,
  enrichment_method = "gsea",
  ranking_measure = "statistic",
  fdr_cutoff = 0.1,
  test_absolute_effects = TRUE,
  additional_grouping_vars = NULL
  ) {
  
  test_mzroll_list(mzroll_list)
  
  checkmate::assertDataFrame(regression_significance)
  checkmate::assertList(pathway_list, names = "unique")
  checkmate::assertChoice(enrichment_method, c("gsea", "fisher"))
  
  if (enrichment_method == "gsea") {
    enrichment_function <- "calculate_pathway_enrichment_gsea"
  } else if (enrichment_method == "fisher") {
    enrichment_function <- "calculate_pathway_enrichment_fisher"
  } else {
    stop(glue::glue("{enrichment_method} is not implemented"))
  }
  
  checkmate::assertChoice(ranking_measure, colnames(regression_significance))
  checkmate::assertNumber(fdr_cutoff, lower = 0, upper = 1)
  checkmate::assertLogical(test_absolute_effects, len = 1)

  regression_significance <- regression_significance %>%
    dplyr::mutate(ranking_measure = !!rlang::sym(ranking_measure))
  
  features_with_significance <- mzroll_list$features %>%
    dplyr::inner_join(regression_significance, by = "groupId") %>%
    dplyr::filter(!is.na(ranking_measure))
  
  if (test_absolute_effects) {
    regression_significance <- regression_significance %>%
      dplyr::mutate(ranking_measure = abs(ranking_measure))
  }

  # format grouping variables as per regression to check
  #   additional_grouping_vars 
  ignore_me <- format_grouping_vars(
    mzroll_list,
    additional_grouping_vars
  )
  
  if (length(additional_grouping_vars) > 0) {
    unnested_vars <- c("term", additional_grouping_vars)
  } else {
    unnested_vars <- "term"
  }
  nesting_vars <- setdiff(colnames(features_with_significance), unnested_vars)
  
  enrichment_by_term <- suppressWarnings(
    features_with_significance %>%
        tidyr::nest(!!!rlang::syms(nesting_vars), .key = "term_data")
    ) %>%
    dplyr::mutate(
      gsea_results = purrr::map(
        term_data,
        !!rlang::sym(enrichment_function),
        pathway_list
    )) %>%
    dplyr::select(-term_data) %>%
    tidyr::unnest_wider(gsea_results)

  enrichment_table <- enrichment_by_term %>%
    dplyr::select(!!!rlang::syms(unnested_vars), pathway_enrichments) %>%
    tidyr::unnest(pathway_enrichments)
  
  if (enrichment_method == "gsea") {
    
    pathway_enrichments_output <- expand_gsea_output(
      enrichment_by_term,
      enrichment_table,
      pathway_list,
      additional_grouping_vars
    )
    
    return(pathway_enrichments_output)
  } else {
    return(enrichment_table)
  }
}

expand_gsea_output <- function(
  enrichment_by_term,
  enrichment_table,
  pathway_list,
  additional_grouping_vars
  ) {
  
  enrichment_table <- enrichment_table %>%
    dplyr::select(-leadingEdge)
  
  # create a table off grobs
  enrichment_table_grobs <- enrichment_by_term %>%
    dplyr::select(-pathway_enrichments) %>%
    tidyr::unnest(gsea_table_grob)
  
  term_plots <- enrichment_table_grobs$gsea_table_grob
  names(term_plots) <- (enrichment_table_grobs %>%
                          dplyr::select(-gsea_table_grob) %>%
                          tidyr::unite(name))$name
  
  
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
      !!!rlang::syms(additional_grouping_vars),
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
  
  return(pathway_enrichments_output)
}

#' Calculate Pathway Enrichments using GSEA
#'
#' @param term_data test statistics for a single term
#' @inheritParams find_pathway_enrichments
#'
#' @return a list containing a summary plot for the term and tibble
#'   containing data and plots for all comparisons performed.
calculate_pathway_enrichment_gsea <- function(
  term_data,
  pathway_list,
  fdr_cutoff = 0.1
  ) {
  
  checkmate::assertDataFrame(term_data)
  checkmate::assertList(pathway_list, names = "unique")
  checkmate::assertNumber(fdr_cutoff, lower = 0, upper = 1)
  
  ranked_coefs <- term_data$ranking_measure
  names(ranked_coefs) <- as.character(term_data$groupId)
  sorted_ranked_coefs <- sort(ranked_coefs, decreasing = TRUE)
  
  n_duplicated_names <- length(ranked_coefs) -
    length(unique(names(ranked_coefs)))
  
  if (n_duplicated_names > 0) {
    stop(glue::glue(
      "There were {n_duplicated_names} duplicated statistics for the same
      - group found. Set the \"additional_grouping_vars\" argument in
      - \"find_pathway_enrichments\" so it matched \"diffex_mzroll\""
      ))
  }
  
  checkmate::assertNumeric(ranked_coefs, finite = TRUE, any.missing = FALSE)
  
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
    dplyr::filter(padj < fdr_cutoff) %>%
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
  
  output <- list(
    pathway_enrichments = fgsea_results,
    gsea_table_grob = list(plot = gsea_table_grob)
  )
  
  return(output)
}

#' Calculate Pathway Enrichments with a Fisher-Exact tests
#'
#' @inheritParams calculate_pathway_enrichment_gsea
#'
#' @return a list containing a tibble of enrichments
calculate_pathway_enrichment_fisher <- function(
  term_data,
  pathway_list,
  fdr_cutoff = 0.1
) {
  
  checkmate::assertDataFrame(term_data)
  checkmate::assertList(pathway_list, names = "unique")
  checkmate::assertNumber(fdr_cutoff, lower = 0, upper = 1)
  
  n_features <- nrow(term_data)
  diffex_features <- term_data %>%
    dplyr::filter(qvalue < fdr_cutoff)
  
  pathways_tibble <- purrr::map2(
    names(pathway_list),
    unname(pathway_list),
    function(x,y){
      tibble::tibble(pathway = x, groupId = y)
      }) %>%
    dplyr::bind_rows() %>%
    tidyr::nest(members = groupId) %>%
    dplyr::mutate(enrichment = purrr::map(
      members,
      fisher_test_enrichment,
      diffex_features,
      n_features
      )) %>%
    dplyr::select(pathway, enrichment) %>%
    tidyr::unnest(enrichment) %>%
    dplyr::mutate(qvalue = p.adjust(p.value, method = "BH"))
  
  output <- list(
    pathway_enrichments = pathways_tibble
  )
  
  return(output)
}


fisher_test_enrichment <- function(members, diffex_features, n_features) {
  
  n_diffex <- nrow(diffex_features)
  n_members <- nrow(members)
  
  n_diffex_members <- sum(members$groupId %in% diffex_features$groupId)
  n_diffex_nonmembers <- n_diffex - n_diffex_members
  n_member_nondiffex <- n_members - n_diffex_members
  n_nonmember_nondiffex <- n_features - n_diffex_members - n_diffex_nonmembers - n_member_nondiffex
  
  contingency_table <- matrix(c(
    n_diffex_members,
    n_diffex_nonmembers,
    n_member_nondiffex,
    n_nonmember_nondiffex
    ),
    byrow = TRUE, ncol = 2)
  
  stopifnot(all(contingency_table) >= 0, sum(contingency_table) == n_features)
  
  fisher.test(contingency_table, alternative = "greater") %>%
    broom::tidy()
}
