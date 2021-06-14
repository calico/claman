library(dplyr)

test_that("All pathway enrichmnet methods work", {
  
  regression_significance <- diffex_mzroll(
    nplug_mzroll_normalized,
    "normalized_log2_abundance",
    "limitation + limitation:DR + 0"
  )
  
  expect_equal(nrow(regression_significance), 1060)
  
  pathway_nest <- nplug_mzroll_normalized$features %>%
    dplyr::select(groupId, pathway) %>%
    tidyr::nest(pathway_members = groupId)
  
  pathway_list <- purrr::map(
    pathway_nest$pathway_members,
    function(x) {
      as.character(x$groupId)
    }
  )
  names(pathway_list) <- pathway_nest$pathway
  
  gsea_enrichments <- find_pathway_enrichments(
    mzroll_list = nplug_mzroll_normalized,
    regression_significance,
    pathway_list,
    test_absolute_effects = FALSE
  )
  
  top_enrichments <- gsea_enrichments$enrichment_table %>%
    dplyr::arrange(padj) %>%
    dplyr::slice(1:2)
  
  expect_equal(
    top_enrichments[["term"]],
    c("limitationNH4:DR", "limitationNH4")
  )
  expect_equal(
    top_enrichments[["pathway"]],
    c("Amino Acids", "Amino Acids")
  )
  expect_equal(
    top_enrichments[["pval"]],
    c(0.000000218, 0.00000892),
    tolerance = 1e-5
  )
  
  fisher_enrichments <- find_pathway_enrichments(
    mzroll_list = nplug_mzroll_normalized,
    regression_significance,
    pathway_list,
    test_absolute_effects = FALSE,
    enrichment_method = "fisher"
  )
  
  top_enrichments <- fisher_enrichments%>%
    dplyr::arrange(qvalue) %>%
    dplyr::slice(1:2)
  
  expect_equal(
    top_enrichments[["term"]],
    c("limitationNH4:DR", "limitationNH4")
  )
  expect_equal(
    top_enrichments[["pathway"]],
    c("Amino Acids", "Amino Acids")
  )
  expect_equal(
    top_enrichments[["p.value"]],
    c(0.0000127, 0.0000516),
    tolerance = 1e-5
    )
  
  
  
})

