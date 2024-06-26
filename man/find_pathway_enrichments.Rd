% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathway_enrichment.R
\name{find_pathway_enrichments}
\alias{find_pathway_enrichments}
\title{Find Pathway Enrichments}
\usage{
find_pathway_enrichments(
  mzroll_list,
  regression_significance,
  pathway_list,
  enrichment_method = "gsea",
  ranking_measure = "statistic",
  fdr_cutoff = 0.1,
  test_absolute_effects = TRUE,
  additional_grouping_vars = NULL
)
}
\arguments{
\item{mzroll_list}{output of \link{process_mzroll} or
\link{process_mzroll_multi}}

\item{regression_significance}{returned by \code{\link{diffex_mzroll}};
a tibble of tests performed.}

\item{pathway_list}{a named list where names are pathways and entries
are groupIds belonging to each pathway.}

\item{enrichment_method}{method used to calculate pathway enrichments.
\describe{
  \item{gsea}{gene set enrichment analysis using \code{fgsea}}
  \item{fisher}{Fisher exact test to look for enrichment of differentially
    abundant metabolite among pathways.}
}}

\item{ranking_measure}{variable in \code{regression_significance} to use
when calculating enrichment.}

\item{fdr_cutoff}{minimum adjusted p-value to use for creating a the
pathway summary plot.}

\item{test_absolute_effects}{If TRUE then only consider the magnitude of
test-statistics when calculating rankings for enrichment; if FALSE then
consider sign as well.}

\item{additional_grouping_vars}{sample-, or measurement-level variables
to groupby when performing regression in addition to groupId.}
}
\value{
a list containing:
\itemize{
  \item{enrichment_table: a tibble containing each term's enrichment for
    every pathway and an enrichment_plot}
  \item{enrichment_plots: pre-generated plots containing the most
    significant enrichments for each term}
  }
}
\description{
Find Pathway Enrichments
}
\examples{

library(dplyr)

regression_significance <- diffex_mzroll(
  nplug_mzroll_normalized,
  "normalized_log2_abundance",
  "limitation + limitation:DR + 0"
)
pathway_nest <- nplug_mzroll_normalized$features \%>\%
  dplyr::select(groupId, pathway) \%>\%
  tidyr::nest(pathway_members = groupId)

pathway_list <- purrr::map(
  pathway_nest$pathway_members,
  function(x) {
    as.character(x$groupId)
  }
)
names(pathway_list) <- pathway_nest$pathway

find_pathway_enrichments(
  mzroll_list = nplug_mzroll_normalized,
  regression_significance,
  pathway_list,
  test_absolute_effects = FALSE,
  enrichment_method = "fisher"
)
}
