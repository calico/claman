% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_mzroll.R
\name{process_mzroll_multi}
\alias{process_mzroll_multi}
\title{Process mzroll multi}
\usage{
process_mzroll_multi(
  mzroll_paths,
  samples_tbl,
  id_strings,
  only_identified = TRUE,
  validate = FALSE,
  exact = FALSE
)
}
\arguments{
\item{mzroll_paths}{a tibble with two variables:
\itemize{
  \item{method_tag: a character vector to tag each dataset with},
  \item{mzroll_db_path: path to a mzrollDB file}
  }}

\item{samples_tbl}{Table of sample metadata}

\item{id_strings}{one or more variables which will be used to match
sample names.}

\item{only_identified}{TRUE/FALSE, filter to only features which were
identified.}

\item{validate}{TRUE/FALSE, use meta-data to only name the subset of
features which were manually validated.}

\item{exact}{if true, an exact match between mzroll names and id_strings
will be found; if false, then substring matches will be used.}
}
\value{
an mzroll_list containing three tibbles:
\itemize{
  \item{peakgroups: one row per unique analyte (defined by a unique
    groupId)},
  \item{samples: one row per unique sample (defined by a unique sampleId)},
  \item{peaks: one row per peak (samples x peakgroups)}
  }
}
\description{
Aggregate multiple mzrollDB datasets which possess the same sample metadata
}
\examples{
mzroll_paths <- tibble::tribble(
  ~method_tag, ~mzroll_db_path,
  "method1", nplug_mzroll(),
  "method2", nplug_mzroll()
)

process_mzroll_multi(
  mzroll_paths,
  nplug_samples,
  "sample_name",
  exact = TRUE
  )
}
