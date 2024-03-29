% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_mzroll.R
\name{process_mzroll_load_peakgroups}
\alias{process_mzroll_load_peakgroups}
\title{Process MzRoll - Load Peakgroups}
\usage{
process_mzroll_load_peakgroups(mzroll_db_con, quant_col = "peakAreaTop")
}
\arguments{
\item{mzroll_db_con}{an SQLite connection to an mzrollDB database}

\item{quant_col}{name of quant column to use for standardization from the
  peaks table.
\code{default = "peakAreaTop"}}
}
\value{
a table of peakgroups
}
\description{
Load peakgroups table and add mz and rt variables to peakgroups based on
  peak-level data.
}
