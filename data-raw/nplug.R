## code to prepare `NPLUG` dataset goes here

library(dplyr)

# save the raw data file to the top-level data directory to future-proof
#   against journal URL changes.

raw_data_url <- "https://www.molbiolcell.org/doi/suppl/10.1091/mbc.e09-07-0597/suppl_file/data_set_2.txt"
local_save <- file.path("assets", "data", "boer2010.Rds")

if (file.exists(local_save)) {
  nplug_raw <- readRDS(local_save)
} else {
  nplug_raw <- readr::read_tsv(raw_data_url) %>%
    dplyr::rename(sample_name = X1)
  
  saveRDS(nplug_raw, local_save)
}

# sample meta-data

nplug_samples <- nplug_raw %>%
  dplyr::select(
    sample_name,
    month = Month,
    replicate = Which,
    DR = Gr,
    exp_ref = `Exp/Ref`,
    extraction = Method
  ) %>%
  dplyr::mutate(
    month = month.abb[month],
    replicate = toupper(replicate)
  )

usethis::use_data(nplug_samples, overwrite = TRUE)

# create mzroll file

nplug_raw <- nplug_raw %>%
  dplyr::mutate(sampleId = 1:n())

# setup measurements
abundances <- nplug_raw %>%
  dplyr::select(-c(sample_name:Method)) %>%
  tidyr::gather(name, peakAreaTop, -sampleId)
  
peakgroups <- abundances %>%
  dplyr::distinct(name) %>%
  dplyr::mutate(groupId = 1:n())

peaks <- abundances %>%
  dplyr::inner_join(peakgroups, by = "name") %>%
  dplyr::select(-name)

mzroll_data <- list(
  samples = nplug_raw %>% dplyr::select(sampleId),
  peakgroups = peakgroups,
  peaks = peaks
) %>%
  clamr::expand_with_mzroll_defaults()

mzroll_db_path <- file.path("inst", "extdata", "nplug.mzrollDB")

mzroll_db_con <- DBI::dbConnect(RSQLite::SQLite(), mzroll_db_path)

purrr::walk2(
  names(mzroll_data),
  mzroll_data,
  DBI::dbCreateTable,
  conn = mzroll_db_con
)

purrr::walk2(
  names(mzroll_data),
  mzroll_data,
  DBI::dbAppendTable,
  conn = mzroll_db_con
)
  
DBI::dbDisconnect(mzroll_db_con)
