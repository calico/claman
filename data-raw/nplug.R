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
    limitation = Nutr,
    exp_ref = `Exp/Ref`,
    extraction = Method
  ) %>%
  dplyr::mutate(
    month = month.abb[month],
    replicate = toupper(replicate)
  ) %>%
  dplyr::mutate(limitation = ifelse(
    limitation %in% c("R", "R2", "R3"),
    "PO4",
    limitation))

# use DR when finding distinct experimental conditions
distinct_exp_conditions <- nplug_samples %>%
  filter(exp_ref == "exp") %>%
  distinct(month, limitation, DR, extraction) %>%
  mutate(condition = 1:n())

# ignore DR when finding distinct reference conditions
distinct_ref_conditions <- nplug_samples %>%
  filter(exp_ref == "ref") %>%
  distinct(month, extraction) %>%
  mutate(condition = max(distinct_exp_conditions$condition) + 1:n())

nplug_samples <- bind_rows(
  nplug_samples %>%
    filter(exp_ref == "exp") %>%
    left_join(
      distinct_exp_conditions,
      by = c("month", "limitation", "DR", "extraction")
    ),
  nplug_samples %>%
    filter(exp_ref == "ref") %>%
    left_join(
      distinct_ref_conditions,
      by = c("month", "extraction")
    )
) %>%
  left_join(
    distinct_ref_conditions %>%
      rename(reference = condition),
    by = c("month", "extraction")
  )

usethis::use_data(nplug_samples, overwrite = TRUE)

# create mzroll file

nplug_raw <- nplug_raw %>%
  dplyr::mutate(sampleId = 1:n())

# setup measurements
abundances <- nplug_raw %>%
  dplyr::select(-c(sample_name:Method)) %>%
  tidyr::gather(compoundName, peakAreaTop, -sampleId)
  
peakgroups <- abundances %>%
  dplyr::distinct(compoundName) %>%
  dplyr::mutate(
    groupId = 1:n(),
    # indicate that the peakgroup has been validated
    label = "g"
  )

peaks <- abundances %>%
  dplyr::inner_join(
    peakgroups %>%
      dplyr::select(groupId, compoundName),
    by = "compoundName"
    ) %>%
  dplyr::select(-compoundName) %>%
  dplyr::mutate(peakId = 1:dplyr::n())

mzroll_data <- list(
  samples = nplug_raw %>% dplyr::select(sampleId, name = sample_name),
  peakgroups = peakgroups,
  peaks = peaks
) %>%
  clamr::expand_with_mzroll_defaults()

mzroll_db_path <- file.path("inst", "extdata", "nplug.mzrollDB")
if (file.exists(mzroll_db_path)) {
  file.remove(mzroll_db_path)
}

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
