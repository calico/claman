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

# save compound metadata

pathways <- list(
  "Pentose Phosphate Pathway" = c("6-phospho-D-gluconate", "Ribose-P",
                                  "D-erythrose-4-phosphate", "D-gluconate",
                                  "D-glucono-delta-lactone-6-phosphate",
                                  "D-sedoheptulose-7-phosphate",
                                  "pyridoxine", "riboflavin"),
  "Glycolysis" = c("1,3-diphopshateglycerate", "3-phosphoglycerate",
                   "D-glyceraldehdye-3-phosphate", "dihydroxy-acetone-phosphate",
                   "fructose-1,6-bisphosphate", "glycerate", "hexose-phosphate",
                   "pyruvate", "trehalose/sucrose", "trehalose-6-phosphate",
                   "Lactate"),
  "Purine Synthesis" = c("adenosine", "cyclic-AMP", "dATP", "deoxyadenosine",
                         "deoxyguanosine", "guanine", "guanosine",
                         "hypoxanthine", "inosine", "thiamine", "xanthine",
                         "dGDP"),
  "Pyrimidines" = c("CTP", "CDP", "UTP", "UDP", "OMP", "cytidine", "cytosine",
                    "dCTP", "dihydroorotate", "dTTP", "N-carbamoyl-L-aspartate",
                    "orotate", "UDP-D-glucose", "UDP-N-acetyl-glucosamine",
                    "uridine"),
  "Amino Acids" = c("alanine", "arginine", "asparagine", "aspartate", "glutamate",
                    "glutamine", "glycine", "histidine", "leucine/isoleucine",
                    "lysine", "methionine", "phenylalanine", "proline",
                    "serine", "threonine", "tryptophan", "tyrosine", "valine"),
  "Amino Acid Synthesis" = c("citrulline", "cystathionine", "histidinol",
                             "ornithine", "phenylpyruvate", "Arginino-succinate",
                             "Carbamyl phosphate", "Homoserine", "Hydroxyphenylpyruvate",
                             "N-acetyl-glutamate", "Prephenic acid", "PRPP"),
  "Amino Acid Derived" = c("glutathione", "glutathione disulfide",
                           "N-acetyl-glucosamine-1-phosphate", "N-acetyl-glutamine",
                           "nicotinate", "pantothenate", "quinolinate",
                           "S-adenosyl-L-methionine", "Dimethylglycine"),
  "TCA cycle" = c("acetyl-CoA", "a-ketoglutarate",  "aconitate", "citrate",
                  "citrate/isocitrate", "fumarate", "isocitrate",
                  "malate", "succinate"),
  "Energetics" = c("ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "FAD", "FMN",
                   "NAD+", "NADH", "NADP+"),
  "Lipid metabolism" = c("choline", "3-hydroxy-3-methylglutaryl-CoA",
                         "trans, trans-farnesyl diphosphate"))

nplug_compounds <- purrr::map(names(pathways),
           function(x){tibble::tibble(pathway = x, compoundName = pathways[[x]])}) %>%
  dplyr::bind_rows() %>%
  dplyr::select(compoundName, pathway)

usethis::use_data(nplug_compounds, overwrite = TRUE)
