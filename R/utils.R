#' \code{calicomics} package
#'
#' Calicomics can read .mzrollDB files created using MAVEN or Quahog into an
#' mzroll_list. These lists can then be modified through normalization,
#' signal flooring, and filtering. Differential abundance analysis of
#' metabolites can then be performed along with a variety of visualizations.
#'
#' @docType package
#' @name calicomics
#' @importFrom dplyr %>%
#' @import ggplot2
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
utils::globalVariables(c(
  ".",
  ":=",
  ".center",
  ".group_center",
  ".loess_fit",
  ".loess_shift",
  ".reference",
  ".reference_diff",
  "abund_mean",
  "abund_se",
  "abundance",
  "anova_signif",
  "centered_log2_abundance",
  "char",
  "compoundId",
  "compoundName",
  "databaseId",
  "datetime",
  "diff_to_median",
  "displayName",
  "double_1",
  "double_2",
  "double_3",
  "double_4",
  "enrichment_plot",
  "ES",
  "FA1",
  "FA2",
  "FA3",
  "FA4",
  "fdr_summary",
  "field",
  "filename",
  "focus_pathway",
  "fgsea_results",
  "groupData",
  "groupId",
  "gsea_results",
  "gsea_table_grob",
  "id",
  "ID string entry",
  "lm_fits",
  "is_discovery",
  "is_focus_pathway",
  "is_tracking_sheet",
  "is_unknown",
  "label",
  "leadingEdge",
  "lipidClass",
  "log2_abundance",
  "median",
  "median_abund",
  "method_tag",
  "MS ID string",
  "MS ID string alternative",
  "mz_delta",
  "mz_median",
  "mzmax",
  "mzmax_adj",
  "mzmin",
  "mzmin_adj",
  "mzroll_db_path",
  "n",
  "n_entries",
  "n_unique_records",
  "name",
  "name_tibble",
  "NES",
  "new_compoundName",
  "newSampleId",
  "nMoreExtreme",
  "num_sn_chains",
  "OH_1",
  "OH_2",
  "OH_3",
  "OH_4",
  "old_sampleId",
  "ordered_groupId",
  "ordered_sampleId",
  "one_peak_data",
  "p.value",
  "p.value.trans",
  "padj",
  "pathway",
  "pathway_entries",
  "pathwayId",
  "pathwayName",
  "peak_label",
  "peakAreaTop",
  "peakMz",
  "peakId",
  "plasmalogen_type",
  "position",
  "processed_mzroll",
  "pval",
  "qvalue",
  "r_position",
  "rt",
  "rt_delta",
  "rt_median",
  "rtmax",
  "rtmax_adj",
  "rtmin",
  "rtmin_adj",
  "sampleId",
  "sample description",
  "samples",
  "scaling_factor",
  "single_1",
  "single_2",
  "single_3",
  "single_4",
  "size",
  "sn1",
  "sn2",
  "sn3",
  "sn4",
  "sn_chains",
  "statistic",
  "sumComposition",
  "systematicCompoundId",
  "term",
  "term_data",
  "total_double",
  "total_OH",
  "total_single",
  "tube",
  "tube label",
  "variable",
  "weights"
))

#' Configure Database Access
#'
#' Setup and test connection to Calico databases
#'
#' @param standards_db type of standards database to connect to:
#' \itemize{
#'   \item{standard: just standards}
#'   \item{extended: standards, external libraries and predicted spectra}
#'   }
#'
#' @return a list containing connections to the standards and systematic
#'   compounds database
#'
#' @export
configure_db_access <- function(standards_db = "standard") {

  # load metabolite MySQL db password
  db_password <- Sys.getenv("metabolite_db")
  if (db_password == "") {
    warning("metabolite_db password missing from .Renviron
              please add it. i.e., metabolite_db=xxxxx")
  }

  standard_dbs <- tibble::tribble(
    ~name, ~db_name,
    "standard", "mass_spec_standards",
    "extended", "mass_spec_standards_extended"
  )

  if (any(class(standards_db) != "character") || length(standards_db) != 1) {
    stop("\"standards_db\" must be a length 1 character vector
           valid entries are: ", paste(standard_dbs$name, collapse = ", "))
  }

  if (!(standards_db %in% standard_dbs$name)) {
    stop(
      standards_db,
      " is not a valid entry for \"standards_db\" valid entries are: ",
      paste(standard_dbs$name, collapse = ", ")
    )
  }
  standard_db_name <- standard_dbs$db_name[standard_dbs$name == standards_db]

  # connect to metabolite MySQL databases
  # (this is behind the firewall)
  mass_spec_standards_con <- DBI::dbConnect(
    RMySQL::MySQL(),
    user = "admin",
    password = db_password,
    dbname = standard_db_name,
    host = "104.196.252.153"
  )

  systematic_compounds_con <- DBI::dbConnect(
    RMySQL::MySQL(),
    user = "admin",
    password = db_password,
    dbname = "systematic_compounds",
    host = "104.196.252.153"
  )

  # present connections
  list(
    mass_spec_standards_con = mass_spec_standards_con,
    systematic_compounds_con = systematic_compounds_con
  )
}

#' Test Mzroll List
#'
#' @param mzroll_list output of \link{process_mzroll} or
#'   \link{process_mzroll_multi}
#'
#' \itemize{
#'   \item{peakgroups: one row per unique analyte (defined by a
#'     unique groupId)},
#'   \item{samples: one row per unique sample (defined by a unique sampleId)},
#'   \item{peaks: one row per peak (samples x peakgroups)}
#'   }
#'
test_mzroll_list <- function(mzroll_list) {
  if (!("list" %in% class(mzroll_list))) {
    stop("\"mzroll_list\" must be a list")
  }

  required_tables <- c("peakgroups", "samples", "peaks")
  provided_tables <- names(mzroll_list)

  missing_required_tables <- setdiff(required_tables, provided_tables)

  if (length(missing_required_tables) != 0) {
    stop(
      "missing required tables: ",
      paste(missing_required_tables, collapse = ", "),
      "; generate mzroll_list with \"process_mzroll\""
    )
  }
  extra_provided_tables <- setdiff(provided_tables, required_tables)
  if (length(missing_required_tables) != 0) {
    warning(
      "extra tables present in mzroll_list: ",
      paste(extra_provided_tables, collapse = ", ")
    )
  }

  # this overlaps with some functions in clamr/clamdb - might want to pull
  # them out as general utils
  required_fields <- tibble::tribble(
    ~tbl, ~variable,
    "peaks", "groupId",
    "peaks", "sampleId",
    "peakgroups", "groupId",
    "samples", "sampleId"
  )

  included_variables <- tibble::tibble(tbl = names(mzroll_list)) %>%
    dplyr::mutate(variable = purrr::map(mzroll_list, colnames)) %>%
    tidyr::unnest(variable)

  absent_required_fields <- required_fields %>%
    dplyr::anti_join(included_variables, by = c("tbl", "variable")) %>%
    dplyr::mutate(message = glue::glue(
      "required variable missing from {tbl}: {variable}"
    ))

  if (nrow(absent_required_fields) != 0) {
    stop(paste(absent_required_fields$message, collapse = "\n"))
  }
}

#' Paths for X0106 examples
#'
#' methods and paths to X0106 example mzroll datasets.
#'
#' @examples
#' \dontrun{
#' X0106_paths <- tibble::tribble(
#'   ~method_tag, ~mzroll_db_path,
#'   "metabolites (neg)", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M001-peakdetector.mzrollDB",
#'   "metabolites (pos) - 1", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M002-peakdetector.mzrollDB",
#'   "metabolites (pos) - 2", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M002-peakdetector2.mzrollDB",
#'   "lipids (neg) - 1", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M004A-1.mzrollDB",
#'   "lipids (neg) - 2", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M004A-2.mzrollDB",
#'   "lipids (neg) - 3", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M004A-2b-AA-BDP-CL.mzrollDB",
#'   "lipids (neg) - 4", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M004A-2b-AA-BDP-CL+Stds.mzrollDB",
#'   "lipids (pos) - 1", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M005A-peakdetector.mzrollDB",
#'   "lipids (pos) - 2", "/tmp/mzkit_assets/X0106_mzrollDBs/X0106-M005A-peakdetectorMGDGTGStds.mzrollDB"
#' )
#'
#' usethis::use_data(X0106_paths, overwrite = TRUE)
#' }
"X0106_paths"

#' Flatten mzroll list to single table
#'
#' Flatten mzroll list to single table using groupId and sampleId
#'
#' @param mzroll_list output of process_mzrollDB functions
#' \itemize{
#'   \item{peakgroups: one row per unique analyte (defined by a unique
#'     groupId)},
#'   \item{samples: one row per unique sample (defined by a unique sampleId)},
#'   \item{peaks: one row per peak (samples x peakgroups)}
#'   }
#'
#' @return mzroll_list (process_mzrollDB functions output) reformatted as a
#'   single table
#'
#' @export
mzroll_table <- function(mzroll_list) {

  # guards
  if (class(mzroll_list) != "list") {
    stop("input \"mzroll_list\" must be a list.")
  }
  if (length(mzroll_list) != 3) {
    stop("input \"mzroll_list\" must contain exactly 3 entries.")
  }
  if (names(mzroll_list[1]) != "peakgroups") {
    stop("input \"mzroll_list\" first element must be named \"peakgroups\"")
  }
  if (names(mzroll_list[2]) != "samples") {
    stop("input \"mzroll_list\" second element must be named \"samples\"")
  }
  if (names(mzroll_list[3]) != "peaks") {
    stop("input \"mzroll_list\" third element must be named \"peaks\"")
  }
  if (!"groupId" %in% colnames(mzroll_list[[1]])) {
    stop("peakgroups table missing required column \"groupId\"")
  }
  if (!"groupId" %in% colnames(mzroll_list[[3]])) {
    stop("peaks table missing required column \"groupId\"")
  }
  if (!"sampleId" %in% colnames(mzroll_list[[2]])) {
    stop("samples table missing required column \"sampleId\"")
  }
  if (!"sampleId" %in% colnames(mzroll_list[[3]])) {
    stop("peaks table missing required column \"sampleId\"")
  }

  mzroll_table <- dplyr::inner_join(
    mzroll_list$peakgroups,
    mzroll_list$peaks,
    by = c("groupId")
  ) %>%
    dplyr::inner_join(mzroll_list$samples, by = c("sampleId"))

  mzroll_table
}

#' Generate relevant lipid components from compound name
#'
#' Produce data frame with lipid components from name
#'
#' @param compound_name vector of compound names
#'
#' @return compound_name vector expanded to large table, with structural
#'   components described as their own columns.
#'
#' @export
lipid_components <- function(compound_name) {
  compound_components <- tibble::tibble("compoundName" = compound_name) %>%
    # general compound information
    dplyr::mutate(
      lipidClass = stringr::str_extract(compoundName, "^.*(?=\\()")
    ) %>%
    dplyr::mutate(
      compound_name_adduct = stringr::str_extract(compoundName, "(?<=\\) ).*$")
    ) %>%
    dplyr::mutate(
      plasmalogen_type = stringr::str_extract(compoundName, "[op]-(?=[0-9]+:)")
    ) %>%
    dplyr::mutate(lipidClass_o_p = ifelse(
      is.na(plasmalogen_type),
      lipidClass,
      paste0(plasmalogen_type, lipidClass)
    )) %>%
    dplyr::mutate(lipidClass_plas = ifelse(
      is.na(plasmalogen_type),
      lipidClass,
      paste0("plas-", lipidClass)
    )) %>%
    # Parse all chains.  Also retrieves number of hydroxyl groups on each
    #   FA (m=1, d=2, t=3), if available
    dplyr::mutate(sn_chains = stringr::str_extract_all(
      compoundName,
      "(?<=[\\(/_])[A-Za-z]?-?[0-9]+:[0-9]+;?O?[0-9]?(?=[\\)/_])"
    )) %>%
    dplyr::mutate(num_sn_chains = sapply(sn_chains, function(x) {
      length(x)
    })) %>%
    dplyr::mutate(num_unique_sn_chains = sapply(sn_chains, function(x) {
      x %>%
        unique() %>%
        length()
    })) %>%
    # sn1
    dplyr::mutate(sn1 = sapply(sn_chains, function(x) {
      x[1]
    })) %>%
    dplyr::mutate(FA1 = stringr::str_extract(sn1, "[0-9]+:[0-9]+")) %>%
    dplyr::mutate(single_1 = as.numeric(
      stringr::str_extract(FA1, ".*(?=:)")
    )) %>%
    dplyr::mutate(single_1 = ifelse(is.na(single_1), 0, single_1)) %>%
    dplyr::mutate(double_1 = as.numeric(
      stringr::str_extract(FA1, "(?<=:).*")
    )) %>%
    dplyr::mutate(double_1 = ifelse(is.na(double_1), 0, double_1)) %>%
    dplyr::mutate(double_1 = ifelse(
      (!is.na(plasmalogen_type) & plasmalogen_type == "p-"),
      double_1 + 1,
      double_1
    )) %>%
    dplyr::mutate(OH_1 = dplyr::case_when(
      grepl("^m", sn1) ~ 1,
      grepl("^d", sn1) ~ 2,
      grepl("^t", sn1) ~ 3,
      grepl(";O4", sn1) ~ 4,
      grepl(";O3", sn1) ~ 3,
      grepl(";O2", sn1) ~ 2,
      grepl(";O1", sn1) ~ 1,
      grepl(";O", sn1) ~ 1,
      TRUE ~ 0
    )) %>%
    # sn2
    dplyr::mutate(sn2 = sapply(sn_chains, function(x) {
      x[2]
    })) %>%
    dplyr::mutate(FA2 = stringr::str_extract(sn2, "[0-9]+:[0-9]+")) %>%
    dplyr::mutate(single_2 = as.numeric(
      stringr::str_extract(FA2, ".*(?=:)")
    )) %>%
    dplyr::mutate(single_2 = ifelse(is.na(single_2), 0, single_2)) %>%
    dplyr::mutate(double_2 = as.numeric(
      stringr::str_extract(FA2, "(?<=:).*")
    )) %>%
    dplyr::mutate(double_2 = ifelse(is.na(double_2), 0, double_2)) %>%
    dplyr::mutate(OH_2 = dplyr::case_when(
      grepl("^m", sn2) ~ 1,
      grepl("^d", sn2) ~ 2,
      grepl("^t", sn2) ~ 3,
      grepl(";O4", sn2) ~ 4,
      grepl(";O3", sn2) ~ 3,
      grepl(";O2", sn2) ~ 2,
      grepl(";O1", sn2) ~ 1,
      grepl(";O", sn2) ~ 1,
      TRUE ~ 0
    )) %>%
    # sn3 (for TGs)
    dplyr::mutate(sn3 = sapply(sn_chains, function(x) {
      x[3]
    })) %>%
    dplyr::mutate(FA3 = stringr::str_extract(sn3, "[0-9]+:[0-9]+")) %>%
    dplyr::mutate(single_3 = as.numeric(
      stringr::str_extract(FA3, ".*(?=:)")
    )) %>%
    dplyr::mutate(single_3 = ifelse(is.na(single_3), 0, single_3)) %>%
    dplyr::mutate(double_3 = as.numeric(
      stringr::str_extract(FA3, "(?<=:).*")
    )) %>%
    dplyr::mutate(double_3 = ifelse(is.na(double_3), 0, double_3)) %>%
    dplyr::mutate(OH_3 = dplyr::case_when(
      grepl("^m", sn3) ~ 1,
      grepl("^d", sn3) ~ 2,
      grepl("^t", sn3) ~ 3,
      grepl(";O4", sn3) ~ 4,
      grepl(";O3", sn3) ~ 3,
      grepl(";O2", sn3) ~ 2,
      grepl(";O1", sn3) ~ 1,
      grepl(";O", sn3) ~ 1,
      TRUE ~ 0
    )) %>%
    # sn4 (for Cardiolipins CLs)
    dplyr::mutate(sn4 = sapply(sn_chains, function(x) {
      x[4]
    })) %>%
    dplyr::mutate(FA4 = stringr::str_extract(sn4, "[0-9]+:[0-9]+")) %>%
    dplyr::mutate(single_4 = as.numeric(
      stringr::str_extract(FA4, ".*(?=:)")
    )) %>%
    dplyr::mutate(single_4 = ifelse(is.na(single_4), 0, single_4)) %>%
    dplyr::mutate(double_4 = as.numeric(
      stringr::str_extract(FA4, "(?<=:).*")
    )) %>%
    dplyr::mutate(double_4 = ifelse(is.na(double_4), 0, double_4)) %>%
    dplyr::mutate(OH_4 = dplyr::case_when(
      grepl("^m", sn4) ~ 1,
      grepl("^d", sn4) ~ 2,
      grepl("^t", sn4) ~ 3,
      grepl(";O4", sn4) ~ 4,
      grepl(";O3", sn4) ~ 3,
      grepl(";O2", sn4) ~ 2,
      grepl(";O1", sn4) ~ 1,
      grepl(";O", sn4) ~ 1,
      TRUE ~ 0
    )) %>%
    # finally, put together all info from chains to describe summed
    #   composition. sometimes, the compound is already reported as a summed
    #   composition. In that case, return the original name.
    dplyr::mutate(total_single = ifelse(
      num_sn_chains == 0,
      as.numeric(stringr::str_extract(
        compoundName,
        "(?<=\\()[A-Za-z]?-?[0-9]+(?=:[0-9]+,)"
      )),
      (single_1 + single_2 + single_3 + single_4)
    )) %>%
    dplyr::mutate(total_single = ifelse(
      is.na(total_single),
      0,
      total_single
    )) %>%
    dplyr::mutate(total_double = ifelse(
      num_sn_chains == 0,
      as.numeric(stringr::str_extract(compoundName, "(?<=:)[0-9]+(?=,)")),
      (double_1 + double_2 + double_3 + double_4)
    )) %>%
    dplyr::mutate(total_double = ifelse(
      is.na(total_double),
      0,
      total_double
    )) %>%
    dplyr::mutate(total_OH = ifelse(
      num_sn_chains == 0,
      as.numeric(stringr::str_extract(compoundName, "(?<=,)[0-9]+(?=-)")),
      OH_1 + OH_2 + OH_3 + OH_4
    )) %>%
    dplyr::mutate(total_OH = ifelse(is.na(total_OH), 0, total_OH)) %>%
    dplyr::mutate(sumComposition = dplyr::case_when(
      num_sn_chains > 0 & total_OH == 0 ~
      glue::glue("{class}({plas}{single}:{double})",
        class = lipidClass,
        plas = ifelse(is.na(plasmalogen_type), "", plasmalogen_type),
        single = total_single,
        double = ifelse(
          (!is.na(plasmalogen_type) & plasmalogen_type == "p-"),
          (total_double - 1),
          total_double
        )
      ),
      num_sn_chains > 0 & total_OH == 1 ~
      glue::glue(
        "{class}({plas}{single}:{double};O)",
        class = lipidClass,
        plas = ifelse(is.na(plasmalogen_type), "", plasmalogen_type),
        single = total_single,
        double = ifelse(
          (!is.na(plasmalogen_type) & plasmalogen_type == "p-"),
          (total_double - 1),
          total_double
        )
      ),
      num_sn_chains > 0 & total_OH > 1 ~
      glue::glue(
        "{class}({plas}{single}:{double};O{num_OH})",
        class = lipidClass,
        plas = ifelse(is.na(plasmalogen_type), "", plasmalogen_type),
        single = total_single,
        double = ifelse(
          (!is.na(plasmalogen_type) & plasmalogen_type == "p-"),
          (total_double - 1),
          total_double
        ),
        num_OH = total_OH
      ),
      num_sn_chains == 0 ~ compoundName
    )) %>%
    dplyr::mutate(etherPlasmalogenCompoundName = ifelse(
      is.na(plasmalogen_type),
      compoundName,
      glue::glue(
        "{class}({single}:{double}e/{sn2_chain})",
        class = lipidClass,
        single = single_1,
        double = double_1,
        sn2_chain = FA2
      )
    )) %>%
    dplyr::mutate(etherPlasmalogenSumComposition = ifelse(
      is.na(plasmalogen_type),
      sumComposition,
      glue::glue(
        "{class}({single}:{double}e)",
        class = lipidClass,
        single = total_single,
        double = total_double
      )
    ))

  compound_components
}
