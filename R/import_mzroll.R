#' Process mzRoll
#'
#' @param mzroll_db_path path to mzroll DB file
#' @param only_identified TRUE/FALSE, filter to only features which were
#'   identified.
#' @param validate TRUE/FALSE, use meta-data to only name the subset of
#'   features which were manually validated.
#' @param method_tag what method was run (for the purpose of matching
#'   meta-data and aggregating).
#'
#' @return a **triple_omic** from **romic** containing three tibbles:
#' \itemize{
#'   \item{features: one row per unique analyte (defined by a unique
#'     groupId)},
#'   \item{samples: one row per unique sample (defined by a unique sampleId)},
#'   \item{measurements: one row per peak (samples x peakgroups)}
#'   }
#'  And, a design list which tracks the variables in each table.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' mzroll_db_path = nplug_mzroll()
#' process_mzroll(mzroll_db_path)
process_mzroll <- function(
  mzroll_db_path,
  only_identified = TRUE,
  validate = FALSE,
  method_tag = ""
  ) {
  
  checkmate::assertFileExists(mzroll_db_path)
  checkmate::assertLogical(only_identified, len = 1)
  checkmate::assertLogical(validate, len = 1)
  checkmate::assertString(method_tag)
  
  mzroll_db_con <- authutils::create_sqlite_con(mzroll_db_path)
  
  # add peakgroup m/z and rt
  peakgroups <- process_mzroll_load_peakgroups(mzroll_db_con)
  # if validate is true, use peakgroup labels to rename unknowns
  peakgroups <- process_mzroll_validate_peakgroups(peakgroups, validate)
  # if only_identified is true, only named compounds are retained, if false
  #   unknowns are labeled by m/z and rt
  peakgroups <- process_mzroll_identify_peakgroups(peakgroups, only_identified)

  # reduce to smaller number of peakgroups features

  possible_peakgroup_reduced_vars <- c(
    "groupId",
    "compoundName",
    "smiles",
    "adductName",
    "tagString",
    "is_unknown",
    "mz",
    "rt",
    "systematicCompoundId",
    "compoundDB",
    "searchTableName",
    "label"
  )

  detected_peakgroup_reduced_vars <- intersect(
    possible_peakgroup_reduced_vars,
    colnames(peakgroups)
  )

  if (nrow(peakgroups) == 0) {
    reduced_peakgroups <- peakgroups %>%
      dplyr::select(!!!rlang::syms(detected_peakgroup_reduced_vars)) %>%
      dplyr::mutate(
        peak_label = "",
        method_tag = method_tag
      )

    warning("No named compounds were found; an empty dataset will be returned")
  } else {
    reduced_peakgroups <- peakgroups %>%
      dplyr::select(!!!rlang::syms(detected_peakgroup_reduced_vars)) %>%
      dplyr::group_by(compoundName) %>%
      dplyr::mutate(peak_label = dplyr::case_when(
        dplyr::n() == 1 ~ compoundName,
        TRUE ~ paste0(compoundName, " (", 1:dplyr::n(), ")")
      )) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(method_tag = method_tag) %>%
      dplyr::arrange(groupId)
  }

  debugr::dwatch(
    msg = "Before summarizing distinct peaks and samples... [calicomics<import_mzroll.R>::process_mzroll]\n"
  )

  # summarize distinct peaks and samples

  peaks <- dplyr::tbl(
    mzroll_db_con,
    dbplyr::sql("SELECT peakId, groupId, sampleId, peakAreaTop FROM peaks")
    ) %>%
    dplyr::collect() %>%
    dplyr::semi_join(peakgroups, by = "groupId") %>%
    dplyr::group_by(groupId) %>%
    dplyr::mutate(
      log2_abundance = log2(peakAreaTop),
      centered_log2_abundance = log2_abundance - mean(log2_abundance)
    ) %>%
    dplyr::select(
      peakId,
      groupId,
      sampleId,
      log2_abundance,
      centered_log2_abundance
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sampleId = as.character(sampleId)) %>%
    dplyr::arrange(peakId)

  debugr::dwatch(
    msg = "Summarized peaks. [calicomics<import_mzroll.R>::process_mzroll]\n"
  )

  samples <- dplyr::tbl(
    mzroll_db_con,
    dbplyr::sql("SELECT sampleId, name, filename FROM samples")
  ) %>%
    dplyr::collect() %>%
    dplyr::arrange(sampleId) %>%
    dplyr::mutate(
      sampleId = as.character(sampleId)
    )
    
  debugr::dwatch(
    msg = "Summarized samples. [calicomics<import_mzroll.R>::process_mzroll]\n"
  )

  # disconnect for sql-lite
  DBI::dbDisconnect(mzroll_db_con)

  debugr::dwatch(
    msg = "About to create list... [calicomics<import_mzroll.R>::process_mzroll]\n"
  )

  mzroll_list <- romic::create_triple_omic(
    measurement_df = peaks,
    feature_df = reduced_peakgroups,
    sample_df = samples,
    feature_pk = "groupId",
    sample_pk = "sampleId",
    omic_type_tag = "mzroll"
  )
  
  debugr::dwatch(
    msg = "Created list. [calicomics<import_mzroll.R>::process_mzroll]\n"
  )

  debugr::dwatch(msg = glue::glue(
    "number of samples in mzroll_list$samples: {nrow(mzroll_list$samples)}"
    ))

  test_mzroll_list(mzroll_list)
  
  return(mzroll_list)
}

#' Process MzRoll - Load Peakgroups
#' 
#' Load peakgroups table and add mz and rt variables to peakgroups based on
#'   peak-level data.
#' 
#' @param mzroll_db_con an SQLite connection to an mzrollDB database
#' 
#' @return a table of peakgroups
process_mzroll_load_peakgroups <- function(mzroll_db_con){
  
  peakgroup_positions <- dplyr::tbl(mzroll_db_con, "peaks") %>%
    dplyr::collect() %>%
    dplyr::group_by(groupId) %>%
    dplyr::summarize(
      mz = sum(peakMz * peakAreaTop) / sum(peakAreaTop),
      rt = sum(rt * peakAreaTop) / sum(peakAreaTop)
    )
  
  peakgroups <- dplyr::tbl(mzroll_db_con, "peakgroups") %>%
    dplyr::collect() %>%
    dplyr::select(-compoundId) %>%
    # peel off smiles if they are present
    tidyr::separate(
      compoundName,
      into = c("compoundName", "smiles"),
      sep = ": ", extra = "drop", fill = "right"
    ) %>%
    dplyr::left_join(peakgroup_positions, by = "groupId")
  
  # if provided, prefer display name over compoundName
  if ("displayName" %in% colnames(peakgroups)) {
    peakgroups <- peakgroups %>%
      dplyr::mutate(compoundName = ifelse(
        is.na(displayName),
        compoundName,
        displayName
      ))
  }
  
  return(peakgroups)
}

#' Process MzRoll - Validate Peakgroups
#' 
#' If validate is True then remove unvalidated compoundNames (based on
#'   peakgroup labels) and identify all unvalidated compounds with an
#'   "is_unknown" variable.
#'   
#' @param peakgroups a table of distinct ions with characteristic m/z and rt
#' @inheritParams process_mzroll
#'
#' @return a table of peakgroups
process_mzroll_validate_peakgroups <- function(peakgroups, validate){
  
  checkmate::assertDataFrame(peakgroups)
  checkmate::assertLogical(validate, len = 1)
  
  if (validate) {
    peakgroups <- peakgroups %>%
      dplyr::mutate(
        compoundName = dplyr::case_when(
          is.na(compoundName) ~ NA_character_,
          searchTableName == "Bookmarks" &
            !stringr::str_detect(string = label, pattern = "b") ~ compoundName,
          # backwards compatibility
          searchTableName == "rumsDB" &
            stringr::str_detect(string = label, pattern = "g") ~ compoundName,
          searchTableName == "clamDB" &
            stringr::str_detect(string = label, pattern = "g") ~ compoundName,
          # directly added EICs
          searchTableName == "EICs" ~ compoundName,
          TRUE ~ NA_character_
        ),
        is_unknown = ifelse(is.na(compoundName), TRUE, FALSE)
      )
  }
  
  return(peakgroups)
  
}

#' Process MzRoll - Identify Peakgroups
#' 
#' Either remove unidentified peakgroups or name unknowns using m/z and rt.
#' 
#' @inheritParams process_mzroll_validate_peakgroups
#' @inheritParams process_mzroll
#' 
#' @return a table of peakgroups
process_mzroll_identify_peakgroups <- function(peakgroups, only_identified){
  
  checkmate::assertDataFrame(peakgroups)
  
  if (only_identified) {
    peakgroups <- peakgroups %>%
      dplyr::filter(!is.na(compoundName))
  } else {
    unknown_groups <- peakgroups %>%
      dplyr::filter(is.na(compoundName)) %>%
      {.$groupId}
    
    # label unknowns based on m/z and rt
    
    unknown_names <- peakgroups %>%
      dplyr::filter(groupId %in% unknown_groups) %>%
      dplyr::mutate(new_compoundName = as.character(glue::glue(
        "unk {round(mz, 3)} @ {round(rt,1)}"
      ))) %>%
      dplyr::select(groupId, new_compoundName)
    
    peakgroups <- peakgroups %>%
      dplyr::left_join(unknown_names, by = "groupId") %>%
      dplyr::mutate(
        is_unknown = ifelse(is.na(compoundName), TRUE, FALSE),
        compoundName = dplyr::case_when(
          !is.na(new_compoundName) ~ new_compoundName,
          TRUE ~ compoundName
        )
      ) %>%
      dplyr::select(-new_compoundName)
  }
  
  return(peakgroups)
}

#' Process mzroll multi
#'
#' Aggregate multiple mzrollDB datasets which possess the same sample metadata
#'
#' @param mzroll_paths a tibble with two variables:
#' \itemize{
#'   \item{method_tag: a character vector to tag each dataset with},
#'   \item{mzroll_db_path: path to a mzrollDB file}
#'   }
#' @inheritParams process_mzroll
#'
#' @return an mzroll_list containing three tibbles:
#' \itemize{
#'   \item{peakgroups: one row per unique analyte (defined by a unique
#'     groupId)},
#'   \item{samples: one row per unique sample (defined by a unique sampleId)},
#'   \item{peaks: one row per peak (samples x peakgroups)}
#'   }
#'
#' @export
process_mzroll_multi <- function(
  mzroll_paths,
  standard_databases,
  sample_sheet_list,
  only_identified = TRUE,
  validate = FALSE
  ) {
  checkmate::assertDataFrame(mzroll_paths)
  if (!all(colnames(mzroll_paths) == c("method_tag", "mzroll_db_path"))) {
    stop("mzroll_paths must contain two columns: method_tag & mnzroll_db_path")
  }
  stopifnot(length(mzroll_paths$method_tag) == length(unique(mzroll_paths$method_tag)))

  checkmate::assertClass(sample_sheet_list, "list")
  stopifnot(names(sample_sheet_list) == c("tracking_sheet_id", "sample_sheet"))

  # call process_mzroll once per dataset
  mzroll_list_nest <- mzroll_paths %>%
    dplyr::mutate(processed_mzroll = purrr::map2(
      mzroll_db_path,
      method_tag,
      process_mzroll,
      standard_databases = standard_databases,
      sample_sheet_list = sample_sheet_list,
      only_identified = only_identified,
      validate = validate
    ))

  aggregate_mzroll_list <- aggregate_mzroll_nest(mzroll_list_nest)

  mzroll_multi_qc(aggregate_mzroll_list)

  return(aggregate_mzroll_list)
}

#' Aggregate mzroll lists
#'
#' Combine mzroll lists into a single mzroll_list
#'
#' @param mzroll_list_nest a nested list of mzroll_lists produced from
#'   \link{process_mzroll}.
#'
#' @return a single mzroll_list
aggregate_mzroll_nest <- function(mzroll_list_nest) {

  # identify samples with shared "tube label"

  mzroll_list_all_samples <- mzroll_list_nest %>%
    dplyr::mutate(samples = purrr::map(processed_mzroll, function(x) {
      x$samples
    })) %>%
    dplyr::select(-method_tag, -processed_mzroll) %>%
    tidyr::unnest(samples)

  consensus_samples <- mzroll_list_all_samples %>%
    dplyr::mutate(sampleId = `tube label`) %>%
    dplyr::select(-method_tag, -mzroll_db_path, -name, -filename) %>%
    dplyr::group_by(sampleId) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # iterate through mzrollDB

  max_peakId <- 0
  max_groupId <- 0

  mzroll_indecies <- 1:nrow(mzroll_list_nest)

  # samples defined upfront as the union of samples from the sample sheet
  # that are observed "tube label" will be the new sampleId
  updated_mzroll_list <- list(
    peakgroups = NULL,
    samples = consensus_samples,
    peaks = NULL
  )

  # peaks and peakgroups are updated to unique values across each method
  # and sampleId is updated to "tube label"

  for (a_mzroll_index in mzroll_indecies) {
    one_mzroll_list <- mzroll_list_nest$processed_mzroll[[a_mzroll_index]]
    stopifnot(
      all(c("peakgroups", "samples", "peaks") %in% names(one_mzroll_list))
    )

    # update groupid and peakId to new unique values

    one_mzroll_list$peakgroups$groupId <-
      one_mzroll_list$peakgroups$groupId + max_groupId
    one_mzroll_list$peaks$peakId <-
      one_mzroll_list$peaks$peakId + max_peakId
    one_mzroll_list$peaks$groupId <-
      one_mzroll_list$peaks$groupId + max_groupId

    # find sample updates sampleId -> tube label

    sampleId_lookup <- one_mzroll_list$samples %>%
      dplyr::select(oldSampleId = sampleId, newSampleId = `tube label`)

    one_mzroll_list$peaks <- one_mzroll_list$peaks %>%
      dplyr::left_join(sampleId_lookup, by = c("sampleId" = "oldSampleId")) %>%
      dplyr::select(
        peakId,
        groupId,
        sampleId = newSampleId,
        log2_abundance,
        centered_log2_abundance
      )

    updated_mzroll_list$peaks <- dplyr::bind_rows(
      updated_mzroll_list$peaks,
      one_mzroll_list$peaks
    )
    updated_mzroll_list$peakgroups <- dplyr::bind_rows(
      updated_mzroll_list$peakgroups,
      one_mzroll_list$peakgroups
    )

    max_peakId <- max(one_mzroll_list$peaks$peakId)
    max_groupId <- max(one_mzroll_list$peakgroups$groupId)
  }

  updated_mzroll_list
}

mzroll_multi_qc <- function(mzroll_list) {
  mzroll_coverage_table <- mzroll_list$peaks %>%
    dplyr::left_join(
      mzroll_list$peakgroups %>%
        dplyr::select(groupId, method_tag),
      by = "groupId"
    ) %>%
    dplyr::left_join(
      mzroll_list$samples %>%
        dplyr::select(sampleId, `tube label`),
      by = "sampleId"
    ) %>%
    dplyr::count(`tube label`, method_tag) %>%
    tidyr::spread(method_tag, n)

  missed_matches <- mzroll_coverage_table %>%
    dplyr::filter_all(dplyr::any_vars(is.na(.)))

  if (nrow(missed_matches) != 0) {
    missed_matches_warning <- missed_matches %>%
      tidyr::gather(-`tube label`, key = "method_tag", value = "n_entries") %>%
      dplyr::filter(is.na(n_entries)) %>%
      dplyr::mutate(out = glue::glue(
        "tube: {`tube label`}, method_tag: {method_tag}"
      )) %>%
      {
        paste(.$out, collapse = "\n")
      }

    stop(
      nrow(missed_matches),
      " tubes only matched a subset of methods; this will cause downstream problems.
        Either update your sample sheet or analyze each method separately.
        Details:\n",
      missed_matches_warning
    )
  }

  return(invisible(0))
}

#' MOVE ME - Add Peakgroup Annotations
#' 
#' TO DO - MOVE TO AUTHUTILS
#' 
#' @export
add_peakgroup_annotations <- function() {
  
  # add standard and systematic standard data if available
  
  if (class(standard_databases) != "NULL") {
    
    # match compounds to standards
    compounds <- dplyr::tbl(
      standard_databases$mass_spec_standards_con,
      "compounds"
    ) %>%
      dplyr::collect()
    
    peakgroup_compounds <- peakgroups %>%
      dplyr::left_join(compounds %>%
                         dplyr::select(compoundName, compoundId, systematicCompoundId),
                       by = "compoundName"
      )
    
    # add pathway of each compound
    peakgroup_compound_pathways <- peakgroup_compounds %>%
      dplyr::left_join(
        query_systematic_compounds(
          unique(peakgroup_compounds$systematicCompoundId) %>%
            .[!is.na(.)],
          standard_databases$systematic_compounds_con
        ) %>%
          summarize_compound_pathways(
            min_pw_size = 5L,
            focus_pathways = c(
              "Glycolysis / Gluconeogenesis",
              "Citrate cycle (TCA cycle)",
              "Pentose phosphate pathway",
              "Biosynthesis of amino acids",
              "Purine metabolism",
              "Pyrimidine metabolism"
            )
          ),
        by = "systematicCompoundId"
      ) %>%
      dplyr::mutate(
        pathway = ifelse(is.na(pathway), "Other", pathway),
        focus_pathway = ifelse(is.na(focus_pathway), "Other", focus_pathway)
      )
    
    if (!only_identified) {
      peakgroup_compound_pathways <- peakgroup_compound_pathways %>%
        dplyr::mutate(
          pathway = ifelse(is_unknown, "Unidentified", pathway),
          focus_pathway = ifelse(is_unknown, "Unidentified", focus_pathway)
        )
    }
    
    annotated_peakgroups <- peakgroup_compound_pathways
  } else {
    annotated_peakgroups <- peakgroups
  }
}

add_sample_metadata <- function() {
  
  if (class(sample_sheet_list) == "list") {
    samples <- augment_samples_with_samplesheet(samples, sample_sheet_list)
    
    # drop data from unmatched samples
    peaks <- peaks %>%
      dplyr::semi_join(samples, by = "sampleId")
    reduced_peakgroups <- reduced_peakgroups %>%
      dplyr::semi_join(peaks, by = "groupId")
  }
  
}
