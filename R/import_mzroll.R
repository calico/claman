#' Process mzRoll
#'
#' @param mzroll_db_path path to mzroll DB file
#' @param only_identified TRUE/FALSE, filter to only features which were
#'   identified.
#' @param validate TRUE/FALSE, use meta-data to only name the subset of
#'   features which were manually validated.
#' @param method_tag what method was run (for the purpose of matching
#'   meta-data and aggregating).
#' @param peakgroup_labels_to_keep Peak groups containing any one of the
#'   characters in this string are retained.
#' \code{default = "*"} (retain all peak groups)
#' @param peakgroup_labels_to_exclude Peak groups containing any one of
#'   the characters in this string are excluded.  Exclusion takes precedence
#'   over inclusion. \code{default = ""} (do not exclude any peak groups).
#' @param quant_col name of quant column to use for standardization from the
#'   peaks table.
#' \code{default = "peakAreaTop"}
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
#' process_mzroll(nplug_mzroll())
process_mzroll <- function(mzroll_db_path,
                           only_identified = TRUE,
                           validate = FALSE,
                           method_tag = "",
                           peakgroup_labels_to_keep = "*",
                           peakgroup_labels_to_exclude = "",
                           quant_col = "peakAreaTop") {
  checkmate::assertFileExists(mzroll_db_path)
  checkmate::assertLogical(only_identified, len = 1)
  checkmate::assertLogical(validate, len = 1)
  checkmate::assertString(method_tag)

  # connect to SQLite .mzrollDB database
  mzroll_db_con <- create_sqlite_con(mzroll_db_path)

  # add peakgroup m/z and rt
  # Issue 14: option to specify quant type other than peakAreaTop
  peakgroups <- process_mzroll_load_peakgroups(mzroll_db_con, quant_col)

  # if validate is true, use peakgroup labels to rename unknowns
  peakgroups <- process_mzroll_validate_peakgroups(peakgroups, validate)

  # if only_identified is true, only named compounds are retained, if false
  #   unknowns are labeled by m/z and rt
  peakgroups <- process_mzroll_identify_peakgroups(peakgroups, only_identified)

  # Issue 13: filter list based on peak group labels
  peakgroups <- peakgroups %>%
    dplyr::rowwise() %>%
    dplyr::filter(is_has_label(label, peakgroup_labels_to_keep, peakgroup_labels_to_exclude)) %>%
    dplyr::ungroup()

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
      dplyr::arrange(groupId) %>%
      dplyr::mutate(groupId = factor(groupId, levels = groupId))
  }

  debugr::dwatch(
    msg = "Before summarizing distinct peaks and samples...[calicomics<import_mzroll.R>::process_mzroll]\n"
  )

  # summarize distinct peaks and samples

  samples <- dplyr::tbl(
    mzroll_db_con,
    dbplyr::sql("SELECT sampleId, name, filename FROM samples")
  ) %>%
    dplyr::collect() %>%
    dplyr::arrange(sampleId) %>%
    dplyr::mutate(
      sampleId = factor(sampleId, levels = sampleId)
    )

  debugr::dwatch(
    msg = "Summarized samples. [calicomics<import_mzroll.R>::process_mzroll]\n"
  )

  peaks <- dplyr::tbl(
    mzroll_db_con,
    dbplyr::sql(paste0("SELECT groupId, sampleId, ", quant_col, " FROM peaks"))
  ) %>%
    dplyr::collect() %>%
    dplyr::mutate(quant_value = !!sym(quant_col)) %>%
    dplyr::mutate(quant_value = as.numeric(quant_value)) %>%
    dplyr::semi_join(peakgroups, by = "groupId") %>%
    dplyr::group_by(groupId) %>%
    dplyr::mutate(
      log2_abundance = log2(quant_value),
      centered_log2_abundance = log2_abundance - mean(log2_abundance)
    ) %>%
    dplyr::select(
      groupId,
      sampleId,
      log2_abundance,
      centered_log2_abundance
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      groupId = factor(
        groupId,
        levels = levels(reduced_peakgroups$groupId)
      ),
      sampleId = factor(
        sampleId,
        levels = levels(samples$sampleId)
      )
    )

  debugr::dwatch(
    msg = "Summarized peaks. [calicomics<import_mzroll.R>::process_mzroll]\n"
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

#' Process Metabolon
#'
#' Import results generated by metabolon into romic structure, using series of CSV
#' output files delivered by Metabolon as inputs.
#'
#' @param peak_areas path to peak areas CSV file from metabolon results,
#' or table containing loaded peak areas
#' @param chem_anno path to chemical annotations CSV file from metabolon results,
#' or table containing loaded chemical annotations
#' @param metadata_info path to metadata CSV file from metabolon results,
#' or table containing loaded metadata
#'
#' @param peak_floor_intensity value to use to replace NA peak measurements (missing values)
#'
#' @details
#' An example dataset is available in the \code{MetabolonR} package,\cr
#' which may be installed via \code{remotes::install_github("miepstei/MetabolonR")}\cr
#'
#' Once installed, execute the function via\cr\cr
#' \code{data <- process_metabolon(}\cr
#' \code{MetabolonR::metaboliteDataWide,}\cr
#' \code{MetabolonR::metaboliteInfo,}\cr
#' \code{MetabolonR::sampleInfo}\cr
#' \code{)}\cr
#' @export
process_metabolon <- function(peak_areas,
                              chem_anno,
                              metadata_info,
                              peak_floor_intensity = 4096) {
  checkmate::assertNumeric(peak_floor_intensity)

  # data tables
  if (inherits(peak_areas, "data.frame") &&
    inherits(chem_anno, "character") &&
    inherits(metadata_info, "character")) {
    peaks <- peak_areas
    anno <- chem_anno
    metadata <- metadata_info

    # rename columns to handle different versions of Metabolon
    if (!"PARENT_SAMPLE_NAME" %in% colnames(peaks) && "SAMPLE_NAME" %in% colnames(peaks)) {
      peaks <- peaks %>% dplyr::rename(PARENT_SAMPLE_NAME = SAMPLE_NAME)
    }

    if (!"PARENT_SAMPLE_NAME" %in% colnames(metadata) && "SAMPLE_NAME" %in% colnames(metadata)) {
      metadata <- metadata %>% dplyr::rename(PARENT_SAMPLE_NAME = SAMPLE_NAME)
    }

    if (!"CLIENT_SAMPLE_NUMBER" %in% colnames(metadata) && "METAB" %in% colnames(metadata)) {
      metadata <- metadata %>% dplyr::rename(CLIENT_SAMPLE_NUMBER = METAB)
    }

    if (!"CHEM_ID" %in% colnames(anno) && "METABOLITE_NAME" %in% colnames(anno)) {
      anno <- anno %>%
        dplyr::mutate(CHEM_ID = as.integer(gsub("[^0-9.-]", "", METABOLITE_NAME)))
    }

    if (!"SMILES" %in% colnames(anno)) {
      anno <- anno %>% dplyr::mutate(SMILES = "")
    }

    if (!"PLATFORM" %in% colnames(anno)) {
      anno <- anno %>% dplyr::mutate(PLATFORM = "")
    }

    if (!"TYPE" %in% colnames(anno)) {
      anno <- anno %>% dplyr::mutate(TYPE = "")
    }
  } else if (inherits(peak_areas, "character") &&
    inherits(chem_anno, "character") &&
    inherits(metadata_info, "character")) {
    checkmate::assertFileExists(peak_areas)
    checkmate::assertFileExists(chem_anno)
    checkmate::assertFileExists(metadata_info)

    # import metabolon results as CSV files
    peaks <- readr::read_csv(peak_areas, show_col_types = FALSE)
    anno <- readr::read_csv(chem_anno, col_types = readr::cols(.default = "c"))
    metadata <- readr::read_csv(metadata_info, col_types = readr::cols(.default = "c"))
  } else {
    stop("Invalid input types - first three arguments must be 3 paths to Metabolon CSV files, or 3 data frames of loaded Metabolon results.\n")
  }

  columns_to_pivot <- colnames(peaks)
  columns_to_pivot <- columns_to_pivot[2:length(columns_to_pivot)]
  samples <- metadata %>%
    dplyr::select(PARENT_SAMPLE_NAME, CLIENT_SAMPLE_NUMBER) %>%
    dplyr::rename(sampleName = PARENT_SAMPLE_NAME, sampleId = CLIENT_SAMPLE_NUMBER) %>%
    dplyr::mutate(sampleId = factor(sampleId, levels = sampleId))

  peaks_updated_noIds <- peaks %>%
    tidyr::pivot_longer(cols = dplyr::all_of(columns_to_pivot)) %>%
    dplyr::mutate(
      groupId = as.integer(gsub("[^0-9.-]", "", name)), # some versions of Metabolon include non-numeric characters for groupId
      sampleName = PARENT_SAMPLE_NAME,
      peakAreaTop = .data$value
    ) %>% # claman is hard-coded to use peakAreaTop as quant method from mzrolldb files
    dplyr::select(-name, -PARENT_SAMPLE_NAME, -.data$value) %>%
    dplyr::left_join(samples, by = "sampleName") %>%
    dplyr::mutate(peakAreaTop = ifelse(is.na(peakAreaTop), peak_floor_intensity, peakAreaTop)) %>%
    dplyr::group_by(groupId) %>%
    dplyr::mutate(
      log2_abundance = log2(peakAreaTop),
      centered_log2_abundance = log2_abundance - mean(log2_abundance)
    ) %>%
    dplyr::ungroup()

  peakId <- seq(length = nrow(peaks_updated_noIds))
  peaks_updated <- cbind(peakId, peaks_updated_noIds)

  peakgroups <- anno %>%
    dplyr::select(CHEM_ID, CHEMICAL_NAME, SMILES, PLATFORM, TYPE) %>%
    dplyr::rename(
      groupId = CHEM_ID,
      compoundName = CHEMICAL_NAME,
      smiles = SMILES,
      method_tag = PLATFORM,
    )

  peakgroups_updated <- peakgroups %>%
    dplyr::mutate(
      tagString = " ",
      mz = NA,
      rt = NA,
      compoundDB = "Metabolon",
      searchTableName = " ",
      label = "g", # assume all Metabolon IDs are good
      groupId = factor(groupId, levels = groupId)
    )

  samples_updated <- samples %>%
    dplyr::rename(name = sampleName)

  peaks_updated$groupId <- factor(
    peaks_updated$groupId,
    levels = levels(peakgroups_updated$groupId)
  )

  peaks_updated$sampleId <- factor(
    peaks_updated$sampleId,
    levels = levels(samples_updated$sampleId)
  )

  # combine metabolon data into mzroll_list with added metadata
  mzroll_list <- romic::create_triple_omic(
    measurement_df = peaks_updated,
    feature_df = peakgroups_updated,
    sample_df = samples_updated,
    feature_pk = "groupId",
    sample_pk = "sampleId",
    omic_type_tag = "mzroll"
  )

  mzroll_list <- merge_samples_tbl(
    mzroll_list = mzroll_list,
    samples_tbl = metadata,
    id_strings = c("PARENT_SAMPLE_NAME"),
    exact = TRUE
  )

  test_mzroll_list(mzroll_list)

  return(mzroll_list)
}

#' Create SQLite Connection
#'
#' Open a connection to an SQLite database
#'
#' @param sqlite_path path to sqlite file
#'
#' @return sqlite connection
create_sqlite_con <- function(sqlite_path) {
  checkmate::assertFileExists(sqlite_path)

  sqlite_con <- DBI::dbConnect(
    RSQLite::SQLite(),
    sqlite_path,
    synchronous = NULL
  )

  return(sqlite_con)
}


#' Process MzRoll - Load Peakgroups
#'
#' Load peakgroups table and add mz and rt variables to peakgroups based on
#'   peak-level data.
#'
#' @param mzroll_db_con an SQLite connection to an mzrollDB database
#' @param quant_col name of quant column to use for standardization from the
#'   peaks table.
#' \code{default = "peakAreaTop"}
#'
#' @return a table of peakgroups
process_mzroll_load_peakgroups <- function(mzroll_db_con, quant_col = "peakAreaTop") {
  peakgroup_positions <- dplyr::tbl(mzroll_db_con, "peaks") %>%
    dplyr::collect() %>%
    dplyr::group_by(groupId) %>%
    dplyr::mutate(quant_value := !!sym(quant_col)) %>%
    dplyr::mutate(quant_value = as.numeric(quant_value)) %>%
    dplyr::summarize(
      mz = sum(peakMz * quant_value) / sum(quant_value),
      rt = sum(rt * quant_value) / sum(quant_value)
    )

  # Issue-19
  # Fixes for R CMD check
  peakgroups <- dplyr::tbl(mzroll_db_con, "peakgroups") %>%
    dplyr::collect() %>%
    dplyr::select(-compoundId)

  # Only separate SMILES if the pattern exists in compoundName
  if (any(stringr::str_detect(peakgroups$compoundName, ": "))) {
    peakgroups <- peakgroups %>%
      tidyr::separate(
        compoundName,
        into = c("compoundName", "smiles"),
        sep = ": ", extra = "drop", fill = "right"
      )
  } else {
    # Add empty smiles column if no separation needed
    peakgroups <- peakgroups %>%
      dplyr::mutate(smiles = NA_character_)
  }

  peakgroups <- peakgroups %>%
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
process_mzroll_validate_peakgroups <- function(peakgroups, validate) {
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
process_mzroll_identify_peakgroups <- function(peakgroups, only_identified) {
  checkmate::assertDataFrame(peakgroups)

  if (only_identified) {
    peakgroups <- peakgroups %>%
      dplyr::filter(!is.na(compoundName))
  } else {
    unknown_groups <- peakgroups %>%
      dplyr::filter(is.na(compoundName)) %>%
      {
        .$groupId
      }

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
#' @inheritParams merge_samples_tbl
#'
#' @return an mzroll_list containing three tibbles:
#' \itemize{
#'   \item{peakgroups: one row per unique analyte (defined by a unique
#'     groupId)},
#'   \item{samples: one row per unique sample (defined by a unique sampleId)},
#'   \item{peaks: one row per peak (samples x peakgroups)}
#'   }
#'
#' @examples
#' mzroll_paths <- tibble::tribble(
#'   ~method_tag, ~mzroll_db_path,
#'   "method1", nplug_mzroll(),
#'   "method2", nplug_mzroll()
#' )
#'
#' process_mzroll_multi(
#'   mzroll_paths,
#'   nplug_samples,
#'   "sample_name",
#'   exact = TRUE
#' )
#' @export

process_mzroll_multi <- function(
    mzroll_paths,
    samples_tbl,
    id_strings,
    only_identified = TRUE,
    validate = FALSE,
    exact = FALSE,
    peakgroup_labels_to_keep = "*",
    peakgroup_labels_to_exclude = "",
    quant_col = "peakAreaTop") {
  checkmate::assertDataFrame(mzroll_paths, min.rows = 2)
  if (!all(colnames(mzroll_paths) == c("method_tag", "mzroll_db_path"))) {
    stop("mzroll_paths must contain two columns: method_tag & mnzroll_db_path")
  }
  stopifnot(length(mzroll_paths$method_tag) == length(unique(mzroll_paths$method_tag)))

  checkmate::assertDataFrame(samples_tbl)
  checkmate::assertLogical(only_identified, len = 1)
  checkmate::assertLogical(validate, len = 1)

  mzroll_list_nest <- mzroll_paths %>%
    dplyr::mutate(
      # read each dataset
      mzroll_list = purrr::map2(
        mzroll_db_path,
        method_tag,
        process_mzroll,
        only_identified = only_identified,
        validate = validate,
        peakgroup_labels_to_keep = peakgroup_labels_to_keep,
        peakgroup_labels_to_exclude = peakgroup_labels_to_exclude,
        quant_col = quant_col
      ),
      # add sample meta-data
      mzroll_list = purrr::map(
        mzroll_list,
        merge_samples_tbl,
        samples_tbl = samples_tbl,
        id_strings = id_strings,
        exact = exact
      )
    )

  aggregate_mzroll_list <- aggregate_mzroll_nest(
    mzroll_list_nest,
    samples_tbl
  )

  mzroll_multi_qc(aggregate_mzroll_list)

  # Issue 7: remove name and samples_tbl_row columns to mimic calicomics output
  aggregate_mzroll_list$samples <- aggregate_mzroll_list$samples %>%
    dplyr::select(-name, -samples_tbl_row)

  aggregate_mzroll_list <- romic::update_tomic(
    aggregate_mzroll_list,
    aggregate_mzroll_list$samples
  )

  return(aggregate_mzroll_list)
}

#' Aggregate mzroll lists
#'
#' Combine mzroll lists into a single mzroll_list
#'
#' @param mzroll_list_nest a nested list of mzroll_lists produced from
#'   \link{process_mzroll}.
#' @inheritParams merge_samples_tbl
#'
#' @return a single mzroll_list
aggregate_mzroll_nest <- function(mzroll_list_nest, samples_tbl) {
  checkmate::assertDataFrame(mzroll_list_nest)
  # check that all mzroll lists have the same design

  designs <- purrr::map(mzroll_list_nest$mzroll_list, function(x) {
    x$design
  })
  if (length(unique(designs)) > 1) {
    stop("mzroll_lists have different designs and cannot be aggregated")
  }

  # identify samples with shared "samples_tbl_row"

  mzroll_list_all_samples <- mzroll_list_nest %>%
    dplyr::mutate(samples = purrr::map(mzroll_list, function(x) {
      x$samples
    })) %>%
    dplyr::select(-method_tag, -mzroll_list) %>%
    tidyr::unnest(samples)

  consensus_samples <- mzroll_list_all_samples %>%
    dplyr::mutate(sampleId = samples_tbl_row) %>%
    dplyr::select(sampleId, name, samples_tbl_row, !!!rlang::syms(colnames(samples_tbl))) %>%
    dplyr::group_by(sampleId) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sampleId = factor(sampleId, levels = sampleId))

  # iterate through mzrollDB

  max_groupId <- 0
  mzroll_indices <- seq_len(nrow(mzroll_list_nest))

  # samples defined upfront as the union of samples from the sample sheet
  # that are observed "samples_tbl_row" will be the new sampleId
  updated_mzroll_list <- list(
    features = NULL,
    samples = consensus_samples,
    measurements = NULL
  )

  # peaks and peakgroups are updated to unique values across each method
  # and sampleId is updated to "samples_tbl_row"

  for (a_mzroll_index in mzroll_indices) {
    one_mzroll_list <- mzroll_list_nest$mzroll_list[[a_mzroll_index]]
    checkmate::assertClass(one_mzroll_list, "tomic")
    checkmate::assertClass(one_mzroll_list, "mzroll")

    # update groupId to new unique values
    # temporarily convert from factors -> integers so they can be added

    one_mzroll_list[["features"]]$groupId <-
      as.integer(one_mzroll_list[["features"]]$groupId) + max_groupId
    one_mzroll_list[["measurements"]]$groupId <-
      as.integer(one_mzroll_list[["measurements"]]$groupId) + max_groupId

    # find sample updates sampleId -> samples_tbl_row

    sampleId_lookup <- one_mzroll_list$samples %>%
      dplyr::select(oldSampleId = sampleId, newSampleId = samples_tbl_row)

    one_mzroll_list$measurements <- one_mzroll_list$measurements %>%
      dplyr::left_join(sampleId_lookup, by = c("sampleId" = "oldSampleId")) %>%
      dplyr::select(-sampleId) %>%
      dplyr::rename(sampleId = newSampleId) %>%
      dplyr::mutate(sampleId = factor(
        sampleId,
        levels = levels(consensus_samples$sampleId)
      )) %>%
      dplyr::select(!!!rlang::syms(colnames(one_mzroll_list$measurements)))

    updated_mzroll_list$measurements <- dplyr::bind_rows(
      updated_mzroll_list$measurements,
      one_mzroll_list$measurements
    )
    updated_mzroll_list$features <- dplyr::bind_rows(
      updated_mzroll_list$features,
      one_mzroll_list$features
    )

    max_groupId <- max(one_mzroll_list[["features"]]$groupId)
  }

  # overwrite one of the original mzrolls to update the romic schema

  one_mzroll_list$features <- updated_mzroll_list$features %>%
    dplyr::mutate(groupId = factor(groupId, levels = groupId))
  one_mzroll_list$measurements <- updated_mzroll_list$measurements %>%
    dplyr::mutate(
      groupId = factor(
        groupId,
        levels = levels(one_mzroll_list$features$groupId)
      )
    )
  # update values and schema
  one_mzroll_list <- romic::update_tomic(
    one_mzroll_list,
    updated_mzroll_list$samples
  )

  return(one_mzroll_list)
}

mzroll_multi_qc <- function(mzroll_list) {
  mzroll_coverage_table <- mzroll_list$measurements %>%
    dplyr::left_join(
      mzroll_list$features %>%
        dplyr::select(groupId, method_tag),
      by = "groupId"
    ) %>%
    dplyr::left_join(
      mzroll_list$samples %>%
        dplyr::select(sampleId, samples_tbl_row),
      by = "sampleId"
    ) %>%
    dplyr::count(samples_tbl_row, method_tag) %>%
    tidyr::spread(method_tag, n)

  missed_matches <- mzroll_coverage_table %>%
    dplyr::filter_all(dplyr::any_vars(is.na(.)))

  if (nrow(missed_matches) != 0) {
    missed_matches_data <- mzroll_list$samples %>%
      dplyr::filter(samples_tbl_row %in% missed_matches$samples_tbl_row)

    missed_matches_w_data <- dplyr::inner_join(missed_matches, missed_matches_data,
      by = c("samples_tbl_row")
    ) %>%
      dplyr::mutate(out = paste0(
        "Sample Description: '",
        `sample description`,
        "', Sample Variable 1: '",
        `sample variable 1`,
        "', Sample Variable 2: '",
        `sample variable 2`,
        "'\n"
      ))

    missed_matches_warning <- missed_matches %>%
      tidyr::gather(
        -samples_tbl_row,
        key = "method_tag",
        value = "n_entries"
      ) %>%
      dplyr::filter(is.na(n_entries)) %>%
      dplyr::mutate(out = glue::glue(
        "row: {samples_tbl_row}, method_tag: {method_tag}"
      )) %>%
      {
        paste(.$out, collapse = "\n")
      }

    stop(
      nrow(missed_matches),
      " rows only matched a subset of methods;
        this will cause downstream problems.
        Either update your sample sheet or analyze each method separately.
        Details:\n",
      missed_matches_warning, "\n",
      "Metadata associated with samples missing matches:\n",
      missed_matches_w_data$out
    )
  }

  return(invisible(0))
}
