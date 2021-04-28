
#' Write flat tidy list
#'
#' @inheritParams test_mzroll_list
#' @param dir_path path to save outputs
#' @param name_preamble start of output file name
#'
#' @return
#' Export three tables:
#' \itemize{
#'     \item{features: one row per features measured (i.e., a metabolite)}
#'     \item{sample: one row per sample}
#'     \item{measurements: one row per measurement (i.e., one metabolite in
#'       one sample)}
#' }
#'
#' @export
write_tidy_list_triple <- function(mzroll_list, dir_path, name_preamble) {
  test_mzroll_list(mzroll_list)

  stopifnot(class(dir_path) == "character", length(dir_path) == 1)
  if (!file.exists(dir_path)) {
    stop("output directory does not exist: ", dir_path)
  }

  stopifnot(class(name_preamble) == "character", length(name_preamble) == 1)

  for (k in names(mzroll_list)) {
    readr::write_tsv(
      mzroll_list[[k]],
      path = file.path(dir_path, paste0(name_preamble, "_", k, ".tsv"))
    )
  }

  invisible(0)
}

#' Write flat tidy list
#'
#' @inheritParams write_tidy_list_triple
#'
#' @return
#' Export one table which is one row per peak, which includes all feature and
#'   sample attributes.
#'
#' @export
write_tidy_list_augmented <- function(mzroll_list, dir_path, name_preamble) {
  test_mzroll_list(mzroll_list)

  stopifnot(class(dir_path) == "character", length(dir_path) == 1)
  if (!file.exists(dir_path)) {
    stop("output directory does not exist: ", dir_path)
  }

  stopifnot(class(name_preamble) == "character", length(name_preamble) == 1)

  mega_table <- mzroll_list$peaks %>%
    dplyr::left_join(mzroll_list$peakgroups, by = "groupId") %>%
    dplyr::left_join(mzroll_list$samples %>%
      dplyr::select(!!!rlang::syms(setdiff(colnames(.), "method_tag"))),
    by = "sampleId"
    )

  readr::write_tsv(
    mega_table,
    path = file.path(dir_path, paste0(name_preamble, "_augmented_table.tsv"))
  )

  invisible(0)
}

#' Write wide output
#'
#' abundances form a matrix with metabolites as rows and samples as columns.
#'   Use transpose to treat samples as rows
#'
#' @inheritParams write_tidy_list_triple
#' @inheritParams diffex_mzroll
#' @param transpose if TRUE then samples will be stored as rows
#'
#' @return
#' Export one table which contains metabolites as rows and samples as columns.
#'
#' @export
write_tidy_list_wide <- function(
  mzroll_list,
  dir_path,
  name_preamble,
  value.var = "log2_abundance",
  transpose = FALSE
  ) {
  test_mzroll_list(mzroll_list)

  stopifnot(class(dir_path) == "character", length(dir_path) == 1)
  if (!file.exists(dir_path)) {
    stop("output directory does not exist: ", dir_path)
  }

  stopifnot(class(name_preamble) == "character", length(name_preamble) == 1)
  stopifnot(class(transpose) == "logical", length(transpose) == 1)

  # structure abundances
  if (transpose) {
    cast_formula <- stats::as.formula("sampleId ~ groupId")
  } else {
    cast_formula <- stats::as.formula("groupId ~ sampleId")
  }

  abundance_matrix <- mzroll_list$peaks %>%
    reshape2::acast(formula = cast_formula, value.var = value.var)

  if (transpose) {
    feature_labels <- rownames(abundance_matrix)
    sample_labels <- colnames(abundance_matrix)
  } else {
    feature_labels <- colnames(abundance_matrix)
    sample_labels <- rownames(abundance_matrix)
  }

  # create a top-left-null section
  n_peakgroup_attr <- ncol(mzroll_list$peakgroups)
  n_sample_attr <- ncol(mzroll_list$samples)

  top_left_void <- matrix(nrow = n_sample_attr, ncol = n_peakgroup_attr)
  if (transpose) {
    top_left_void <- t(top_left_void)
  }
  # remove one row to leave space for row attribute labels
  top_left_void <- top_left_void[-1, , drop = FALSE]

  # setup sample and peakgroup attributes

  ordered_samples <- mzroll_list$samples %>%
    dplyr::mutate(sampleId = factor(sampleId, levels = feature_labels)) %>%
    dplyr::arrange(sampleId) %>%
    dplyr::mutate_all(as.character)

  ordered_groups <- mzroll_list$peakgroups %>%
    dplyr::mutate(groupId = factor(groupId, levels = sample_labels)) %>%
    dplyr::arrange(groupId) %>%
    dplyr::mutate_if(is.numeric, round, 3) %>%
    dplyr::mutate_all(as.character)

  if (transpose) {
    left_matrix <- rbind(
      top_left_void,
      matrix(colnames(ordered_samples), nrow = 1),
      unname(as.matrix(ordered_samples))
    )


    top_right_matrix <- rbind(
      matrix(colnames(ordered_groups), nrow = 1),
      unname(as.matrix(ordered_groups))
    ) %>%
      t()

    bottom_right_matrix <- cbind(
      matrix(rownames(abundance_matrix), ncol = 1),
      abundance_matrix
    )

    right_matrix <- rbind(
      top_right_matrix,
      bottom_right_matrix
    )

    output <- cbind(left_matrix, right_matrix)
  } else {
    left_matrix <- rbind(
      top_left_void,
      matrix(colnames(ordered_groups), nrow = 1),
      unname(as.matrix(ordered_groups))
    )


    top_right_matrix <- rbind(
      matrix(colnames(ordered_samples), nrow = 1),
      unname(as.matrix(ordered_samples))
    ) %>%
      t()

    bottom_right_matrix <- cbind(
      matrix(rownames(abundance_matrix), ncol = 1),
      abundance_matrix
    )

    right_matrix <- rbind(
      top_right_matrix,
      bottom_right_matrix
    )

    output <- cbind(left_matrix, right_matrix)
  }

  output %>%
    as.data.frame() %>%
    readr::write_tsv(
      path = file.path(dir_path, paste0(name_preamble, "_", "wide.tsv")),
      col_names = FALSE
    )

  invisible(0)
}
