#' Omics Heatmap
#'
#' Generate a basic heatmap of all of your metabolomic and lipidomic data
#'
#' @inheritParams test_mzroll_list
#' @param feature.var variable from "peakgroups" to use as a unique feature
#'   label.
#' @param sample.var variable from "samples" to use as a unique sample label.
#' @param value.var which variable in "peaks" to use for quantification.
#' @param cluster_dim dimensions to cluster (row, column, or both)
#' @param change_threshold values with a more extreme absolute change will be
#'   thresholded to this value.
#' @param plot_type plotly (for interactivity) or grob (for a static ggplot)
#'
#' @examples
#' library(dplyr)
#' mzroll_list <- authutils::get_clamr_assets("X0083-M001A.mzrollDB") %>%
#'   process_mzroll(
#'     standard_databases = NULL,
#'     method_tag = "M001A",
#'     only_identified = TRUE
#'   ) %>%
#'   floor_peaks(12)
#'
#' plot_heatmap(
#'   mzroll_list,
#'   sample.var = "sampleId",
#'   feature.var = "peak_label",
#'   change_threshold = 5,
#'   cluster_dim = "rows",
#'   plot_type = "grob"
#' )
#' \dontrun{
#' plot_heatmap(
#'   mzroll_list,
#'   sample.var = "sampleId",
#'   feature.var = "peak_label",
#'   change_threshold = 5,
#'   cluster_dim = "rows"
#' )
#' }
#'
#' @export
plot_heatmap <- function(mzroll_list,
                         feature.var = "groupId",
                         sample.var = "name",
                         value.var = "centered_log2_abundance",
                         cluster_dim = "both",
                         change_threshold = 3,
                         plot_type = "plotly") {
  test_mzroll_list(mzroll_list)

  checkmate::assertString(feature.var)
  checkmate::assertString(sample.var)
  checkmate::assertString(value.var)

  stopifnot(feature.var %in% colnames(mzroll_list$peakgroups))
  stopifnot(sample.var %in% colnames(mzroll_list$samples))
  stopifnot(value.var %in% colnames(mzroll_list$peaks))

  stopifnot(
    class(cluster_dim) == "character",
    cluster_dim %in% c("rows", "columns", "both")
  )
  stopifnot(
    class(change_threshold) %in% c("numeric", "integer"),
    length(change_threshold) == 1,
    change_threshold > 0
  )
  stopifnot(
    class(plot_type) == "character",
    length(plot_type) == 1,
    plot_type %in% c("plotly", "grob")
  )

  # add additional labels

  augmented_peaks <- mzroll_list$peaks %>%
    # add labels if needed
    dplyr::left_join(mzroll_list$samples %>%
      dplyr::select(!!quo(c("sampleId", sample.var))),
    by = "sampleId"
    ) %>%
    dplyr::left_join(mzroll_list$peakgroups %>%
      dplyr::select(!!quo(c("groupId", feature.var))),
    by = "groupId"
    )

  # format list to matrix

  # convert groupId and sampleId to factors so they are ordered appropriately

  cluster_orders <- hclust_order(
    mzroll_list$peaks, "groupId", "sampleId",
    value.var, cluster_dim
  )

  # save classes of sampleId and groupId so appropriate class coercion occurs
  class(cluster_orders$rows) <- class(mzroll_list$peakgroups$groupId)
  class(cluster_orders$columns) <- class(mzroll_list$samples$sampleId)

  # order rows and columns

  distinct_features <- augmented_peaks %>%
    dplyr::distinct(groupId, !!rlang::sym(feature.var))

  if (cluster_dim == "columns") {

    # order by factor or alpha-numerically

    if (class(distinct_features[[feature.var]]) %in% c("factor", "ordered")) {
      # retain previous ordering

      ordered_distinct_features <- distinct_features %>%
        dplyr::arrange(!!rlang::sym(feature.var))
    } else {
      ordered_distinct_features <- distinct_features %>%
        dplyr::arrange(!!rlang::sym(feature.var))
    }

    ordered_distinct_features <- ordered_distinct_features %>%
      dplyr::mutate(ordered_groupId = factor(groupId, levels = groupId))
  } else {

    # order with hclust

    ordered_distinct_features <- distinct_features %>%
      dplyr::left_join(tibble::tibble(groupId = cluster_orders$rows) %>%
        dplyr::mutate(order = 1:dplyr::n()),
      by = "groupId"
      ) %>%
      dplyr::arrange(order) %>%
      dplyr::mutate(ordered_groupId = factor(groupId, levels = groupId))
  }

  distinct_samples <- augmented_peaks %>%
    dplyr::distinct(sampleId, !!rlang::sym(sample.var))

  if (cluster_dim == "rows") {

    # order by factor or alpha-numerically

    if (class(distinct_samples[[sample.var]]) %in% c("factor", "ordered")) {
      # retain previous ordering

      ordered_distinct_samples <- distinct_samples %>%
        dplyr::arrange(!!rlang::sym(sample.var))
    } else {
      ordered_distinct_samples <- distinct_samples %>%
        dplyr::arrange(!!rlang::sym(sample.var))
    }

    ordered_distinct_samples <- ordered_distinct_samples %>%
      dplyr::mutate(ordered_sampleId = factor(sampleId, levels = sampleId))
  } else {

    # order with hclust

    ordered_distinct_samples <- distinct_samples %>%
      dplyr::left_join(tibble::tibble(sampleId = cluster_orders$columns) %>%
        dplyr::mutate(order = 1:dplyr::n()),
      by = "sampleId"
      ) %>%
      dplyr::arrange(order) %>%
      dplyr::mutate(ordered_sampleId = factor(sampleId, levels = sampleId))
  }

  # update labels
  ordered_distinct_features <- ordered_distinct_features %>%
    dplyr::mutate(label = stringr::str_wrap(
      as.character(!!rlang::sym(feature.var)),
      width = 40,
      exdent = 2
    ))

  ordered_distinct_samples <- ordered_distinct_samples %>%
    dplyr::mutate(label = stringr::str_wrap(
      as.character(!!rlang::sym(sample.var)),
      width = 40,
      exdent = 2
    ))

  # setup abundance values
  ordered_augmented_peaks <- augmented_peaks %>%
    # order all rows and columns
    dplyr::left_join(ordered_distinct_features %>%
      dplyr::select(groupId, ordered_groupId),
    by = "groupId"
    ) %>%
    dplyr::left_join(ordered_distinct_samples %>%
      dplyr::select(sampleId, ordered_sampleId),
    by = "sampleId"
    ) %>%
    # threshold max
    dplyr::mutate(abundance = pmax(
      pmin(
        !!rlang::sym(value.var),
        change_threshold
      ),
      -1 * change_threshold
    ))

  heatmap_theme <- theme_minimal() +
    theme(
      text = element_text(size = 16, color = "black"),
      title = element_text(size = 20, color = "black"),
      axis.text.x = element_text(angle = 60, hjust = 1),
      axis.text.y = element_text(size = 4),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 18),
      legend.position = "top",
      strip.background = element_rect(fill = "gray80")
    )

  heatmap_plot <- ggplot(
    ordered_augmented_peaks,
    aes(x = ordered_sampleId, y = ordered_groupId, fill = abundance)
  ) +
    geom_tile() +
    scale_fill_gradient2(
      expression(log[2] ~ abundance),
      low = "steelblue1",
      mid = "black",
      high = "yellow",
      midpoint = 0
    ) +
    scale_x_discrete(sample.var,
      breaks = ordered_distinct_samples$ordered_sampleId,
      labels = ordered_distinct_samples$label
    ) +
    scale_y_discrete(
      breaks = ordered_distinct_features$ordered_groupId,
      labels = ordered_distinct_features$label,
      position = "right"
    ) +
    expand_limits(fill = c(-1 * change_threshold, change_threshold)) +
    heatmap_theme

  if (plot_type == "grob") {
    return(heatmap_plot)
  } else if (plot_type == "plotly") {
    suppressWarnings(
      plotly::ggplotly(heatmap_plot) %>%
        plotly::layout(margin = 0)
    )
  } else {
    stop("undefined plotting type logic")
  }
}

#' Hierarchical clustering order
#'
#' Format and hierarchically cluster a data.frame. If hclust could not
#'   normally be produced (usually because no samples are in common for a
#'   feature) pad the matrix with zeros and still calculate the distance
#'
#' @param df data.frame to cluster
#' @param feature_pk variable uniquely defining a row
#' @param sample_pk variable uniquely definining a sample
#' @param measurement_var measurement variables
#' @param cluster_dim rows, columns, or both
#'
#' @return a list containing a hierarchically clustered set of rows
#'   and/or columns
#'
#' @export
hclust_order <- function(df,
                         feature_pk,
                         sample_pk,
                         measurement_var,
                         cluster_dim) {
  stopifnot("data.frame" %in% class(df))
  stopifnot(length(feature_pk) == 1, feature_pk %in% colnames(df))
  stopifnot(length(sample_pk) == 1, sample_pk %in% colnames(df))
  stopifnot(length(measurement_var) == 1, measurement_var %in% colnames(df))
  stopifnot(
    length(cluster_dim) == 1,
    cluster_dim %in% c("rows", "columns", "both")
  )

  quant_matrix <- df %>%
    dplyr::select_(.dots = stats::setNames(
      list(feature_pk, sample_pk, measurement_var),
      c("feature_id", "sample_id", "quant")
    )) %>%
    reshape2::acast(feature_id ~ sample_id, value.var = "quant")

  output <- list()

  if (cluster_dim %in% c("rows", "both")) {
    cluster_rows <- try(
      stats::hclust(stats::dist(quant_matrix), method = "ward.D2"),
      silent = TRUE
    )

    # if distance cannot be computed (because of missing values) pad with
    # zeros and recalculate
    if (class(cluster_rows) == "try-error") {
      pad_matrix <- matrix(0, ncol = 2, nrow = nrow(quant_matrix))
      colnames(pad_matrix) <- c("pad1", "pad2")
      quant_matrix_pad <- cbind(quant_matrix, pad_matrix)

      cluster_rows <- stats::hclust(
        stats::dist(quant_matrix_pad),
        method = "ward.D2"
        )
    }

    output$rows <- rownames(quant_matrix)[cluster_rows$order]
  }

  if (cluster_dim %in% c("columns", "both")) {
    cluster_cols <- try(
      stats::hclust(stats::dist(t(quant_matrix)), method = "ward.D2"),
      silent = TRUE
    )

    # if distance cannot be computed (because of missing values) pad with
    # zeros and recalculate
    if (class(cluster_cols) == "try-error") {
      pad_matrix <- matrix(0, ncol = 2, nrow = ncol(quant_matrix))
      colnames(pad_matrix) <- c("pad1", "pad2")
      quant_matrix_pad <- cbind(t(quant_matrix), pad_matrix)

      cluster_cols <- stats::hclust(
        stats::dist(quant_matrix_pad),
        method = "ward.D2"
        )
    }

    output$columns <- colnames(quant_matrix)[cluster_cols$order]
  }

  output
}


#' Plot Barplots
#'
#' @inheritParams test_mzroll_list
#' @param groupIds groupIds for compounds to plot
#'
#' @export
plot_barplot <- function(mzroll_list, groupIds) {
  stopifnot(class(groupIds) %in% c("numeric", "integer"))

  reduced_groups <- mzroll_list$peakgroups %>%
    dplyr::filter(groupId %in% groupIds)

  reduced_peaks <- mzroll_list$peaks %>%
    dplyr::inner_join(reduced_groups, by = "groupId") %>%
    dplyr::left_join(mzroll_list$samples, by = "sampleId")

  standard_errors <- reduced_peaks %>%
    dplyr::group_by(groupId, peak_label, `sample description`) %>%
    dplyr::summarize(
      abund_mean = mean(log2_abundance),
      abund_se = stats::sd(log2_abundance) / sqrt(dplyr::n())
    ) %>%
    dplyr::ungroup()

  ggplot(standard_errors, aes(x = `sample description`)) +
    geom_col(aes(y = 2^abund_mean)) +
    geom_errorbar(aes(
      ymin = 2^(abund_mean - abund_se),
      ymax = 2^(abund_mean + abund_se)
    )) +
    facet_wrap(~peak_label, scales = "free_y") +
    scale_x_discrete("Sample description") +
    scale_y_continuous(expression(IC %+-% "1SE")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}
