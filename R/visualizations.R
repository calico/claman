#' Plot Barplots
#'
#' @inheritParams test_mzroll_list
#' @param groupIds groupIds for compounds to plot
#' @param grouping_vars variables to use for ordering samples and calculating
#'   standard errors
#' @inheritParams diffex_mzroll
#' @param fill_var variable to use for filling bars
#'
#' @examples
#'
#' plot_barplot(
#'   mzroll_list = nplug_mzroll_normalized,
#'   groupIds = 1:4,
#'   grouping_vars = c("limitation", "DR"),
#'   value_var = "normalized_log2_abundance",
#'   fill_var = "limitation"
#' )
#' @export
plot_barplot <- function(mzroll_list,
                         groupIds,
                         grouping_vars,
                         value_var,
                         fill_var = NULL) {
  checkmate::assertClass(mzroll_list, "tomic")
  checkmate::assertClass(mzroll_list, "mzroll")
  stopifnot(class(groupIds) %in% c("numeric", "integer", "factor", "ordered"))
  checkmate::assertCharacter(grouping_vars)
  purrr::walk(
    grouping_vars,
    checkmate::assertChoice,
    mzroll_list$design$samples$variable
  )
  checkmate::assertChoice(value_var, mzroll_list$design$measurements$variable)
  checkmate::assertClass(mzroll_list$measurements[[value_var]], "numeric")
  if (!is.null(fill_var)) {
    checkmate::assertChoice(fill_var, mzroll_list$design$samples$variable)
  }

  reduced_groups <- mzroll_list$features %>%
    dplyr::filter(groupId %in% groupIds)

  ordered_samples <- mzroll_list$samples %>%
    dplyr::arrange(!!!rlang::syms(grouping_vars)) %>%
    dplyr::mutate(
      sample_label = paste(!!!rlang::syms(grouping_vars)),
      sample_label = factor(
        sample_label,
        levels = unique(sample_label)
      )
    )

  reduced_peaks <- mzroll_list$measurements %>%
    dplyr::inner_join(reduced_groups, by = "groupId") %>%
    dplyr::left_join(ordered_samples, by = "sampleId")

  standard_errors <- reduced_peaks %>%
    dplyr::group_by(groupId, peak_label, sample_label) %>%
    dplyr::summarize(
      abund_mean = mean(log2_abundance),
      abund_se = stats::sd(log2_abundance) / sqrt(dplyr::n()),
      .groups = "drop"
    )

  if (!is.null(fill_var)) {
    standard_errors <- standard_errors %>%
      dplyr::left_join(
        ordered_samples %>%
          dplyr::select(sample_label, !!rlang::sym(fill_var)),
        by = "sample_label"
      )
  }

  if (!is.null(fill_var)) {
    grob <- ggplot(standard_errors, aes(
      x = sample_label,
      fill = !!rlang::sym(fill_var)
    ))
  } else {
    grob <- ggplot(standard_errors, aes(x = sample_label))
  }

  grob +
    geom_col(aes(y = 2^abund_mean)) +
    geom_errorbar(
      data = standard_errors %>%
        dplyr::filter(!is.na(abund_se)),
      aes(ymin = 2^ (abund_mean - abund_se), ymax = 2^ (abund_mean + abund_se))
    ) +
    facet_wrap(~peak_label, scales = "free_y") +
    scale_x_discrete("Sample description") +
    scale_y_continuous(expression(IC %+-% "1SE")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}
