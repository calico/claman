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
