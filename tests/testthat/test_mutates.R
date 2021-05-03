context("Test Mutates")

test_that("Remove constant name", {
  name_set <- c("aaaxyzbbb", "aaaklbbb", "aaaxyzbbb")
  truncated_name_set <- c("xyz", "kl", "xyz")

  constant_names <- claman::remove_constant_name(name_set)

  expect_equal(constant_names, truncated_name_set)
})

test_that("Test batch centering", {
  
  centered_batches <- normalize_peaks(
    claman::nplug_mzroll_augmented,
    "center batches",
    quant_peak_varname = "log2_abundance",
    norm_peak_varname = "normalized_log2_abundance",
    batch_varnames = "month",
    centering_fxn = median
    ) %>%
    romic::triple_to_tidy()
  
  differences <- centered_batches$data %>%
    dplyr::group_by(groupId, month) %>%
    dplyr::summarize(
      group_median = median(normalized_log2_abundance, na.rm = TRUE),
      .groups = "drop"
      ) %>%
    dplyr::group_by(groupId) %>%
    dplyr::summarize(diff = diff(range(group_median)), .groups = "drop") %>%
    {.$diff}
  
  expect_equal(differences, rep(0, length(differences)))
})

test_that("Test reference samples and reference conditions", {
  
  reference_condition_normalized <- claman::normalize_peaks(
    claman::nplug_mzroll_augmented,
    normalization_method = "reference condition",
    quant_peak_varname = "log2_abundance",
    norm_peak_varname = "normalized_log2_abundance",
    condition_varname = "condition",
    reference_varname = "reference"
  )
  
  reference_sample_normalized <- claman::normalize_peaks(
    claman::nplug_mzroll_augmented,
    normalization_method = "reference sample",
    quant_peak_varname = "log2_abundance",
    norm_peak_varname = "normalized_log2_abundance",
    batch_varnames = c("month", "extraction"),
    reference_varname = "exp_ref",
    reference_values = "ref"
  )
  
  # reference condition and reference samples should just be off by
  # an intercept
  
  compare_normalizations <- reference_condition_normalized$measurements %>%
    dplyr::select(
      groupId,
      sampleId,
      condition_normalized = normalized_log2_abundance
      ) %>%
    dplyr::left_join(
      reference_sample_normalized$measurements %>%
        dplyr::select(
          groupId,
          sampleId,
          sample_normalized = normalized_log2_abundance,
        ),
      by = c("groupId", "sampleId")
    ) %>%
    dplyr::mutate(diff = sample_normalized - condition_normalized)
  
  expect_equal(
    compare_normalizations$diff,
    rep(0, nrow(compare_normalizations))
    )
  })

test_that("Flooring works and is maintained", {
  
  floored_peaks <- claman::floor_peaks(nplug_mzroll_augmented, 12)
  
  expect_equal(
    floored_peaks$measurements$log2_abundance >= 12,
    rep(TRUE, nrow(floored_peaks$measurements))
  )
  
  median_polished <- claman::normalize_peaks(
    claman::nplug_mzroll_augmented,
    "median polish",
    quant_peak_varname = "log2_abundance",
    norm_peak_varname = "normalized_log2_abundance",
    log2_floor_value = 12
  )

  expect_equal(
    median_polished$measurements$normalized_log2_abundance >= 12,
    rep(TRUE, nrow(median_polished$measurements))
  )

  reference_condition_normalized <- claman::normalize_peaks(
    claman::nplug_mzroll_augmented,
    normalization_method = "reference condition",
    quant_peak_varname = "log2_abundance",
    norm_peak_varname = "normalized_log2_abundance",
    condition_varname = "condition",
    reference_varname = "reference",
    log2_floor_value = 12
  )
  
  # is normalization actually doing something
  expect_false(isTRUE(all.equal(
    reference_condition_normalized$measurements$log2_abundance,
    reference_condition_normalized$measurements$normalized_log2_abundance
  )))
  
  # is flooring maintained
  expect_equal(
    reference_condition_normalized$measurements$normalized_log2_abundance >= 12,
    rep(TRUE, nrow(reference_condition_normalized$measurements))
  )
  
  reference_sample_normalized <- claman::normalize_peaks(
    claman::nplug_mzroll_augmented,
    normalization_method = "reference sample",
    quant_peak_varname = "log2_abundance",
    norm_peak_varname = "normalized_log2_abundance",
    batch_varnames = c("month", "extraction"),
    reference_varname = "exp_ref",
    reference_values = "ref",
    log2_floor_value = 12
  )
  
  # is normalization doing something
  expect_false(isTRUE(all.equal(
    reference_condition_normalized$measurements$log2_abundance,
    reference_condition_normalized$measurements$normalized_log2_abundance
  )))
  
  # is flooring maintained
  expect_equal(
    reference_sample_normalized$measurements$normalized_log2_abundance >= 12,
    rep(TRUE, nrow(reference_sample_normalized$measurements))
  )
  
  reference_sample_nofloor <- claman::normalize_peaks(
    claman::nplug_mzroll_augmented,
    normalization_method = "reference sample",
    quant_peak_varname = "log2_abundance",
    norm_peak_varname = "normalized_log2_abundance",
    batch_varnames = c("month", "extraction"),
    reference_varname = "exp_ref",
    reference_values = "ref"
  )
  # this should be identical reference_sample_normalized aside from
  # group-level intercept
  
  inconsistencies <- reference_sample_normalized$measurements %>%
    dplyr::select(
      groupId,
      sampleId,
      abund_unfloored = normalized_log2_abundance
      ) %>%
    dplyr::left_join(
      reference_sample_nofloor$measurements %>%
        dplyr::select(
          groupId,
          sampleId,
          abund_floored = normalized_log2_abundance
          ),
      by = c("groupId", "sampleId")
    ) %>%
    dplyr::mutate(diff = abund_unfloored - abund_floored) %>%
    dplyr::group_by(groupId) %>%
    dplyr::filter(all(abund_unfloored >= 12.5)) %>%
    dplyr::summarize(val = diff(range(diff)))
    
  expect_equal(
    inconsistencies$val,
    rep(0, nrow(inconsistencies))
  )
})