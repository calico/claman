test_that("Filtering by groupId works", {
  filtered_mzroll <- filter_groupIds(
    claman::nplug_mzroll_normalized,
    groupIds = factor(1:10),
    invert = FALSE
  )

  expect_equal(as.numeric(filtered_mzroll$features$groupId), 1:10)

  filtered_mzroll <- filter_groupIds(
    claman::nplug_mzroll_normalized,
    groupIds = factor(1:10),
    invert = TRUE
  )

  expect_equal(sum(c(1:11) %in% filtered_mzroll$features$groupId), 1)
})
