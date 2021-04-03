context("Test Mutates")

test_that("Remove constant name", {
  name_set <- c("aaaxyzbbb", "aaaklbbb", "aaaxyzbbb")
  truncated_name_set <- c("xyz", "kl", "xyz")

  constant_names <- calicomics::remove_constant_name(name_set)

  expect_equal(constant_names, truncated_name_set)
})
