context("basic")


test_that("u_count"), {
  vector <- c("c", "d", "v", "d", "t", "c", "a")
  expect_identical(u_count(vector), 5)
})


test_that("mean_na"), {
  vector <- c(2, 2, 4, NA, 2, 6, 10, 12, 14)
  expect_identical(mean_na(vector), 6.5)
})


test_that("sd_na"), {
  vector <- c(2, 2, 4, NA, 2, 6, 10, 12, 14)
  expect_identical(signif(sd_na(vector), 6), 4.86973)
})


test_that("cv_na"), {
  vector <- c(2, 2, 4, NA, 2, 6, 10, 12, 14)
  expect_identical(signif(cv_na(vector), 6), 0.749189)
})


test_that("median_na"), {
  vector <- c(2, 2, 4, NA, 2, 6, 10, 12, 14)
  expect_identical(signif(median_na(vector), 6), 5)
})


test_that("sum_na"), {
  vector <- c(2, 2, 4, NA, 2, 6, 10, 12, 14)
  expect_identical(signif(sum_na(vector), 6), 52)
})


