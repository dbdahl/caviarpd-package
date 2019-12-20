context("dummy-test")

# skip("dummy-test")

test_that("Dummy test.", {
  requireLevel(1)
  testthat::expect_gt(0.5, 0.01)
})

