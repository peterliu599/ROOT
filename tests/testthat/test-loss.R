test_that("loss() computes SE-like objective and returns Inf on degenerate denom", {
  D <- data.frame(vsq = c(1, 4, 9), w = c(1, 1, 1))
  # No change
  expect_equal(loss(1, integer(0), D), sqrt(sum(D$vsq) / (sum(D$w)^2)))
  # Set one index to 0 -> smaller denominator effect
  out <- loss(0, 1L, D)
  expect_true(is.finite(out))
  # Force denominator zero
  D2 <- data.frame(vsq = c(1, 4), w = c(0, 0))
  expect_equal(loss(1, integer(0), D2), Inf)
})
