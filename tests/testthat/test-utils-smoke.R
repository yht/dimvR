test_that("mse computes numeric mean squared error", {
  expect_equal(dimvR:::mse(c(1, 2, 3), c(1, 2, 4)), 1 / 3)
  expect_equal(dimvR:::mse(c(1, NA, 3), c(1, 2, 5)), 2)
})

test_that("cache helpers store and retrieve objects consistently", {
  X <- data.frame(a = 1:3, b = c(2, 4, 6))
  obj <- list(ok = TRUE, value = 42)

  key <- dimvR:::cache_key("demo", X)
  expect_true(is.character(key))
  expect_match(key, "^demo_")

  # cache miss
  expect_null(dimvR:::get_cached("demo", X))

  # cache hit
  dimvR:::set_cached("demo", X, obj)
  got <- dimvR:::get_cached("demo", X)
  expect_equal(got, obj)
})
