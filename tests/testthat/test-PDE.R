test_that("PDE produces NA when NA in inputs", {
  mods1 <- c(29, 25, 21, 22)
  mods2 <- c(3, 23, 21, NA)
  expect_true(is.na(PDE(mods1, mods2)))
})

test_that("PDE produces a single value", {
  mods1 <- c(29, 25, 21, 22)
  mods2 <- c(3, 23, 21, 4)
  expect_length(PDE(mods1, mods2), 1)
})

test_that("PDE is 1 when mod and obs are the same", {
  mods1 <- c(29, 25, 21, 22)
  mods2 <- c(29, 25, 21, 22)
  expect_equal(PDE(mods1, mods2), 1)
})
