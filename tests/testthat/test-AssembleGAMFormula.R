
test_that("AssembleGAMFormula creates an object of class formula", {
  test.gam.table <- data.frame(
    type = c("smooth", "smooth"),
    dims = c(2, 2),
    term = c("lon", "bcurrentU"),
    term2 = c("lat", "bcurrentV"),
    bs = c("ds", "ds"),
    k = c(10, 10),
    m = c(1, 1),
    m2 = c(0.5, 0.5)
  )
  x <- AssembleGAMFormula(gam.table = test.gam.table,yvar = "dogfish", hgam = FALSE)
  expect_s3_class(object = x, class = "formula")
})

test_that("basic smooth term works", {
  gam.table <- data.frame(
    term = "x1",
    dims = c(2),
    type = "smooth",
    term2 = NA,
    bs = NA,
    k = NA,
    m = NA,
    m2 = NA
  )

  formula_out <- AssembleGAMFormula("y", gam.table, hgam = FALSE)
  # have to use 'deparse' to avoid comparing environments
  expect_equal(deparse(formula_out), deparse(as.formula("y ~ s(x1)")))
})

test_that("handles smooth and factor terms", {
  gam.table <- data.frame(
    term = c("x1", "x2"),
    type = c("smooth", "factor"),
    term2 = c(NA, NA),
    bs = c(NA, NA),
    k = c(NA, NA),
    m = c(NA, NA),
    m2 = c(NA, NA)
  )

  formula_out <- AssembleGAMFormula(yvar = "y", gam.table, hgam = FALSE)
  expect_equal(deparse(formula_out), deparse(as.formula("y ~ s(x1) + as.factor(x2)")))
})


test_that("handles m and m2 together", {
  gam.table <- data.frame(
    term = "x3",
    type = "smooth",
    term2 = NA,
    bs = "tp",
    k = 4,
    m = 1,
    m2 = 2
  )

  formula_out <- AssembleGAMFormula("y", gam.table)
  expect_equal(deparse(formula_out), deparse(as.formula("y ~ s(x3,bs='tp',k=4,m=c(1,2))")))
})
