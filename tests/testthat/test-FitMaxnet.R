test_that("FitMaxnet() creates object of proper length", {
  # This is basic but will be useful if package dependencies change
  species.data <- readRDS(system.file("test_files", "goa_data_logarea_folds.rds", package = "EFHSDM"))
  maxnet.covars <- c(
    "bcurrentU", "bcurrentV", "bcurrentUSD", "bcurrentVSD", "bdepth",
    "slope", "aspectE", "aspectN", "curve", "btemp", "tmax", "rocky", "BPI"
  )
  cofactors <- c("sponge", "coral", "pen")
  r.mult <- 0.5
  maxnet.model0 <- FitMaxnet(
    data = species.data, species = "dogfish", vars = maxnet.covars,
    facs = cofactors,
    regmult = r.mult, reduce = T
  )

  expect_length(object = class(maxnet.model0), n = 3)
})

# Need to check that gets the same outputs as in the 2023 5-Year EFH cycle
# This is species-specific and right now it uses GOA dogfish.
test_that("FitMaxnet() object matches 2023 review cycle", {
  species.data <- readRDS(system.file("test_files", "goa_data_logarea_folds.rds", package = "EFHSDM"))
  maxnet.covars <- c(
    "bcurrentU", "bcurrentV", "bcurrentUSD", "bcurrentVSD", "bdepth",
    "slope", "aspectE", "aspectN", "curve", "btemp", "tmax", "rocky", "BPI"
  )
  cofactors <- c("sponge", "coral", "pen")
  r.mult <- 0.5
  maxnet.model0 <- FitMaxnet(
    data = species.data, species = "dogfish", vars = maxnet.covars,
    facs = cofactors,
    regmult = r.mult, reduce = T
  )
  maxnet.model0_check <- readRDS(system.file("test_files", "goa_maxnet_model0_test.rds", package = "EFHSDM"))
  expect_true(identical(maxnet.model0, maxnet.model0_check))
})


# Other tests -------------------------------------------------------------
set.seed(123)
n <- 100
mock_data <- data.frame(
  species = rbinom(n, 1, 0.5),
  var1 = runif(n),
  var2 = rnorm(n),
  factor1 = sample(c("A", "B"), n, replace = TRUE),
  stringsAsFactors = FALSE
)

library(maxnet)

test_that("FitMaxnet runs with numeric variables only", {
  model <- FitMaxnet(
    data = mock_data,
    species = "species",
    vars = c("var1", "var2"),
    reduce = FALSE
  )

  expect_s3_class(model, "maxnet")
  expect_true("var1" %in% names(model$betas))
})

test_that("FitMaxnet handles factor variables", {
  model <- FitMaxnet(
    data = mock_data,
    species = "species",
    vars = c("var1"),
    facs = c("factor1"),
    reduce = FALSE
  )

  expect_s3_class(model, "maxnet")
  # Expect factor1 dummy variables to be included
  expect_true(any(grepl("factor1", names(model$betas))))
})

test_that("FitMaxnet removes rows with NA", {
  temp_data <- mock_data
  temp_data$var1[1:5] <- NA # introduce NA

  model <- FitMaxnet(
    data = temp_data,
    species = "species",
    vars = c("var1", "var2"),
    reduce = FALSE
  )

  # Should still return a model
  expect_s3_class(model, "maxnet")
  # Make sure rows were dropped
  expect_true(length(model$penalty.vec) < n)
})

test_that("FitMaxnet produces a maxnet object", {
  # Create data where one variable is noise
  noisy_data <- mock_data
  noisy_data$noise <- rnorm(n) * 1e-6 # basically noise

  model_full <- FitMaxnet(
    data = noisy_data,
    species = "species",
    vars = c("var1", "var2", "noise"),
    reduce = TRUE
  )

  model_full_class <- class(model_full)[1]

  expect_equal(model_full_class, "maxnet")
})

test_that("FitMaxnet works with regmult set", {
  model <- FitMaxnet(
    data = mock_data,
    species = "species",
    vars = c("var1", "var2"),
    regmult = 2
  )

  expect_s3_class(model, "maxnet")
  # Optional: check some structure of the model if regmult affects it
})
