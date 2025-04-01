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
