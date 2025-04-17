library(terra)
library(maxnet)

# testing data -------------------------------------------------------------
# Create some presence/absence data to fit a model
species.data <- readRDS(system.file("test_files", "goa_data_logarea_folds.rds", package = "EFHSDM"))
data(raster_stack) # not sure if this will work in context. This is GOA data only.
raster.stack <- terra::rast(raster_stack)
names(raster.stack) <- c("lon", "lat", "bdepth", "btemp", "slope", "sponge")

maxnet.covars <- c("bdepth", "btemp", "slope")
cofactors <- c("sponge")
r.mult <- 0.5
maxnet.model0 <- FitMaxnet(
  data = species.data,
  species = "dogfish",
  vars = maxnet.covars,
  facs = cofactors,
  regmult = r.mult,
  reduce = TRUE
)

# tests -------------------------------------------------------------------
test_that("cloglog prediction works and returns a SpatRaster", {
  result <- MakeMaxEntAbundance(model = maxnet.model0, maxent.stack = raster.stack,
                                scale.fac = 1.619, type = "cloglog")
  expect_s4_class(result, "SpatRaster")
  expect_equal(ncell(result), ncell(raster.stack[[1]]))
})

test_that("maxnet prediction type works", {
  result <- MakeMaxEntAbundance(model = maxnet.model0, maxent.stack = raster.stack,
                                scale.fac = 1.619, type = "maxnet")
  expect_s4_class(result, "SpatRaster")
})

test_that("does not error with NULL or NA filename", {
  expect_no_error(MakeMaxEntAbundance(maxnet.model0, raster.stack, filename = NULL))
  expect_no_error(MakeMaxEntAbundance(maxnet.model0, raster.stack, filename = NA))
})
