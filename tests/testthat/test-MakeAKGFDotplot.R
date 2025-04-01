test_that("MakeAKGFDotPlot() creates a ggplot object", {
  # Do I need this system.file() part?
  species.data <- readRDS(system.file("test_files", "goa_data_logarea_folds.rds", package = "EFHSDM"))
  hd <- quantile(species.data[species.data[, "dogfish"] > 0, "dogfish"], probs = .9)
  testfig <- MakeAKGFDotplot(
    presence = species.data[species.data[, "dogfish"] > 0, ],
    absence = species.data[species.data[, "dogfish"] == 0, ],
    highdensity = species.data[species.data[, "dogfish"] > hd, ],
    region = tolower("GOA"),
    dataCRS = raster::crs(raster.stack, asText = TRUE),
    title.name = figure.name
  )
  expect_true(is(testfig, "ggplot"))
})
