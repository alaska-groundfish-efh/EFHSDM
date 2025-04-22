test_that("MakeAKGFDotPlot() creates a ggplot object", {
  # use the system.file() part if you are using data that are automatically included in the package
  species.data <- readRDS(system.file("test_files", "goa_data_logarea_folds.rds", package = "EFHSDM"))
  data("raster_stack", package = "EFHSDM")
  hd <- quantile(species.data[species.data[, "dogfish"] > 0, "dogfish"], probs = .9)

  testfig <- suppressMessages( # suppress akgfmaps message about EPSG:3338
    MakeAKGFDotplot(
    presence = species.data[species.data[, "dogfish"] > 0, ],
    absence = species.data[species.data[, "dogfish"] == 0, ],
    highdensity = species.data[species.data[, "dogfish"] > hd, ],
    region = tolower("GOA"),
    dataCRS = raster::crs(raster_stack, asText = TRUE),
    title.name = "This figure is a test"
  )
  )
  expect_true(is(testfig, "ggplot"))
})
