#' Plot EFH comparison maps
#' @description This function makes a canned version of the EFH comparison plots. Could use some more work, and is a little less flexible than other functions at current.
#' @param old raster of EFH values
#' @param new new raster of EFH values
#' @param main character; title for the figure
#' @param background raster used to define the survey area, often a covariate raster
#' @param leg.name character; name for the legend
#' @param leg.labels character vector; names to appear in the legend
#' @param region character; which region is being used, must correpond to an akgfmaps baselayer
#' @param col.vec vector of colors to be used for the categories (usually length 4)
#' @param label.pos coordinates for the title box
#' @param leg.pos decimal values for placing the legend
#' @param ext.adjust.x vector length 2 to adjust the horizontal extent
#' @param ext.adjust.y vector length 2 to adjust the vertical extent
#' @param nonEFH integer; the value corresponding to non-EFH in the old and new rasters
#' @importFrom akgfmaps get_base_layers
#' @importFrom magrittr %>%
#' @importFrom terra classify
#' @importFrom terra values
#'
#' @return ggplot object with the comparison map
#' @export
#'
#' @examples
PlotEFHComparison <- function(old = NA, new = NA, main = "", background, leg.name = NULL, leg.labels, region, col.vec = NA,
                              label.pos = NA, leg.pos = NA, ext.adjust.x = NA, ext.adjust.y = NA, nonEFH = 1) {

  # start by getting the map
  region <- tolower(region)
  if (region %in% c(
    "ebs", "bs.all", "sebs", "bs.south", "ecs", "ebs.ecs", "ai",
    "ai.west", "ai.central", "ai.east", "goa", "goa.west", "goa.east"
  )) {
    MAP <- akgfmaps::get_base_layers(select.region = region, set.crs = "auto")
  } else {
    stop("region not recognized")
  }

  # set up some plot parameters
  if (region %in% c("bs.all", "ebs")) {
    if (is.na(label.pos)) {
      label.pos <- c(-550000, 810000)
    }
    if (is.na(leg.pos)) {
      leg.pos <- c(0.03, 0.21)
    }
    if (is.na(ext.adjust.x)) {
      ext.adjust.x <- c(0, 0)
    }
    if (is.na(ext.adjust.y)) {
      ext.adjust.y <- c(0, 0)
    }
  }
  if (region == "goa") {
    if (is.na(label.pos)) {
      label.pos <- c(-1350000, 660000)
    }
    if (is.na(leg.pos)) {
      leg.pos <- c(0.01, 0.7)
    }
    if (is.na(ext.adjust.x)) {
      ext.adjust.x <- c(200000, -100000)
    }
    if (is.na(ext.adjust.y)) {
      ext.adjust.y <- c(0, 0)
    }
  }
  if (region == "ai") {
    if (is.na(label.pos)) {
      label.pos <- c(-900000, 490000)
    }
    if (is.na(leg.pos)) {
      leg.pos <- c(0.01, .38)
    }
    if (is.na(ext.adjust.x)) {
      ext.adjust.x <- c(0, 0)
    }
    if (is.na(ext.adjust.y)) {
      ext.adjust.y <- c(0, 0)
    }
  }

  if (is.na(col.vec)) {
    col.vec <- c("dodgerblue", "orange", "orchid3")
  }

  dummy.sf0 <- stars::st_as_stars(is.na(background) == F)
  dummy.sf<- sf::st_cast(sf::st_as_sf(dummy.sf0,merge = TRUE),"POLYGON") # cast the polygons to polylines
  dummy.sf2 <- sf::st_transform(dummy.sf, sf::st_crs(MAP$akland))
  dummy.sf3 <- dummy.sf2[dummy.sf2$layer == 1, ]


  old.col <- NULL
  new.col <- NULL
  mix.col <- NULL

  # Check which rasters are available
  try(old.present <- is.na(old@crs) == F)
  try(new.present <- is.na(new@crs) == F)

  if (exists("old.present") & exists("new.present")) {
    comp.raster <- terra::classify(background, breaks = c(-Inf, Inf))

    # find which areas are EFH in each version
    old2 <- terra::classify(x = old, breaks = c(0, nonEFH + .5, Inf))
    new2 <- terra::classify(x = new, breaks = c(0, nonEFH + .5, Inf))

    vals <- terra::values(comp.raster)

    both <- which(terra::values(old2) == 2 & terra::values(new2) == 2)
    justold <- which(terra::values(old2) == 2 &
                       ((terra::values(new2) == 1) | is.na(terra::values(new2))))
    justnew <- which(((terra::values(old2) == 1) | is.na(terra::values(old2))) &
                       terra::values(new2) == 2)

    # assign new values to each category
    vals[justold] <- 2
    vals[justnew] <- 3
    vals[both] <- 4

    comp.raster <- terra::setValues(x = comp.raster, values = vals)

    efhpoly0 <- stars::st_as_stars(comp.raster)
    efhpoly<-  sf::st_as_sf(efhpoly0,merge = TRUE)
    efhpoly2 <- efhpoly[efhpoly$layer > 1, ]

    # this is a kludge to make sure there are always all 3 types
    layer.present <- unique(efhpoly2$layer)
    n <- 0
    if (2 %in% layer.present == F) {
      efhpoly2$layer[which.min(sf::st_area(efhpoly2))] <- 2
      n <- 1
    }
    if (3 %in% layer.present == F) {
      efhpoly2$layer[order(sf::st_area(efhpoly2))[1 + n]] <- 3
      n <- n + 1
    }
    if (4 %in% layer.present == F) {
      efhpoly2$layer[order(sf::st_area(efhpoly2))[1 + n]] <- 4
    }

    area1 <- sum(terra::values(old) > nonEFH, na.rm = T)
    area2 <- sum(terra::values(new) > nonEFH, na.rm = T)
    main <- paste0(main, "\nChange = ", round((area2 / area1 - 1) * 100, 1), "%")

    old.col <- col.vec[1]
    new.col <- col.vec[2]
    mix.col <- col.vec[3]
  } else {
    if (exists("old.present")) {
      efhpoly0 <- stars::st_as_stars(old > nonEFH)
      efhpoly <- sf::st_as_sf(efhpoly0,merge = TRUE)
      efhpoly2 <- efhpoly[efhpoly$layer != 0, ]
      old.col <- col.vec[1]
    }
    if (exists("new.present")) {
      efhpoly0 <- stars::st_as_stars(new > nonEFH)
      efhpoly <- sf::st_as_sf(efhpoly0,merge = TRUE)
      efhpoly2 <- efhpoly[efhpoly$layer != 0, ]
      new.col <- col.vec[2]
    }
  }

  out.plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = dummy.sf3, fill = "grey95") +
    ggplot2::geom_sf(data = MAP$survey.area, fill = "grey95") +
    ggplot2::geom_sf(data = efhpoly2, ggplot2::aes(fill = factor(layer)), col = 1, linewidth = .05) +
    ggplot2::geom_sf(data = MAP$akland, fill = "grey40") +
    ggplot2::geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) +
    ggplot2::geom_sf(data = MAP$bathymetry, col = "grey60", size = .25) +
    ggplot2::coord_sf(xlim = MAP$plot.boundary$x + ext.adjust.x, ylim = MAP$plot.boundary$y + ext.adjust.y) +
    ggplot2::geom_label(
      data = data.frame(x = label.pos[1], y = label.pos[2], label = main),
      ggplot2::aes(x = x, y = y, label = main, hjust = 0, vjust = 1), size = 5
    ) +
    ggplot2::scale_fill_manual(values = c(old.col, new.col, mix.col), labels = leg.labels, name = leg.name) +
    ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) +
    ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      legend.position.inside = leg.pos, legend.justification = c(0, 1),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = NA, color = NA)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, shape = 15)))

  return(out.plot)
}
