# this script will make some more generalized functions for mapping things with the akgfmaps package

# rpackages <- c("viridis", "stars", "sf", "gridExtra", "patchwork", "MASS", "scales", "labeling", "ggplot2")
#
# which_not_installed <- which(rpackages %in% rownames(installed.packages()) == FALSE)
#
# if (length(which_not_installed) > 1) {
#   install.packages(rpackages[which_not_installed], dep = TRUE)
# }
#
# if ("akgfmaps" %in% rownames(installed.packages()) == F) {
#   devtools::install_github("sean-rohan-noaa/akgfmaps", build_vignettes = TRUE)
# }
# rm(rpackages, which_not_installed)

# require(ggplot2)
# require(akgfmaps)
# require(viridis)
# require(stars)
# require(sf)
# require(gridExtra)
# require(patchwork)
# require(MASS)
# require(scales)
# require(labeling)
#
#
#' Make dotplot with akfgmaps
#'
#' @description makes a pretty good dotplot; needs additional testing with areas other than "bs.all", "goa", and "ai".
#' @details Settings are set to work with 8x8 in 300 res output, so results may vary when making different size figures.
#' @param presence data frame with coordinates (lon/lat) for presence locations
#' @param absence data frame with coordinates (lon/lat) for absence locations
#' @param highdensity data frame with coordinates (lon/lat) for high density locations
#' @param dataCRS CRS for the data, if different from the maps
#' @param region character; which region is being used, must correpond to an akgfmaps baselayer
#' @param survey.area sf or raster layer to use, if the desired survey area is different from akgfmaps
#' @param ext.adjust vector of length four, representing xmin,xmax, ymin,ymax, to adjust the extent
#' @param legend.pos vector of length two with a decimal location for the legend, use NA to suppress
#' @param title.name character; a title for the figure, or NA to suppress
#' @param title.count logical; should a count of the total number of presence locations be added to the title
#' @param title.pos vector of length two with coordinates for the title box
#' @param abs.name character; a name for the absences in the legend
#' @param pres.name character; a name for the presences in the legend
#' @param hd.name characer; a name for the high density presences in the legend
#' @param abs.col color for absence dots
#' @param pres.col color for presence dots
#' @param hd.col color for high density dots
#' @param abs.shape shape for presence dots
#' @param pres.shape shape for absence dots
#' @param hd.shape shape for high density dots
#' @param abs.size size for absence dots
#' @param pres.size size for presence dots
#' @param hd.size size for high density dots
#' @importFrom akgfmaps get_base_layers
#' @importFrom magrittr %>%
#'
#' @return a ggplot with the desired figure
#' @export
#'
#' @examples
MakeAKGFDotplot <- function(presence,
                            absence = NA,
                            highdensity = NA,
                            dataCRS = NA, #
                            region,
                            survey.area = "default",
                            ext.adjust = "default", #
                            legend.pos = "default", #
                            title.name = NA,
                            title.count = T,
                            title.pos = "default",
                            abs.name = "absent",
                            pres.name = "present",
                            hd.name = "top 10%",
                            abs.col = "steelblue4",
                            pres.col = "orange",
                            hd.col = "red",
                            abs.shape = 16,
                            pres.shape = 1,
                            hd.shape = 16,
                            abs.size = .1,
                            pres.size = 1,
                            hd.size = 1) {
  # start by getting the map
  region <- tolower(region)
  #browser()
  if (region %in% c(
    "ebs", "bs.all", "sebs", "bs.south", "ecs", "ebs.ecs", "ai",
    "ai.west", "ai.central", "ai.east", "goa", "goa.west", "goa.east"
  )) {
    MAP <- akgfmaps::get_base_layers(select.region = region, set.crs = "auto",use.survey.bathymetry = TRUE)
  } else {
    stop("region not recognized")
  }

  # set up some essentials for later
  if (ext.adjust == "default") {
    if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("goa", "goa.west", "goa.east")) {
      ext.adjust <- c(200000, -100000, 0, 0)
    }
  }
  if (length(legend.pos)<2){
    if(legend.pos == "default") {
      if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
        legend.pos <- c(0.12, 0.20)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
        legend.pos <- c(0.12, 0.18)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        legend.pos <- c(0.08, 0.41)
      }
    }else{
      legend.pos<-NA
    }
  }
  if (length(title.pos)<2){
    if(title.pos == "default") {
      if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
        title.pos <- c(-900000, 490000)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
        title.pos <- c(-550000, 800000)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        title.pos <- c(-1300000, 630000)
      }
    }else{
      title.pos<-NA
    }
  }

  # detect and set up the outline, for now it is required
  if (is.character(survey.area) && survey.area[1] == "default") {
    survey.sf <- MAP$survey.area
  } else {
    if (class(survey.area)[1] == "sf") {
      survey.sf <- sf::st_transform(survey.area, sf::st_crs(MAP$akland))
    }
    if (class(survey.area)[1] == "RasterLayer") {
      survey.sf1<-stars::st_as_stars(is.na(survey.area))
      survey.sf2<-sf::st_as_sf(survey.sf1,merge = TRUE)
      survey.sf3<-sf::st_cast(survey.sf2,"POLYGON") # cast the polygons to polylines

      survey.sf <- sf::st_transform(survey.sf3, sf::st_crs(MAP$akland))[1:(nrow(survey.sf3) - 1), ]
    }
  }

  leg.col <- NULL
  leg.shape <- NULL
  leg.size <- NULL
  leg.name <- NULL

  # Now go through and set up the dot locations and add them to legend
  if (is.data.frame(absence)) {
    if (is.na(dataCRS)) {
      abs.sf <- sf::st_as_sf(x = absence, coords = c("lon", "lat"), crs = sf::st_crs(MAP$akland))
    } else {
      abs.sf0 <- sf::st_as_sf(x = absence, coords = c("lon", "lat"), crs = dataCRS)
      abs.sf <- sf::st_transform(abs.sf0, sf::st_crs(MAP$akland))
    }
    leg.name <- abs.name
    leg.col <- abs.col
    leg.shape <- abs.shape
    leg.size <- abs.size
    abs.fac <- 1
  } else {
    abs.fac <- 0
  }

  if (is.na(dataCRS)) {
    pres.sf <- sf::st_as_sf(x = presence, coords = c("lon", "lat"), crs = sf::st_crs(MAP$akland))
  } else {
    pres.sf0 <- sf::st_as_sf(x = presence, coords = c("lon", "lat"), crs = dataCRS)
    pres.sf <- sf::st_transform(pres.sf0, sf::st_crs(MAP$akland))
  }
  leg.name <- c(leg.name, pres.name)
  leg.col <- c(leg.col, pres.col)
  leg.shape <- c(leg.shape, pres.shape)
  leg.size <- c(leg.size, pres.size)
  pres.fac <- abs.fac + 1

  if (is.data.frame(highdensity)) {
    if (is.na(dataCRS)) {
      high.sf <- sf::st_as_sf(x = highdensity, coords = c("lon", "lat"), crs = sf::st_crs(MAP$akland))
    } else {
      high.sf0 <- sf::st_as_sf(x = highdensity, coords = c("lon", "lat"), crs = dataCRS)
      high.sf <- sf::st_transform(high.sf0, sf::st_crs(MAP$akland))
    }
    leg.name <- c(leg.name, hd.name)
    leg.col <- c(leg.col, hd.col)
    leg.shape <- c(leg.shape, hd.shape)
    leg.size <- c(leg.size, hd.size)
    hd.fac <- pres.fac + 1
  }

  # set up the basic map, will add more customization later
  dotplot <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = survey.sf, fill = "grey95") +
    ggplot2::geom_sf(data = MAP$akland, fill = "grey40") +
    ggplot2::geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) +
    ggplot2::geom_sf(data = MAP$bathymetry, color = "grey60")

  # add the dots
  if (is.data.frame(absence)) {
    dotplot <- dotplot +
      ggplot2::geom_sf(data = abs.sf, alpha = .25, size = abs.size, shape = abs.shape, ggplot2::aes(color = factor(abs.fac)))
  }
  dotplot <- dotplot +
    ggplot2::geom_sf(data = pres.sf, size = pres.size, ggplot2::aes(color = factor(pres.fac)), shape = pres.shape, stroke = .8)
  if (is.data.frame(highdensity)) {
    dotplot <- dotplot +
      ggplot2::geom_sf(data = high.sf, size = hd.size, shape = hd.shape, ggplot2::aes(color = factor(hd.fac)))
  }

  # add the themes
  dotplot <- dotplot +
    ggplot2::coord_sf(xlim = MAP$plot.boundary$x + ext.adjust[1:2], ylim = MAP$plot.boundary$y + ext.adjust[3:4]) +
    ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) +
    ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      legend.position = legend.pos,
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = NA, color = NA)
    )

  # add a title
  if (is.na(title.name) == F && length(title.pos) == 2) {
    if (title.count) {
      title.name <- paste0(title.name, "\nN = ", format(nrow(pres.sf), big.mark = ","))
    }
    dotplot <- dotplot +
      ggplot2::geom_label(
        data = data.frame(x = title.pos[1], y = title.pos[2], label = title.name),
        ggplot2::aes(x = x, y = y, label = label, hjust = 0, vjust = 1), size = 5
      )
  }

  # add a legend
  if (length(legend.pos)==2) {
    dotplot <- dotplot +
      ggplot2::scale_color_manual(name = NULL, values = leg.col, labels = leg.name) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = leg.shape, size = leg.size)))
  }
  return(dotplot)
}


#' Map a continuous variable
#'
#' @description Makes a pretty good map of any continuous variable
#' @param region character; which region is being used, must correpond to an akgfmaps baselayer
#' @param density.map raster or sf; a map to be plotted
#' @param buffer numeric; a value between 0 and 1 that can be used to remove outlying value that mess up the color scale
#' @param survey.area sf or raster layer to use, if the desired survey area is different from akgfmaps
#' @param ext.adjust vector of length four, representing xmin,xmax, ymin,ymax, to adjust the extent
#' @param legend.pos vector of length two with decimal position of legend box
#' @param legend.title character; a title for the legend
#' @param col.palette character a color palette for
#' @param col.palette.limits vector of length two giving start and end points for the color palette
#' @param title.name character; a title for the figure, or NA to suppress
#' @param title.pos vector of length two with coordinates for for the legend, use NA to suppress
#' @param barheight numeric; height of the bar in the legend
#' @importFrom akgfmaps get_base_layers
#' @importFrom magrittr %>%
#'
#' @return a ggplot object of the map
#' @export
#'
#' @examples
MakeAKGFDensityplot <- function(region,
                                density.map,
                                buffer = 1,
                                survey.area = "default",
                                ext.adjust = "default", # vector of length four, representing xmin,xmax, ymin,ymax, to adjust the extent
                                legend.pos = "default", # a vector of length two with a location for the legend, use NA to suppress
                                legend.title = NA,
                                col.palette = "plasma",
                                col.palette.limits = c(0, 1),
                                title.name = NA,
                                title.pos = "default",
                                barheight = 5) {

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

  # set up some essentials for later
  if (ext.adjust == "default") {
    if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("goa", "goa.west", "goa.east")) {
      ext.adjust <- c(200000, -100000, 0, 0)
    }
  }
  if (length(legend.pos)<2){
    if(legend.pos == "default") {
      if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
        legend.pos <- c(0.08, 0.28)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
        legend.pos <- c(0.12, 0.18)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        legend.pos <- c(0.08, 0.51)
      }
    }else{
      legend.pos<-NA
    }
  }
  if (length(title.pos)<2){
    if(title.pos == "default") {
      if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
        title.pos <- c(-900000, 490000)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
        title.pos <- c(-550000, 800000)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        title.pos <- c(-1300000, 630000)
      }
    }else{
      title.pos<-NA
    }
  }

  # detect and set up the outline, for now it is required
  if (is.character(survey.area) && survey.area[1] == "default") {
    survey.sf <- MAP$survey.area
  } else {
    if (class(survey.area)[1] == "sf") {
      survey.sf <- sf::st_transform(survey.area, sf::st_crs(MAP$akland))[1:(nrow(survey.area) - 1), ]
    }
    if (class(survey.area)[1] == "RasterLayer") {
      survey.sf1<-stars::st_as_stars(is.na(survey.area))
      survey.sf2<-sf::st_as_sf(survey.sf1,merge = TRUE)
      survey.sf3<-sf::st_cast(survey.sf2,"POLYGON") # cast the polygons to polylines

      survey.sf <- sf::st_transform(survey.sf3, sf::st_crs(MAP$akland))[1:(nrow(survey.sf3) - 1), ]
    }
  }

  # set up the density component
  if (class(density.map)[1] == "sf") {
    density.sf <- sf::st_transform(density.map, sf::st_crs(MAP$akland))
    names(density.sf[1]) <- "density"
  }
  if (class(density.map)[1] == "RasterLayer") {
    density.dat0 <- as.data.frame(raster::xyFromCell(density.map, 1:raster::ncell(density.map)))
    vals <- raster::getValues(density.map)

    # convert the raster to ggplot format
    density.dat <- data.frame(lat = density.dat0$y, lon = density.dat0$x, density = vals)

    density.sf0 <- sf::st_as_sf(x = subset(density.dat, is.na(density) == F), coords = c("lon", "lat"), crs = density.map@crs)
    density.sf <- sf::st_transform(density.sf0, sf::st_crs(MAP$akland))
  }

  # often helps to remove some of the high points
  upper<-stats::quantile(density.sf$density, buffer, na.rm = T)
  density.sf$density[density.sf$density > upper] <- upper

  # set up the basic map, will add more customization later
  densityplot <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = density.sf, ggplot2::aes(col = density), size = .05)+
    ggplot2::geom_sf(data = survey.sf, fill = NA) +
    ggplot2::geom_sf(data = MAP$akland, fill = "grey40") +
    ggplot2::geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) +
    ggplot2::geom_sf(data = MAP$bathymetry, color = "grey60",size=.25)+
    ggplot2::coord_sf(xlim = MAP$plot.boundary$x + ext.adjust[1:2], ylim = MAP$plot.boundary$y + ext.adjust[3:4]) +
    ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) +
    ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) +
    viridis::scale_color_viridis(
      option = col.palette, begin = col.palette.limits[1], end = col.palette.limits[2],
      na.value = NA, name = legend.title, labels = scales::comma_format(big.mark = ",")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      legend.position = legend.pos, legend.margin = ggplot2::margin(0, 0, 0, 0),
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = NA, color = NA)
    )

  # add a title
  if (is.na(title.name) == F && length(title.pos) == 2) {
    densityplot <- densityplot +
      ggplot2::geom_label(
        data = data.frame(x = title.pos[1], y = title.pos[2], label = title.name),
        ggplot2::aes(x = x, y = y, label = label, hjust = 0, vjust = 1), size = 5
      )
  }
  # add a legend
  if (length(legend.pos) == 2) {
    densityplot <- densityplot +
      ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top", title.hjust = .5, barheight = barheight))
  }
  return(densityplot)
}


#' Make EFH map with akgfmaps
#' @description Create EFH map for a species/lifestage combination
#' @param region character; which region is being used, must correspond to an akgfmaps baselayer
#' @param efh.map raster of factors to be interpretted as EFH
#' @param drop a buffer that prevents outlining of small isolated pixels, makes it look nicer
#' @param survey.area sf or raster layer to use, if the desired survey area is different from akgfmaps
#' @param ext.adjust vector of length four, representing xmin,xmax, ymin,ymax, to adjust the extent
#' @param legend.pos vector of length two with decimal position of legend box
#' @param legend.labels vector of labels for the legend, should match number of factors
#' @param legend.title character; a title for the legend
#' @param col.palette character a color palette for
#' @param col.palette.limits vector of length two giving start and end points for the color palette
#' @param title.name character; a title for the figure, or NA to suppress
#' @param title.pos vector of length two with coordinates for for the legend, use NA to suppress
#' @importFrom akgfmaps get_base_layers
#' @importFrom magrittr %>%
#'
#' @return raster of EFH quantiles (subareas)
#' @export
#'
#' @examples
MakeAKGFEFHplot <- function(region,
                            efh.map,
                            drop = 10^8,
                            survey.area = "default",
                            ext.adjust = "default", # vector of length four, representing xmin,xmax, ymin,ymax, to adjust the extent
                            legend.pos = "default", # a vector of length two with a location for the legend, use NA to suppress
                            legend.labels = c("95%", "75%", "50%", "25%"),
                            legend.title = NA,
                            col.palette = "viridis",
                            col.palette.limits = c(0, 1),
                            title.name = NA,
                            title.pos = "default") {

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

  # set up some essentials for later
  if (ext.adjust == "default") {
    if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
      ext.adjust <- c(0, 0, 0, 0)
    }
    if (region %in% c("goa", "goa.west", "goa.east")) {
      ext.adjust <- c(200000, -100000, 0, 0)
    }
  }
  if (length(legend.pos)<2){
    if(legend.pos == "default") {

      if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
        legend.pos <- c(.07, .28)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
        legend.pos <- c(0.12, 0.18)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        legend.pos <- c(0.08, 0.51)
      }
    }else{
      legend.pos<-NA
    }
  }
  if (length(title.pos)<2){
    if(title.pos == "default") {
      if (region %in% c("ai", "ai.west", "ai.central", "ai.east")) {
        title.pos <- c(-900000, 490000)
      }
      if (region %in% c("ebs", "bs.all", "bs.south", "sebs", "ecs", "ebs.ecs")) {
        title.pos <- c(-550000, 800000)
      }
      if (region %in% c("goa", "goa.west", "goa.east")) {
        title.pos <- c(-1300000, 630000)
      }
    }else{
      title.pos<-NA
    }
  }
  # detect and set up the outline, for now it is required
  if (is.character(survey.area) && survey.area[1] == "default") {
    survey.sf <- MAP$survey.area
  } else {
    if (class(survey.area)[1] == "sf") {
      survey.sf <- sf::st_transform(survey.area, sf::st_crs(MAP$akland))[1:(nrow(survey.area) - 1), ]
    }
    if (class(survey.area)[1] == "RasterLayer") {
      survey.sf1<-stars::st_as_stars(is.na(survey.area))
      survey.sf2<-sf::st_as_sf(survey.sf1,merge = TRUE)
      survey.sf3<-sf::st_cast(survey.sf2,"POLYGON") # cast the polygons to polylines

      survey.sf <- sf::st_transform(survey.sf3, sf::st_crs(MAP$akland))[1:(nrow(survey.sf3) - 1), ]
    }
  }

  # set up the factor maps
  efh.vals <- raster::getValues(efh.map)
  efh.vals[efh.vals == 1] <- NA

  # convert the raster to polygons
  efhpoly0 <- stars::st_as_stars(efh.map)
  efhpoly <- sf::st_as_sf(efhpoly0,merge = TRUE)
  efhpoly2 <- efhpoly[efhpoly$layer != 1, ]

  # we'll need a new outline
  efh.dummy.raster <- raster::raster(efh.map)
  efh.vals2 <- is.na(efh.vals) == F
  efh.dummy.raster <- raster::setValues(efh.dummy.raster, values = efh.vals2)

  efhdummy0 <- stars::st_as_stars(efh.dummy.raster)
  efhdummy <- sf::st_cast(sf::st_as_sf(efhdummy0,merge = TRUE))
  efhdummy2 <- sf::st_transform(efhdummy, sf::st_crs(MAP$akland))

  # Now we need to get rid of a lot of the tiniest bits, which we'll do by dropping the smallest areas
  efhdummy.poly <- sf::st_cast(efhdummy2, "POLYGON")
  areas <- sf::st_area(efhdummy.poly)

  outside <- order(areas, decreasing = T)[1]
  toosmall <- which(as.numeric(areas) < drop)

  efh.x <- efhdummy2$layer[-c(outside, toosmall)]
  efh.y <- efhdummy2$geometry[-c(outside, toosmall)]
  efhdummy3 <- sf::st_sf(efh.x, efh.y)


  # set up the basic map, will add more customization later
  efhplot <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = survey.sf, fill = "grey95")+
    ggplot2::geom_sf(data = efhpoly2, ggplot2::aes(fill = as.factor(layer)), col = NA) +
    ggplot2::geom_sf(data = efhdummy3,fill=NA, size = .3) +
    ggplot2::geom_sf(data = MAP$akland, fill = "grey40") +
    ggplot2::geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) +
    ggplot2::geom_sf(data = MAP$bathymetry, color = "grey60",size=.25)+

  # add the themes
  efhplot <- efhplot +
    ggplot2::coord_sf(xlim = MAP$plot.boundary$x + ext.adjust[1:2], ylim = MAP$plot.boundary$y + ext.adjust[3:4]) +
    ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) +
    ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) +
    viridis::scale_fill_viridis(discrete = T, name = legend.title, labels = legend.labels) +
    ggplot2::theme_bw() +
    ggplot2:: theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      legend.position = legend.pos,
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = NA, color = NA))

  # add a title
  if (is.na(title.name) == F && length(title.pos) == 2) {
    efhplot <- efhplot +
      ggplot2::geom_label(
        data = data.frame(x = title.pos[1], y = title.pos[2], label = title.name),
        ggplot2::aes(x = x, y = y, label = label, hjust = 0, vjust = 1), size = 5
      )
  }

  # add a legend
  if (length(legend.pos) == 2) {
    efhplot <- efhplot +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, shape = 15)))
  }
  return(efhplot)
}

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
    comp.raster <- raster::cut(background, breaks = c(-Inf, Inf))

    # find which areas are EFH in each version
    old2 <- raster::cut(x = old, breaks = c(0, nonEFH + .5, Inf))
    new2 <- raster::cut(x = new, breaks = c(0, nonEFH + .5, Inf))

    vals <- raster::getValues(comp.raster)

    both <- which(raster::getValues(old2) == 2 & raster::getValues(new2) == 2)
    justold <- which(raster::getValues(old2) == 2 &
                       ((raster::getValues(new2) == 1) | is.na(raster::getValues(new2))))
    justnew <- which(((raster::getValues(old2) == 1) | is.na(raster::getValues(old2))) &
                       raster::getValues(new2) == 2)

    # assign new values to each category
    vals[justold] <- 2
    vals[justnew] <- 3
    vals[both] <- 4

    comp.raster <- raster::setValues(x = comp.raster, values = vals)

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

    area1 <- sum(raster::getValues(old) > nonEFH, na.rm = T)
    area2 <- sum(raster::getValues(new) > nonEFH, na.rm = T)
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
    ggplot2::geom_sf(data = efhpoly2, ggplot2::aes(fill = factor(layer)), col = 1, size = .05) +
    ggplot2::geom_sf(data = MAP$akland, fill = "grey40") +
    ggplot2::geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) +
    ggplot2::geom_sf(data = MAP$bathymetry, col = "grey60", size = .25) +
    ggplot2::coord_sf(xlim = MAP$plot.boundary$x + ext.adjust.x, ylim = MAP$plot.boundary$y + ext.adjust.y) +
    ggplot2::geom_label(
      data = data.frame(x = label.pos[1], y = label.pos[2], label = main),
      aes(x = x, y = y, label = main, hjust = 0, vjust = 1), size = 5
    ) +
    ggplot2::scale_fill_manual(values = c(old.col, new.col, mix.col), labels = leg.labels, name = leg.name) +
    ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) +
    ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      legend.position = leg.pos, legend.justification = c(0, 1),
      axis.title = ggplot2::element_blank(), axis.text = element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = NA, color = NA)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, shape = 15)))

  return(out.plot)
}



#' Effects plot
#'
#' @description Improved version of the effects plot that uses ggplot and AKGF to keep the style consistent.
#' @param effects.list list of data frames describing the effects for each term in a model or ensemble
#' @param region character; required region to be used with lon/lat map, must correpond to an akgfmaps baselayer,
#' @param crs CRS for lat/lon coordinates if different from akgf
#' @param nice.names data frame linking names to nicer version for publication figures
#' @param vars character vector with names or list indices to be plotted
#' @importFrom akgfmaps get_base_layers
#' @importFrom magrittr %>%
#'
#' @return list of ggplot objects containing the individual panels and effects
#' @export
#'
#' @examples
Effectsplot <- function(effects.list, region = NA, crs = NA, nice.names = NULL, vars = "all") {

  # check the variable names and restrict things to those requested
  if (length(vars)>1 &&vars != "all" & is.character(vars)) {
    vars2 <- which(names(effects.list) %in% vars)
  } else {
    if (vars == "all") {
      vars2 <- 1:length(effects.list)
    } else {
      vars2 <- vars
    }
  }

  if (length(vars2) < length(vars)) {
    missing <- vars[which(vars %in% names(effects.list) == F)]
    stop("Variables [", paste0(missing, collapse = ", "), "] not found in effects list")
  }

  # trim out any undesired effects
  effects.list2 <- list()
  for (i in 1:length(vars2)) {
    effects.list2[[i]] <- effects.list[[vars2[i]]]
    names(effects.list2)[i] <- names(effects.list)[vars2[i]]
  }

  # if names aren't supplied, use the terms from the model
  if (is.null(nice.names)) {
    nice.names <- data.frame(
      var = unlist(strsplit(names(effects.list), "[*]")),
      name = unlist(strsplit(names(effects.list), "[*]"))
    )
  }

  out.list <- list()
  for (j in 1:length(effects.list2)) {
    e.data <- effects.list2[[j]]

    # if necessary, transform the coordinates
    if (nrow(e.data) == 1600) {
      effect.names <- strsplit(names(effects.list2)[j], split = "[*]")[[1]]

      xname <- nice.names$name[which(nice.names$var == effect.names[1])]
      yname <- nice.names$name[which(nice.names$var == effect.names[2])]

      if (xname %in% c("lon", "Longitude", "long")) {
        con.data <- grDevices::contourLines(x = sort(unique(e.data$x)), y = sort(unique(e.data$y)), z = matrix(nrow = 40, ncol = 40, data = e.data$effect))

        con.data2 <- data.frame(lon = con.data[[1]]$x, lat = con.data[[1]]$y, effect = con.data[[1]]$level, group = 1)
        label.spots <- data.frame(con.data2[round(nrow(con.data2) / 2), 1:2], tag = con.data[[1]]$level)
        for (i in 2:length(con.data)) {
          c.dat <- data.frame(lon = con.data[[i]]$x, lat = con.data[[i]]$y, effect = con.data[[i]]$level, group = i)
          label.spots <- rbind(label.spots, data.frame(c.dat[round(nrow(c.dat) / 2), 1:2], tag = con.data[[i]]$level))
          con.data2 <- rbind(con.data2, c.dat)
        }
        if (nrow(label.spots) > 10) {
          picks <- seq(from = 1, to = nrow(label.spots), by = floor(nrow(label.spots) / 10))
          label.spots <- label.spots[picks, ]
        }

        MAP <- akgfmaps::get_base_layers(select.region = tolower(region), set.crs = "auto")

        if (is.na(crs)) {
          crs <- MAP$crs
        }
        # ok, now you can transform it to the right projection
        e.ef <- sf::st_as_sf(x = con.data2, coords = c("lon", "lat"), crs = crs)
        e.ef2 <- sf::st_transform(e.ef, sf::st_crs(MAP$akland))

        spots.sf <- sf::st_as_sf(x = label.spots, coords = c("lon", "lat"), crs = crs)
        spots.sf2 <- sf::st_transform(spots.sf, sf::st_crs(MAP$akland))
        spot.data2 <- data.frame(sf::st_coordinates(spots.sf2), label = spots.sf2$tag)

        e.data2 <- data.frame(sf::st_coordinates(e.ef2), e.ef2$effect, group = e.ef2$group)
        names(e.data2) <- c("lon", "lat", "effect", "group")

        ext.adjust.x <- c(0, 0)
        ext.adjust.y <- c(0, 0)
        if (tolower(region) == "goa") {
          ext.adjust.x <- c(400000, -200000)
          ext.adjust.y <- c(-490000, 600000)
          MAP$graticule$degree_label[c(1,3,5,7,9)]<-""
          MAP$lon.breaks<-c(-170,-160,-150,-140,-130)
        }
        if (tolower(region) == "ai") {
          ext.adjust.x <- c(150000, -400000)
          ext.adjust.y <- c(-390000, 500000)
        }

        var.plot <- ggplot2::ggplot() +
          ggplot2::geom_sf(data = MAP$akland, fill = "grey70") +
          ggplot2::geom_sf(data = MAP$bathymetry, col = "grey60") +
          ggplot2::geom_path(data = e.data2, ggplot2::aes(x = lon, y = lat, group = group), size = 1) +
          ggplot2::geom_label(
            data = spot.data2, ggplot2::aes(x = X, y = Y, label = label), fill = grDevices::rgb(1, 1, 1, .9),
            label.size = NA, size = 4, label.padding = ggplot2::unit(.10, "lines"), nudge_x = -1000
          ) +
          ggplot2::coord_sf(xlim = MAP$plot.boundary$x + ext.adjust.x, ylim = MAP$plot.boundary$y + ext.adjust.y) +
          ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) +
          ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.border = ggplot2::element_rect(color = "black", fill = NA),
            panel.background = ggplot2::element_rect(fill = NA, color = "black"),
            axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 12),
            plot.background = ggplot2::element_rect(fill = NA, color = NA)
          )
      } else {
        # for options other than lat.lon

        con.data <- grDevices::contourLines(
          x = sort(unique(e.data[, 1])), y = sort(unique(e.data[, 2])),
          z = matrix(nrow = 40, ncol = 40, data = e.data$effect), nlevels = 10
        )
        con.data2 <- data.frame(x = con.data[[1]]$x, y = con.data[[1]]$y, effect = con.data[[1]]$level, group = 1)
        label.spots <- data.frame(con.data2[round(nrow(con.data2) / 2), 1:2], tag = con.data[[1]]$level)

        for (i in 2:length(con.data)) {
          c.dat <- data.frame(x = con.data[[i]]$x, y = con.data[[i]]$y, effect = con.data[[i]]$level, group = i)
          label.spots <- rbind(label.spots, data.frame(c.dat[round(nrow(c.dat) / 2), 1:2], tag = con.data[[i]]$level))
          con.data2 <- rbind(con.data2, c.dat)
        }
        names(e.data)[1:2] <- c("x", "y")

        var.plot <- ggplot2::ggplot() +
          ggplot2::geom_path(data = con.data2, ggplot2::aes(x = x, y = y, group = group)) +
          ggplot2::geom_label(data = label.spots, ggplot2::aes(x = x, y = y, label = round(tag, 1)), fill = "white", label.size = NA) +
          ggplot2::xlab(xname) +
          ggplot2::ylab(yname) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            panel.border = ggplot2::element_rect(color = "black", fill = NA),
            panel.background = ggplot2::element_rect(fill = NA, color = "black"),
            axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 12),
            plot.background = ggplot2::element_rect(fill = NA, color = NA)
          )
        if (xname == "Current Velocity East (m/s)") {
          var.plot <- var.plot + ggplot2::geom_vline(xintercept = 0, linetype = 3) + ggplot2::geom_hline(yintercept = 0, linetype = 3)
        }
        if (xname == "Current Velocity East SD (m/s)") {
          var.plot <- var.plot + ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 3)
        }
      }
    }
    if (nrow(e.data) == 100) {
      xname <- nice.names$name[which(nice.names$var == names(effects.list2)[j])]

      span<-max(e.data$effect,na.rm=T)-min(e.data$effect,na.rm=T)
      upper.lim<-max(e.data$effect,na.rm=T)+3*span
      lower.lim<-min(e.data$effect,na.rm=T)-3*span

      y.lim<-c(ifelse(min(e.data$lower,na.rm=T)<lower.lim,lower.lim,NA),
               ifelse(max(e.data$upper,na.rm=T)>upper.lim,upper.lim,NA))

      # now the single dimension smoothed terms
      var.plot <- ggplot2::ggplot() +
        ggplot2::geom_line(data = e.data, ggplot2::aes(x = x, y = effect)) +
        ggplot2::geom_line(data = e.data, ggplot2::aes(x = x, y = upper), linetype = 2) +
        ggplot2::geom_line(data = e.data, ggplot2::aes(x = x, y = lower), linetype = 2)
      if(any(!is.na(y.lim))){var.plot<-var.plot+ylim(y.lim)}
      var.plot<-var.plot +ggplot2::xlab(xname) +
        ggplot2::ylab("Variable Effect") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          panel.background = ggplot2::element_rect(fill = NA, color = "black"),
          axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 12),
          plot.background = ggplot2::element_rect(fill = NA, color = NA)
        )
    }
    if (nrow(e.data) < 10) {
      # now the factors
      xname <- nice.names$name[which(nice.names$var == names(effects.list2)[j])]
      e.data$x <- as.numeric(as.character(e.data$x))

      var.plot <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = e.data, ggplot2::aes(y = effect, yend = effect, x = x - .35, xend = x + .35), size = 2) +
        ggplot2::geom_segment(data = e.data, ggplot2::aes(y = lower, yend = lower, x = x - .35, xend = x + .35), size = 1, linetype = 2) +
        ggplot2::geom_segment(data = e.data, ggplot2::aes(y = upper, yend = upper, x = x - .35, xend = x + .35), size = 1, linetype = 2) +
        ggplot2::xlab(xname) +
        ggplot2::ylab("Variable Effect") +
        ggplot2::scale_x_continuous(breaks = (e.data$x)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          panel.background = ggplot2::element_rect(fill = NA, color = "black"),
          axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 12),
          plot.background = ggplot2::element_rect(fill = NA, color = NA)
        )
    }

    # save the plot for later
    out.list[[j]] <- var.plot
    names(out.list)[j] <- names(effects.list2)[j]
  }
  return(out.list)
}
