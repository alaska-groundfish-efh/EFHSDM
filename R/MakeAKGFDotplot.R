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
#' data("region_data_all")
#' region.data <- region_data_all
#' region.data$sponge <- as.integer(region.data$sponge > 0)
#' region.data$coral <- as.integer(region.data$coral > 0)
#' region.data$pen <- as.integer(region.data$pen > 0)
#' region.data$logarea <- log(region.data$area)

#' data("raster_stack")
#' raster.stack <- terra::rast(raster_stack)

#' species.data <- subset(region.data, year >= 2012)
#' species <- "a_atf"
#' hd <- stats::quantile(species.data[species.data[, species] > 0, species], .9)

#' MakeAKGFDotplot(
#'   presence = species.data[species.data[, species] > 0, ],
#'   absence = species.data[species.data[, species] == 0, ],
#'   highdensity = species.data[species.data[, species] >= hd, ],
#'   dataCRS = terra::crs(raster.stack), region = "goa",
#'   title.name = "Adult arrowtooth flounder"
#' )

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
      legend.position.inside = legend.pos,
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
