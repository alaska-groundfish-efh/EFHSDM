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
    survey.sf <- terra::vect(MAP$survey.area)
    survey.sf <- terra::project(survey.sf, "epsg:3338")
  } else {
    if (class(survey.area)[1] == "sf") {
      survey.sf <- terra::project(survey.area, "epsg:3338")
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
  if (class(density.map)[1] == "SpatRaster") {
    density.dat0 <- as.data.frame(terra::xyFromCell(density.map, 1:terra::ncell(density.map)))
    vals <- terra::values(density.map)

    # convert the raster to ggplot format
    density.dat <- data.frame(lat = density.dat0$y, lon = density.dat0$x, density=vals)
    colnames(density.dat) <- c("lat", "lon", "density")
    density.sf0 <- terra::vect(x = subset(density.dat, is.na(density) == F),
                               geom = c("lon", "lat"), crs = "epsg:3338") #formerly crs = density.map@crs
    #density.sf <- sf::st_transform(density.sf0, sf::st_crs(MAP$akland))
    density.sf <- density.sf0 #don't need to transform
  }

  # often helps to remove some of the high points
  upper<-stats::quantile(density.sf$density, buffer, na.rm = T)
  density.sf$density[density.sf$density > upper] <- upper

  # project the bounding box for plotting terra spatvectors
  plot.boundary <- terra::ext(density.sf)
  plot.boundary.df2 <- data.frame(x=c(plot.boundary[1], plot.boundary[2]),
                                  y=c(plot.boundary[3], plot.boundary[4]))
  plot.boundary <- plot.boundary.df2

  # set up the basic map, will add more customization later
  #geom_sf changed to geom_spatvector from the tidyterra package
  densityplot <- ggplot2::ggplot() +
    tidyterra::geom_spatvector(data = survey.sf, fill = NA) +
    tidyterra::geom_spatvector(data = terra::project(terra::vect(MAP$akland), "epsg:3338"), fill = "grey40") +
    tidyterra::geom_spatvector(data = terra::project(terra::vect(MAP$bathymetry), "epsg:3338"), color = "grey60",linewidth=.25) +
    tidyterra::geom_spatvector(data = terra::project(terra::vect(MAP$graticule), "epsg:3338"), color = "grey70", alpha = 0.5) +
    tidyterra::geom_spatvector(data = density.sf, ggplot2::aes(col = density), size = .05) +
    ggplot2::coord_sf(xlim = plot.boundary$x , ylim = plot.boundary$y ) + #linked to plot.boundary above
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
      legend.position.inside = legend.pos, legend.margin = ggplot2::margin(0, 0, 0, 0),
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
