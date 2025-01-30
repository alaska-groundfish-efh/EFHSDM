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
    survey.sf <- terra::vect(MAP$survey.area)
    survey.sf <- terra::project(survey.sf, "EPSG:3338")
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
  efh.vals <- terra::values(efh.map)
  efh.vals[efh.vals == 1] <- NA

  # convert the raster to polygons
  efhpoly <- terra::as.polygons(efh.map)
  names(efhpoly) <- "lyr1"
  efhpoly2 <- efhpoly[efhpoly$lyr1 != 1, ]
  efhpoly2 <- terra::project(efhpoly2, "EPSG:3338")

  # we'll need a new outline
  efhdummy2 <- efhpoly2
  efhdummy3 <- terra::aggregate(efhdummy2)


  # set up the basic map, will add more customization later
  efhplot <- ggplot2::ggplot() +
    tidyterra::geom_spatvector(data = terra::project(terra::vect(MAP$akland), "EPSG:3338"), fill = "grey40") +
    tidyterra::geom_spatvector(data = terra::project(terra::vect(MAP$graticule), "EPSG:3338"), color = "grey70", alpha = 0.5) +
    tidyterra::geom_spatvector(data = terra::project(terra::vect(MAP$bathymetry), "EPSG:3338"), color = "grey60",linewidth=.25) +
    tidyterra::geom_spatvector(data = survey.sf, fill = "grey95") +
    tidyterra::geom_spatvector(data = efhpoly2, ggplot2::aes(fill = as.factor(lyr1)), col = NA) +
    tidyterra::geom_spatvector(data = efhdummy3,fill=NA, linewidth = .3)
  # project the bounding box for plotting terra spatvectors
  plot.boundary <- terra::ext(survey.sf)
  plot.boundary.df2 <- data.frame(x=c(plot.boundary[1], plot.boundary[2]),
                                  y=c(plot.boundary[3], plot.boundary[4]))
  plot.boundary <- plot.boundary.df2

  # add the themes
  efhplot <- efhplot +
    ggplot2::coord_sf(xlim = plot.boundary$x , ylim = plot.boundary$y ) + #linked to plot.boundary above
    ggplot2::scale_x_continuous(name = "Longitude", breaks = MAP$lon.breaks) +
    ggplot2::scale_y_continuous(name = "Latitude", breaks = MAP$lat.breaks) +
    viridis::scale_fill_viridis(discrete = T, name = legend.title, labels = legend.labels) +
    ggplot2::theme_bw() +
    ggplot2:: theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      legend.position.inside = legend.pos,
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
