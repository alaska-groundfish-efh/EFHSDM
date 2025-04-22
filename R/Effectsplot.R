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

        con.data2 <- data.frame("lon" = con.data[[1]]$x,
                                "lat" = con.data[[1]]$y,
                                "effect" = con.data[[1]]$level,
                                "group" = 1)
        label.spots <- data.frame(con.data2[round(nrow(con.data2) / 2), 1:2], "tag" = con.data[[1]]$level)
        for (i in 2:length(con.data)) {
          c.dat <- data.frame("lon" = con.data[[i]]$x, "lat" = con.data[[i]]$y, "effect" = con.data[[i]]$level, group = i)
          label.spots <- rbind(label.spots, data.frame(c.dat[round(nrow(c.dat) / 2), 1:2], "tag" = con.data[[i]]$level))
          con.data2 <- rbind(con.data2, c.dat)
        }
        if (nrow(label.spots) > 10) {
          picks <- seq(from = 1, to = nrow(label.spots), by = floor(nrow(label.spots) / 10))
          label.spots <- label.spots[picks, ]
        }

        MAP <- akgfmaps::get_base_layers(select.region = tolower(region), set.crs = "auto")

        if (is.na(crs)) {
          crs <- "epsg:3338"
        }
        # ok, now you can transform it to the right projection
        e.ef <- terra::vect(x = con.data2, geom = c("lon", "lat"), crs = crs)
        e.ef2 <- terra::project(e.ef, "epsg:3338")

        spots.sf <- terra::vect(x = label.spots, geom = c("lon", "lat"), crs = crs)
        spots.sf2 <- terra::project(spots.sf, "epsg:3338")
        spot.data2 <- data.frame(terra::geom(spots.sf2), label = spots.sf2$tag)

        e.data2 <- data.frame(terra::geom(e.ef2), e.ef2$effect, group = e.ef2$group)
        e.data2 <- e.data2[, c(3,4,6,7)] #select columns
        names(e.data2) <- c("lon", "lat", "effect", "group")

        ext.adjust.x <- c(0, 0)
        ext.adjust.y <- c(0, 0)
        if (tolower(region) == "goa") {
          ext.adjust.x <- c(400000, 400000)
          ext.adjust.y <- c(500000, 900000)
          MAP$graticule$degree_label[c(1,3,5,7,9)]<-""
          MAP$lon.breaks<-c(-170,-160,-150,-140,-130)
        }
        if (tolower(region) == "ai") {
          ext.adjust.x <- c(150000, -400000)
          ext.adjust.y <- c(-390000, 500000)
        }

        var.plot <- ggplot2::ggplot() +
          tidyterra::geom_spatvector(data = terra::project(terra::vect(MAP$akland), "epsg:3338"), fill = "grey70") +
          tidyterra::geom_spatvector(data = terra::project(terra::vect(MAP$bathymetry), "epsg:3338"), col = "grey60") +
          ggplot2::geom_path(data = e.data2, ggplot2::aes(x = lon, y = lat, group = group), linewidth = 1) +
          ggplot2::geom_label(
            data = spot.data2, ggplot2::aes(x = x, y = y, label = label), fill = grDevices::rgb(1, 1, 1, .9),
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
        con.data2 <- data.frame("x" = con.data[[1]]$x, "y" = con.data[[1]]$y, "effect" = con.data[[1]]$level, group = 1)
        label.spots <- data.frame(con.data2[round(nrow(con.data2) / 2), 1:2], "tag" = con.data[[1]]$level)

        for (i in 2:length(con.data)) {
          c.dat <- data.frame("x" = con.data[[i]]$x, "y" = con.data[[i]]$y, "effect" = con.data[[i]]$level, group = i)
          label.spots <- rbind(label.spots, data.frame(c.dat[round(nrow(c.dat) / 2), 1:2], "tag" = con.data[[i]]$level))
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

      # start the plot
      var.plot <- ggplot2::ggplot() +
        ggplot2::geom_line(data = e.data, ggplot2::aes(x = x, y = effect))

      # check if there are CIs, and plot the CIs if present
      if("lower"%in%names(e.data)){
        y.lim<-c(ifelse(min(e.data$lower,na.rm=T)<lower.lim,lower.lim,NA),
                 ifelse(max(e.data$upper,na.rm=T)>upper.lim,upper.lim,NA))
        var.plot<-var.plot+
          ggplot2::geom_line(data = e.data, ggplot2::aes(x = x, y = upper), linetype = 2) +
          ggplot2::geom_line(data = e.data, ggplot2::aes(x = x, y = lower), linetype = 2)
      }else{
        # If no CIs, just continue on
        y.lim<-c(lower.lim,upper.lim)
      }

      if(any(!is.na(y.lim))){var.plot<-var.plot+ggplot2::ylim(y.lim)}
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
        ggplot2::geom_segment(data = e.data, ggplot2::aes(y = effect, yend = effect, x = x - .35, xend = x + .35), linewidth = 2)

      if("lower"%in%names(e.data)){
        var.plot<-var.plot+
          ggplot2::geom_segment(data = e.data, ggplot2::aes(y = lower, yend = lower, x = x - .35, xend = x + .35), linewidth = 1, linetype = 2) +
          ggplot2::geom_segment(data = e.data, ggplot2::aes(y = upper, yend = upper, x = x - .35, xend = x + .35), linewidth = 1, linetype = 2)
      }
      var.plot+
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
