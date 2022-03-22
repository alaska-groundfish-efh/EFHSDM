# This script contains various functions mostly related to either mapping or plotting the output from other analyses
# In general, they should be capatible with all types of models, though be careful to set the appropriate options
#
# rpackages <- c("rgdal", "sp","gstat","viridis","raster")
#
# which_not_installed <- which(rpackages %in% rownames(installed.packages()) == FALSE)
#
# if(length(which_not_installed) > 1){
#   install.packages(rpackages[which_not_installed], dep = TRUE)
# }
# rm(rpackages,which_not_installed)
#
# require(rgdal)
# require(sp)
# require(gstat)
# require(viridis)
# require(raster)


#' Load EFH data
#'
#' @description Quickly set up the covariate rasters.
#' @details  Should be called before the Load Map function. It would be a good idea to integrate the two, but I haven't gotten to it. Also, it is unlikely to useful except for this project. Deprecated: Likely to be deleted soon.
#'
#' @param region character; "AI", "GOA", or "EBS"
#' @param raster.path character; path to the main folder for the EFH project
#'
#' @return does not return anything, but adds various things to the global environment
#' @export
#'
#' @examples
LoadEFHData<-function(region,raster.path="//akc0ss-n086/SEA_Programs/RACE_EFH_variables/Variables"){

  region<-toupper(region)
  region2<-ifelse(region=="EBS","BS",region)

  # (network location is \\akc0ss-n086/SEA_Programs/RACE_EFH_variables)
  bathy <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Bathy"))
  slope <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Slope"))
  tmax <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Tmax"))
  btemp <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Btemp"))
  btemp<-raster::crop(x = btemp,y=bathy)
  BPI <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/BPI"))
  BPI<-raster::crop(x = BPI,y=bathy)

  AspectE <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Aspect_East"))
  AspectN <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Aspect_North"))

  Curve <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Curve_Mean"))

  lat <- raster::init(bathy, v ='y')
  lat <- raster::mask(lat, bathy,overwrite = F)
  lon <- raster::init(bathy, v ='x')
  lon <- raster::mask(lon, bathy,overwrite = F)
  coral <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Coralfactor"))
  sponge <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Spongefactor"))
  whips <- raster::raster(paste0(raster.path,"/Variables_",region,"_1km/Whipsfactor"))

  east<-raster::raster(paste0(raster.path,"/Variables_",region,"_1km/ROMSbcurrentEastings"))
  north<-raster::raster(paste0(raster.path,"/Variables_",region,"_1km/ROMSbcurrentNorthings"))
  eastSD<-raster::raster(paste0(raster.path,"/Variables_",region,"_1km/ROMSbEastingsSD"))
  northSD<-raster::raster(paste0(raster.path,"/Variables_",region,"_1km/ROMSbNorthingsSD"))
  if(region=="EBS"){
    phi <- raster::raster(paste0(raster.path,"/Variables_EBS_1km/phi"))

    raster_stack <- raster::stack(lon,lat,bathy,slope,AspectE,AspectN,Curve,btemp,east,north,eastSD,northSD,tmax,phi,BPI, sponge, coral, whips)
    names(raster_stack) <- c("lon","lat","bdepth","slope","aspectE","aspectN","curve","btemp","bcurrentU","bcurrentV",
                             "bcurrentUSD","bcurrentVSD", "tmax","phi","BPI","sponge","coral","pen")
  }
  # GOA and AI don't have sediment grabs to calculate phi, so there is a "rockiness" variable instead
  if(region%in%c("AI","GOA")){
    rocky<-raster::raster(paste0(raster.path,"/Variables_",region,"_1km/rocky"))

    raster_stack <- raster::stack(lon,lat,bathy,slope,AspectE,AspectN,Curve,btemp,east,north,eastSD,northSD,tmax,rocky,BPI, sponge, coral, whips)
    names(raster_stack) <- c("lon","lat","bdepth","slope","aspectE","aspectN","curve","btemp","bcurrentU","bcurrentV",
                             "bcurrentUSD","bcurrentVSD", "tmax","rocky","BPI","sponge","coral","pen")
  }
  # rather than return a value, the finished raster.stack is written direct to the global environment
  # that way other function can use it as a default value
  assign(x = "raster.stack", value = raster_stack, envir = .GlobalEnv)
  assign(x="covars.vec",value=names(raster_stack), envir=.GlobalEnv)
}

# Deprecated: Likely to be deleted soon
# This function loads the various map parameters used to make plots of different regions.
# The function assigns various values directly to the global environment, rather than
# returning them.
# Note that loading the map requires an appropriate parameters file, described elsewhere.
# NOte, most of this is somewhat unnecessary as we are now using the akgfmaps package for the
# publication quality figures, but it might still be useful for quick and dirty plotting.

#' Title
#'
#' @param region character; "AI", "GOA", or "EBS"
#' @param parameter.file a csv file containing various mapping paramters
#' @param covariate.raster the raster stack produced by LoadEFHData
#' @param coast.file file path for a shapefile of the Alaska landmass, very slow to load
#' @param GOA.mask file path an additional shapefile more specific to the GOA survey area
#'
#' @return does not return anything, but adds various things to the global environment
#' @export
#'
#' @examples
LoadMap<-function(region,parameter.file="G:/Harris/EFH_copy/Map_Settings.csv",covariate.raster=raster.stack,
                  coast.file="//akc0ss-n086/SEA_Programs/RACE_EFH_variables/shapefiles",
                  GOA.mask="//akc0ss-n086/SEA_Programs/RACE_EFH_variables/shapefiles"){
  region<-toupper(region)

  aea.proj <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

  # this is all for setting up the mapping parameters
  all.region.data<-read.csv(parameter.file)[-1]
  region.data<-subset(all.region.data,Region==region)

  yaxis<-subset(region.data,Var=="yaxis")[,3]
  xaxis<-subset(region.data,Var=="xaxis")[,3]

  yaxis.ticks <- rgdal::project(as.matrix(subset(region.data,Var=="yaxis.ticks")[,3:4]), aea.proj)
  xaxis.ticks <- rgdal::project(as.matrix(subset(region.data,Var=="xaxis.ticks")[,3:4]), aea.proj)
  yaxis.ticks <- yaxis.ticks[,2]
  xaxis.ticks <- xaxis.ticks[,1]

  ext.pol <- sp::Polygon(cbind(subset(region.data,Var=="ext.x")$Value1,
                               subset(region.data,Var=="ext.y")$Value1))
  ext.pol <- sp::Polygons(list(ext.pol), region)

  # Note, GOA and EBS are already transformed
  #if(region%in%c("EBS","GOA")){
  if(region%in%c("GOA","EBS")){
    ext.pol.proj<-sp::SpatialPolygons(list(ext.pol), proj4string = CRS(aea.proj))
  }else{
    ext.pol <- sp::SpatialPolygons(list(ext.pol), proj4string = CRS("+proj=longlat +datum=WGS84"))
    ext.pol.proj <- sp::spTransform(ext.pol, CRS(aea.proj))
  }

  wrld_p <- maps::map("world", interior = FALSE, plot = FALSE)
  llCRS <- sp::CRS("+proj=longlat +ellps=WGS84")
  wrld_sp <- maptools::map2SpatialLines(wrld_p, proj4string = llCRS)
  prj_new <- sp::CRS("+proj=moll")
  wrld_proj <- sp::spTransform(wrld_sp, prj_new)
  wrld_grd <- sp::gridlines(wrld_sp, easts = c(-178, seq(-176,178, 2), 180), norths = seq(-75, 75, 2), ndiscr = 100)
  wrld_grd_proj2 <- sp::spTransform(wrld_grd, CRS(aea.proj))

  assign(x = "x.axis", value = xaxis, envir =.GlobalEnv)
  assign(x = "y.axis", value = yaxis, envir =.GlobalEnv)

  assign(x = "x.ticks", value = xaxis.ticks, envir =.GlobalEnv)
  assign(x = "y.ticks", value = yaxis.ticks, envir =.GlobalEnv)

  assign(x="map.extent", value=ext.pol.proj, envir =.GlobalEnv)

  assign(x = "world.grid", value = wrld_grd_proj2, envir =.GlobalEnv)

  alaska.coast <- rgdal::readOGR(dsn = coast.file, layer = "namerica_dcw", verbose = F)
  alaska.raster<-raster::raster(covariate.raster)
  alaska.raster2<-raster::rasterize(alaska.coast,alaska.raster)

  if(region=="GOA"){
    mask1<-rgdal::readOGR(dsn = GOA.mask, layer = "completeGOAgrid", verbose = F)
    mask2<-raster::raster(covariate.raster)
    mask3<-raster::rasterize(mask1,mask2)
    assign("GOA.mask",value=mask3,envir=.GlobalEnv)
  }

  assign(x= "ak.coast", value = alaska.coast, envir = .GlobalEnv)
  assign(x= "ak.raster", value = alaska.raster2, envir = .GlobalEnv)
}

#' Calculate RMSE
#'
#' @description A quick function to calculate the RMSE.
#' @param pred vector of predictions
#' @param obs vector of observations
#'
#' @return returns RMSE of obs and preds
#' @export
#'
#' @examples
RMSE<-function(pred,obs){
  keep<-which(is.na(pred)==F & is.na(obs)==F)
  return(sqrt(sum((pred[keep]-obs[keep])^2)/length(pred[keep])))
}

#' Calculate PDE
#'
#' @description A quick function to calculate the PDE
#' @param pred vector of predictions
#' @param obs vector of observations
#'
#' @return returns the estimate percent deviance explained, assuming a Poisson distribution
#' @export
#'
#' @examples
PDE<-function(obs,pred){
  term1<-obs*log(obs/mean(obs))
  term1[is.nan(term1)]<-0
  term2<-obs-mean(obs)

  nulldev<-2*sum(term1-term2)

  pred[pred==0]<-.00001
  term1<-obs*log(obs/pred)
  term2<-obs-pred

  pdev<-2*sum(term1-term2)
  return(1-(pdev/nulldev))
}



#' Add grid
#'
#' @description Add a map of Alaska and other features to an existing map
#' @details This function adds a map of AK and other features to an existing map, including a legend. It has many options for customizing the appearance. The whole "AddGrid" system is old and ought to be removed in favor of the newer akgfmaps plotting, when I get the time. Deprecated: Likely to be deleted soon.
#' @param legPosition
#' @param horiz
#' @param legName
#' @param legVals
#' @param legend.size
#' @param col.vec
#' @param depth
#' @param dlevels
#' @param dlabels
#' @param dlabel.size
#' @param bathy
#' @param grid
#' @param ext
#' @param map.grid
#' @param grid.col
#' @param xticks
#' @param yticks
#' @param xaxis
#' @param yaxis
#' @param axistext.size
#' @param coast
#' @param land.col
#' @param border.col
#' @param plot.yaxis
#' @param plot.xaxis
#'
#' @return
#' @export
#'
#' @examples
AddGrid<-function(legPosition="bottomleft",                # where should the legend be drawn, NA for no legend
                  horiz=F,                                 # should the legend be in horizontal format
                  legName="Percentiles",                   # Name of legend, if applicable
                  legVals=c("95%","75%","50%","25%"),      # legend labels for the EFH cuts
                  legend.size=1,                           # multiplier for legend size
                  col.vec=viridis(4),                      # colors for the EFH categories
                  depth=F,                                 # Should depth contours be drawn
                  dlevels=c(100,200,500),                  # depth contours to be shown, if depth=T
                  dlabels=F,                               # should depth contours be labelled
                  dlabel.size=.5,                          # multiplier for contour labels
                  bathy=raster.stack$bdepth,               # raster used to draw the depth contours
                  grid=F,                                  # should a lat/lon grid be added
                  ext=map.extent,                          # the extent to plot,
                  map.grid=world.grid,                     # the lat/lon grid to plot
                  grid.col=rgb(.8,.8,.8,.5),               # a color for the grid lines
                  xticks=x.ticks,yticks=y.ticks,           # locations for the ticks on the grid
                  xaxis=x.axis,yaxis=y.axis,               # labels for the ticks
                  axistext.size=1,                         # multiplier for axis text size
                  coast=ak.coast,                          # a raster or sp object to plot the land
                  land.col="grey30",                        # color for the land
                  border.col="grey20",                      # color for the border of land polygons
                  plot.yaxis=T,plot.xaxis=T){              # should the x or y axis be plotted
  suppressWarnings(plot(coast, col = land.col, add = TRUE,legend=F,border=border.col))
  if(suppressWarnings(is.na(xaxis[1])==F)){
    if(xaxis[1]==""){xaxis<-rep("",times=length(xticks))}
    axis(1, at = xticks, labels = xaxis,cex.axis=axistext.size,gap.axis = .25)
  }
  if(suppressWarnings(is.na(yaxis[1])==F)){
    if(yaxis[1]==""){yaxis<-rep("",times=length(yticks))}
    axis(2, at = yticks, labels = yaxis,cex.axis=axistext.size,gap.axis = .25)
  }
  if(grid){ plot(world.grid, lty = 3, col =grid.col, add = T,lwd=.1,ext=map.extent)}
  legend(legPosition, legend = legVals, pch = 15, col = col.vec, bty = "n", pt.cex = 2 ,
         title = legName, cex = legend.size,horiz=horiz)
  if(depth==T){
    raster::contour(raster.stack[["bdepth"]], levels = dlevels, col = rgb(0,0,0,.6),
            lwd = 0.05, drawlabels = dlabels, add = TRUE,labcex=dlabel.size)
  }
}


#' Plot EFH (soon to be deprecated)
#' @description A function that provides a quick method of making EFH maps. Currently replaced with akgfmaps function.
#' @param map a raster with factor values equivalent to EFH
#' @param map.ext an extent object if one wishes to narrow the plot
#' @param outline logical; should an outline be drawn around the EFH areas
#' @param outline.lwd numeric; the thickness of the outline
#' @param col.vec vector of colors to be used in the map
#' @param title character; a main title
#' @param background a color for non-EFH areas
#' @param ylab character; a y axis label
#' @param xlab character; a x axis label
#' @param zlim a vector of length two, which limits the factors
#' @param label.size numeric; the size of the map labels
#'
#' @return does not return anything, but makes a plot
#' @export
#'
#' @examples
plotEFH<-function(map,                                      # raster with factor values equivalent to EFH
                  map.ext=NULL,                             # an extent object to narrow the plot
                  outline=F,                                # Should an outline be drawn around EFH spots
                  outline.lwd=1,                            # size of outline line
                  col.vec=viridis(4),                       # a set of colors or a color scale
                  title="",                                 # main title
                  background=NULL,                          # background color for areas that are not EFH
                  ylab = "Latitude", xlab = "Longitude",    # labels for the lat/lon
                  zlim=NA,                                  # a vector of length two, limiting the factors
                  label.size=1){                            # multiplier for axis labels

  if(is.null(map.ext)){map.ext<-raster::extent(map)}

  suppressWarnings(if(is.na(zlim)==T){
    zlim = c(2-as.integer(is.null(background)==F), raster::maxValue(map))
  })

  raster::plot(map, main = title, xaxt = "n", yaxt = "n", box = F, col = c(background,col.vec),
       ext = map.ext, legend = FALSE, ylab = ylab, xlab = xlab, horiz = TRUE,
       zlim = zlim,cex.lab=label.size)
  if(outline==T){
    dummy<-raster::setValues(map,values = as.integer(raster::getValues(map)>=2))
    dummy2<-raster::setValues(map,values = as.integer(is.na(raster::getValues(map))))
    raster::contour(x = dummy2,add=T,drawlabels=F,lwd=outline.lwd,levels=.25,method="simple")
  }
}


#' Plot abundance (deprecated)
#'
#' @description This function primarily provides abundance or CPUE-type maps. However, it does a decent job with any raster that is on a continuous scale, and is highly customizable. Deprecated: Likely to be deleted soon in favor of the AKGF version
#' @param map
#' @param map.ext
#' @param outline
#' @param outline.lwd
#' @param col.vec
#' @param legend.name
#' @param legend.text
#' @param legend
#' @param title
#' @param zmin
#' @param zquant
#' @param zmax
#' @param center.scale
#' @param back.col
#' @param ylab
#' @param xlab
#' @param label.size
#'
#' @return
#' @export
#'
#' @examples
plotAbundance<-function(map,                                # A raster of abundances, or a similar metric
                        map.ext=NULL,                       # an optional extent object to narrow the plot
                        outline=F,                          # should an outline be drawn around the study area
                        outline.lwd=.5,                     # size of the outline
                        col.vec=plasma(255),                # a scale of colors to be used
                        legend.name="Abundance (count)",    # a name for the legend
                        legend.text=.65,                    # mulitplier for legend text
                        legend=T,                           # should the legend be drawn
                        title="",                           # a main title for the plot
                        zmin=0,                             # a minimum value for the plot data, lower values are changed to this
                        zquant=.98,                         # if no zmax is supplied, a quantile to compute one
                        zmax=NA,                            # a maximum value for the plot data, higher values are reduced to this
                        center.scale=F,                     # make the scale centered on zero (useful for some plots)
                        back.col=NA,                        # a background color for NA values
                        ylab="Latitude",xlab="Longitude",   # labels for the axes
                        label.size=1){                      # multiplier for axis label size

  abund.raster<-map
  if(is.null(map.ext)){map.ext<-raster::extent(map)}

  # if no zmax supplied, sample the raster and choose the 95th percentile, so this rules out the rare
  # extreme predictions that tend to screw up the scale
  if(is.na(zmax)){

    sample <- raster::sampleRandom(map,min(10000,ncell(map)), na.rm = TRUE)
    zmax <- stats::quantile(sample[is.finite(sample)], probs = zquant, na.rm = TRUE, names = FALSE)
  }
  if(is.na(zmin)){
    sample <- raster::sampleRandom(map,min(10000,ncell(map)), na.rm = TRUE)
    zmin <- stats::quantile(sample[is.finite(sample)], probs = 1-zquant, na.rm = TRUE, names = FALSE)
  }
  if(center.scale==T){
    zmax<-max(abs(zmin),abs(zmax))
    zmin<--zmax
  }

  abund.raster[raster::getValues(abund.raster)>zmax]<-zmax
  abund.raster[raster::getValues(abund.raster)<zmin]<-zmin

  raster::plot(abund.raster, main = title, xaxt = "n", yaxt = "n", box = F, col = col.vec, ext = map.ext,
       legend.shrink = 0.5, axis.args = list(cex.axis = legend.text*.8),cex.lab=label.size,legend=legend,
       legend.args = list(text = legend.name, cex = legend.text,cex.lab = legend.text, side = 1, line = 2),
       horiz = TRUE, ylab = ylab, xlab = xlab,zlim = c(zmin, zmax),colNA=back.col)
  if(outline==T){
    dummy<-raster::setValues(map,values = as.integer(is.na(raster::getValues(abund.raster))))
    raster::contour(x = dummy,add=T,drawlabels=F,lwd=outline.lwd,levels=.25,method="simple")
  }
}

#' Dot plot of sample locations (soon to be deprecated)
#'
#' @description This function plots the sample locations for a species because you need to interleave different pieces, all the parameters for the AddGrid function are needed. Will expand to add background
#' @param train.data
#' @param test.data
#' @param species
#' @param figure.name
#' @param legPosition
#' @param background
#' @param background.file
#' @param coast
#' @param back.col
#' @param train.col
#' @param test.col
#' @param point.size
#' @param text.size
#' @param legendtext.size
#' @param axistext.size
#' @param horiz
#' @param title
#' @param dlabels
#' @param dlevels
#' @param dlabel.size
#'
#' @return
#' @export
#'
#' @examples
plotDots<-function(train.data,                     # data set with lat and lon, used for training
                   test.data=NA,                   # data set with lat and lon used for tests (optional)
                   species,                        # species name or index number for the data
                   figure.name=NULL,               # a figure name to be displayed
                   legPosition="bottomleft",       # position for the legend
                   background=raster.stack,        # background raster to be plotted before points
                   background.file=NA,             # filename to load for background
                   coast=ak.coast,                 # coast shapefile for the background
                   back.col=NA,                    # background color, NA will plot a blank background
                   train.col=plasma(4)[3],         # color for main points
                   test.col=plasma(4)[1],          # second color for validation points
                   point.size=1,                   # multiplier for the point sizes
                   text.size=1,                    # multiplier for axis labels and text
                   legendtext.size=1,              # multiplier for legend text
                   axistext.size=1,                # multiplier for axis text
                   horiz=F,                        # should the legend be in horizontal format
                   title=NA,                       # a title to go above the figure
                   dlabels=T,                      # should depth contours be labeled
                   dlevels=c(100,200,500),         # levels for the depth contour, NA will suppress plotting
                   dlabel.size=1){                 # multiplier for depth label size
  train<-train.data[train.data[,species]>0,]

  if(is.na(background.file)==F){
    background=raster(background.file)
  }

  n.presence<-sum(train.data[,species]>0)+ifelse(is.na(test.data),0,sum(test.data[,species]>0))
  if(n.presence>1000){
    train<-train[sample(x = 1:nrow(train),size = 1000,replace = F),]
  }

  if(is.na(test.data)==F){
    leg.text=c("Training","Testing")
    leg.col=c(train.col,test.col)
  }else{
    leg.text=""
    leg.col=NA
  }

  # needs to plot background, then the land, then the dots last
  raster::plot(background$bdepth,col=back.col,xaxt = "n", yaxt = "n",legend = FALSE, ylab = "Latitude",
       xlab = "Longitude",main=title,cex.lab=text.size)

  AddGrid(legPosition =NA,depth=T,land.col="grey30",dlevels=dlevels,dlabels=dlabels,
          axistext.size = axistext.size,dlabel.size = dlabel.size,coast = coast)
  if(is.null(figure.name)){
    legend(legPosition,legend=leg.text, pch = c(16,1), col = leg.col,
           bty = "n", cex = legendtext.size,horiz=horiz,title = paste0("n = ",n.presence))
  }else{
    legend(legPosition,title=figure.name, pch = c(16,1), col = leg.col,
           bty = "n", cex = legendtext.size,horiz=horiz,legend = paste0("n = ",n.presence))
  }
  points(x=train$lon,y=train$lat,col=train.col,pch=1,cex=.6*point.size)
  if(is.na(test.data)==F){
    test<-test.data[test.data[,species]>0,]
    points(x=test$lon,y=test$lat,col=test.col,pch=1,cex=1*point.size)
  }
}



#' EFH map comparison
#' @description This is a function to quickly make a raster that compares the differences between two EFH rasters. The EFH rasters are usually composed of factor values, so you need to know the meaning of said values.  Usually, 1 is <5%, 2 is 5-25%, 3 is 25-50%, 4 is 50-75%, and 5 is >75%.
#' @param old raster; map of the original EFH area
#' @param new raster; map of a new EFH area to be compared to old
#' @param nonEFH integer; the value corresponding to non-EFH
#' @param background raster; any raster that can be used to mask the relevant area
#'
#' @return returns a new raster, with coding 1 = non-EFH, 2 = EFH in old but not new;
#' 3 = EFH in new but not old, and 4 = EFH in both
#' @export
#'
#' @examples
EFHComparison<-function(old,                    # a EFH raster (with discrete, ordered values)
                        new,                    # a second EFH raster to be compared to the first
                        nonEFH=1,               # a value equal to or less than is considered not EFH
                        background){            # a raster that can be used as a mask for the relevant area

  # this assigns a default value of 1 to anywhere not an NA
  out.raster<-raster::cut(background,breaks=c(-Inf,Inf))

  # find which areas are EFH in each version
  old2<-raster::cut(x = old,breaks=c(0,nonEFH+.5,Inf))
  new2<-raster::cut(x = new,breaks=c(0,nonEFH+.5,Inf))

  vals<-raster::getValues(out.raster)

  both<-which(raster::getValues(old2)==2 & raster::getValues(new2)==2)
  justold<-which(raster::getValues(old2)==2 &
                   ((raster::getValues(new2)==1)| is.na(raster::getValues(new2))))
  justnew<-which(((raster::getValues(old2)==1)| is.na(raster::getValues(old2))) &
                   raster::getValues(new2)==2)

  #assign new values to each category
  vals[justold]<-2
  vals[justnew]<-3
  vals[both]<-4

  out.raster<-raster::setValues(x = out.raster,values = vals)
  return(out.raster)
}


#' Plot EFH comparison
#' @description Deprecated: Likely to be deleted soon in favor of AKGF version - a function to plot the comparison map produced by the above function
#'
#' @param map raster; map of factor values, usually produced by EFHComparison
#' @param map.ext an extent that can limit the plot area
#' @param col.vec vector of colors to be used in the map
#' @param title character; a main title
#' @param ylab character; a label for the y axis
#' @param xlab character; a label for the x axis
#'
#' @return nothing; but creates a plot
#' @export
#'
#' @examples
plotComparison<-function(map,                                                       # the EFH comparison raster
                         map.ext=NULL,                                              # an extent to limit the plotting
                         col.vec = c("grey85","#D6556DFF","#43DF71FF","#440154C1"), # colors for the categories
                         title="",                                                  # a main title
                         ylab = "Latitude", xlab = "Longitude"){                    # labels for the axes

  if(is.null(map.ext)){map.ext<-extent(map)}
  raster::plot(map, main = title, xaxt = "n", yaxt = "n", box = F, col = col.vec,
       ext = map.ext, legend = FALSE, ylab = ylab, xlab = xlab)
}



#' Find EFH breaks
#' @description Find the break points used for EFH for a given abundance raster, given settings.
#' @details This has undergone some recent changes. With the additional of the "sanity check" elsewhere, we no longer recommend supplying the data set. The default threshold of .0513 is the poisson abundance equivalent to a 5% encounter prob.
# and the project has vacillated between the percentile and cumulative methods quite a bit.
#' @param abund.raster raster; map of predicted abundance
#' @param method character; "cumulative" or "percentile" for method of EFH calculation
#' @param threshold numeric; a threshold to use with the "percentile" method, default is equivalent to 5% prob in Poisson
#' @param quantiles vector of quantiles (other than 0 and 1, that are desired)
#' @param data optional data frame with columns "lat" & "lon" to draw the sample from, instead of entire raster
#'
#' @return vector of breaks for the specified quantiles
#' @export
#'
#' @examples
FindEFHbreaks<-function(abund.raster,                  # an abundance raster
                        method="cumulative",           # a method, currently "cumulative" or "percentile"
                        threshold=.0513,               # a threshold to use with the "percentile" method, default is equivalent to 5% prob
                        quantiles=c(.05,.25,.5,.75),   #
                        data=NULL){                    #

  # format quantiles and correct for any errors
  quants<-sort(unique(c(0,quantiles,1)))

  # choose an EFH method
  if(method=="percentile"){
    sample <- na.omit(raster::getValues(abund.raster))
    sample[sample <= threshold] <- NA
    breaks <- quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
    breaks[1]<-0
    breaks[length(breaks)]<-Inf
  }
  if(method=="cumulative"){
    # Decide whether to sample at given locations or use the whole thing
    if(is.null(data)){
      vals<-na.omit(sort(raster::getValues(abund.raster)))
      vals2<-cumsum(vals)/sum(vals)
    }else{
      vals<-na.omit(sort(raster::extract(abund.raster,data.frame(data$lon,data$lat))))
      vals2<-cumsum(vals)/sum(vals)
    }

    # Loop to calculate the breaks
    breaks<-c(0,rep(NA,length(quants)-2),Inf)
    while(length(unique(na.omit(breaks)))!=length(quants)){
      for(j in 2:(length(quants)-1)){
        breaks[j]<-vals[which(vals2>quants[j])[1]]
      }
      vals<-vals[-length(vals)]
      vals2<-cumsum(vals)/sum(vals)
    }
  }
  return(breaks)
}



#' Cross-validate the model
#' @description Implement n-fold cross-validation, based on either pre-defined folds or randomly
#' @details Works for all models currently in use, including maxnet. For quasi-poisson or negbinom, just use type="gam" and it'll work. Outputs a list with 2 or more elements, so you may need to fish through it to find what you want.
#'
#' @param model a fitted model from either mgcv::gam or maxnet::maxnet
#' @param regmult a regularization for maxnet models
#' @param model.type character; the type of model being made ("maxnet","cloglog","hgam","gam")
#' @param scale.preds should the predictions be scaled (assumes count data for now)
#' @param data data frame, presumably the same one used to fit the model
#' @param key character; an identifier for each record, such as hauljoin
#' @param species character; the species or a column in the data set, optional for gams
#' @param folds integer; number of folds if random CV is to be used, otherwise ignored
#' @param group character; a column or variable name to be used for pre-defined folds in CV, or "random"
#'
#' @return a list of 4 elements; 1- data frame with observations, predictions, and CV out-of-bag predictions
#'                               2- list of models generated for each CV fold
#'                               3- scale factor for the original model
#'                               4- vector of scale factors used for each CV model
#' @export
#'
#' @examples
CrossValidateModel<-function(model,
                             regmult=1,
                             model.type,
                             scale.preds=F,
                             data,
                             key=NA,
                             species=NA,
                             folds=10,
                             group="random"){

  if(model.type!="maxnet"){
    species<-ifelse(model$family$family=="ziplss",as.character(formula(model)[[1]])[[2]],as.character(formula(model))[[2]])
  }

  if(model.type%in%c("maxnet","cloglog","hgam","gam")==F){
    stop("Model type not recognized")
  }
  model.list<-list()
  scale.factor=1

  # first, check if we need to do randomization, going to apportion things so that even distribution of
  # presences/absences is evenish
  if(group=="random"){
    # a bunch of checks figure out how many presence and absences
    n.tot<-nrow(data)
    n.pres<-sum(data[,species]>0)
    pres<-which(data[,species]>0)
    n.abs<-sum(data[,species]==0)
    abs<-which(data[,species]==0)

    pres.base<-floor(n.pres/folds)
    pres.add<-n.pres%%folds
    abs.base<-floor(n.abs/folds)
    abs.add<-n.abs%%folds

    # basically, assign the presences and absences separately, so they are even
    pres.scheme<-c(rep(pres.base+1,times=pres.add),rep(pres.base,times=folds-pres.add))
    abs.scheme<-c(rep(abs.base,times=folds-abs.add),rep(abs.base+1,times=abs.add))

    pres.randos<-sample(1:n.pres,size=n.pres,replace=F)
    abs.randos<-sample(1:n.abs,size=n.abs,replace=F)

    pres.group<-rep(LETTERS[1],times=pres.scheme[1])
    abs.group<-rep(LETTERS[1],times=abs.scheme[1])
    for(i in 2:folds){
      pres.group<-c(pres.group,rep(LETTERS[i],times=pres.scheme[i]))
      abs.group<-c(abs.group,rep(LETTERS[i],times=abs.scheme[i]))
    }
    pres.data<-data.frame(data[pres[pres.randos],],group=pres.group)
    abs.data<-data.frame(data[abs[abs.randos],],group=abs.group)

    data<-rbind(pres.data,abs.data)
    group<-"group"
  }

  # now, we can get down to business, first set up some basic information
  n.folds<-length(unique(data[,group]))
  fold.vec<-sort(unique(data[,group]))
  fold.table<-table(data[,group])
  start.vec<-cumsum(c(1,fold.table[1:(length(fold.table))-1]))
  names(start.vec)<-names(fold.table)
  end.vec<-cumsum(fold.table)

  out.names<-NULL
  #Set up the output table
  if(is.na(key)==F){
    out.names<-key
  }
  if(sum(c("lon","lat")%in%names(data))==2){
    out.names<-c(out.names,"lon","lat")
  }

  out.names<-c(out.names,group,"abund","pred","prob","cvpred","cvprob","error","cverror")
  error.data<-as.data.frame(matrix(data=NA,nrow = nrow(data),ncol=length(out.names)))
  colnames(error.data)<-out.names

  if(scale.preds){scale.factors<-rep(NA,length=n.folds)}

  #progress bar
  pb <- txtProgressBar(min = 0, max = n.folds, style = 3)

  # For each fold, first set up the data, then use the appropriate method to generate predictions
  for(i in 1:n.folds){
    train.data<-data[data[,group]!=fold.vec[i],]
    test.data<-data[data[,group]==fold.vec[i],]


    error.data[start.vec[i]:end.vec[i],group]<-fold.vec[i]
    if(is.na(key)==F){
      error.data[start.vec[i]:end.vec[i],1]<-test.data[,key]
    }
    if("lon"%in%out.names){
      error.data$lon[start.vec[i]:end.vec[i]]<-test.data$lon
      error.data$lat[start.vec[i]:end.vec[i]]<-test.data$lat
    }
    error.data[start.vec[i]:end.vec[i],group]<-fold.vec[i]
    error.data$abund[start.vec[i]:end.vec[i]]<-test.data[,species]

    if(model.type=="maxnet"){
      preds<-exp(predict(object = model,newdata=test.data,response="link")+model$entropy)
      probs<-predict(object = model,newdata=test.data,type="cloglog")
      # then on to the cv model
      vars0<-names(model$samplemeans)
      facs<-vars0[vars0%in%names(model$varmax)==F]

      try(cv.model<-FitMaxnet(data = train.data,species = species,vars = names(model$varmax),facs = facs,regmult = regmult))
      if(exists("cv.model")){
        cvpreds<-exp(predict(object = cv.model,newdata=test.data,response="link")+cv.model$entropy)
        cvprobs<-predict(object = cv.model,newdata=test.data,type="cloglog")
      }else{
        cvpreds<-rep(NA,times=nrow(test.data))
        cvprobs<-rep(NA,times=nrow(test.data))
        cv.model<-NA
      }
    }
    if(model.type=="cloglog"){
      preds<-exp(mgcv::predict.gam(object = model,newdata=test.data,type="link"))
      probs<-mgcv::predict.gam(object = model,newdata=test.data,type="response")

      try(cv.model<-FitGAM(data = train.data,reduce=F,family.gam = "binomial",select=F,
                           link.fx = "cloglog",gam.formula = formula(model),verbose = F))
      if(exists("cv.model")){
        cvpreds<-exp(mgcv::predict.gam(object = cv.model,newdata=test.data,type="link"))
        cvprobs<-mgcv::predict.gam(object = cv.model,newdata=test.data,type="response")
      }else{
        cvpreds<-rep(NA,times=nrow(test.data))
        cvprobs<-rep(NA,times=nrow(test.data))
        cv.model<-NA
      }
    }
    if(model.type=="hgam"){
      preds<-mgcv::predict.gam(object = model,newdata=test.data,type="response")
      probs<-1-exp(-exp(mgcv::predict.gam(model,newdata=test.data)[,2]))

      try(cv.model<-FitHurdleGAM(density.formula = formula(model)[[1]],prob.formula = formula(model)[[2]],
                                 data = train.data,reduce = F,verbose = F,select = F))
      if(exists("cv.model")){
        cvpreds<-mgcv::predict.gam(object = cv.model,newdata=test.data,type="response")
        cvprobs<-1-exp(-exp(mgcv::predict.gam(cv.model,newdata=test.data)[,2]))

      }else{
        cvpreds<-rep(NA,times=nrow(test.data))
        cvprobs<-rep(NA,times=nrow(test.data))
        cv.model<-NA
      }
    }
    if(model.type=="gam"){
      if(strsplit(model$family$family,split="[()]")[[1]][1]=="Negative Binomial"){
        gamfam<-"nb"
        theta<-as.numeric(strsplit(model$family[[1]],split="[()]")[[1]][2])
        probs<-1-dnbinom(0,mu = mgcv::predict.gam(model,newdata=test.data,type="response"),size = theta)
      }else{
        gamfam<-model$family$family
        probs<-(1-dpois(0,mgcv::predict.gam(object = model,newdata=test.data,type="response")))
      }
      preds<-mgcv::predict.gam(object = model,newdata=test.data,type="response")

      try(cv.model<-FitGAM(data = train.data,reduce=F,family.gam = gamfam,select=F,
                           link.fx = model$family$link,gam.formula = formula(model),verbose = F))
      if(exists("cv.model")){
        cvpreds<-mgcv::predict.gam(object = cv.model,newdata=test.data,type="response")
        if(strsplit(model$family$family,split="[()]")[[1]][1]=="Negative Binomial"){
          cvtheta<-as.numeric(strsplit(cv.model$family[[1]],split="[()]")[[1]][2])
          cvprobs<-1-dnbinom(0,mu = cvpreds,size = cvtheta)
        }else{
          cvprobs<-1-dpois(0,cvpreds)
        }
      }else{
        cvpreds<-rep(NA,times=nrow(test.data))
        cvprobs<-rep(NA,times=nrow(test.data))
        cv.model<-NA
      }
    }

    # save a few things and get ready for the next cycle
    if(scale.preds){
      scale.factors[i]<-mean(train.data[,species])/mean(cvpreds)
      cvpreds<-cvpreds*scale.factors[i]

    }
    error.data$pred[start.vec[i]:end.vec[i]]<-preds
    error.data$prob[start.vec[i]:end.vec[i]]<-probs
    error.data$cvpred[start.vec[i]:end.vec[i]]<-cvpreds
    error.data$cvprob[start.vec[i]:end.vec[i]]<-cvprobs
    model.list[[i]]<-cv.model

    # need to remove a few old objects
    suppressWarnings(rm(cv.model,cv.preds))
    setTxtProgressBar(pb, i)
  }

  if(scale.preds){
    scale.factor<-mean(data[,species])/mean(preds)
    error.data$pred<-error.data$pred*scale.factor
  }

  close(pb)

  error.data$error<-error.data$pred-error.data$abund
  error.data$cverror<-error.data$cvpred-error.data$abund

  if(scale.preds){
    return(list(data=error.data,models=model.list,scale.factor=scale.factor,scale.factors=scale.factors))
  }else{
    return(list(data=error.data,models=model.list))
  }
}


#' Make cross-validation plots
#'
#' @description Use the model and the cross validated errors to make some residual plots
#' @details This is low priority, but could be much better.
#' @param error.data data frame containing observation, predictions, and CV predictions
#' @param method character; a method to be passed to the cor function
#' @param make.hist # should the histograms be plotted, they sometimes fail and may need to be turned off
#'
#' @return nothing, but creates some plots
#' @export
#'
#' @examples
MakeCrossValidationPlots<-function(error.data,           # a data frame, typically from the CrossvalidateModel function
                                   method="pearson",     # a method for the residuals, accepts pearson and spearman
                                   make.hist=T){         # should the histograms be plotted, they sometimes fail and may need to be turned off

  #remove any bad data
  keepers<-unique(which(is.infinite(error.data$pred)==F & is.infinite(error.data$cvpred)==F))
  if(length(keepers)<nrow(error.data)){
    warning(paste0(nrow(error.data)-length(keepers)," Infinite values detected and removed; estimates of error may be too low"))
  }

  error.data2<-error.data
  # compute the summary stats for the output
  if(method=="spearman"){
    error.data2$abund<-rank(error.data2$abund)
    error.data2$pred<-rank(error.data2$pred)
    error.data2$cvpred<-rank(error.data2$cvpred)
  }

  main.regr<-lm(error.data2$pred[keepers]~error.data2$abund[keepers])
  main.r2<-summary(main.regr)$r.squared
  main.rmse<-sqrt(sum((na.omit(error.data$abund[keepers]-error.data$pred[keepers]))^2)/nrow(na.omit(error.data[keepers,])))

  #now need to do the cv tests, which should already be in a nice format from the CV function
  cv.regr<-lm(error.data2$cvpred[keepers]~error.data2$abund[keepers])
  cv.r2<-summary(cv.regr)$r.squared
  cv.rmse<-sqrt(sum((na.omit(error.data$abund[keepers]-error.data$cvpred[keepers]))^2)/nrow(na.omit(error.data[keepers,])))

  print(paste("Full model",method,"Rsq =",round(main.r2,2)))
  print(paste("CV",method,"Rsq =",round(cv.r2,2)))
  print(paste("Full model RMSE =",round(main.rmse,3)))
  print(paste("CV RMSE =",round(cv.rmse,3)))

  # make all the plots
  old.par<-par()[c("mfcol","family","mar","xaxs","yaxs")]
  par(mfcol = c(ifelse(make.hist==T,3,2),2), family = "sans", mar = c(4,4,3,1))

  qqnorm((error.data$pred[keepers] - error.data$abund[keepers]), main = "Model Predictions")
  qqline((error.data$pred[keepers] - error.data$abund[keepers]))
  if(make.hist==T){hist((error.data$pred[keepers] - error.data$abund[keepers]), xlab = "Residuals", main = "")}
  pred.max <- ifelse(method=="pearson",quantile(error.data2$pred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  abund.max <- ifelse(method=="pearson",quantile(error.data2$pred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  plot.max<-max(pred.max,abund.max)*1.1

  plot(y=error.data2$pred[keepers], x=error.data2$abund[keepers], ylim = c(0,plot.max), xlim = c(0,plot.max),
       ylab = ifelse(method=="pearson","Predicted","Predicted Ranks"),
       xlab = ifelse(method=="pearson","Observed","Observed Ranks"),main = "", pch = 20)
  abline(coef = c(0,1), lty = 2)
  abline(main.regr,col=2)
  text(1, plot.max*.9, paste(method,"R-squared = ", signif(main.r2,2)), pos = 4)


  #Plots for test/CV data
  pred.max <- ifelse(method=="pearson",quantile(error.data2$cvpred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  abund.max <- ifelse(method=="pearson",quantile(error.data2$cvpred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  plot.max<-max(pred.max,abund.max)*1.1

  qqnorm((error.data$cvpred[keepers] - error.data$abund[keepers]), main = "Test data")
  qqline((error.data$cvpred[keepers] - error.data$abund[keepers]))
  if(make.hist==T){hist((error.data$cvpred[keepers] - error.data$abund[keepers]), xlab = "Residuals", main = "")}

  plot(y=error.data2$cvpred[keepers], x=error.data2$abund[keepers], ylim = c(0,plot.max), xlim = c(0,plot.max),
       ylab = ifelse(method=="pearson","Predicted","Predicted Ranks"),
       xlab = ifelse(method=="pearson","Observed","Observed Ranks"),main = "", pch = 20)
  abline(coef = c(0,1), lty = 2)
  abline(cv.regr,col=2)
  text(1, plot.max*.9, paste(method,"R-squared = ", signif(cv.r2,2)), pos = 4)
  suppressWarnings(par(old.par))
}


#' Make variance rasters
#' @description This function makes a non-parametric estimate of spatial variance based on the cross validation models. In order to do this, it needs to hold a lot of data in memory at once, so this can take awhile.
#' @param model.list list of models produced by the CV folds
#' @param raster.stack raster stack of the covariates used for the model
#' @param model.type character; the type of model ("maxnet","cloglog","hgam","gam")
#' @param scale.factor numeric; a scale factor to be applied
#' @param efh.break numeric; the EFH breakpoint for the full model, optionally creates an extra map
#'
#' @return raster of estimate non-parametric variance in model predictions
#' @export
#'
#' @examples
MakeVarianceRasters<-function(model.list,            # a list of models for each cv fold
                              raster.stack,          # a stack of covariates for the model
                              model.type,            # the type of model
                              scale.factor = 1,      # should the output be scaled
                              efh.break=NA){         # the efh breakpoint from the full model


  # to speed things up and preserve memory, we're going to separate out the NAs
  data.spots<-which(is.na(raster::getValues(raster.stack[[1]]))==F)
  for(r in 2:raster::nlayers(raster.stack)){
    spots<-which(is.na(raster::getValues(raster.stack[[r]]))==F)
    data.spots<-data.spots[data.spots%in%spots]
  }
  data<-extract(raster.stack,data.spots)

  #check for empty models
  list.index<-1
  model.list2<-list()
  for(m in 1:length(model.list)){
    if(is.na(model.list[[m]])==F){
      model.list2[[list.index]]<-model.list[[m]]
      list.index<-list.index+1
    }
  }

  out.data<-matrix(nrow=nrow(data),ncol = length(model.list2))

  # kind of complicated method of detecting an offset
  if(model.type!="maxnet"){
    terms<-AutodetectGAMTerms(model.list[[1]],hgam = "d")
    has.offset<-"offset"%in%terms$type
    facs<-terms$term[terms$type=="factor"]
    offset.name<-terms$term[terms$type=="offset"]

    # check for an offset, kind of kludgy solution for now
    if(model.type=="hgam" & has.offset){
      data<-data.frame(data,mean(model.list2[[1]]$offset[[1]]))
      pterms<-AutodetectGAMTerms(model.list[[1]],hgam = "p")
      facs<-unique(c(facs,pterms$term[pterms$type=="factor"]))
    }
    if(model.type%in%c("cloglog","gam")){
      data<-data.frame(data,mean(model.list2[[1]]$offset))
    }
    names(data)[ncol(data)]<-offset.name

    for(f in 1:length(facs)){
      data[,facs[f]]<-as.factor(data[,facs[f]])
    }
  }


  #progress bar
  pb <- txtProgressBar(min = 0, max = length(model.list2), style = 3)

  # loop through and get predictions for each model
  for(m in 1:length(model.list2)){
    if(model.type=="maxnet"){
      out.data[,m]<-exp(predict(model.list2[[m]],newdata=data,type="link")+model.list2[[m]]$entropy)
    }
    if(model.type=="cloglog"){
      out.data[,m]<-exp(mgcv::predict.gam(object = model.list2[[m]],newdata=data,type="link"))
    }
    if(model.type=="hgam"){
      out.data[,m]<-mgcv::predict.gam(object = model.list2[[m]],newdata=data,type="response")
    }
    if(model.type=="gam"){
      out.data[,m]<-mgcv::predict.gam(object = model.list2[[m]],newdata=data,type="response")
    }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  # now tally things up

  raster.template<-raster::raster(raster.stack)

  variances<-apply(X = out.data*scale.factor,MARGIN = 1,FUN = var)
  var.vec<-rep(NA,times=raster::ncell(raster.stack))
  var.vec[data.spots]<-variances
  var.raster<-raster::setValues(raster.template,values = var.vec)

  if(is.na(efh.break)==F){
    percents<-apply(X = out.data>efh.break,MARGIN = 1,FUN = sum)/ncol(out.data)
    per.vec<-rep(NA,times=raster::ncell(raster.stack))
    per.vec[data.spots]<-percents
    per.raster<-raster::setValues(raster.template,values = per.vec)
    return(list(var.raster,per.raster))
  }else{
    return(var.raster)
  }
}


#' Plot effects
#'
#' @description Plot effects from an effects list
#' @details If the cv.models are provided, it will calculate confidence intervals based on that.
#' @param effects list of data frames with the effects, such as from the GetGAMEffects function
#' @param nice.names data frame linking abbreviated names to nicer versions to be used in plottingv(optional)
#' @param vars vector of which variables to plot, or "all"
#' @param land shapefile used to insert landmasses into map of lat/lon effect
#'
#' @return nothing, but creates plots
#' @export
#'
#' @examples
plotEffects<-function(effects,                  #
                      nice.names=NULL,          #
                      vars="all",               #
                      land=NULL){

  # check the variable names and restrict things to those requested
  if(vars!="all" & is.character(vars)){
    vars2<-which(names(effects)%in%vars)
  }else{
    if(vars=="all"){
      vars2<-1:length(effects)
    }else{
      vars2<-vars
    }
  }

  if(length(vars2)<length(vars)){
    missing<-vars[which(vars%in%names(effects)==F)]
    stop("Variables [",paste0(missing,collapse = ", "),"] not found in effects list")
  }

  # if names aren't supplied, use the terms from the model
  if(is.null(nice.names)){
    nice.names<-data.frame(var=unlist(strsplit(names(effects),"[*]")),
                           name=unlist(strsplit(names(effects),"[*]")))
  }

  # set up parameters for the plots
  oldpar<-par()
  n.col<-ifelse(length(vars2)>4,3,ifelse(length(vars2)>1,2,1))
  n.rows<-ceiling(length(vars2)/n.col)
  par(mfrow = c(n.rows, n.col), mar = c(4,4,1,0.01),oma=c(.5,.5,.5,.5))

  for(i in 1:length(vars2)){
    var.dat<-effects[[vars2[i]]]

    # 2d variables
    if("x"%in%names(var.dat) & "y"%in%names(var.dat)){

      x.name<-strsplit(names(effects)[vars2[i]],split="[*]")[[1]][1]
      y.name<-strsplit(names(effects)[vars2[i]],split="[*]")[[1]][2]

      x.lab<-nice.names$name[nice.names$var==x.name]
      y.lab<-nice.names$name[nice.names$var==y.name]

      xseq<-unique(var.dat$x)
      yseq<-unique(var.dat$y)

      # There is some special handling for certain types of plots
      if(x.lab%in%c("Longitude","longitude","lon","long")){
        plot(NA,xlab=x.lab,ylab=y.lab,xlim=c(min(xseq),max(xseq)),ylim=c(min(yseq),max(yseq)),
             xaxt = "n",yaxt = "n")
        if(is.null(land)==F){
          suppressWarnings(plot(land, col = "grey90", add = TRUE,legend=F,border="grey80"))
        }
        box(which="plot")
        contour(x = xseq,y=yseq,z=matrix(nrow=length(yseq),ncol=length(xseq),data = matrix(var.dat$effect),byrow=F),
                add=T,nlevels=10)
      }else{
        plot(NA,xlab=x.lab,ylab=y.lab,xlim=c(min(xseq),max(xseq)),ylim=c(min(yseq),max(yseq)))
        contour(x = xseq,y=yseq,z=matrix(data=var.dat$effect,byrow=T,nrow=length(xseq),ncol=length(yseq)),add=T,nlevels=10)
      }

      if(x.lab%in%c("Current Velocity East (m/s)","bcurrentU")){
        abline(h=0,v=0,lty=3)
      }
      if(x.lab%in%c("Current Velocity East SD (m/s)","bcurrentUSD")){
        abline(a=0,b=1,lty=3)
      }
    }

    # plots for 1 d variables
    if("x"%in%names(var.dat) & "y"%in%names(var.dat)==F & length(var.dat$effect)>10){

      y.lab<-ifelse(i%%n.col==1| n.col==1,"Variable effect" ,NA)
      x.lab<-nice.names$name[nice.names$var==names(effects)[vars2[i]]]

      if("upper"%in%names(var.dat)){
        y.lim<-c(max(min(var.dat$lower[is.finite(var.dat$lower)]),min(var.dat$effect-5,na.rm = T)),
                 min(max(var.dat$upper[is.finite(var.dat$upper)]),max(var.dat$effect+5,na.rm = T)))
        plot(x=var.dat$x,y=var.dat$effect,type="l",ylab = y.lab, xlab =x.lab,ylim=y.lim)
        lines(x=var.dat$x,y=var.dat$upper,lty=3)
        lines(x=var.dat$x,y=var.dat$lower,lty=3)
      }else{
        y.lim<-c(min(var.dat$effect[is.finite(var.dat$effect)])-.2,max(var.dat$effect[is.finite(var.dat$effect)])+.2)
        plot(x=var.dat$x,y=var.dat$effect,type="l",ylab = y.lab, xlab =x.lab,ylim=y.lim)
      }
    }

    # plots for factors
    if("x"%in%names(var.dat) & "y"%in%names(var.dat)==F & length(var.dat$effect)<=10){

      y.lab<-ifelse(i%%n.col==1| n.col==1,"Variable effect" ,NA)
      x.lab<-nice.names$name[nice.names$var==names(effects)[vars2[i]]]

      if("upper"%in%names(var.dat)){
        y.lim<-c(max(min(var.dat$lower[is.finite(var.dat$lower)]),min(var.dat$effect-5,na.rm = T)),
                 min(max(var.dat$upper[is.finite(var.dat$upper)]),max(var.dat$effect+5,na.rm = T)))
        plot(x=as.factor(var.dat$x),y=var.dat$effect,type="l",ylab = y.lab, xlab =x.lab,ylim=y.lim)
        plot(x=as.factor(var.dat$x),y=var.dat$upper,lty=3,lwd=.4,add=T)
        plot(x=as.factor(var.dat$x),y=var.dat$lower,lty=3,lwd=.4,add=T)
      }else{
        y.lim<-c(min(var.dat$effect[is.finite(var.dat$effect)])-.2,max(var.dat$effect[is.finite(var.dat$effect)])+.2)
        plot(x=as.factor(var.dat$x),y=var.dat$effect,type="l",ylab = y.lab, xlab =x.lab,ylim=y.lim)
      }
    }
  }
  suppressWarnings(par(oldpar))
}




