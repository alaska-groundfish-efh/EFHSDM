# this script will make some more generalized functions for mapping things with the akgfmaps package

rpackages <- c("viridis", "stars","sf","gridExtra","patchwork","MASS","scales","labeling","ggplot2")

which_not_installed <- which(rpackages %in% rownames(installed.packages()) == FALSE)

if(length(which_not_installed) > 1){
  install.packages(rpackages[which_not_installed], dep = TRUE)
}

if("akgfmaps" %in% rownames(installed.packages()) == F){
  devtools::install_github("sean-rohan-noaa/akgfmaps", build_vignettes = TRUE)
}
rm(rpackages,which_not_installed)

require(ggplot2)
require(akgfmaps)
require(viridis)
require(stars)
require(sf)
require(gridExtra)
require(patchwork)
require(MASS)
require(scales)
require(labeling)
# 
# example.raster<-raster("Y:/RACE_EFH_variables/Trawl_Models/AI/Adult_arrowtooth_flounder/ensemble_abundance")
# dataCRS<-example.raster@crs
# 
# ai.catch<-read.csv("Y:/RACE_EFH_variables/Trawl_Models/AI/all_AI_data_2021.csv")
# 
# presence<-ai.catch[ai.catch$a_rex>0,c("lon","lat")]
# absence<-ai.catch[ai.catch$a_rex==0,c("lon","lat")]
# highdensity<-ai.catch[ai.catch$a_rex>=quantile(ai.catch$a_rex[ai.catch$a_rex>0],.9),c("lon","lat")]
# 

# makes a pretty good dotplot; needs additional testing with areas other than "bs.all", "goa", and "ai".
# also, settings are set to work with 8x8 300 res output, so results may vary when making different size figures
#' Title
#'
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
#'
#' @return a ggplot with the desired figure
#' @export
#'
#' @examples
MakeAKGFDotplot<-function(presence,                  
                          absence=NA,                
                          highdensity=NA,            
                          dataCRS=NA,                # 
                          region,
                          survey.area="default",
                          ext.adjust="default",      # 
                          legend.pos="default",      # 
                          title.name=NA,
                          title.count=T,
                          title.pos="default",
                          abs.name="absent",
                          pres.name="present",
                          hd.name="top 10%",
                          abs.col="steelblue4",
                          pres.col="orange",
                          hd.col="red",
                          abs.shape=16,
                          pres.shape=1,
                          hd.shape=16,
                          abs.size=.1,
                          pres.size=1,
                          hd.size=1){
  # start by getting the map
  region<-tolower(region)
  if(region%in%c("ebs","bs.all","sebs","bs.south","ecs","ebs.ecs","ai",
                 "ai.west","ai.central","ai.east","goa","goa.west","goa.east")){
    MAP<-get_base_layers(select.region = region,set.crs = "auto")
  }else{
    stop("region not recognized")
  }
  
  #set up some essentials for later
  if(ext.adjust=="default"){
    if(region%in%c("ai","ai.west","ai.central","ai.east")){ext.adjust<-c(0,0,0,0)}
    if(region%in%c("ebs","bs.all","bs.south","sebs","ecs","ebs.ecs")){ext.adjust<-c(0,0,0,0)}
    if(region%in%c("goa","goa.west","goa.east")){ext.adjust<-c(200000,-100000,0,0)}
  }
  if(is.na(legend.pos)==F && legend.pos=="default"){
    if(region%in%c("ai","ai.west","ai.central","ai.east")){legend.pos<-c(0.12,0.20)}
    if(region%in%c("ebs","bs.all","bs.south","sebs","ecs","ebs.ecs")){legend.pos<-c(0.12,0.18)}
    if(region%in%c("goa","goa.west","goa.east")){legend.pos<-c(0.08,0.41)}
  }
  if(is.na(title.pos)==F && title.pos=="default"){
    if(region%in%c("ai","ai.west","ai.central","ai.east")){title.pos<-c(-900000,490000)}
    if(region%in%c("ebs","bs.all","bs.south","sebs","ecs","ebs.ecs")){title.pos<-c(-550000,800000)}
    if(region%in%c("goa","goa.west","goa.east")){title.pos<-c(-1300000,630000)}
  }
  
  # detect and set up the outline, for now it is required
  if(is.character(survey.area) && survey.area[1]=="default"){
    survey.sf<-MAP$survey.area
  }else{
    if(class(survey.area)[1]=="sf"){
      survey.sf<-st_transform(survey.area,st_crs(MAP$akland))
    }
    if(class(survey.area)[1]=="RasterLayer"){
      survey.sf0 <- st_as_stars(is.na(survey.area)) %>%
        st_as_sf(merge = TRUE) %>% # this is the raster to polygons part
        st_cast("POLYGON") # cast the polygons to polylines
      
      survey.sf<-st_transform(survey.sf0,st_crs(MAP$akland))[1:(nrow(survey.sf0)-1),]
      
    }
  }
  
  leg.col<-NULL
  leg.shape<-NULL
  leg.size<-NULL
  leg.name<-NULL
  
  # Now go through and set up the dot locations and add them to legend
  if(is.na(absence)==F){
    if(is.na(dataCRS)){
      abs.sf<-st_as_sf(x=absence,coords = c("lon","lat"),crs=st_crs(MAP$akland))
    }else{
      abs.sf0<-st_as_sf(x=absence,coords = c("lon","lat"),crs=dataCRS)
      abs.sf<-st_transform(abs.sf0,st_crs(MAP$akland))
    }
    leg.name<-abs.name
    leg.col<-abs.col
    leg.shape<-abs.shape
    leg.size<-abs.size
    abs.fac<-1
  }else{
    abs.fac<-0
  }
  
  if(is.na(dataCRS)){
    pres.sf<-st_as_sf(x=presence,coords = c("lon","lat"),crs=st_crs(MAP$akland))
  }else{
    pres.sf0<-st_as_sf(x=presence,coords = c("lon","lat"),crs=dataCRS)
    pres.sf<-st_transform(pres.sf0,st_crs(MAP$akland))
  }
  leg.name<-c(leg.name,pres.name)
  leg.col<-c(leg.col,pres.col)
  leg.shape<-c(leg.shape,pres.shape)
  leg.size<-c(leg.size,pres.size)
  pres.fac<-abs.fac+1
  
  if(is.na(highdensity)==F){
    if(is.na(dataCRS)){
      high.sf<-st_as_sf(x=highdensity,coords = c("lon","lat"),crs=st_crs(MAP$akland))
    }else{
      high.sf0<-st_as_sf(x=highdensity,coords = c("lon","lat"),crs=dataCRS)
      high.sf<-st_transform(high.sf0,st_crs(MAP$akland))
    }
    leg.name<-c(leg.name,hd.name)
    leg.col<-c(leg.col,hd.col)
    leg.shape<-c(leg.shape,hd.shape)
    leg.size<-c(leg.size,hd.size)
    hd.fac<-pres.fac+1
  }
  
  # set up the basic map, will add more customization later
  dotplot<-ggplot()+
    geom_sf(data=survey.sf,fill="grey95")+
    geom_sf(data=MAP$akland,fill="grey40")+
    geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) +
    geom_sf(data = MAP$bathymetry,color="grey60")
  
  # add the dots
  if(is.na(absence)==F){dotplot<-dotplot+geom_sf(data=abs.sf,alpha=.25,size=abs.size,shape=abs.shape,aes(color=factor(abs.fac)))}
  dotplot<-dotplot+geom_sf(data=pres.sf,size=pres.size,aes(color=factor(pres.fac)),shape=pres.shape,stroke=.8)
  if(is.na(highdensity)==F){dotplot<-dotplot+geom_sf(data=high.sf,size=hd.size,shape=hd.shape,aes(color=factor(hd.fac)))}  
  
  # add the themes
  dotplot<-dotplot+
    coord_sf(xlim = MAP$plot.boundary$x+ext.adjust[1:2], ylim = MAP$plot.boundary$y+ext.adjust[3:4])+
    scale_x_continuous(name = "Longitude",breaks = MAP$lon.breaks) +
    scale_y_continuous(name = "Latitude",breaks = MAP$lat.breaks) +
    theme_bw()+
    theme(panel.border = element_rect(color = "black",fill = NA),
          panel.background = element_rect(fill = NA,color = "black"),
          legend.key = element_rect(fill = NA,color = "grey30"),
          legend.position = legend.pos,
          axis.title = element_blank(), axis.text = element_text(size = 12),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12),
          plot.background = element_rect(fill = NA, color = NA))
  
  # add a title
  if(is.na(title.name)==F && is.na(title.pos)==F){
    if(title.count){
      title.name<-paste0(title.name,"\nN = ",format(nrow(pres.sf),big.mark = ","))
    }
    dotplot<-dotplot+
      geom_label(data=data.frame(x=title.pos[1],y=title.pos[2],label=title.name),
                 aes(x=x,y=y,label=label,hjust=0,vjust=1),size=5)
  }
  
  # add a legend
  if(is.na(legend.pos)==F){
    dotplot<-dotplot+
      scale_color_manual(name=NULL,values=leg.col,labels=leg.name)+
      guides(color=guide_legend(override.aes = list(shape=leg.shape,size=leg.size)))
  }
  return(dotplot)
}


# Makes a pretty good map of any continuous variable

#' Title
#'
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
#'
#' @return a ggplot object of the map
#' @export
#'
#' @examples
MakeAKGFDensityplot<-function(region,
                              density.map,
                              buffer=1,
                              survey.area="default",
                              ext.adjust="default",      # vector of length four, representing xmin,xmax, ymin,ymax, to adjust the extent
                              legend.pos="default",      # a vector of length two with a location for the legend, use NA to suppress
                              legend.title=NA,
                              col.palette="plasma",
                              col.palette.limits=c(0,1),
                              title.name=NA,
                              title.pos="default"){
  
  # start by getting the map
  region<-tolower(region)
  if(region%in%c("ebs","bs.all","sebs","bs.south","ecs","ebs.ecs","ai",
                 "ai.west","ai.central","ai.east","goa","goa.west","goa.east")){
    MAP<-get_base_layers(select.region = region,set.crs = "auto")
  }else{
    stop("region not recognized")
  }
  
  #set up some essentials for later
  if(ext.adjust=="default"){
    if(region%in%c("ai","ai.west","ai.central","ai.east")){ext.adjust<-c(0,0,0,0)}
    if(region%in%c("ebs","bs.all","bs.south","sebs","ecs","ebs.ecs")){ext.adjust<-c(0,0,0,0)}
    if(region%in%c("goa","goa.west","goa.east")){ext.adjust<-c(200000,-100000,0,0)}
  }
  if(is.na(legend.pos)==F && legend.pos=="default"){
    if(region%in%c("ai","ai.west","ai.central","ai.east")){legend.pos<-c(.07,.28)}
    if(region%in%c("ebs","bs.all","bs.south","sebs","ecs","ebs.ecs")){legend.pos<-c(0.12,0.18)}
    if(region%in%c("goa","goa.west","goa.east")){legend.pos<-c(0.08,0.51)}
  }
  if(is.na(title.pos)==F && title.pos=="default"){
    if(region%in%c("ai","ai.west","ai.central","ai.east")){title.pos<-c(-900000,490000)}
    if(region%in%c("ebs","bs.all","bs.south","sebs","ecs","ebs.ecs")){title.pos<-c(-550000,800000)}
    if(region%in%c("goa","goa.west","goa.east")){title.pos<-c(-1300000,630000)}
  }
  
  # detect and set up the outline, for now it is required
  if(is.character(survey.area) && survey.area[1]=="default"){
    survey.sf<-MAP$survey.area
  }else{
    if(class(survey.area)[1]=="sf"){
      survey.sf<-st_transform(survey.area,st_crs(MAP$akland))[1:(nrow(survey.area)-1),]
    }
    if(class(survey.area)[1]=="RasterLayer"){
      survey.sf0 <- st_as_stars(is.na(survey.area)) %>%
        st_as_sf(merge = TRUE) %>% # this is the raster to polygons part
        st_cast("POLYGON") # cast the polygons to polylines
      
      survey.sf<-st_transform(survey.sf0,st_crs(MAP$akland))[1:(nrow(survey.sf0)-1),]
      
    }
  }
  
  #set up the density component
  if(class(density.map)[1]=="sf"){
    density.sf<-st_transform(density.map,st_crs(MAP$akland))
    names(density.sf[1])<-"density"
  }
  if(class(density.map)[1]=="RasterLayer"){
    density.dat0<-as.data.frame(xyFromCell(density.map,1:ncell(density.map)))
    vals<-getValues(density.map)
    
    # convert the raster to ggplot format
    density.dat<-data.frame(lat=density.dat0$y,lon=density.dat0$x,density=vals)
    
    density.sf0<-st_as_sf(x=subset(density.dat,is.na(density)==F),coords = c("lon","lat"),crs=density.map@crs)
    density.sf<-st_transform(density.sf0,st_crs(MAP$akland))
  }
  
  # often helps to remove some of the high points
  density.sf$density[density.sf$abund>quantile(density.sf$density,buffer,na.rm=T)]<-quantile(density.sf$density,buffer,na.rm=T)
  
  # set up the basic map, will add more customization later
  densityplot<-ggplot()+
    geom_sf(data=survey.sf,fill="grey95")+
    geom_sf(data=MAP$akland,fill="grey40")+
    geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) +
    geom_sf(data = MAP$bathymetry,color="grey60")
  
  # add the density plot
  densityplot<-densityplot+
    geom_sf(data=density.sf,aes(col=density),size=.05)
  
  # add the themes
  densityplot<-densityplot+
    coord_sf(xlim = MAP$plot.boundary$x+ext.adjust[1:2], ylim = MAP$plot.boundary$y+ext.adjust[3:4])+
    scale_x_continuous(name = "Longitude",breaks = MAP$lon.breaks) +
    scale_y_continuous(name = "Latitude",breaks = MAP$lat.breaks) +
    scale_color_viridis(option=col.palette,begin = col.palette.limits[1],end = col.palette.limits[2],
                        na.value=NA,name=legend.title,labels=comma_format(big.mark = ","))+
    theme_bw()+
    theme(panel.border = element_rect(color = "black",fill = NA),
          panel.background = element_rect(fill = NA,color = "black"),
          legend.key = element_rect(fill = NA,color = "grey30"),
          legend.position = legend.pos,
          axis.title = element_blank(), axis.text = element_text(size = 12),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12),
          plot.background = element_rect(fill = NA, color = NA))
  
  # add a title
  if(is.na(title.name)==F && is.na(title.pos)==F){
    densityplot<-densityplot+
      geom_label(data=data.frame(x=title.pos[1],y=title.pos[2],label=title.name),
                 aes(x=x,y=y,label=label,hjust=0,vjust=1),size=5)
  }
  # add a legend
  if(is.na(legend.pos[1])==F){
    densityplot<-densityplot+
      guides(color=guide_colorbar(title.position = "top",title.hjust = .5,barheight = 5))
  }
  return(densityplot)
}


#
#' Title
#'
#' @param region character; which region is being used, must correpond to an akgfmaps baselayer
#' @param efh.map raster of factors to be interpretted as EFH
#' @drop a buffer that prevents outlining of small isolated pixels, makes it look nicer
#' @param survey.area sf or raster layer to use, if the desired survey area is different from akgfmaps
#' @param ext.adjust vector of length four, representing xmin,xmax, ymin,ymax, to adjust the extent
#' @param legend.pos vector of length two with decimal position of legend box
#' @param legend.labels vector of labels for the legend, should match number of factors
#' @param legend.title character; a title for the legend
#' @param col.palette character a color palette for 
#' @param col.palette.limits vector of length two giving start and end points for the color palette
#' @param title.name character; a title for the figure, or NA to suppress
#' @param title.pos vector of length two with coordinates for for the legend, use NA to suppress
#'
#' @return raster of EFH quantiles (subareas)
#' @export
#'
#' @examples
MakeAKGFEFHplot<-function(region,
                         efh.map,
                         drop=10^8,
                         survey.area="default",
                         ext.adjust="default",      # vector of length four, representing xmin,xmax, ymin,ymax, to adjust the extent
                         legend.pos="default",      # a vector of length two with a location for the legend, use NA to suppress
                         legend.labels=c("95%","75%","50%","25%"),
                         legend.title=NA,
                         col.palette="viridis",
                         col.palette.limits=c(0,1),
                         title.name=NA,
                         title.pos="default"){
  
  # start by getting the map
  region<-tolower(region)
  if(region%in%c("ebs","bs.all","sebs","bs.south","ecs","ebs.ecs","ai",
                 "ai.west","ai.central","ai.east","goa","goa.west","goa.east")){
    MAP<-get_base_layers(select.region = region,set.crs = "auto")
  }else{
    stop("region not recognized")
  }
  
  #set up some essentials for later
  if(ext.adjust=="default"){
    if(region%in%c("ai","ai.west","ai.central","ai.east")){ext.adjust<-c(0,0,0,0)}
    if(region%in%c("ebs","bs.all","bs.south","sebs","ecs","ebs.ecs")){ext.adjust<-c(0,0,0,0)}
    if(region%in%c("goa","goa.west","goa.east")){ext.adjust<-c(200000,-100000,0,0)}
  }
  if(is.na(legend.pos)==F && legend.pos=="default"){
    if(region%in%c("ai","ai.west","ai.central","ai.east")){legend.pos<-c(.07,.28)}
    if(region%in%c("ebs","bs.all","bs.south","sebs","ecs","ebs.ecs")){legend.pos<-c(0.12,0.18)}
    if(region%in%c("goa","goa.west","goa.east")){legend.pos<-c(0.08,0.51)}
  }
  if(is.na(title.pos)==F && title.pos=="default"){
    if(region%in%c("ai","ai.west","ai.central","ai.east")){title.pos<-c(-900000,490000)}
    if(region%in%c("ebs","bs.all","bs.south","sebs","ecs","ebs.ecs")){title.pos<-c(-550000,800000)}
    if(region%in%c("goa","goa.west","goa.east")){title.pos<-c(-1300000,630000)}
  }
  
  # detect and set up the outline, for now it is required
  if(is.character(survey.area) && survey.area[1]=="default"){
    survey.sf<-MAP$survey.area
  }else{
    if(class(survey.area)[1]=="sf"){
      survey.sf<-st_transform(survey.area,st_crs(MAP$akland))[1:(nrow(survey.area)-1),]
    }
    if(class(survey.area)[1]=="RasterLayer"){
      survey.sf0 <- st_as_stars(is.na(survey.area)) %>%
        st_as_sf(merge = TRUE) %>% # this is the raster to polygons part
        st_cast("POLYGON") # cast the polygons to polylines
      
      survey.sf<-st_transform(survey.sf0,st_crs(MAP$akland))[1:(nrow(survey.sf0)-1),]
      
    }
  }
  
  # set up the factor maps
  efh.vals<-getValues(efh.map)
  efh.vals[efh.vals==1]<-NA
  
  # convert the raster to polygons
  efhpoly <- st_as_stars(efh.map) %>%
    st_as_sf(merge = TRUE) 
  efhpoly2<-efhpoly[efhpoly$layer!=1,]
  
  #we'll need a new outline
  efh.dummy.raster<-raster(efh.map)
  efh.vals2<-is.na(efh.vals)==F
  efh.dummy.raster<-setValues(efh.dummy.raster,values = efh.vals2)
  
  efhdummy <- st_as_stars(efh.dummy.raster) %>%
    st_as_sf(merge = TRUE) %>% # this is the raster to polygons part
    st_cast("MULTILINESTRING") # cast the polygons to polylines
  
  efhdummy2<-st_transform(efhdummy,st_crs(MAP$akland))
  
  # Now we need to get rid of a lot of the tiniest bits, which we'll do by dropping the smallest areas
  efhdummy.poly<-st_cast(efhdummy2,"POLYGON")
  areas<-st_area(efhdummy.poly)
  
  outside<-order(areas,decreasing = T)[1]
  toosmall<-which(as.numeric(areas)<drop)
  
  efh.x<-efhdummy2$layer[-c(outside,toosmall)]
  efh.y<-efhdummy2$geometry[-c(outside,toosmall)]
  efhdummy3<-st_sf(efh.x,efh.y)
  
  # set up the basic map, will add more customization later
  efhplot<-ggplot()+
    geom_sf(data=survey.sf,fill="grey95")+
    geom_sf(data=MAP$akland,fill="grey40")+
    geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) +
    geom_sf(data = MAP$bathymetry,color="grey60")
  
  # add the efh polys
  efhplot<-efhplot+
    geom_sf(data=efhpoly2,aes(fill=as.factor(layer)),col=NA)+
    geom_sf(data=efhdummy3,size=.3)
  
  # add the themes
  efhplot<-efhplot+
    coord_sf(xlim = MAP$plot.boundary$x+ext.adjust[1:2], ylim = MAP$plot.boundary$y+ext.adjust[3:4])+
    scale_x_continuous(name = "Longitude",breaks = MAP$lon.breaks) +
    scale_y_continuous(name = "Latitude",breaks = MAP$lat.breaks) +
    scale_fill_viridis(discrete=T,name=legend.title,labels=legend.labels)+
    theme_bw()+
    theme(panel.border = element_rect(color = "black",fill = NA),
          panel.background = element_rect(fill = NA,color = "black"),
          legend.key = element_rect(fill = NA,color = "grey30"),
          legend.position = legend.pos,
          axis.title = element_blank(), axis.text = element_text(size = 12),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12),
          plot.background = element_rect(fill = NA, color = NA))
  
  # add a title
  if(is.na(title.name)==F && is.na(title.pos)==F){
    efhplot<-efhplot+
      geom_label(data=data.frame(x=title.pos[1],y=title.pos[2],label=title.name),
                 aes(x=x,y=y,label=label,hjust=0,vjust=1),size=5)
  }
  
  # add a legend
  if(is.na(legend.pos[1])==F){
    efhplot<-efhplot+
      guides(color=guide_legend(override.aes=list(size=4,shape=15)))
  }
  return(efhplot)
}

#######################################################################
# this function makes a canned version of the EFH comparison plots. Could use some more work, and is a little less
# flexible than other functions at current.
#' Title
#'
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
#'
#' @return ggplot object with the comparison map
#' @export
#'
#' @examples
PlotEFHComparison<-function(old=NA,new=NA,main="",background,leg.name=NULL,leg.labels,region,col.vec=NA,
                            label.pos=NA,leg.pos=NA,ext.adjust.x=NA,ext.adjust.y=NA,nonEFH=1){
  
  # start by getting the map
  region<-tolower(region)
  if(region%in%c("ebs","bs.all","sebs","bs.south","ecs","ebs.ecs","ai",
                 "ai.west","ai.central","ai.east","goa","goa.west","goa.east")){
    MAP<-get_base_layers(select.region = region,set.crs = "auto")
  }else{
    stop("region not recognized")
  }
  
  # set up some plot parameters
  if(region%in%c("bs.all","ebs")){
    if(is.na(label.pos)){label.pos<-c(-550000,810000)}
    if(is.na(leg.pos)){leg.pos<-c(0.03,0.21)}
    if(is.na(ext.adjust.x)){ext.adjust.x<-c(0,0)}
    if(is.na(ext.adjust.y)){ext.adjust.y<-c(0,0)}
  }
  if(region=="goa"){
    if(is.na(label.pos)){label.pos<-c(-1350000,660000)}
    if(is.na(leg.pos)){leg.pos<-c(0.01,0.7)}
    if(is.na(ext.adjust.x)){ext.adjust.x<-c(200000,-100000)}
    if(is.na(ext.adjust.y)){ext.adjust.y<-c(0,0)}
  }
  if(region=="ai"){
    if(is.na(label.pos)){label.pos<-c(-900000,490000)}
    if(is.na(leg.pos)){leg.pos<-c(0.01,.38)}
    if(is.na(ext.adjust.x)){ext.adjust.x<-c(0,0)}
    if(is.na(ext.adjust.y)){ext.adjust.y<-c(0,0)}
  }
  
  if(is.na(col.vec)){col.vec<-c("dodgerblue","orange","orchid3")}
  
  dummy.sf <- st_as_stars(is.na(background)==F) %>%
    st_as_sf(merge = TRUE) %>% # this is the raster to polygons part
    st_cast("POLYGON") # cast the polygons to polylines
  dummy.sf2<-st_transform(dummy.sf,st_crs(MAP$akland))
  dummy.sf3<-dummy.sf2[dummy.sf2$layer==1,]
  
  
  old.col=NULL
  new.col=NULL
  mix.col=NULL
  
  # Check which rasters are available
  try(old.present<-is.na(old@crs)==F)
  try(new.present<-is.na(new@crs)==F)
  
  if(exists("old.present") & exists("new.present")){
    comp.raster<-cut(background,breaks=c(-Inf,Inf))
    
    # find which areas are EFH in each version
    old2<-cut(x = old,breaks=c(0,nonEFH+.5,Inf))
    new2<-cut(x = new,breaks=c(0,nonEFH+.5,Inf))
    
    vals<-getValues(comp.raster)
    
    both<-which(getValues(old2)==2 & getValues(new2)==2)
    justold<-which(getValues(old2)==2 &
                     ((getValues(new2)==1)| is.na(getValues(new2))))
    justnew<-which(((getValues(old2)==1)| is.na(getValues(old2))) &
                     getValues(new2)==2)
    
    #assign new values to each category
    vals[justold]<-2
    vals[justnew]<-3
    vals[both]<-4
    
    comp.raster<-setValues(x = comp.raster,values = vals)
    
    efhpoly <- st_as_stars(comp.raster) %>%
      st_as_sf(merge = TRUE) 
    efhpoly2<-efhpoly[efhpoly$layer>1,]
    
    # this is a kludge to make sure there are always all 3 types
    layer.present<-unique(efhpoly2$layer)
    n=0
    if(2%in%layer.present==F){
      efhpoly2$layer[which.min(st_area(efhpoly2))]<-2
      n<-1
    }
    if(3%in%layer.present==F){
      efhpoly2$layer[order(st_area(efhpoly2))[1+n]]<-3
      n<-n+1
    }
    if(4%in%layer.present==F){
      efhpoly2$layer[order(st_area(efhpoly2))[1+n]]<-4
    }
    
    area1<-sum(getValues(old)>nonEFH,na.rm=T)
    area2<-sum(getValues(new)>nonEFH,na.rm=T)
    main<-paste0(main,"\nChange = ",round((area2/area1-1)*100,1),"%")
    
    old.col<-col.vec[1]
    new.col<-col.vec[2]
    mix.col<-col.vec[3]
  }else{
    if(exists("old.present")){
      efhpoly <- st_as_stars(old>nonEFH) %>%
        st_as_sf(merge = TRUE) 
      efhpoly2<-efhpoly[efhpoly$layer!=0,]
      old.col<-col.vec[1]
    }
    if(exists("new.present")){
      efhpoly <- st_as_stars(new>nonEFH) %>%
        st_as_sf(merge = TRUE) 
      efhpoly2<-efhpoly[efhpoly$layer!=0,]
      new.col<-col.vec[2]
    }
  }
  
  out.plot<-ggplot()+
    geom_sf(data=dummy.sf3,fill="grey95")+
    geom_sf(data=MAP$survey.area,fill="grey95")+
    geom_sf(data=efhpoly2,aes(fill=factor(layer)),col=1,size=.05)+
    geom_sf(data=MAP$akland,fill="grey40")+
    geom_sf(data = MAP$graticule, color = "grey70", alpha = 0.5) +
    geom_sf(data = MAP$bathymetry,col="grey60",size=.25) +
    coord_sf(xlim = MAP$plot.boundary$x+ext.adjust.x, ylim = MAP$plot.boundary$y+ext.adjust.y) +
    geom_label(data=data.frame(x=label.pos[1],y=label.pos[2],label=main),
               aes(x=x,y=y,label=main,hjust=0,vjust=1),size=5)+
    scale_fill_manual(values=c(old.col,new.col,mix.col),labels=leg.labels,name=leg.name)+
    scale_x_continuous(name = "Longitude",breaks = MAP$lon.breaks) +
    scale_y_continuous(name = "Latitude",breaks = MAP$lat.breaks) +
    theme_bw()+
    theme(panel.border = element_rect(color = "black",fill = NA),
          panel.background = element_rect(fill = NA,color = "black"),
          legend.key = element_rect(fill = NA,color = NA),
          legend.position = leg.pos,legend.justification = c(0,1),
          axis.title = element_blank(), axis.text = element_text(size = 12),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12),
          plot.background = element_rect(fill = NA, color = NA))+
    guides(color=guide_legend(override.aes=list(size=4,shape=15)))
  
  return(out.plot)
}
