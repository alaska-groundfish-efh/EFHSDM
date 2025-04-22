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
    sample <- stats::na.omit(terra::values(abund.raster))
    sample[sample <= threshold] <- NA
    breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
    breaks[1]<-0
    breaks[length(breaks)]<-Inf
  }
  if(method=="cumulative"){
    # Decide whether to sample at given locations or use the whole thing
    if(is.null(data)){
      vals<-stats::na.omit(sort(terra::values(abund.raster)))
      vals2<-cumsum(vals)/sum(vals)
    }else{
      vals<-stats::na.omit(sort(terra::extract(abund.raster,data.frame(data$lon,data$lat))))
      vals2<-cumsum(vals)/sum(vals)
    }

    # Loop to calculate the breaks
    breaks<-c(0,rep(NA,length(quants)-2),Inf)
    while(length(unique(stats::na.omit(breaks)))!=length(quants)){
      for(j in 2:(length(quants)-1)){
        breaks[j]<-vals[which(vals2>quants[j])[1]]
      }
      vals<-vals[-length(vals)]
      vals2<-cumsum(vals)/sum(vals)
    }
  }

  #create matrix for terra classification
  breaks.matrix <- c(breaks[1], breaks[2], 1,
                     breaks[2], breaks[3], 2,
                     breaks[3], breaks[4], 3,
                     breaks[4], breaks[5], 4,
                     breaks[5], breaks[6], 5)
  breaks <- matrix(breaks.matrix, ncol=3, byrow=TRUE)

  return(breaks)
}
