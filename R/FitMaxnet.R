#' Fit MaxEnt model
#' @description Fit a MaxEnt model for presence/absence
#'
#' @details
#' This function fits a maxnet model. Maxnet is somewhat more finicky about NA values and
#' covariate names, so code has been added to make sure those match
#' @param data a data frame containing the covariates and presence/absence data
#' @param species character; the name of the column containing the dependent variable
#' @param vars a vector of names for columns that contain the covariates
#' @param reduce Logical; removes covariates with 0 influence from the model
#' @param regmult Numeric; a  regularization multiplier value
#' @param facs a vector of names for columns for covariates that should be treated as factors
#'
#' @return a fitted maxnet model
#' @export
#'
#' @examples
#' species.data <- readRDS(system.file("test_files", "goa_data_logarea_folds.rds", package = "EFHSDM"))
#' maxnet.covars <- c(
#'   "bcurrentU", "bcurrentV", "bcurrentUSD", "bcurrentVSD", "bdepth",
#'   "slope", "aspectE", "aspectN", "curve", "btemp", "tmax", "rocky", "BPI"
#' )
#' cofactors <- c("sponge", "coral", "pen")
#' r.mult <- 0.5
#' maxnet.model0 <- FitMaxnet(
#'   data = species.data, species = "dogfish", vars = maxnet.covars,
#'   facs = cofactors,
#'   regmult = r.mult, reduce = T
#' )

FitMaxnet<-function(data,
                    species,
                    vars,
                    reduce=F,
                    regmult=1,
                    facs=NULL){

  presence.vec<-as.integer(data[,species]>0)
  maxnet.data<-data[,c(vars,facs)]

  if(length(facs)>0){
    for(f in facs){
      maxnet.data[,f]<-as.factor(maxnet.data[,f])
    }
  }
  # unlike maxent, any NAs will crash maxnet
  drops<-NULL
  #Need to filter all NAs out
  for(i in 1:ncol(maxnet.data)){
    drops<-c(drops,which(is.na(maxnet.data[,i])))
  }
  if(length(drops)>0){
    presence.vec<-presence.vec[-unique(drops)]
    maxnet.data<-maxnet.data[-unique(drops),]
  }

  maxnet.model<-maxnet::maxnet(p = presence.vec,data = maxnet.data,regmult = regmult)

  # remove variables that aren't contributing
  if(reduce){
    m.coefs<-MaxnetCoefs(maxnet.model)
    badvars<-names(m.coefs)[which(m.coefs==0)]

    if(length(badvars)>0){
      badcols<-which(names(maxnet.data)%in%badvars)
      maxnet.data2<-maxnet.data[,-badcols]
      maxnet.model<-maxnet::maxnet(p = presence.vec,data = maxnet.data2,regmult = regmult)
    }
  }
  return(maxnet.model)
}
