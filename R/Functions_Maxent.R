# This one will institute functions for running maxent and maxnet
# Maxent is on its way out, and most functions are specifically designed for maxnet, and will require modification
# to function with maxent

# Setting lists of required packages & installing it
# rpackages <- c("raster", "PresenceAbsence", "maxnet","ENMeval")
#
# # rJava, rgeos, maps
# which_not_installed <- which(rpackages %in% rownames(installed.packages()) == FALSE)
#
# if(length(which_not_installed) > 1){
#   install.packages(rpackages[which_not_installed], dep = TRUE)
# }
# rm(rpackages,which_not_installed)
#
# require(raster)
# require(PresenceAbsence)
# require(maxnet)
# require(ENMeval)

####################################################################################################################################
##################################################MAXENT MODEL######################################################################
####################################################################################################################################
# Pipeline:
# Step 1: Load data and maps using the LoadData and LoadMap functions.
# Step 2: Fit the maxent model or maxnet model using FitMaxent or FitMaxnet
# Step 3: Use the model from step 2 to make an abundance prediction using MakeMaxnetAbundance
# Step 4: Find the break points using the FindEFHBreaks function from the load maps script
# Step 5: Use the "cut" function with those break to convert the abundance estimate to an EFH map

# A variety of additional evaluation functions are available after this point

# This function fits a maxnet model. Maxnet is somewhat more finicky about NA values and
# covariate names, so code has been added to make sure those match

#' Fit MaxEnt model
#' @description Fit a MaxEnt model for presence/absence
#'
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


#' Make abundance/distribution raster from MaxEnt model
#'
#' @description Make an abundance/distribution raster from any dismo or maxnet model. Depending on the settings, it will apply various masks and the cloglog transformation.
#' @param model a fitted maxnet model
#' @param maxent.stack a raster stack containing all covariates
#' @param scale.fac numeric; scale factor that will multiply the model predictions
#' @param land raster; a mask to be applied to the output usually landmasses
#' @param mask raster; an additional mask to apply to the output, if needed
#' @param type character; use "maxnet" for prob suitable habitat or "cloglog" for approximate abundance
#' @param clamp Logical; just leave this as F; shoudl the covariates be restricted to the range observe in fitting the model?
#' @param filename a filename to save results, "" writes to active memory
#'
#' @return a raster map with the desired prediction
#' @export
#'
#' @examples
MakeMaxEntAbundance<-function(model,
                              maxent.stack,
                              scale.fac=1,
                              land=NULL,
                              mask=NULL,
                              type="cloglog",
                              clamp=F,
                              filename=""){

  #correct a common mistake
  if(is.null(filename)||is.na(filename)){filename<-""}

  # Identify the type and make the main prediction
  type=tolower(type)
  # Somewhat counterintuitive, but this is the type if using the cloglog link to make an abundance estimate
  if(type=="cloglog"){
    # since ENMeval 2.0, they got rid of the useful function and I need to make the predictions the long way
    dat<-raster::getValues(maxent.stack)

    # using predict with maxnet will quietly remove the NAs, so need to track them manually
    na.spots<-which(apply(X = dat,MARGIN = 1,FUN = function(x){return(any(is.na(x)))}))
    dat.spots<-which(seq(1:nrow(dat))%in%na.spots==F)

    preds<-stats::predict(model,dat[dat.spots,],type="link")
    preds2<-exp(preds+model$ent)*scale.fac
    new.vals<-vector(length=nrow(dat))
    new.vals[na.spots]<-NA
    new.vals[dat.spots]<-preds2
    habitat.prediction<-raster::setValues(x = raster::raster(maxent.stack),values = new.vals)
  }
  # this makes a habitat suitability map from a maxnet model
  if(type=="maxnet"){
    dat<-raster::getValues(maxent.stack)

    # using predict with maxnet will quietly remove the NAs, so need to track them manually
    na.spots<-which(apply(X = dat,MARGIN = 1,FUN = function(x){return(any(is.na(x)))}))
    dat.spots<-which(seq(1:nrow(dat))%in%na.spots==F)

    preds<-stats::predict(model,dat[dat.spots,],type="cloglog")
    new.vals<-vector(length=nrow(dat))
    new.vals[na.spots]<-NA
    new.vals[dat.spots]<-preds
    habitat.prediction<-raster::setValues(x = raster::raster(maxent.stack),values = new.vals)
  }
  # need to add a check to see about the strange problems with the EBS
  if(is.null(land)==F){
    habitat.prediction<-raster::extend(x=habitat.prediction,y=land)
  }

  # For some reason, the crs info isn't always carrying over
  habitat.prediction@crs<-maxent.stack@crs
  if(filename!=""){raster::writeRaster(x = habitat.prediction,filename = filename, overwrite = TRUE)}

  # Apply additional masks if necessary
  if(is.null(land)==F){
    habitat.prediction<-raster::mask(habitat.prediction, land, inverse = T, overwrite = TRUE,
                                     filename = filename)
  }
  # Apply additional masks if necessary
  if(is.null(mask)==F){
    habitat.prediction<-raster::mask(habitat.prediction, mask, overwrite = TRUE,
                                     filename = filename)
  }
  return(habitat.prediction)
}


#' Get MaxEnt effects
#'
#' @description Grab the estimated covariate effects and variances for all the effects.
#' @param model a maxnet model
#' @param data data frame; typically the same data used to fit the model
#' @param cv.models list;a list of models from cross validation used for confidence intervals
#' @param vars a vector of the names of model terms to be plotted, or "all"
#' @param maxnet2d a list of vectors of variables to be treated as 2D; used for consistency with the GAMs
#' @param add.entropy logical; should be model entropy be added back to the effects, choose F if not estimating abundance
#' @param scale.factor numeric; a scale factor to multiply abundance predictions
#' @param scale character; should results be on the scale of the linear predictor"log" or abundance "abund"
#'
#' @return a list of data frames containing the estimated covariate effects and variance measures for those effects
#' @export
#'
#' @examples
GetMaxnetEffects<-function(model,
                           data,
                           cv.models=NULL,                    #
                           vars="all",                        #
                           maxnet2d=NULL,                     #
                           add.entropy=T,                     #
                           scale.factor=1,                    #
                           scale="log"){                      #

  # check the variable names and restrict things to those requested
  xvars<-names(model$varmax)
  xfacs<-names(model$samplemeans)[names(model$samplemeans)%in%xvars==F]

  xvars2d<-NULL
  #now, remove the ones that overlap with the 2d variables
  if(length(maxnet2d)>0){
    #check that specified vars are actually in the model
    xvars2d0<-xvars[xvars%in%unlist(maxnet2d)]

    for(j in 1:length(maxnet2d)){
      if(sum(maxnet2d[[j]]%in%xvars2d0)>0){
        xvars2d<-c(xvars2d,paste(maxnet2d[[j]],collapse="*"))
      }
    }
    xvars<-xvars[xvars%in%unlist(maxnet2d)==F]
  }

  if(vars!="all"){
    if(length(maxnet2d)>0){
      xvars2d<-xvars2d[xvars2d%in%vars]
    }else{
      xvars2d<-NULL
    }
    xvars<-xvars[xvars%in%vars]
    xfacs<-xfacs[xfacs%in%vars]
  }
  allvars<-c(xvars2d,xvars,xfacs)
  ntot<-length(allvars)

  if(ntot==0){
    warning("Warning: Variables ",vars," not found in model")
    return(list(data.frame(NA,effect=ifelse(scale=="log",log(.01),0))))
  }

  # if names aren't supplied, use the terms from the model
  if(is.null(nice.names)){
    xvar2dnames<-unlist(strsplit(xvars2d,"[*]"))
    nice.names<-data.frame(var=c(xvar2dnames,xvars,xfacs),name=c(xvar2dnames,xvars,xfacs),stringsAsFactors = F)
  }
  # will try to fix it so that the order is consistent with the gams
  vars2<-vector(length=ntot)
  varnames<-vector(length=ntot)
  v=1
  for(i in 1:nrow(nice.names)){
    if(nice.names[i,1]%in%allvars){
      vars2[v]<-which(allvars==nice.names[i,1])
      varnames[v]<-nice.names[i,2]
      v<-v+1
    }
  }

  # Set up the parameters for the plots
  ent<-model$entropy*as.integer(add.entropy)
  names.vec<-NULL
  out.list<-list()
  list.index<-1

  if(length(xvars2d)>0){
    for(i in 1:length(xvars2d)){
      x.name<-strsplit(xvars2d[i],split="[*]")[[1]][1]
      y.name<-strsplit(xvars2d[i],split="[*]")[[1]][2]

      xseq<-seq(from=min(data[,x.name]),to=max(data[,x.name]),length.out=40)
      yseq<-seq(from=min(data[,y.name]),to=max(data[,y.name]),length.out=40)

      if(x.name%in%names(model$varmax)){
        effect.x<-ent+as.vector(maxnet::response.plot(model,v=x.name,plot=F,type="link")$pred)
      }else{
        effect.x<-log(.01)
      }
      if(y.name%in%names(model$varmax)){
        effect.y<-ent+as.vector(maxnet::response.plot(model,v=y.name,plot=F,type="link")$pred)
      }else{
        effect.y<-log(.01)
      }

      # this is too big though, so shrink it to 40x40 for compatibility with the GAMS
      if(length(effect.x)==1){
        effect.x<-rep(effect.x,40)
      }else{
        effect.x<-effect.x[round(1:40*2.5)]
      }
      if(length(effect.y)==1){
        effect.y<-rep(effect.y,40)
      }else{
        effect.y<-effect.y[round(1:40*2.5)]
      }
      dat<-data.frame(x=rep(xseq,40),y=rep(yseq,each=40),effect=rep(effect.x,40)+rep(effect.y,each=40))
      if(scale=="abund"){
        dat$effect<-exp(dat$effect)*scale.factor
      }else{
        dat$effect<-dat$effect+log(scale.factor)
      }
      out.list[[list.index]]<-dat
      list.index<-list.index+1
    }
  }

  if(length(xvars)>0){
    for(i in 1:length(xvars)){
      # calculate the main effects
      dat<-data.frame(x=seq(from=model$varmin[xvars[i]],to=model$varmax[xvars[i]],length.out = 100))
      suppressWarnings(dat$effect<-ent+as.vector(maxnet::response.plot(model,v=xvars[i],plot=F,type="link")$pred))

      if(scale=="abund"){
        dat$effect<-exp(dat$effect)*scale.factor
      }else{
        dat$effect<-dat$effect+log(scale.factor)
      }

      # if supplied, calculate the cv effects in a loop
      if(is.null(cv.models)==F){
        cv.dat<-matrix(data=NA,nrow=100,ncol=length(cv.models))
        for(f in 1:length(cv.models)){
          if(is.list(cv.models[[f]])){
            suppressWarnings(cv.dat[,f]<-maxnet::response.plot(cv.models[[f]],type = "link",v = xvars[i],plot=F)$pred
                             +cv.models[[f]]$entropy*as.integer(add.entropy))
          }else{
            cv.dat[,f]<-NA
          }
        }
        if(scale=="abund"){
          cv.dat<-exp(cv.dat)*scale.factor
        }else{
          cv.dat<-cv.dat+log(scale.factor)
        }
        colnames(cv.dat)<-paste0("CV",1:length(cv.models))

        uppers<-apply(X = cv.dat,MARGIN = 1,FUN = "quantile",probs=.95,na.rm=T)
        lowers<-apply(X = cv.dat,MARGIN = 1,FUN = "quantile",probs=.05,na.rm=T)
        cv.var<-apply(X = cv.dat,MARGIN = 1,FUN = var,na.rm=T)

        out.list[[list.index]]<-data.frame(dat,var=cv.var,upper=uppers,lower=lowers,cv.dat)
        list.index<-list.index+1
      }else{
        out.list[[list.index]]<-dat
        list.index<-list.index+1
      }
    }
  }
  if(length(xfacs)>0){
    for(i in 1:length(xfacs)){
      # calculate the main effects
      dat<-data.frame(model$levels[xfacs[i]])
      names(dat)<-"x"
      suppressWarnings(dat$effect<-ent+as.vector(maxnet::response.plot(model,v=xfacs[i],plot=F,type="link")$pred))

      if(scale=="abund"){
        dat$effect<-exp(dat$effect)*scale.factor
      }else{
        dat$effect<-dat$effect+log(scale.factor)
      }
      # if supplied, calculate the cv effects in a loop
      if(is.null(cv.models)==F){
        cv.dat<-matrix(data=NA,nrow=nrow(dat),ncol=length(cv.models))
        for(f in 1:length(cv.models)){
          if(is.list(cv.models[[f]])){
            cv.dat[,f]<-maxnet::response.plot(cv.models[[f]],type = "link",
                                      v = xfacs[i],plot=F,levels = dat$x)$pred+ent
          }else{
            cv.dat[,f]<-NA
          }
        }
        if(scale=="abund"){
          cv.dat<-exp(cv.dat)*scale.factor
        }else{
          cv.dat<-cv.dat+log(scale.factor)
        }
        colnames(cv.dat)<-paste0("CV",1:length(cv.models))

        uppers<-apply(X = cv.dat,MARGIN = 1,FUN = "quantile",probs=.95,na.rm=T)
        lowers<-apply(X = cv.dat,MARGIN = 1,FUN = "quantile",probs=.05,na.rm=T)
        cv.var<-apply(X = cv.dat,MARGIN = 1,FUN = var,na.rm=T)

        out.list[[list.index]]<-data.frame(dat,var=cv.var,upper=uppers,lower=lowers,cv.dat)
        list.index<-list.index+1
      }else{
        out.list[[list.index]]<-dat
        list.index<-list.index+1
      }
    }
  }
  names(out.list)<-c(xvars2d,xvars,xfacs)
  return(out.list)
}


#' Get MaxEnt coefficients
#'
#' @description Get the number of coefficients for each variable in a maxnet model
#' @details If you want things to match the gams, supply a list for maxnet2d
#' @param model a fitted maxnet model
#' @param maxnet2d a list of vectors of variables to be treated as 2D; used for consistency with the GAMs
#'
#' @return a named vector listing the number of coefficients that apply to each term
#' @export
#'
#' @examples
MaxnetCoefs<-function(model,
                      maxnet2d=NULL){
  # This first section that is to find how many coefficients each covariate uses
  covar.names<-names(model$samplemeans)
  beta.names<-strsplit(names(model$betas),split=c("\\(|\\:|\\^|\\)"))

  if(length(maxnet2d)>0){
    coef.vec2<-vector(length=length(maxnet2d))
    vars1d<-covar.names[covar.names%in%unlist(maxnet2d)==F]
    maxnet2d.name<-unlist(lapply(X = maxnet2d,FUN = function(x){paste(x,collapse="*")}))

    for(i in 1:length(maxnet2d)){
      coef2d1<-unlist(lapply(X=beta.names,FUN=function(x,y){y%in%x},y=maxnet2d[[i]][1]))
      coef2d2<-unlist(lapply(X=beta.names,FUN=function(x,y){y%in%x},y=maxnet2d[[i]][2]))

      coef.vec2[i]<-sum((coef2d1+coef2d2)>0)
      names(coef.vec2)<-maxnet2d.name
    }
  }else{
    coef.vec2<-NULL
    vars1d<-covar.names
  }

  coef.vec<-vector(length=length(vars1d))
  for(i in 1:length(vars1d)){
    coef.vec[i]<-sum(unlist(lapply(X=beta.names,FUN=function(x,y){y%in%x},y=vars1d[i])))
  }
  names(coef.vec)<-vars1d
  return(c(coef.vec2,coef.vec))
}


#' Get MaxEnt stats
#' @description This function institutes a jackknife estimate of the importance of each covariate for a maxnet model. Can sometimes run slowly if there are convergence issues.
#' @param model a fitted maxnet model
#' @param maxnet2d a list of vectors of variables to be treated as 2D; used for consistency with the GAMs
#' @param regmult numeric; a regularization penalty multiplier
#' @param data data frame; typically the same data used to fit the model
#' @param species character; the name of the column for the dependent variable in the model
#'
#' @return a named vector with estimates of the relative deviance explained by each term
#' @export
#'
#' @examples
MaxnetStats<-function(model,                 # a maxnet model
                      maxnet2d=NULL,         # a list of terms that are to be considered jointly
                      regmult=1,             # a multiplier for regularization
                      data,                  # the data used to fit that model
                      species){              # the species or column of data for the dependent variable

  # This first section that is to find how many coefficients each covariate uses
  covar.names<-names(model$samplemeans)

  if(length(maxnet2d)>0){
    coef.vec2<-vector(length=length(maxnet2d))
    vars1d<-covar.names[covar.names%in%unlist(maxnet2d)==F]
    maxnet2d.name<-unlist(lapply(X = maxnet2d,FUN = function(x){paste(x,collapse="*")}))
  }else{
    vars1d<-covar.names
  }

  #detect factors, and remove them from the vars1d vector
  min.names<-names(model$varmin)
  facs<-vars1d[vars1d%in%min.names==F]
  vars1d<-vars1d[vars1d%in%min.names]

  pb <- utils::txtProgressBar(min = 0, max = length(c(maxnet2d,vars1d,facs)), style = 3)
  pb.i<-1
  # This part makes the estimate of deviance explained for each covariate
  # each category of covariate needs to be handled separately
  if(length(maxnet2d)>0){
    dev.vec2d<-rep(NA,length=length(maxnet2d))
    for(i in 1:length(maxnet2d)){
      #print(paste0("testing deviance for ",maxnet2d.name[i]))
      d.covars<-covar.names[covar.names%in%maxnet2d[[i]]==F]

      try(test.model<-FitMaxnet(data=data,species = species,vars = d.covars,facs = facs,regmult = regmult))
      if(exists("test.model")){
        dev.vec2d[i]<-test.model$dev.ratio[length(test.model$dev.ratio)]
        rm(test.model)
        utils::setTxtProgressBar(pb, pb.i)
        pb.i<-pb.i+1
      }else{
        close(pb)
        break
      }
    }
    names(dev.vec2d)<-maxnet2d.name
  }else{
    dev.vec2d<-NULL
    maxnet2d.name<-NULL
  }

  # 1D covariates
  dev.vec<-rep(NA,length=length(vars1d))
  for(i in 1:length(vars1d)){
    #print(paste0("testing deviance for ",vars1d[i]))
    try(test.model<-FitMaxnet(data=data,species = species,vars = c(unlist(maxnet2d),vars1d[-i]),facs = facs,regmult = regmult))
    if(exists("test.model")){
      dev.vec[i]<-test.model$dev.ratio[length(test.model$dev.ratio)]
      rm(test.model)
      utils::setTxtProgressBar(pb, pb.i)
      pb.i<-pb.i+1
    }else{
      close(pb)
      break
    }
  }
  names(dev.vec)<-vars1d

  # factors
  fac.vec<-rep(NA,length=length(facs))
  for(i in 1:length(facs)){
    #print(paste0("testing deviance for ",facs[i]))
    try(test.model<-FitMaxnet(data=data,species = species,vars = c(unlist(maxnet2d),vars1d),facs = facs[-i],regmult = regmult))
    if(exists("test.model")){
      fac.vec[i]<-test.model$dev.ratio[length(test.model$dev.ratio)]
      rm(test.model)
      utils::setTxtProgressBar(pb, pb.i)
      pb.i<-pb.i+1
    }else{
      close(pb)
      break
    }
  }
  close(pb)

  names(fac.vec)<-facs

  out.dev.vec<-c(dev.vec2d,dev.vec,fac.vec)
  if(sum(is.na(out.dev.vec))==0){
    dev.lost<-1-out.dev.vec/model$dev.ratio[length(model$dev.ratio)]
    dev.exp<-dev.lost/sum(dev.lost,na.rm=T)*100

    #sometimes you end up with negative deviance, so correct that
    if(min(dev.exp)<0){
      dev.exp<-dev.exp-min(dev.exp)
      dev.exp<-dev.exp/sum(dev.exp)*100
    }
  }else{
    dev.exp<-rep(NA,times=length(c(maxnet2d,vars1d,facs)))
  }
  names(dev.exp)<-c(maxnet2d.name,vars1d,facs)
  return(dev.exp)
}
