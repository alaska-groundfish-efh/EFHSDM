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

  pred[pred<.00001]<-.00001
  term1<-obs*log(obs/pred)
  term1[is.nan(term1)]<-0
  term2<-obs-pred

  pdev<-2*sum(term1-term2)
  return(1-(pdev/nulldev))
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
    sample <- stats::na.omit(raster::getValues(abund.raster))
    sample[sample <= threshold] <- NA
    breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
    breaks[1]<-0
    breaks[length(breaks)]<-Inf
  }
  if(method=="cumulative"){
    # Decide whether to sample at given locations or use the whole thing
    if(is.null(data)){
      vals<-stats::na.omit(sort(raster::getValues(abund.raster)))
      vals2<-cumsum(vals)/sum(vals)
    }else{
      vals<-stats::na.omit(sort(raster::extract(abund.raster,data.frame(data$lon,data$lat))))
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
    species<-ifelse(model$family$family=="ziplss",as.character(stats::formula(model)[[1]])[[2]],as.character(stats::formula(model))[[2]])
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
  pb <- utils::txtProgressBar(min = 0, max = n.folds, style = 3)

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
                           link.fx = "cloglog",gam.formula = stats::formula(model),verbose = F))
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

      try(cv.model<-FitHurdleGAM(density.formula = stats::formula(model)[[1]],prob.formula = stats::formula(model)[[2]],
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
        probs<-1-stats::dnbinom(0,mu = mgcv::predict.gam(model,newdata=test.data,type="response"),size = theta)
      }else{
        gamfam<-model$family$family
        probs<-(1-stats::dpois(0,mgcv::predict.gam(object = model,newdata=test.data,type="response")))
      }
      preds<-mgcv::predict.gam(object = model,newdata=test.data,type="response")

      try(cv.model<-FitGAM(data = train.data,reduce=F,family.gam = gamfam,select=F,
                           link.fx = model$family$link,gam.formula = stats::formula(model),verbose = F))
      if(exists("cv.model")){
        cvpreds<-mgcv::predict.gam(object = cv.model,newdata=test.data,type="response")
        if(strsplit(model$family$family,split="[()]")[[1]][1]=="Negative Binomial"){
          cvtheta<-as.numeric(strsplit(cv.model$family[[1]],split="[()]")[[1]][2])
          cvprobs<-1-stats::dnbinom(0,mu = cvpreds,size = cvtheta)
        }else{
          cvprobs<-1-stats::dpois(0,cvpreds)
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
    suppressWarnings(rm(cv.model,cvpreds))
    utils::setTxtProgressBar(pb, i)
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
#' @param make.hist should the histograms be plotted, they sometimes fail and may need to be turned off
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

  main.regr<-stats::lm(error.data2$pred[keepers]~error.data2$abund[keepers])
  main.r2<-summary(main.regr)$r.squared
  main.rmse<-sqrt(sum((stats::na.omit(error.data$abund[keepers]-error.data$pred[keepers]))^2)/nrow(stats::na.omit(error.data[keepers,])))

  #now need to do the cv tests, which should already be in a nice format from the CV function
  cv.regr<-stats::lm(error.data2$cvpred[keepers]~error.data2$abund[keepers])
  cv.r2<-summary(cv.regr)$r.squared
  cv.rmse<-sqrt(sum((stats::na.omit(error.data$abund[keepers]-error.data$cvpred[keepers]))^2)/nrow(stats::na.omit(error.data[keepers,])))

  print(paste("Full model",method,"Rsq =",round(main.r2,2)))
  print(paste("CV",method,"Rsq =",round(cv.r2,2)))
  print(paste("Full model RMSE =",round(main.rmse,3)))
  print(paste("CV RMSE =",round(cv.rmse,3)))

  # make all the plots
  old.par<-graphics::par()[c("mfcol","family","mar","xaxs","yaxs")]
  graphics::par(mfcol = c(ifelse(make.hist==T,3,2),2), family = "sans", mar = c(4,4,3,1))

  stats::qqnorm((error.data$pred[keepers] - error.data$abund[keepers]), main = "Model Predictions")
  stats::qqline((error.data$pred[keepers] - error.data$abund[keepers]))
  if(make.hist==T){graphics::hist((error.data$pred[keepers] - error.data$abund[keepers]), xlab = "Residuals", main = "")}
  pred.max <- ifelse(method=="pearson",stats::quantile(error.data2$pred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  abund.max <- ifelse(method=="pearson",stats::quantile(error.data2$pred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  plot.max<-max(pred.max,abund.max)*1.1

  plot(y=error.data2$pred[keepers], x=error.data2$abund[keepers], ylim = c(0,plot.max), xlim = c(0,plot.max),
       ylab = ifelse(method=="pearson","Predicted","Predicted Ranks"),
       xlab = ifelse(method=="pearson","Observed","Observed Ranks"),main = "", pch = 20)
  graphics::abline(coef = c(0,1), lty = 2)
  graphics::abline(main.regr,col=2)
  graphics::text(1, plot.max*.9, paste(method,"R-squared = ", signif(main.r2,2)), pos = 4)


  #Plots for test/CV data
  pred.max <- ifelse(method=="pearson",stats::quantile(error.data2$cvpred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  abund.max <- ifelse(method=="pearson",stats::quantile(error.data2$cvpred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  plot.max<-max(pred.max,abund.max)*1.1

  stats::qqnorm((error.data$cvpred[keepers] - error.data$abund[keepers]), main = "Test data")
  stats::qqline((error.data$cvpred[keepers] - error.data$abund[keepers]))
  if(make.hist==T){graphics::hist((error.data$cvpred[keepers] - error.data$abund[keepers]), xlab = "Residuals", main = "")}

  plot(y=error.data2$cvpred[keepers], x=error.data2$abund[keepers], ylim = c(0,plot.max), xlim = c(0,plot.max),
       ylab = ifelse(method=="pearson","Predicted","Predicted Ranks"),
       xlab = ifelse(method=="pearson","Observed","Observed Ranks"),main = "", pch = 20)
  graphics::abline(coef = c(0,1), lty = 2)
  graphics::abline(cv.regr,col=2)
  graphics::text(1, plot.max*.9, paste(method,"R-squared = ", signif(cv.r2,2)), pos = 4)
  suppressWarnings(graphics::par(old.par))
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
    if(is.list(model.list[[m]])){
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
  pb <- utils::txtProgressBar(min = 0, max = length(model.list2), style = 3)

  # loop through and get predictions for each model
  for(m in 1:length(model.list2)){
    if(model.type=="maxnet"){
      out.data[,m]<-exp(stats::predict(model.list2[[m]],newdata=data,type="link")+model.list2[[m]]$entropy)
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
    utils::setTxtProgressBar(pb, m)
  }
  close(pb)
  # now tally things up

  raster.template<-raster::raster(raster.stack)

  variances<-apply(X = out.data*scale.factor,MARGIN = 1,FUN = stats::var)
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





