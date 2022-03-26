# This script includes functions useful for making an ensemble from other models
# It assumes some data formats from other functions in this "package"
#
#' Make model ensemble
#' @description Make an ensemble model vector with weights for each model
#' @details This function calculates the weights for an ensemble model using a very simple rmse method. Future work may allow for more complex methods. One can supply a list of dataset (old way), or one can more sensibly just supply a named vector of RMSE values (better way)
#' @param ensemble.list list of data frames from which rmse can be calculated, not recommended
#' @param rmse named vector of rmse values to calculate weights
#' @param names.vec vector of names corresponding to the list or rmse vector
#' @param minimum numeric; drops any model given less weight than this value
#' @param model.types vector of names of models, only one model of a given type will advance
#'
#' @return a named vector with weights for each model
#' @export
#'
#' @examples
MakeEnsemble<-function(ensemble.list=NA,
                       rmse=NA,
                       names.vec=NA,
                       minimum=NA,
                       model.types=NULL ){


  # if the rmse vector isn't supplied, compute it from the ensemble.list
  if(any(is.na(rmse))){
    if(is.na(names.vec)){
      names.vec<-names(ensemble.list)
    }

    rmse<-vector(length=length(ensemble.list))

    # find the rmse for each set
    for(m in 1:length(ensemble.list)){
      ensemble.dat<-stats::na.omit(ensemble.list[[m]])
      if(nrow(ensemble.dat)>0){
        rmse[m]<-sqrt(sum((ensemble.dat$abund-ensemble.dat$pred)^2)/nrow(ensemble.dat))
      }else{
        rmse[m]<-NA
      }
    }
  }else{
    if(is.na(names.vec)){
      names.vec<-names(rmse)
    }
  }

  # these days we are using the square of RMSE to weight things
  rmse2<-rmse^2

  # weight it, by the inverse
  weights<-(1/(rmse2))/sum(1/(rmse2),na.rm=T)
  names(weights)<-names.vec
  weights[is.na(weights)]<-0

  # apply any adjustments
  if(length(model.types)>0){
    for(m in unique(model.types)){
      index<-which(model.types==m)
      weights2<-weights[index]
      keep<-index[which.max(weights2)]
      drop<-index[index!=keep]
      weights[drop]<-0
    }
    weights<-weights/sum(weights)
  }

  # drop anything that's below the cutoff
  if(is.na(minimum)==F){
    included<-which(weights>=minimum)
    set.to.zero<-which(weights<minimum)
    weights2<-weights[included]
    weights2<-weights2/sum(weights2)
    weights[set.to.zero]<-0
    weights[included]<-weights2
  }
  return(weights)
}

#' Validate ensemble
#' @description This function makes residual plots for the ensemble model and returns a helpful data frame of ensemble predictions
#' @param pred.list list of data frames containing observations and predictions from constituent models
#' @param model.weights named vector of weights for each model
#' @param method character; method for the cor() function
#' @param output logical; should the data frame be returned, or just plots
#' @param make.plots logical; should the plots be output, or just the data frame
#' @param key character; a column name present in pred.list that works as a case identifier; usually "hauljoin"
#' @param latlon logical; should lat and lon be included in the output table, must be present in pred.list
#' @param group character; name of a column to be included in the output table representing CV folds
#' @param histogram logical; should histograms be plotted (they often don't look nice)
#'
#' @return data frame of ensemble observations vs predictions (pred and cvpred are the same for the ensemble)
#' @export
#'
#' @examples
ValidateEnsemble<-function(pred.list,
                           model.weights,
                           method="pearson",
                           output=T,
                           make.plots=T,
                           key=NA,
                           latlon=T,
                           group=NA,
                           histogram=F){

  ensemble.preds<-data.frame(abund=pred.list[[1]]$abund,pred=0)

  if(latlon){
    ensemble.preds<-cbind(pred.list[[1]][,c("lon","lat")],ensemble.preds)
  }
  if(is.na(group)==F){
    ensemble.preds<-cbind(pred.list[[1]][,group],ensemble.preds)
    names(ensemble.preds)[1]<-group
  }
  if(is.na(key)==F){
    ensemble.preds<-cbind(pred.list[[1]][,key],ensemble.preds)
    names(ensemble.preds)[1]<-key
  }

  # Now use the weights to find out the predictions
  for(m in 1:length(pred.list)){
    if(model.weights[m]>0){
      ensemble.preds$pred<-ensemble.preds$pred+pred.list[[m]]$pred*model.weights[m]
    }
  }
  ensemble.preds$prob<-1-stats::dpois(0,ensemble.preds$pred)
  ensemble.preds$cvprob<-ensemble.preds$prob

  keepers<-which(is.na(ensemble.preds$pred)==F & is.infinite(ensemble.preds$pred)==F)

  # compute the summary stats for the output
  ensemble.preds2<-ensemble.preds
  if(method=="spearman"){
    ensemble.preds2$abund<-rank(ensemble.preds$abund)
    ensemble.preds2$pred<-rank(ensemble.preds$pred)
  }

  ensemble.rmse<-sqrt(sum((ensemble.preds$abund-ensemble.preds$pred)^2,na.rm=T)/nrow(ensemble.preds))
  regr <- stats::lm(ensemble.preds2$pred[keepers]~ensemble.preds2$abund[keepers])
  rsqr <- summary(regr)$r.squared

  # ploting portion
  if(make.plots==T){
    old.par<-graphics::par()[c("mfcol","family","mar","xaxs","yaxs")]
    graphics::par(mfcol = c(2+as.integer(histogram),1), family = "sans", mar = c(4,4,3,1))

    stats::qqnorm((ensemble.preds2$abund - ensemble.preds2$pred))
    stats::qqline((ensemble.preds2$abund - ensemble.preds2$pred))

    if(histogram){graphics::hist((ensemble.preds2$abund - ensemble.preds2$pred), xlab = "Residuals", main = "")}

    pred.max <- ifelse(method=="pearson",stats::quantile(ensemble.preds2$pred[keepers],probs=.99,na.rm=T),nrow(ensemble.preds2))
    abund.max <- ifelse(method=="pearson",stats::quantile(ensemble.preds2$pred[keepers],probs=.99,na.rm=T),nrow(ensemble.preds2))
    plot.max<-base::max(pred.max,abund.max)*1.1

    plot(y=ensemble.preds2$pred[keepers], x=ensemble.preds2$abund[keepers], ylim = c(0,plot.max), xlim = c(0,plot.max),
         ylab = ifelse(method=="pearson","Predicted","Predicted Ranks"),
         xlab = ifelse(method=="pearson","Observed","Observed Ranks"),main = "", pch = 20)
    graphics::abline(coef = c(0,1), lty = 2)
    graphics::abline(regr,col=2)
    graphics::text(1, plot.max*.9, paste(method,"R-squared = ", signif(rsqr,2)), pos = 4)

    suppressWarnings(graphics::par(old.par))
  }

  print(paste("Ensemble",method,"Rsq =",round(rsqr,3)))
  print(paste("Ensemble RMSE =",round(ensemble.rmse,2)))

  ensemble.preds$error<-ensemble.preds$pred-ensemble.preds$abund

  if(output==T){return(ensemble.preds)}
}


#' Make ensemble abundance
#' @description Make the ensemble into a raster.
#' @details The abundance rasters should be a list, and must match the ordering of the weights. After this, use the FindEFHbreaks function and others to generate an EFH map as normal.
#' @param model.weights a vector of model weights
#' @param abund.list list of abundance rasters corresponding to model weights
#' @param filename character; a file name for the abundance raster
#'
#' @return raster formed from weighted average of abund.list
#' @export
#'
#' @examples
MakeEnsembleAbundance<-function(model.weights,
                                abund.list,
                                filename=""){

  # intialize a raster with appropriate properties
  good.abund<-which(is.na(abund.list)==F)
  new.abund<-raster::raster(abund.list[[good.abund[1]]])
  new.abund<-raster::setValues(new.abund,values = 0)

  #To make the raster, we'll need another loop
  for(i in 1:length(model.weights)){
    if(model.weights[i]>0){
      new.abund<-new.abund+abund.list[[i]]*model.weights[i]
    }
  }
  if(filename!=""){
    raster::writeRaster(x = new.abund,filename = filename,overwrite=T)
  }
  return(new.abund)
}


#' Get ensemble variance
#'
#' @description This function computes the variance in abundance based on a formula from Burnham and Anderson
#' @param model.weights a vector of model weights
#' @param variance.list list of rasters of the variance in model predictions
#' @param abund.list list of rasters of the predicted abundance from each model
#' @param ensemble.abund raster of ensemble predicted abundance
#'
#' @return raster of ensemble predicted variance in predictions
#' @export
#'
#' @examples
GetEnsembleVariance<-function(model.weights,variance.list,abund.list,ensemble.abund){

  keepers<-which(model.weights>0)
  weights2<-model.weights[keepers]
  variance.list2<-list()
  abund.list2<-list()

  for(k in 1:length(keepers)){
    variance.list2[[k]]<-variance.list[[keepers[k]]]
    abund.list2[[k]]<-abund.list[[keepers[k]]]
  }

  data.spots<-which(is.na(raster::getValues(abund.list2[[1]]))==F)
  if(length(abund.list2)>1){
    for(r in 2:length(abund.list2)){
      spots<-which(is.na(raster::getValues(abund.list2[[r]]))==F)
      data.spots<-data.spots[data.spots%in%spots]
    }
  }
  e.abund<-raster::extract(ensemble.abund,data.spots)

  dat<-matrix(nrow=length(data.spots),ncol=length(keepers))

  for(m in 1:length(keepers)){
    a.dat<-raster::extract(abund.list2[[m]],data.spots)
    v.dat<-raster::extract(variance.list2[[m]],data.spots)

    dat[,m]<-weights2[m]*sqrt(v.dat+(a.dat-e.abund)^2)
  }

  std.error<-apply(dat,MARGIN = 1,FUN = sum)

  out.raster<-raster::raster(abund.list2[[1]])
  val.vec<-raster::getValues(out.raster)
  val.vec[data.spots]<-std.error

  out.raster<-raster::setValues(out.raster,val.vec)
  return(out.raster)
}

#' Get ensemble effects
#'
#' @description This function takes a list of effects from individual models, and the model weights, and computes an ensemble effect estimate. If CV models are available, then this is used for the CI.
#' @param effects.list list of lists of data frame, such as those produced from the GetGAMEffects function
#' @param model.weights vector of numeric weights for each model, must match order of effects.list
#' @param vars character; vector of desired term names or "all"
#' @param scale character; should effects be in "log" or "abund" scale
#'
#' @return list of data frames containing the estimated covariate effects in the ensemble
#' @export
#'
#' @examples
GetEnsembleEffects<-function(effects.list,
                             model.weights,
                             vars="all",
                             scale="log"){

  # detect and remove any models that are NA or have zero weight, and the accompanying types
  effects.list1<-list()
  model.weights1<-NULL
  list.index<-1

  for(m in 1:length(effects.list)){
    if(is.list(effects.list[[m]]) & model.weights[m]>0){
      effects.list1[[list.index]]<-effects.list[[m]]
      names(effects.list1)[list.index]<-names(effects.list)[m]
      model.weights1<-c(model.weights1,model.weights[m])
      list.index<-list.index+1
    }
  }


  #get a list of the actual terms in the ensemble, which is potentially a lot in this case
  var.table<-data.frame(type=NULL,dims=NULL,term=NULL)
  for(m in 1:length(effects.list1)){
    for(n in 1:length(effects.list1[[m]])){
      var.table0<-data.frame(type=ifelse(nrow(effects.list1[[m]][[n]])<10,"factor","smooth"),
                             dims=ifelse("y"%in%names(effects.list1[[m]][[n]]),2,1),
                             term=names(effects.list1[[m]])[n])
      var.table<-rbind(var.table,var.table0)
    }
  }
  var.table<-unique(var.table)

  # now limit to only the vars actually in the models
  if(vars[1]=="all"){
    var.table1<-var.table
  }else{
    var.table1<-var.table[var.table$term%in%vars]
  }

  vars2d<-var.table$term[var.table$dims==2]
  vars1d<-var.table$term[var.table$dims==1 & var.table$type=="smooth"]
  varsf<-var.table$term[var.table$type=="factor"]

  #little cosmetic fix, so that lat*lon is always first
  if("lon*lat"%in%vars2d){vars2d<-unique(c("lon*lat",vars2d))}

  out.list<-list()
  list.index<-1
  # First, the 2D variables
  if(length(vars2d)>0){
    for(v in 1:length(vars2d)){

      temp<-effects.list1[[which(unlist(lapply(effects.list1,FUN=function(x){return(vars2d[v]%in%names(x))})))[1]]]

      # Assumes that the models were fit using the same data
      x.seq<-temp[[which(names(temp)==vars2d[v])]]$x
      y.seq<-temp[[which(names(temp)==vars2d[v])]]$y

      # loop through, grab the effects, and take the weighted average
      v2.preds<-as.data.frame(matrix(nrow=length(x.seq),ncol=length(effects.list1)))
      colnames(v2.preds)<-names(model.weights1)
      for(m in 1:length(effects.list1)){

        if(vars2d[v]%in%names(effects.list1[[m]])){
          v2.preds[,m]<-exp(effects.list1[[m]][[vars2d[v]]]$effect)*model.weights1[m]
        }else{
          v2.preds[,m]<-0
        }
      }

      v2.dat<-data.frame(x=x.seq,y=y.seq,effect=apply(v2.preds,MARGIN = 1,FUN = sum))

      # For the ensemble, default to abund, and change scale to log after the fact
      if(scale=="log"){
        v2.dat$effect[v2.dat$effect==0]<-.01
        v2.dat$effect<-log(v2.dat$effect)
      }

      out.list[[list.index]]<-v2.dat
      names(out.list)[list.index]<-vars2d[v]
      list.index<-list.index+1
    }
  }

  # now the one dimensional variables
  if(length(vars1d)>0){
    for(v in 1:length(vars1d)){

      temp<-effects.list1[[which(unlist(lapply(effects.list1,FUN=function(x){return(vars1d[v]%in%names(x))})))[1]]]

      x.seq<-temp[[which(names(temp)==vars1d[v])]]$x

      v.preds<-as.data.frame(matrix(nrow = length(x.seq),ncol=length(model.weights1)))
      colnames(v.preds)<-names(model.weights1)
      m.vec<-NULL

      # first get the main effect
      for(m in 1:length(effects.list1)){

        if(vars1d[v]%in%names(effects.list1[[m]])){
          v.preds[,m]<-exp(effects.list1[[m]][[vars1d[v]]]$effect)*model.weights1[m]
          m.vec<-c(m.vec,m)
        }else{
          v.preds[,m]<-0
        }
      }

      main.pred<-apply(X = v.preds,MARGIN = 1,FUN = sum)

      # For the ensemble, default to abund, and change scale to log after the fact
      if(scale=="log"){
        main.pred[main.pred==0]<-.01
        main.pred<-log(main.pred)
      }

      n.folds<-sum(unlist(lapply(strsplit(names(effects.list1[[m.vec[1]]][[vars1d[v]]]),split=""),
                                 FUN=function(x){return(tolower(paste0(x[1:2],collapse=""))=="cv")})))

      # now get the cv effects
      if(n.folds>0){

        cv.effects.dat<-as.data.frame(matrix(nrow = length(x.seq),ncol=n.folds))

        for(f in 1:n.folds){
          cv.preds<-as.data.frame(matrix(nrow = length(x.seq),ncol=length(effects.list1)))

          for(m in 1:length(effects.list1)){

            if(m%in%m.vec){
              fold.match<-which(unlist(lapply(strsplit(names(effects.list1[[m]][[vars1d[v]]]),split=""),
                                              FUN=function(x){return(tolower(paste0(x[1:2],collapse=""))=="cv")})))

              cv.preds[,m]<-exp(effects.list1[[m]][[vars1d[v]]][,fold.match[f]])*model.weights1[m]
            }else{
              cv.preds[,m]<-0
            }
          }
          cv.effects.dat[,f]<-apply(X = cv.preds,MARGIN = 1,FUN = sum,na.rm=T)

          if(scale=="log"){
            cv.effects.dat[,f][cv.effects.dat[,f]==0]<-.01
            cv.effects.dat[,f]<-log(cv.effects.dat[,f])
          }
        }
        names(cv.effects.dat)<-paste0("CV",1:ncol(cv.effects.dat))

        uppers<-apply(X = cv.effects.dat,MARGIN = 1,FUN = "quantile",probs=.95,na.rm=T)
        lowers<-apply(X = cv.effects.dat,MARGIN = 1,FUN = "quantile",probs=.05,na.rm=T)


        out.list[[list.index]]<-data.frame(x=x.seq,effect=main.pred,upper=uppers,lower=lowers,cv.effects.dat)
        names(out.list)[list.index]<-vars1d[v]
      }else{
        out.list[[list.index]]<-data.frame(x.seq,effect=main.pred)
        names(out.list)[list.index]<-vars1d[v]
      }
      list.index<-list.index+1
    }
  }

  # Now do the factors
  if(length(varsf)>0){
    for(v in 1:length(varsf)){

      temp<-effects.list1[[which(unlist(lapply(effects.list1,FUN=function(x){return(varsf[v]%in%names(x))})))[1]]]
      f.seq<-factor(temp[[which(names(temp)==varsf[v])]]$x)

      f.preds<-as.data.frame(matrix(nrow = length(f.seq),ncol=length(model.weights1)))
      m.vec<-NULL

      # first get the main effect
      for(m in 1:length(model.weights1)){
        if(varsf[v]%in%names(effects.list1[[m]])){
          f.preds[,m]<-exp(effects.list1[[m]][[varsf[v]]]$effect)*model.weights1[m]
          m.vec<-c(m.vec,m)
        }else{
          f.preds[,m]<-0
        }
      }
      main.pred<-apply(X = f.preds,MARGIN = 1,FUN = sum)

      # For the ensemble, default to abund, and change scale to log after the fact
      if(scale=="log"){
        main.pred[main.pred==0]<-.01
        main.pred<-log(main.pred)
      }

      n.folds<-sum(unlist(lapply(strsplit(names(effects.list1[[m.vec[1]]][[varsf[v]]]),split=""),
                                 FUN=function(x){return(tolower(paste0(x[1:2],collapse=""))=="cv")})))

      # now get the cv effects
      if(n.folds>0){

        cv.effects.dat<-as.data.frame(matrix(nrow = length(f.seq),ncol=n.folds))


        for(f in 1:n.folds){
          cv.preds<-as.data.frame(matrix(nrow = length(f.seq),ncol=length(effects.list1)))

          for(m in 1:length(effects.list1)){
            if(m%in%m.vec){
              fold.match<-which(unlist(lapply(strsplit(names(effects.list1[[m]][[varsf[v]]]),split=""),
                                              FUN=function(x){return(tolower(paste0(x[1:2],collapse=""))=="cv")})))

              cv.preds[,m]<-exp(effects.list1[[m]][[varsf[v]]][,fold.match[f]])*model.weights1[m]
            }else{
              cv.preds[,m]<-0
            }
          }
          cv.effects.dat[,f]<-apply(X = cv.preds,MARGIN = 1,FUN = sum,na.rm=T)

          if(scale=="log"){
            cv.effects.dat[,f][cv.effects.dat[,f]==0]<-.01
            cv.effects.dat[,f]<-log(cv.effects.dat[,f])
          }
        }

        names(cv.effects.dat)<-paste0("CV",1:ncol(cv.effects.dat))

        uppers<-apply(X = cv.effects.dat,MARGIN = 1,FUN = "quantile",probs=.95,na.rm=T)
        lowers<-apply(X = cv.effects.dat,MARGIN = 1,FUN = "quantile",probs=.05,na.rm=T)

        out.list[[list.index]]<-data.frame(x=f.seq,effect=main.pred,upper=uppers,lower=lowers,cv.effects.dat)
        names(out.list)[list.index]<-varsf[v]
      }else{
        out.list[[list.index]]<-data.frame(x=f.seq,effect=main.pred)
        names(out.list)[list.index]<-varsf[v]
      }


      list.index<-list.index+1
    }
  }
  return(out.list)
}

