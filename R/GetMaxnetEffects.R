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
