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
    pres.data<-data.frame(data[pres[pres.randos],],"group"=pres.group)
    abs.data<-data.frame(data[abs[abs.randos],],"group"=abs.group)

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
      preds<-exp(stats::predict(object = model,newdata=test.data,response="link")+model$entropy)
      probs<-stats::predict(object = model,newdata=test.data,type="cloglog")
      # then on to the cv model
      vars0<-names(model$samplemeans)
      facs<-vars0[vars0%in%names(model$varmax)==F]

      try(cv.model<-FitMaxnet(data = train.data,species = species,vars = names(model$varmax),facs = facs,regmult = regmult))
      if(exists("cv.model")){
        cvpreds<-exp(stats::predict(object = cv.model,newdata=test.data,response="link")+cv.model$entropy)
        cvprobs<-stats::predict(object = cv.model,newdata=test.data,type="cloglog")
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
