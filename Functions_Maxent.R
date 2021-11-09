# This one will institute functions for running maxent and maxnet
# Maxent is on its way out, and most functions are specifically designed for maxnet, and will require modification
# to function with maxent

require(dismo)
require(raster)
require(PresenceAbsence)
require(maxnet)
require(ENMeval)

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
FitMaxnet<-function(data,                        # A data frame containing both the covariates and presence data
                    species,                     # The column from the data frame that has presence data, does not need to be binary
                    vars,                        # a vector of names for columns that contain the covariates
                    reduce=F,                    # removes variables with 0 influence from the model
                    regmult=1,                   # a  regularization multiplier value
                    facs=NULL){                  # a vector of names for columns for covariates that should be treated as factors
        
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
        
        maxnet.model<-maxnet(p = presence.vec,data = maxnet.data,regmult = regmult)
        
        # remove variables that aren't contributing
        if(reduce){
                m.coefs<-MaxnetCoefs(maxnet.model)
                badvars<-names(m.coefs)[which(m.coefs==0)]
                
                if(length(badvars)>0){
                        badcols<-which(names(maxnet.data)%in%badvars)
                        maxnet.data2<-maxnet.data[,-badcols]
                        maxnet.model<-maxnet(p = presence.vec,data = maxnet.data2,regmult = regmult)
                }
        }
        return(maxnet.model)
}

# This function makes an abundance/distribution raster from any dismo or maxnet model,
# Depending on the settings, it will apply various masks and the cloglog transformation
MakeMaxEntAbundance<-function(model,              # a model fit using either dismo::maxent or maxnet
                              maxent.stack,       # a stack of covariate rasters, dismo::maxent has special requirements
                              scale.fac=1,        # a scale factor to mutiple the results, if estimating abundance
                              land=NULL,          # a land raster to use as a mask
                              mask=NULL,          # any additional masks that should be applied
                              type="cloglog",     # the type of output desired (options: "maxnet","cloglog")
                              clamp=F,            # shoudl the covariates be restricted to the range observe in fitting the model
                              filename=""){       # a filename to save results, "" writes to active memory
        
        #correct a common mistake
        if(is.null(filename)||is.na(filename)){filename<-""}
        
        # Identify the type and make the main prediction
        type=tolower(type)
        # Somewhat counterintuitive, but this is the type if using the cloglog link to make an abundance estimate
        if(type=="cloglog"){
                # since ENMeval 2.0, they got rid of the useful function and I need to make the predictions the long way
                dat<-getValues(maxent.stack)
                
                # using predict with maxnet will quietly remove the NAs, so need to track them manually
                na.spots<-which(apply(X = dat,MARGIN = 1,FUN = function(x){return(any(is.na(x)))}))
                dat.spots<-which(seq(1:nrow(dat))%in%na.spots==F)
                
                preds<-predict(model,dat[dat.spots,],type="link")
                preds2<-exp(preds+model$ent)*scale.fac
                new.vals<-vector(length=nrow(dat))
                new.vals[na.spots]<-NA
                new.vals[dat.spots]<-preds2
                habitat.prediction<-setValues(x = raster(maxent.stack),values = new.vals)
        }
        # this makes a habitat suitability map from a maxnet model
        if(type=="maxnet"){
                dat<-getValues(maxent.stack)
                
                # using predict with maxnet will quietly remove the NAs, so need to track them manually
                na.spots<-which(apply(X = dat,MARGIN = 1,FUN = function(x){return(any(is.na(x)))}))
                dat.spots<-which(seq(1:nrow(dat))%in%na.spots==F)
                
                preds<-predict(model,dat[dat.spots,],type="cloglog")
                new.vals<-vector(length=nrow(dat))
                new.vals[na.spots]<-NA
                new.vals[dat.spots]<-preds
                habitat.prediction<-setValues(x = raster(maxent.stack),values = new.vals)
        }
        # need to add a check to see about the strange problems with the EBS
        if(compareRaster(x = habitat.prediction,land,stopiffalse = F)==F){
                habitat.prediction<-extend(x=habitat.prediction,y=land)
        }
        
        # For some reason, the crs info isn't always carrying over
        habitat.prediction@crs<-maxent.stack@crs
        if(filename!=""){writeRaster(x = habitat.prediction,filename = filename, overwrite = TRUE)}
        
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

# Function to plot effects from a maxnet model
# if the cv.models are provided, it will calculate confidence intervals based on that
plotMaxnetEffects<-function(model,                             # a maxnet model
                            data=NULL,                         # data set with the covariates, only used for the rug
                            cv.models=NULL,                    # a list of models from cross validation used for confidence intervals
                            nice.names=NULL,                   # a dataframe with two columns matching variable names to nicer names for plotting
                            output=F,                          # should a list of effects be returned
                            make.plot=T,                       # should it actually make the plots
                            vars="all",                        # a vector of model terms to be plotted, or "all"
                            maxnet2d=NULL,                     # a list of 2d variables
                            add.entropy=T,                     # should be model entropy be added back to the effects
                            scale.factor=1,                    # a scale factor to multiply predictions
                            scale="log"){                      # should results be on the scale of the linear predictor(log) or abundance (abund)
        
        # check the variable names and restrict things to those requested
        xvars<-names(model$varmax)
        xfacs<-names(model$samplemeans)[names(model$samplemeans)%in%xvars==F]
        
        #now, remove the ones that overlap with the 2d variables
        if(length(maxnet2d)>0){
                #check that specified vars are actually in the model
                xvars2d0<-xvars[xvars%in%unlist(maxnet2d)]
                
                #need a tricksy check to make sure the whole of2d terms made it through
                for(j in 1:length(maxnet2d)){
                        if(sum(maxnet2d[[j]]%in%xvars2d0)==1){
                                spot<-which(maxnet2d[[j]]%in%xvars2d0)
                                if(spot==1){
                                        before<-NULL
                                }else{
                                        before<-1:(spot-1)
                                }
                                if(spot==length(xvars2d0)){
                                        after<-NULL
                                }else{
                                        after<-(spot+1):length(xvars2d0)
                                }
                                xvars2d0<-c(xvars2d0[before],maxnet2d[[j]],xvars2d0[after])
                        }
                }
                
                
                xvars<-xvars[xvars%in%unlist(maxnet2d)==F]
                
                unlist(lapply(maxnet2d,FUN=function(x,y){sum(y%in%x)},y=xvars2d0))==2
                xvars2d<-NULL
                for(i in 1:length(maxnet2d)){
                        if(sum(maxnet2d[[i]]%in%xvars2d0)==2){
                                xvars2d<-c(xvars2d,paste(maxnet2d[[i]],collapse="*"))
                        }
                }
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
        if(make.plot){
                oldpar<-par()
                n.col<-ifelse(ntot>4,3,ifelse(ntot>1,2,1))
                n.rows<-ceiling(ntot/n.col)
                par(mfrow = c(n.rows, n.col), mar = c(4,4,1,0.01),oma=c(.5,.5,.5,.5))
        }
        if(length(xvars2d)>0){
                for(i in 1:length(xvars2d)){
                        x.name<-strsplit(xvars2d[i],split="[*]")[[1]][1]
                        y.name<-strsplit(xvars2d[i],split="[*]")[[1]][2]
                        
                        xseq<-seq(from=min(data[,x.name]),to=max(data[,x.name]),length.out = 40)
                        yseq<-seq(from=min(data[,y.name]),to=max(data[,y.name]),length.out = 40)
                        
                        if(x.name%in%names(model$varmax)){
                                effect.x<-ent+as.vector(response.plot(model,v=x.name,plot=F,type="link"))
                        }else{
                                effect.x<-log(.01)
                        }
                        if(y.name%in%names(model$varmax)){
                                effect.y<-ent+as.vector(response.plot(model,v=y.name,plot=F,type="link"))
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
                        
                        
                        if(make.plot){
                                xlab<-nice.names$name[nice.names$var==x.name]
                                ylab<-nice.names$name[nice.names$var==y.name]
                                plot(NA,xlab=xlab,ylab=ylab,xlim=c(min(xseq),max(xseq)),ylim=c(min(yseq),max(yseq)))
                                contour(x = xseq,y=yseq,z=matrix(data=dat$effect,byrow=T,nrow=40,ncol=40),add=T,nlevels=10)
                                if(xlab%in%c("Current Velocity East (m/s)","bcurrentU")){
                                        abline(h=0,v=0,lty=3)
                                }
                                if(xlab%in%c("Current Velocity East SD (m/s)","bcurrentUSD")){
                                        abline(a=0,b=1,lty=3)
                                }
                        }
                        out.list[[list.index]]<-dat
                        list.index<-list.index+1
                }
        }
        
        if(length(xvars)>0){
                for(i in 1:length(xvars)){
                        # calculate the main effects
                        dat<-data.frame(x=seq(from=model$varmin[xvars[i]],to=model$varmax[xvars[i]],length.out = 100))
                        suppressWarnings(dat$effect<-ent+as.vector(response.plot(model,v=xvars[i],plot=F,type="link")))
                        
                        if(scale=="abund"){
                                dat$effect<-exp(dat$effect)*scale.factor
                        }else{
                                dat$effect<-dat$effect+log(scale.factor)
                        }
                        
                        # if supplied, calculate the cv effects in a loop
                        if(is.null(cv.models)==F){
                                cv.dat<-matrix(data=NA,nrow=100,ncol=length(cv.models))
                                for(f in 1:length(cv.models)){
                                        if(is.na(cv.models[[f]])==F){
                                                suppressWarnings(cv.dat[,f]<-response.plot(cv.models[[f]],type = "link",v = xvars[i],plot=F)
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
                                y.lim<-c(min(lowers),max(uppers))
                                
                                out.list[[list.index]]<-data.frame(dat,var=cv.var,upper=uppers,lower=lowers,cv.dat)
                                list.index<-list.index+1
                        }else{
                                y.lim<-c(min(dat$effect),max(dat$effect))
                                out.list[[list.index]]<-dat
                                list.index<-list.index+1
                        }
                        
                        # actually make the plots
                        if(make.plot){
                                y.lab<-ifelse(i%%n.col==1| n.col==1,"Variable effect" ,NA)
                                x.lab<-nice.names$name[nice.names$var==xvars[i]]
                                
                                plot(x=dat$x,y=dat$effect,type="l",ylab = y.lab, xlab =x.lab,ylim=y.lim)
                                
                                if(is.null(data)==F){rug(x=data[,xvars[i]])}
                                if(is.null(cv.models)==F){
                                        lines(x=dat[,1],y=uppers,lty=3)
                                        lines(x=dat[,1],y=lowers,lty=3)
                                }
                        }
                }
        }
        if(length(xfacs)>0){
                for(i in 1:length(xfacs)){
                        # calculate the main effects
                        dat<-data.frame(model$levels[xfacs[i]])
                        names(dat)<-"x"
                        suppressWarnings(dat$effect<-ent+as.vector(response.plot(model,v=xfacs[i],plot=F,type="link")))
                        
                        if(scale=="abund"){
                                dat$effect<-exp(dat$effect)*scale.factor
                        }else{
                                dat$effect<-dat$effect+log(scale.factor)
                        }
                        # if supplied, calculate the cv effects in a loop
                        if(is.null(cv.models)==F){
                                cv.dat<-matrix(data=NA,nrow=nrow(dat),ncol=length(cv.models))
                                for(f in 1:length(cv.models)){
                                        if(is.na(cv.models[[f]])==F){
                                                cv.dat[,f]<-response.plot(cv.models[[f]],type = "link",
                                                                          v = xfacs[i],plot=F,levels = dat$x)+ent
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
                                y.lim<-c(min(lowers),max(uppers))
                                
                                out.list[[list.index]]<-data.frame(dat,var=cv.var,upper=uppers,lower=lowers,cv.dat)
                                list.index<-list.index+1
                        }else{
                                y.lim<-c(min(dat$effect),max(dat$effect))
                                out.list[[list.index]]<-dat
                                list.index<-list.index+1
                        }
                        
                        # actually make the plots
                        if(make.plot){
                                y.lab<-ifelse(i%%n.col==1| n.col==1,"Variable effect" ,NA)
                                x.lab<-nice.names$name[nice.names$var==xfacs[i]]
                                
                                plot(x=as.factor(dat$x),y=dat$effect,ylab = y.lab, xlab =x.lab,ylim=y.lim)
                                
                                if(is.null(data)==F){rug(x=data[,xfacs[i]])}
                                if(is.null(cv.models)==F){
                                        plot(x=as.factor(dat[,1]),y=uppers,lty=3,lwd=.4,add=T)
                                        plot(x=as.factor(dat[,1]),y=lowers,lty=3,lwd=.4,add=T)
                                }
                        }
                        
                }
                
        }
        if(make.plot){suppressWarnings(par(oldpar))}
        names(out.list)<-c(xvars2d,xvars,xfacs)
        
        if(output){return(out.list)}
}


# This one just returns the number of coefficients for each variable in a maxnet model
# if you want things to match the gams, supply a list for maxnet2d
MaxnetCoefs<-function(model,                           # the maxnet model
                      maxnet2d=NULL){                  # a list of vectors containing variables that are considered together
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


# This function institutes a jackknife estimate of the importance of each covariate for a maxnet model
# Can sometimes run slowly if there are convergence issues
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
        
        pb <- txtProgressBar(min = 0, max = length(c(maxnet2d,vars1d,facs)), style = 3)
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
                                setTxtProgressBar(pb, pb.i)
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
                        setTxtProgressBar(pb, pb.i)
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
                        setTxtProgressBar(pb, pb.i)
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
