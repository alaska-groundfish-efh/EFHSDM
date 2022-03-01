# EFHSDM (in progress)

Authors:
@jeremyharris7
@James-Thorson-NOAA
@Ned-Laman-NOAA
@MargaretSiple-NOAA
@jodipirtle


# Statement of purpose
This package is designed to produce SDMs and SDM visualizations as part of the 2022 EFH 5-year review. It is designed to be moderately flexible and can be expanded in the future, but for now assumes that abundance prediction is the end goal. 

## Dependencies
This project requires the following packages. Bear in mind that the `maxnet` package is still under development and changes occasionally. The `akgfmaps` package must be installed from Sean Rohan's GitHub.

```{r eval=FALSE}
# Packages
xtable, XML, raster, rgdal, gstat, sp, sf, stars, akgfmaps, ggplot, viridis, gridExtra, patchwork, MASS, scales, labeling, maxnet, ENMeval, PresenceAbsence, mgcv
```

## Roadmap
Sections of the analysis are included as separate scripts. The general strategy is to use the functions provided in `Functions_Maxent.R` and `Functions_GamModel.R` to produce models, abundance rasters, effects estimates, and other outputs in a standard-ish format. Then, the scripts `Functions_LoadMap.R` and `Functions_Ensemble.R` provide more general methods for mapping or plotting model outputs and combining inferences from multiple models.

The `Meatgrinder.R` script provides an example of a workflow that combines these functions to make an ensemble SDM and the accompanying maps. It contains control logic designed to accommodate the needs of the 2022 EFH 5-year Review process, which involved running and keeping track of over 200 individual species/lifestages.


1) `Functions_Maxent.R` - this script provides the functions for quickly using maxnet models. 
2) `Functions_GamModel.R` - this script provides functions for conducting several operations with GAMs. The code has been tested with binomial, poisson, negative binomial, and ziplss hurdle models. 
3) `Functions_LoadMap.R` - this script provides versatile functions for mapping or plotting the model types mentioned above. It also contains some functions used to set up data specific to the EFH program that are unlikely to be broadly useful. 
4) `Functions_Ensemble.R` - this script provides functions designed to combine models and model products produced in the previous scripts
5) `Functions_Xtable.R` - this script produces a few standard output tables to summarize model or ensemble performance.
6) Functions_akgfmaps.R - this script provides standard functions to plot maps using the akgfmaps package that are of higher quality than those included elsewhere

Functions may be general or specific to a type of model

### GAM Specific Functions

| Function               |Use                                                                |
|:-----------------------|:------------------------------------------------------------------|
| `FitGAM()`             | Fit and term selection for binomial, poisson, and negbin GAMs |
| `FitHurdleGAM()`       | Fit and term selection for ziplss hurdle GAMs                 |
| `MakeGAMAbundance()`   | Produce a prediction raster from a GAM model and covariate stack |
| `GetGAMEffects()`      | Produce covariate effect estimates and CIs for GAM models        |
| `GAMStats()`           | Produce jacknife estimates of relative deviance explained (slow) |

* The functions `AutodetectGAMTerms()` and `AssembleGAMFormula()` are used internally by other
functions, but are not strictly necessary for all users.

### Maxnet Specific Functions
| Function               |Use                                                                |
|:-----------------------|:------------------------------------------------------------------|
| `FitMaxnet()`              | Fitting and term selection for maxnet models                      |
| `GetMaxnetEffects()`       | Produces covariate effect estimates and CIs for maxnet models     |
| `MakeMaxEntAbundance()`    | Produces a prediction raster from a GAM model and covariate stack |
| `MaxnetCoefs()`            | Returns the number of coefficients pertaining to each covariate   |
| `MaxnetStats()`            | Produces jackknife estimates of relative deviance explained       |


### General Functions
| Function               |Use                                                                |
|:-----------------------|:------------------------------------------------------------------|
| `RMSE()`                   | Calculated the RMSE for a set of observations and predictions     |
| `PDE()`                    | Calculated the deviance explained (Poisson) from obs and preds    |
| `CrossValidateModel()`     | Conducts k-fold cross-validatation of GAM or maxnet models        |
| `FindEFHbreaks()`          | Returns the EFH breaks for a given abundance map                  |
| `MakeVarianceRasters()`    | Produces raster of predicted variance based on CV folds (slow)    |

### Ensemble Functions
| Function       |Use                                                            |
|:-----------------------|:------------------------------------------------------------------|
| `MakeEnsemble()`           | Calculates the weights for a set a models based on RMSE           |
| `ValidateEnsemble()`       | Makes some basic validation plots and returns prediction data     |
| `MakeEnsembleAbundance()`  | Combines multiple rasters into a weighted average                 |
| `GetEnsembleVariance()`    | Estimates the variance in ensemble predictions                    |
| `GetEnsembleEffects()`     | Produces covariate effect estimates for an ensemble               |

There are numerous plotting functions. The ones based on `akgfmaps` are recommended.

### Work Flow

First one will need to prepare any data and covariate rasters. Data should be organized in a data frame with columns for species and any covariates, offsets, etc. Rasters should be combined into a stack. 

The functions are typically called top to bottom. Begin by fitting a model using FitGAM, FitHurdleGAM, or FitMaxnet. The resulting model is used with a MakeGAMAbundance (or MakeMaxEntAbundance) to create an abundance raster, GetGAMEffects to estimate covariate effects, and GAMStats to obtain covariate contributions. Then one follows with the generic CrossValidateMdoel to produce CV models and output a useful data frame of both in-bag and out-of-bag predictions. This data frame can be used to calculate PDE and RMSE. MakeVarianceRasters produces a variance map. FindEFHbreaks gives the abundance thresholds that define each EFH quantile. The final EFH map is made by passing the EFHbreaks

### Simplified Vignette

First, begin by loading the data and the covariate rasters. For example purposes, we will used only the last  years of data and only a few covariates. We will also reduce the resolution of the rasters.


```{r include=F}
source("C:/Users/jeremy.harris/Work/EFH_SDM_code/EFH_SDM/Functions_LoadMap.R")
source("C:/Users/jeremy.harris/Work/EFH_SDM_code/EFH_SDM/Functions_Maxent.R")
source("C:/Users/jeremy.harris/Work/EFH_SDM_code/EFH_SDM/Functions_GamModel.R")
source("C:/Users/jeremy.harris/Work/EFH_SDM_code/EFH_SDM/Functions_Xtable.R")
source("C:/Users/jeremy.harris/Work/EFH_SDM_code/EFH_SDM/Functions_Ensemble.R")
source("C:/Users/jeremy.harris/Work/EFH_SDM_code/EFH_SDM/Functions_akgfmaps.R")
```
``` r
region.data<-read.csv("Y:/RACE_EFH_Variables/Trawl_Models/GOA/all_GOA_data_2021.csv")
region.data<-subset(region.data,year>=2012)
region.data$sponge<-as.integer(region.data$sponge>0)
region.data$logarea<-log(region.data$area)

bathy <- raster("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Bathy")
btemp <- raster("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Btemp")
btemp<-crop(x = btemp,y=bathy)
lat <- raster::init(bathy, v ='y')
lat <- raster::mask(lat, bathy,overwrite = F)
lon <- raster::init(bathy, v ='x')
lon <- raster::mask(lon, bathy,overwrite = F)
slope <- raster("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Slope")
sponge <- raster("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Spongefactor")

raster.stack <-raster::stack(lon,lat,bathy,btemp,slope,sponge)
names(raster.stack) <-c("lon","lat","bdepth","btemp","slope","sponge")
```

Next we will fit a basic Poisson model and generate an abundance map
``` r
gam.form<-formula("a_atf ~ s(lon,lat,bs = 'ds',m=c(1,.5), k=10) + s(bdepth, bs='tp',m=1,k=4) + s(btemp, bs='tp',m=1,k=4) + s(slope, bs='tp',m=1,k=4) + offset(logarea)")
poisson.model<-FitGAM(gam.formula=gam.form,data=region.data,family.gam="poisson")
```

Now we can make a map
``` r
poisson.abundance<-MakeGAMAbundance(poisson.model,raster.stack)
abundance.plot<-MakeAKGFDensityplot(region="goa",density.map=poisson.abundance,buffer=.98,title.name="Adult ATF",legend.title="Abundance")
print(abundance.plot)
```

Get the covariate effects
``` r
poisson.effects<-GetGAMEffects(poisson.model,data=region.data)
plotEffects(poisson.effects)
```

Crossvalidate the model
``` r
poisson.cv<-CrossValidateModel(model=poisson.model,data=region.data,folds=10,model.type="gam",key="hauljoin")
poisson.preds<-poisson.cv[[1]]
head(poisson.preds)
RMSE(obs=poisson.preds$abund,pred=poisson.preds$cvpred)
PDE(obs=poisson.preds$abund,pred=poisson.preds$cvpred)

```
And find the EFH
``` r
poisson.breaks<-FindEFHbreaks(poisson.abundance,method="percentile")
poisson.efh<-raster::cut(poisson.abundance,poisson.breaks)
efh.plot<-MakeAKGFEFHplot(region="goa",efh.map=poisson.efh,title.name="Adult ATF",legend.title="Percentiles")
print(efh.plot)
```

## Legal disclaimer
>This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.
