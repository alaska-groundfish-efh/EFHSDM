# EFHSDM (in progress)

A package is designed to produce SDMs and SDM visualizations as part of the 2023 EFH 5-year review. It is designed to be moderately flexible and can be expanded in the future, but for now assumes that abundance prediction is the end goal. 

The most recent version of this package was built in R 4.4.1.

Authors:

@jeremyharris7  
@James-Thorson-NOAA  
@MasonSmith-NOAA
@MargaretSiple-NOAA  
@jodipirtle
@Ned-Laman-NOAA  


# Installation
`EFHSDM` can be installed using the following code:
```r
devtools::install_github("afsc-gap-products/akgfmaps", build_vignettes=TRUE)
devtools::install_github("alaska-groundfish-efh/EFHSDM@main", dependencies = TRUE, build_vignettes = FALSE)

# Development version - not currently recommended unless you are an active developer
# devtools::install_github("alaska-groundfish-efh/EFHSDM@dev", dependencies = TRUE, build_vignettes = FALSE)
```


## Dependencies
This project requires the following packages. Bear in mind that the `maxnet` package is still under development and changes occasionally. The `akgfmaps` package must be installed from the afsc-gap-products org page.

```r
# Packages
xtable, XML, raster, gstat, sf, stars, akgfmaps, ggplot, viridis, gridExtra, MASS, scales, labeling, maxnet, ENMeval, PresenceAbsence, mgcv
```

# Roadmap
Sections of the analysis are included as separate scripts. The general strategy is to use the functions provided in `Functions_Maxent.R` and `Functions_GamModel.R` to produce models, abundance rasters, effects estimates, and other outputs in a standard-ish format. Then, the scripts `Functions_LoadMap.R` and `Functions_Ensemble.R` provide more general methods for mapping or plotting model outputs and combining inferences from multiple models.

1) `Functions_Maxent.R` - this script provides the functions for quickly using maxnet models. 
2) `Functions_GamModel.R` - this script provides functions for conducting several operations with GAMs. The code has been tested with binomial, poisson, negative binomial, and ziplss hurdle models. 
3) `Functions_LoadMap.R` - this script provides versatile functions for mapping or plotting the model types mentioned above. It also contains some functions used to set up data specific to the EFH program that are unlikely to be broadly useful. 
4) `Functions_Ensemble.R` - this script provides functions designed to combine models and model products produced in the previous scripts
5) `Functions_Xtable.R` - this script produces a few standard output tables to summarize model or ensemble performance.
6) Functions_akgfmaps.R - this script provides standard functions to plot maps using the akgfmaps package that are of higher quality than those included elsewhere

Functions may be general or specific to a type of model

## GAM Specific Functions

| Function               |Use                                                                |
|:-----------------------|:------------------------------------------------------------------|
| `FitGAM()`             | Fit and term selection for binomial, poisson, and negbin GAMs |
| `FitHurdleGAM()`       | Fit and term selection for ziplss hurdle GAMs                 |
| `MakeGAMAbundance()`   | Produce a prediction raster from a GAM model and covariate stack |
| `GetGAMEffects()`      | Produce covariate effect estimates and CIs for GAM models        |
| `GAMStats()`           | Produce jacknife estimates of relative deviance explained (slow) |

* The functions `AutodetectGAMTerms()` and `AssembleGAMFormula()` are used internally by other
functions, but are not strictly necessary for all users.

## Maxnet Specific Functions
| Function               |Use                                                                |
|:-----------------------|:------------------------------------------------------------------|
| `FitMaxnet()`              | Fitting and term selection for maxnet models                      |
| `GetMaxnetEffects()`       | Produces covariate effect estimates and CIs for maxnet models     |
| `MakeMaxEntAbundance()`    | Produces a prediction raster from a GAM model and covariate stack |
| `MaxnetCoefs()`            | Returns the number of coefficients pertaining to each covariate   |
| `MaxnetStats()`            | Produces jackknife estimates of relative deviance explained       |


## General Functions
| Function               |Use                                                                |
|:-----------------------|:------------------------------------------------------------------|
| `RMSE()`                   | Calculated the RMSE for a set of observations and predictions     |
| `PDE()`                    | Calculated the deviance explained (Poisson) from obs and preds    |
| `CrossValidateModel()`     | Conducts k-fold cross-validatation of GAM or maxnet models        |
| `FindEFHbreaks()`          | Returns the EFH breaks for a given abundance map                  |
| `MakeVarianceRasters()`    | Produces raster of predicted variance based on CV folds (slow)    |

## Ensemble Functions
| Function       |Use                                                            |
|:-----------------------|:------------------------------------------------------------------|
| `MakeEnsemble()`           | Calculates the weights for a set a models based on RMSE           |
| `ValidateEnsemble()`       | Makes some basic validation plots and returns prediction data     |
| `MakeEnsembleAbundance()`  | Combines multiple rasters into a weighted average                 |
| `GetEnsembleVariance()`    | Estimates the variance in ensemble predictions                    |
| `GetEnsembleEffects()`     | Produces covariate effect estimates for an ensemble               |

There are numerous plotting functions. The ones based on `akgfmaps` are recommended.

## Workflow

First one will need to prepare any data and covariate rasters. Data should be organized in a data frame with columns for species and any covariates, offsets, etc. Rasters should be combined into a stack. 

The functions are typically called top to bottom. Begin by fitting a model using `FitGAM()`, `FitHurdleGAM()`, or `FitMaxnet()`. The resulting model is used with `MakeGAMAbundance()` (or `MakeMaxEntAbundance()`) to create an abundance raster, `GetGAMEffects()` to estimate covariate effects, and `GAMStats()` to obtain covariate contributions. Then one follows with `CrossValidateModel()` to output a useful dataframe of both in-bag and out-of-bag predictions. This dataframe can be used to calculate PDE and RMSE. `MakeVarianceRasters()` produces a variance map. `FindEFHbreaks()` gives the abundance thresholds that define each EFH quantile. The final EFH map is made by passing the result of `EFHbreaks()` to the `cut()` function from the terra package. 

## Simple example

All of the rasters used for the EFH 2023 Five-Year Review are stored on a shared drive at NOAA. For the purposes of the following example, the datasets are embedded in the package. For additional rasters and datasets, contact the package developers or submit a data [product request](https://github.com/alaska-groundfish-efh/product-requests/issues).

For example purposes, we will used only the last  years of data and only a few covariates. *Note that this means that the map you produce will look different from the final map produced in the 2023 EFH 5-year Review.*

### Load EFHSDM
```r
library(akgfmaps)
library(EFHSDM)
```

### Load the rasters

Currently, the raster for this example is stored with the package, and is just called `raster_stack`. You may eventually want to make your own raster stack and use it here. 
``` r
region.data <- subset(region_data_all, year >= 2012)
region.data$sponge <- as.integer(region.data$sponge > 0)
region.data$logarea <- log(region.data$area)


raster.stack <- terra::rast(terra::unwrap(raster_stack))
#raster.stack <- raster::raster(raster.stack) # option to turn this spatraster stack into a regular raster stack (only works if you don't crop )
```

### Next we will fit a basic Poisson model and generate an abundance map
``` r
gam.form <- formula("a_atf ~ s(lon,lat,bs = 'ds',m=c(1,.5), k=10) + s(bdepth, bs='tp',m=1,k=4) + s(btemp, bs='tp',m=1,k=4) + s(slope, bs='tp',m=1,k=4) + offset(logarea)")

# Note: to change the species/lifestage in this example, you can use any of the names in the columns of the region.data object:
head(names(region.data[-c(1:29,185)]))

poisson.model <- FitGAM(gam.formula = gam.form, data = region.data, family.gam = "poisson")
```

### Make the abundance map
``` r
poisson.abundance <- MakeGAMAbundance(model = poisson.model, r.stack = raster.stack)
abundance.plot <- MakeAKGFDensityplot(region = "goa", density.map = poisson.abundance, buffer = .98, title.name = "Adult ATF", legend.title = "Abundance")

# Display abundance plot. Note: this may take a minute to render!
print(abundance.plot)

png("AbundancePlot.png",width = 8,height = 3.5,units = 'in',res = 120)
abundance.plot
dev.off()
```
![Raster of ATF abundance](https://github.com/alaska-groundfish-efh/EFHSDM/blob/main/inst/readme_images/AbundancePlot.png)


### Get the covariate effects
``` r
poisson.effects <- GetGAMEffects(poisson.model, data = region.data)
Effectsplot(poisson.effects,region="goa",crs=raster.stack@crs)
```

### Crossvalidate the model
``` r
poisson.cv <- CrossValidateModel(model = poisson.model, data = region.data, folds = 10, model.type = "gam", key = "hauljoin")
poisson.preds <- poisson.cv[[1]]
head(poisson.preds)
RMSE(obs = poisson.preds$abund, pred = poisson.preds$cvpred)
PDE(obs = poisson.preds$abund, pred = poisson.preds$cvpred)

```

### And find the EFH
``` r
poisson.breaks <- FindEFHbreaks(poisson.abundance, method = "percentile")
poisson.efh <- terra::classify(poisson.abundance, poisson.breaks)
efh.plot <- MakeAKGFEFHplot(region = "goa", efh.map = poisson.efh, title.name = "Adult ATF", legend.title = "Percentiles")
```

### Make the EFH map 
``` r
# See note below; sometimes rendering this plot takes up too much memory in R and you have to save it to a file to see it.
print(efh.plot)

# Note: If the map does not appear in a Plot window in R, you can see it by writing the file to a .png:
png("EFHMap.png",width = 8,height = 3.5,units = 'in',res = 120)
efh.plot
dev.off()
```
![Raster of ATF EFH](https://github.com/alaska-groundfish-efh/EFHSDM/blob/main/inst/readme_images/EFHMap.png)

### Ensembles
An ensemble is not presented here, as it would take too long and be too complex to organize in this document. However, here is a pseudo-code illustration of how an ensemble can be constructed from previously constructed poisson and maxnet models.

First one performs all the steps shown above from both models. Then use the RMSE values to get the weights.

`MakeEnsemble(rmse = c(maxnet.rmse, poisson.rmse))`

Fit estimates of the ensemble can also be calculated, though the ensemble itself can not be easily crossvalidated:

`ValidateEnsemble(pred.list = list(maxnet.preds, poisson.preds), model.weights = ensemble.weights)`

Then make a new abundance raster, which is just the weighted average of the consistuent abundance rasters:

`MakeEnsembleAbundance(model.weights = ensemble.weights, abund.list = list(maxnet.abundance, poisson.abundance))`

EFH can be made in exactly the same way as for the consistuent models:

```r
FindEFHbreaks(ensemble.abundance, method = "percentile")
terra::cut(ensemble.abundance, ensemble.breaks)
```

## Legal disclaimer
>This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.
