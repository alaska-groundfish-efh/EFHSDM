# EFH_SDM

# Statement of purpose
This package is designed to produce SDMs and SDM visualizations. It is designed to be moderately flexible, but assumes that abundance prediction is the end goal. 

## Roadmap
The package is broken down in several sections indicated by different files. The general strategy is to use the functions provided in Functions_Maxent and Functions_GamModel to produce models, abundance rasters, effects estimates, and other outputs in a standard-ish format. Then, the scripts Functions_LoadMap and Functions_Ensemble provide more general methods for mapping or plotting model outputs and combining inferences from multiple models.

The Meatgrinder script provides an example of a workflow that can be used to combine the many functions to make a SDM ensemble and the accompanying maps. However, it contains a great deal of control logic designed to accommodate the needs of the 2022 EFH 5-year Review process which involved running and keeping track of over 200 individual species/lifestages.

In general, the workflow for making an SDM begins by setting up the data, which will vary by project. Then one use a "Fit" function to make the model. These functions are wrappers that automate some common tasks such as term selection related to fitting said model. The generic CrossValidateModel function provides estimates of model performance. "GetEffects" functions estimate the covariate effects of the model. "MakeAbundance" functions map predictions.

1) Functions_Maxent.R - this script provides the functions for quickly using maxnet models. 
2) Functions_GamModel.R - this script provides functions for conducting several operations with GAMs. The code has been tested with binomial, poisson, negative binomial, and ziplss hurdle models. 
3) Functions_LoadMap.R - this script provides versatile functions for mapping or plotting the model types mentioned above. It also contains some functions used to set up data specific to the EFH program that are unlikely to be broadly useful. 
4) Functions_Ensemble.R - this script provides functions designed to combine models and model products produced in the previous scripts
5) Functions_Xtable.R - this script produces a few standard output tables to summarize model or ensemble performance.
6) Functions_akgfmaps.R - (under construction) this script provides standard functions to plot maps using the akgfmaps package that are of higher quality than functions included elsewhere




## Legal disclaimer
This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.
