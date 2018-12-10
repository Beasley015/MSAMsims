# MSAMsims
Data simulations and testing multi-species abundance models

The R scripts and associated files herein contain simulated abundance data for testing the accuracy and some of the assumptions of multi-species abundance models (MSAMs). 

Currently, the following scripts are in the repository:
* MySim.SansCovs.R - The base model. No covariates or data augmentation.
  - abundsanscovs.txt - BUGS script for base model
  - modsanscovs3.rds - output of base model

* SansCovsAugmented.R - No covariates, but includes data augmentation to account for 2 nondetected species.
  - augmentsanscovs.txt - BUGS script
  - augsanscovs.rds - output of model with augmented dataset
  
* SimWithCovs.R - Currently non-functioning script for a model with covariates but without data augmentation
