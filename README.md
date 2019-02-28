# MSAMsims
Data simulations and testing multi-species abundance models

The R scripts and associated files herein contain simulated abundance data for testing the accuracy and some of the assumptions of multi-species abundance models (MSAMs). 

Currently, the following scripts are in the repository:
* MySim.SansCovs.R - The base model. No covariates or data augmentation.
  - abundsanscovs.txt - BUGS script for base model (Also works in JAGS)
  - modsanscovs3.rds - output of base model

* SansCovsAugmented.R - No covariates, but includes data augmentation to account for 1 nondetected species.
  - augmentsanscovs.txt - BUGS script (Also works in JAGS)
  - augsanscovs.rds - output of model with augmented dataset
  
* SimWithCovs.R - Model with covariates but without data augmentation
  - abundmitcovs.txt - BUGS script
  - covmodel.RDS - output of model with covariates

* WithCovsAugmented.R - Model with covariates and data augmentation for 1 nondetected species.
  - augmitcovs.txt - BUGS script
  - modaugcovs.RDS - output of model with covariates and data augmentation
