############################################################
# Simulated data for testing multi-species abundance model #
# With occupancy and detection covariate                   #
# Includes data augmentation for nondetected species       #
############################################################

#Load packages
library(vcdExtra)
library(TeachingDemos)
library(ggplot2)
library(R2OpenBUGS)
library(abind)

setwd("c:/users/beasley/dropbox/MSAMsims")

#Prelim data: sites, survey, seed -----------------------------------------------
set.seed(15)

J <- 30 #sites
K <- 3 #surveys per site
specs<-11 #Number of species

Ks<-rep(K, J) 
#Ks is a vector of length J indicationg # of sampling periods per site

#Site covariate (abundance)
cov1<-rnorm(n = J, mean = 10, sd = 3)
cov1scale <- as.vector(scale(cov1))

#Survey covariate (detection)
detcov <- matrix(runif(J*K, 0, 10), nrow = J, ncol = K)

#Simulating abundance data --------------------------------------------------
mean.lambdas <- rlogseries(specs, 0.75) #Draw lambdas from a logseries distribution
mean.lambdas[11] <- 0.25 #nondetected species should have low abundance for realism


#log-scale intercept
alpha0 <- log(mean.lambdas)
#response to site covariate: all species have positive response
alpha1 <- alpha1 <- rep(1, specs)

#Get new lambdas
log.lambdas <- matrix(NA, nrow = specs, ncol = J)
for(i in 1:specs){
  log.lambdas[i,] <- alpha0[i] + alpha1[i]*cov1scale #log link function
}

#inverse link transformation
lambdas <- exp(log.lambdas)

#create list of abundance vectors
ns <- matrix(NA, nrow = specs, ncol = J)
for(a in 1:specs){
  ns[a,] <- rpois(n = J, lambda = lambdas[a,])
}

rowSums(ns) #total abundances

#Simulated observation process --------------------------------------------
some.det <- runif(n = specs - 1, min = 0.2, max = 0.6) #simulate detection probs
#These are low to mid detection values
no.det <- 0 #one species was not detected

mean.det <- c(some.det, no.det)

#Detection intercept and cov responses
beta0<-qlogis(mean.det) #put it on logit scale

#Detection cov does not affect detectability
beta1 <- rnorm(n = specs, mean = 0, sd = 0.01)

#Logit link function
logit.p <- array(NA, dim = c(J, K, specs))
for(i in 1:specs){
  logit.p[,,i] <- beta0[i] + beta1[i]*detcov
}

p <- plogis(logit.p)

#Simulate observation data
L<-list()

for(b in 1:specs){
  y<-matrix(NA, ncol = K, nrow = J)
  for(a in 1:K){
    y[,a]<-rbinom(n = J, size = ns[b,], prob = p[,,b])
  }
  L[[b]]<-y
}

#Smash it into array
obsdata<-array(as.numeric(unlist(L)), dim=c(J, K, specs-1))

#Sanity check: convert obs data into abundance matrix
maxobs<-apply(obsdata, c(1,3), max)
print(maxobs)
#looks fine

#Number of observed species
n <- length(some.det)

#Augment data with all-zero matrices
n.aug <- 2
augmats <- array(0, dim = c(J, K, n.aug))

augdata <- abind(obsdata, augmats, along = 3)

maxobs <- apply(augdata, c(1,3), max)
