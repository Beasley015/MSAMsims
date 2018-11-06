###############################
# Simulated abundance data for evaluating multi-species abundance model.
# No covariates.
# Includes data augmentation to model species that were never detected during sampling.
###############################

#install and load packages ####
library(vcdExtra)
library(vegan)

#Prelim data: sites, survey, seed -----------------------------------------------
J <- 30 #sites
K <- 3 #surveys per site
specs<-12 #Number of species

Ks<-rep(K, J) #Ks is a vector of length J indicationg # of sampling periods per site

#simulated values for covs would go here

set.seed(15) #ensures sim is same each time

#Simulating abundance data --------------------------------------------------
mean.lambdas <- rlogseries(specs, 0.75) 
#Draw lambdas from a logseries distribution

alpha0 <- log(mean.lambdas) #log-scale intercept
#abundance responses to any covariates would go here as alpha1, alpha2, etc

log.lambdas <- alpha0  #this is your log link function. add covs here
lambdas <- exp(log.lambdas)  #inverse link transformation

#most of these steps won't matter until covs are added

#create list of abundance vectors
nlist<-list()
for(a in 1:specs){
  nlist[[a]] <- rpois(n = J, lambda = lambdas[a])
}

ns<-do.call(rbind, nlist) #turn abundance vectors into abundance matrix
rowSums(ns) #total abundances

rotate.ns<-t(ns) #I might need an inverted matrix later

#Simulated observation process ------------------------------------------
some.det <- runif(n = specs-2, min = 0.4, max = 0.8)#simulate mean detection probs
#These are mid to high detection probabilities
no.det <- rep(0, 2)
#Two species will not be detected
mean.det <- c(some.det, no.det)

beta0<-qlogis(mean.det) #put it on logit scale

#responses to detection covs would go here

logit.p<-beta0 #logit link function. Det covs go here
p <- plogis(logit.p) #Transform it back

#Simulate observation data
L<-list()

for(b in 1:specs){
  y<-matrix(NA, ncol = K, nrow = J)
  for(a in 1:K){
    y[,a]<-rbinom(n = J, size = ns[b,], prob = p[b])
  }
  L[[b]]<-y
}
#I suppose I could make that a function but I'm lazy

#Smash it into array
obsdata<-array(as.numeric(unlist(L)), dim=c(J, K, specs-2)) 
#Nondetected species were removed from observation data
