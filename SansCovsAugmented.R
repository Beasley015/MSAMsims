###############################
# Simulated abundance data for evaluating multi-species abundance model.
# No covariates.
# Includes data augmentation to model species that were never detected 
# during sampling.
###############################

#Setup -------------------------------------------------------------------------
library(vcdExtra)
library(vegan)
library(R2OpenBUGS)
library(abind)
library(tidyverse)

setwd("c:/users/beasley/dropbox/MSAMsims")

set.seed(15) #ensures sim is same each time

#Prelim data: sites, survey, seed -----------------------------------------------
J <- 30 #sites
K <- 3 #surveys per site
specs<-11 #Number of species

Ks<-rep(K, J) 
#Ks is a vector of length J indicationg # of sampling periods per site

#simulated values for covs would go here

#Simulating abundance data --------------------------------------------------
mean.lambdas <- rlogseries(specs, 0.75) 
mean.lambdas[11] <- 0.5 #nondetected species should have low abundance for realism
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
some.det <- runif(n = specs-1, min = 0.4, max = 0.8)#simulate mean detection probs
#These are mid to high detection probs
no.det <- 0
#One species will not be detected
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
obsdata<-array(as.numeric(unlist(L)), dim=c(J, K, specs-1)) 
#Nondetected species were removed from observation data

#Number of observed species
n <- length(some.det)

#Augment data with all-zero matrices
n.aug <- 2
augmats <- array(0, dim = c(J, K, n.aug))

augdata <- abind(obsdata, augmats, along = 3)

maxobs <- apply(augdata, c(1,3), max)

#Write model and send to Gibbs sampler ------------------------------------------
cat("
    model{
    
    #Define hyperprior distributions

    omega ~ dunif(0,1)
    
    a0.mean ~ dunif(0,1)
    mu.a0 <- log(a0.mean)-log(1-a0.mean)
    tau.a0 ~ dgamma(0.1, 0.1)
    
    b0.mean ~ dunif(0,1)
    mu.b0 <- log(b0.mean)-log(1-b0.mean)
    tau.b0 ~ dgamma(0.1, 0.1)
    
    for(i in 1:(n+n.aug)){
    #create priors from distributions above
    w[i] ~ dbern(omega)
    #w[i] indicates whether or not species is exposed to sampling

    a0[i] ~ dnorm(mu.a0, tau.a0)
    
    b0[i] ~ dnorm(mu.b0, tau.b0)
    
    #Loop within a loop to estimate abund of spec i at site j
    for(j in 1:J){
    lambda[j,i] <- exp(a0[i])
    mu.lambda[j,i] <- lambda[j,i]*w[i]
    Z[j,i] ~ dpois(mu.lambda[j,i])
    #Z is the estimated abundance matrix
    
    #Loop within loops for estimating det of spec i at site j at time k
    for(k in 1:K[j]){
    p[j,k,i] <- b0[i]
    logit.p[j,k,i] <- 1 / (1 + exp(-p[j,k,i]))
    obsdata[j,k,i] ~ dbin(logit.p[j,k,i], Z[j,i])
    }
    }
    }

    #Estimate total richness (N) by adding observed (n) and unobserved (n0) species
    n0<-sum(w[(n+1):(n+n.aug)])
    N<-n+n0

    }
    ", file = "augmentsanscovs.txt")

#Compile data
datalist<-list(n = n, n.aug = n.aug, J=J, K=Ks, obsdata=augdata)

#Specify parameters to return to R
params<-list('Z','lambda','a0','b0', 'mu.a0', 'mu.b0', 'tau.a0', 'tau.b0', 'N')

#Generate initial values
init.values<-function(){
  omega.guess <- runif(1,0,1)
  lambda.guess <- runif(1,0,5)
  list(omega = omega.guess, 
       w=c(rep(1,n), rbinom(n = n.aug,size=1,prob=omega.guess)),
       a0 = rnorm(n = (n+n.aug), mean = mean(alpha0)),
       b0 = rnorm(n = (n+n.aug), mean = runif(1,0,1)),
       Z = maxobs
  )
}

# augmodel <- bugs(model.file = "augmentsanscovs.txt", data = datalist, 
#                  n.chains = 3, parameters.to.save = params, inits = init.values,
#                  n.burnin = 5000, n.iter = 8000, debug = T)
# saveRDS(augmodel, file = "augsanscovs.RDS")

augmodel <- readRDS(file = "augsanscovs.RDS")

# Model evaluation: regional and site-level richness ----------------------------
Ns <- augmodel$sims.list$N
mean(Ns); quantile(Ns, c(0.025, 0.25, 0.75, 0.975))

ggplot()+
  geom_histogram(aes(x = Ns), binwidth = 1)

# Model evaluation: posterior distributions -------------------------------------
post.abund <- augmodel$sims.list$mu.a0
hist(post.abund)

post.det <- augmodel$sims.list$mu.b0
hist(post.det)
#detection isn't bimodal, so that's good

# Compare observed, estimated, and true abundance ------------------------------
Zs <- augmodel$sims.list$Z

#Species
spec.abunds <- apply(Zs, c(1,3), sum)
spec.mean <- apply(spec.abunds[,-12], 2, mean)
spec.obs <- apply(t(maxobs[,-12]), 1, sum)
spec.true <- rowSums(ns)

spec.comp <- data.frame(Rank = rank(spec.true), Estimated = spec.mean, 
                        Observed = spec.obs, True = spec.true)

ggplot(data = spec.comp, aes(x = Rank))+
  geom_point(aes(y = Estimated, color = "Estimated"))+
  geom_smooth(aes(y = Estimated, fill = "Estimated", color = "Estimated"),
              show.legend = F, alpha = 0.25)+
  geom_point(aes(y = Observed, color = "Observed"))+
  geom_smooth(aes(y = Observed, fill = "Observed", color = "Observed"), alpha = 0.25)+
  geom_point(aes(y = True, color = "True"))+
  geom_smooth(aes(y = True, fill = "True", color = "True"), alpha = 0.25)+
  scale_fill_manual(breaks = c("Estimated", "Observed", "True"),
                    values = c("blue", "black", "red"), guide = F)+
  scale_color_manual(breaks = c("Estimated", "Observed", "True"),
                     values = c("blue", "black", "red"))+
  theme_bw()

#Sites
site.abund <- apply(Zs, c(1,2), sum)
site.mean <- apply(site.abund, 2, mean)
site.true <- colSums(ns)
site.obs <- apply(t(maxobs), 2, sum)

site.comp <- data.frame(Rank = rank(site.true), True = site.true, 
                        Observed = site.obs, Estimated = site.mean)

ggplot(data = site.comp, aes(x = Rank))+
  geom_point(aes(y = True, color = "True"))+
  geom_smooth(aes(y = True, color = "True", fill = "True"))+
  geom_point(aes(y = Observed, color = "Observed"))+
  geom_smooth(aes(y = Observed, color = "Observed", fill = "Observed"))+
  geom_point(aes(y = Estimated, color = "Estimated"))+
  geom_smooth(aes(y = Estimated, color = "Estimated", fill = "Estimated"))+
  scale_color_manual(breaks = c("True", "Estimated", "Observed"), values = 
                       c("blue", "black", "red"))+
  scale_fill_manual(breaks = c("True", "Estimated", "Observed"), values = 
                      c("blue", "black", "red"), guide = F)+
  theme_bw()

#Make sure estimates aren't biased at high or low abundances --------------------
#By species
specmod <- lm(data = spec.comp, Estimated~True)

ggplot(data = spec.comp, aes(x = True, y = Estimated))+
  geom_point()+
  geom_smooth()
#Not bad

#By site
sitemod <- lm(data = site.comp, Estimated~True)

ggplot(data = site.comp, aes(x = True, y = Estimated))+
  geom_point()+
  geom_smooth()

#It looks like I have a tradeoff when it comes to det probabilities
#If they're low, the model overestimates abundance (still closer than obs)
#If they're high, the model won't include the nondetected species