############################################################
# Simulated data for testing multi-species abundance model #
# With covariates, no data augmentation                    #
############################################################

#Load packages
library(vcdExtra)
library(tidyverse)
library(R2OpenBUGS)
library(agricolae)

setwd("c:/users/beasley/documents/MSAMsims")

#Global variables: sites, survey, seed ----------------------------------------
J <- 30 #sites
K <- 3 #surveys per site
specs<-10 #Number of species

Ks<-rep(K, J) #Ks is a vector of length J indicationg # of sampling periods per site

set.seed(15) #insures sim is same each time

#create abundance covariates
cov1<-rnorm(n = J, mean = 10, sd = 5)
cov1scale <- as.vector(scale(cov1))

cov2 <- rnorm(n = J, mean = 5, sd = 1)
cov2scale <- as.vector(scale(cov2))

#detection
detcov <- matrix(runif(J*K, 0, 10), nrow = J, ncol = K)

#Simulating abundance data ---------------------------------------------------
mean.lambdas <- runif(n = 10, min = 0.25, max = 2.5)

#create responses to covariates
alpha0 <- log(mean.lambdas) #log-scale intercept

#All species have negative response
alpha1 <- rep(-1, specs)

#Only two species respond
alpha2 <- c(rep(1, 2), rnorm(n = specs-2, mean = 0, sd = 0.1))

#Put covs and responses in log link function
log.lambdas <- matrix(NA, nrow = specs, ncol = J)

for(i in 1:specs){
  log.lambdas[i,] <- alpha0[i] + alpha1[i]*cov1scale + alpha2[i]*cov2scale
}

#inverse link transformation
lambdas <- exp(log.lambdas)

#create list of abundance vectors
ns <- matrix(NA, nrow = specs, ncol = J)
for(a in 1:specs){
  ns[a,] <- rpois(n = J, lambda = lambdas[a,])
}

rowSums(ns) #total abundances

rotate.ns<-t(ns)

#Simulated observation process --------------------------------------------

#Wide range of detection values
# mean.det <- runif(n = specs, min = 0.1, max = 0.9) #simulate mean detection probs

#Low detection values
# mean.det <- runif(n = specs, min = 0.1, max = 0.5)

#Mid detection values 
mean.det <- runif(n = specs, min = 0.3, max = 0.7)

#High detection values
# mean.det <- runif(n = specs, min = 0.5, max = 0.9)


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
obsdata<-array(as.numeric(unlist(L)), dim=c(J, K, specs))

#Convert obsdata into abundance matrix
maxobs<-apply(obsdata, c(1,3), max)

#Run the model ------------------------------------------------------------
#Write model
cat("
    model{
    
    #Define prior distributions
    
    #Intercepts
    a0.mean ~ dnorm(0,0.001)
    sigma.a0 ~ dunif(0,10)
    tau.a0 <- 1/(sigma.a0*sigma.a0)
    
    b0.mean ~ dnorm(0,0.001)
    sigma.b0 ~ dunif(0,10)
    tau.b0 <- 1/(sigma.b0*sigma.b0)

    #Covs
    a1.mean ~ dnorm(0, 0.1)
    sigma.a1 ~ dunif(0,10)
    tau.a1 <- 1/(sigma.a1*sigma.a1)

    a2.mean ~ dnorm(0, 0.1)
    sigma.a2 ~ dunif(0,10)
    tau.a2 <- 1/(sigma.a2*sigma.a2)

    b1.mean ~ dnorm(0, 0.1)
    sigma.b1 ~ dunif(0,10)
    tau.b1 <- 1/(sigma.b1*sigma.b1)
    
    for(i in 1:specs){
    #create priors from distributions above
    a0[i] ~ dnorm(a0.mean, tau.a0)
    a1[i] ~ dnorm(a1.mean, tau.a1)
    a2[i] ~ dnorm(a2.mean, tau.a2)
    
    b0[i] ~ dnorm(b0.mean, tau.b0)
    b1[i] ~ dnorm(b1.mean, tau.b1)
    
    #Loop within a loop to estimate abund of spec i at site j
    for(j in 1:J){
    log(lambda[j,i]) <- a0[i] + a1[i]*cov1[j] + a2[i]*cov2[j]
    Z[j,i] ~ dpois(lambda[j,i])
    
    #Loop within loops for estimating det of spec i at site j at time k
    for(k in 1:K[j]){
    p[j,k,i] <- b0[i] + b1[i]*detcov[j,k]
    logit.p[j,k,i] <- 1 / (1 + exp(-p[j,k,i]))
    obsdata[j,k,i] ~ dbin(logit.p[j,k,i], Z[j,i])
    }
    }
    }
    }
    ", file = "abundmitcovs.txt")

#Compile data into list
datalist<-list(specs=specs, J=J, K=Ks, obsdata=obsdata, cov1 = cov1scale, 
               cov2 = cov2scale, detcov = detcov)

#Specify parameters
params<-list('Z','lambda','a0','b0', 'a0.mean', 'b0.mean', 'a1', 'a2', 'b1')

#Specify initial values
init.values<-function(){
  list(a0 = rnorm(n = specs, mean = mean(alpha0)),
       b0 = rnorm(n = specs, mean = mean(beta0)),
       a1 = rnorm(n = specs), a2 = rnorm(n = specs), b1 = rnorm(n = specs),
       Z = matrix(maxobs, nrow = J, ncol = specs)
  )
}

#Send model to Gibbs sampler 
# mod1<-bugs(model.file = "abundmitcovs.txt", data = datalist, inits = init.values,
#             parameters.to.save = params, n.chains = 3, n.burnin = 6000,
#             n.iter = 12000, debug = T)

# saveRDS(mod1, file = "ranged_det.RDS")
# saveRDS(mod1, file = "lo_det.RDS")
# saveRDS(mod1, file = "mid_det.RDS")
# saveRDS(mod1, file = "hi_det.RDS")

rangedmod <- readRDS(file = "ranged_det.RDS")
lomod <- readRDS(file = "lo_det.RDS")
midmod <- readRDS(file = "mid_det.RDS")
himod <- readRDS(file = "hi_det.RDS")

# Model evaluation: abundance estimates --------------------------------------
#Get Z matrices from model objects
get.zs <- function(x){
  y <- x$sims.list$Z
  return(y)
}

rangeZ <- get.zs(x = rangedmod)
loZ <- get.zs(x = lomod)
midZ <- get.zs(x = midmod)
hiZ <- get.zs(x = himod)

#Get estimated site-level abundances
site.abunds <- function(x){
  y <- apply(x, c(2,3), mean)
  y2 <- rowSums(y)
  return(y2)
}

range.abund <- site.abunds(rangeZ)
lo.abund <- site.abunds(loZ)
mid.abund <- site.abunds(midZ)
hi.abund <- site.abunds(hiZ)

#Put it in a data frame
abunds <- data.frame(True = colSums(ns), Rank = rank(colSums(ns)), Varied = range.abund,
                     Low = lo.abund, Mid = mid.abund, High = hi.abund, 
                     Observed = rowSums(maxobs))

#Plot it
ggplot(data = abunds, aes(x = Rank))+
  geom_point(aes(y = True, color = "True"))+
  geom_smooth(aes(y = True, color = "True", fill = "True"), alpha = 0.2)+
  geom_point(aes(y = Varied, color = "Varied"))+
  geom_smooth(aes(y = Varied, color = "Varied", fill = "Varied"), alpha = 0.2)+
  geom_point(aes(y = Low, color = "Low"))+
  geom_smooth(aes(y = Low, color = "Low", fill = "Low"), alpha = 0.2)+
  geom_point(aes(y = Mid, color = "Mid"))+
  geom_smooth(aes(y = Mid, color = "Mid", fill = "Mid"), alpha = 0.2)+
  geom_point(aes(y = High, color = "High"))+
  geom_smooth(aes(y = High, color = "High", fill = "High"), alpha = 0.2)+
  geom_point(aes(y = Observed, color = "Observed"))+
  geom_smooth(aes(y = Observed, color = "Observed"), alpha = 0.2)+
  scale_color_viridis_d()+
  scale_fill_viridis_d(guide = F)+
  labs(x = "Sites (Ranked)", y = "Abundance")+
  theme_bw()+
  theme(legend.title = element_blank())

#Are there differences between estimates?
abund.melt <- gather(abunds, c(True, Varied:Observed), key = "Source", value = "Abund")

abund.aov <- aov(data = abund.melt, Abund ~ Source + Rank)
summary(abund.aov) #Mostly varies between sites but det plays a role

#Which estimates are the same?
abund.diffs <- HSD.test(y = abund.aov, trt = "Source")
abund.diffs$groups #Observed is different than true, all others not sig. different

# Model evaluation: covariate significance -----------------------------------
get.a1 <- function(x){
  y <- x$sims.list$a1
  y2 <- as.data.frame(y)
  colnames(y2) <- c(1:10)
  return(y2)
}

range.a1 <- get.a1(rangedmod)
lo.a1 <- get.a1(lomod)
mid.a1 <- get.a1(midmod)
hi.a1 <- get.a1(himod)

long.a1 <- function(x, source){
  y <- gather(x, key = "Species", value = "a1")
  y$Source = source
  return(y)
}

range.a1.long <- long.a1(x = range.a1, source = "Varied")
lo.a1.long <- long.a1(x = lo.a1, source = "Low")
mid.a1.long <- long.a1(x = mid.a1, source = "Mid")
hi.a1.long <- long.a1(x = hi.a1, source = "High")

all.a1 <- as.data.frame(rbind(range.a1.long, lo.a1.long, mid.a1.long, hi.a1.long))

all.a1 %>%
  group_by(Source, Species)%>%
  summarise(mean = mean(a1), quant025 = quantile(a1, 0.025), 
                    quant975 = quantile(a1, 0.975)) %>%
                    {. ->> sum.a1}

ggplot(data = sum.a1, aes(x = Species, color = Source))+
  geom_point(aes(y = mean))+
  geom_errorbar(aes(ymin = quant025, ymax = quant975))+
  scale_color_viridis_d()+
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  geom_hline(aes(yintercept = -1), linetype = "dashed", color = "red")+
  labs(y = "Alpha 1")+
  theme_bw()+
  theme(axis.title.x = element_blank())
  
# Model evaluation: relationship between true and estimated data -------------


