############################################################
# Simulated data for testing multi-species abundance model #
# With covariates, no data augmentation                    #
############################################################

#Load packages
library(vcdExtra)
library(TeachingDemos)
library(ggplot2)
library(R2OpenBUGS)

setwd("c:/users/beasley/dropbox/MSAMsims")

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
mean.lambdas <- rlogseries(n = specs, prob = 0.75)

#create responses to covariates
alpha0 <- log(mean.lambdas) #log-scale intercept

#All species have negative response
alpha1 <- rep(-1, specs)

#Only two species respond
alpha2 <- c(rep(2, 2), rnorm(n = specs-2, mean = 0, sd = 0.1))

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
mean.det <- runif(n = specs, min = 0.4, max = 0.9) #simulate mean detection probs
#These are mid to high detection values

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
mod1<-bugs(model.file = "abundmitcovs.txt", data = datalist, inits = init.values,
            parameters.to.save = params, n.chains = 3, n.burnin = 2000,
            n.iter = 6000, debug = T)

saveRDS(mod1, file = "covmodel.RDS")

mod1 <- readRDS(file = "covmodel.RDS")

#Test model accuracy -------------------------------------------------------
a1s <- mod1$sims.list$a1
mean.a1 <- apply(a1s, 2, mean)
lo.a1 <- apply(a1s, 2, quantile, 0.025)
hi.a1 <- apply(a1s, 2, quantile, 0.975)

ggplot()+
  geom_point(aes(x = 1:10, y = alpha1), color = "red")+
  geom_point(aes(x = 1:10, y = mean.a1))+
  geom_errorbar(aes(x = 1:10, ymin = lo.a1, ymax = hi.a1))

a2s <- mod1$sims.list$a2
mean.a2 <- apply(a2s, 2, mean)
lo.a2 <- apply(a2s, 2, quantile, 0.025)
hi.a2 <- apply(a2s, 2, quantile, 0.975)

ggplot()+
  geom_point(aes(x = 1:10, y = alpha2), color = "red")+
  geom_point(aes(x = 1:10, y = mean.a2))+
  geom_errorbar(aes(x = 1:10, ymin = lo.a2, ymax = hi.a2))

b1s <- mod1$sims.list$b1
mean.b1 <- apply(b1s, 2, mean)
lo.b1 <- apply(b1s, 2, quantile, 0.025)
hi.b1 <- apply(b1s, 2, quantile, 0.975)

ggplot()+
  geom_point(aes(x = 1:10, y = beta1), color = "red")+
  geom_point(aes(x = 1:10, y = mean.b1))+
  geom_errorbar(aes(x = 1:10, ymin = lo.b1, ymax = hi.b1))

#Test model assumptions ----------------------------------------------------
#Is estimated abundance closer to the true value than observed?
Z <- mod1$sims.list$Z
mean.Z <- apply(Z, c(2,3), mean)

site.mean <- apply(mean.Z, 1, sum)
site.true <- colSums(ns)
site.obs <- apply(maxobs, 1, sum)

site.comp <- data.frame(rank(site.true), site.mean, site.obs, site.true)
colnames(site.comp) <- c("sites.ranked", "site.mean", "site.obs", "site.true")

ggplot(data = site.comp, aes(x = sites.ranked))+
  geom_point(aes(y = site.mean, fill = "Estimated", color = "Estimated"))+
  geom_smooth(aes(y = site.mean, fill = "Estimated", color = "Estimated"),
              show.legend = F)+
  geom_point(aes(y = site.obs, fill = "Observed", color = "Observed"))+
  geom_smooth(aes(y = site.obs, fill = "Observed", color = "Observed"),
              show.legend = F)+
  geom_point(aes(y = site.true, fill = "True", color = "True"))+
  geom_smooth(aes(y = site.true, fill = "True", color = "True"), show.legend= F)+
  scale_fill_manual(breaks = c("Estimated", "Observed", "True"),
                    values = c("blue", "black", "red"), guide = F)+
  scale_color_manual(breaks = c("Estimated", "Observed", "True"),
                     values = c("blue", "black", "red"))+
  labs(x = "Sites (Ranked)", y = "Abundance")+
  theme_bw()+
  theme(legend.title = element_blank())

#Is there a 1:1 relationship between estimated and true abundance?
linears <- data.frame(as.vector(t(ns)),as.vector(apply(Z, c(2,3), mean)))
colnames(linears) <- c("True", "Estimated")

ggplot(data = linears, aes(x = True, y = Estimated))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_bw()

summary(lm(data = linears, Estimated~True))
