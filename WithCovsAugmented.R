############################################################
# Simulated data for testing multi-species abundance model #
# With occupancy and detection covariate                   #
# Includes data augmentation for nondetected species       #
############################################################

#Load packages
library(vcdExtra)
library(R2OpenBUGS)
library(abind)
library(tidyverse)
library(vegan)

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
alpha1 <- rnorm(n = specs, mean = 0, sd = 0.25) + 1

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

nbyspec <- rowSums(ns) #total abundances

#Simulated observation process --------------------------------------------
some.det <- runif(n = specs, min = 0, max = 0.7) #simulate detection probs
#Wide range in detection
mean.det <- some.det
mean.det[min(nbyspec)] <- 0 #Least abundant species will not be detected

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

#Sanity check: convert obs data into abundance matrix
maxobs<-apply(obsdata, c(1,3), max)
print(maxobs)
#looks fine

#Remove non-detected species
nondet <- which(apply(obsdata, 3, sum) == 0)

new.obs <- obsdata[,,-nondet]

#Number of observed species
n <- specs - 1

#Augment data with all-zero matrices
n.aug <- 1
augmats <- array(0, dim = c(J, K, n.aug))

augdata <- abind(new.obs, augmats, along = 3)

maxobs <- apply(augdata, c(1,3), max)

# Run the model -----------------------------------------------------------
#Write model
cat("
    model{
    
    #Define hyperprior distributions

    omega ~ dunif(0,1)
    
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

    b1.mean ~ dnorm(0, 0.1)
    sigma.b1 ~ dunif(0,10)
    tau.b1 <- 1/(sigma.b1*sigma.b1)
    
    for(i in 1:(n + n.aug)){
    #create priors from distributions above

    w[i] ~ dbern(omega)
    #w[i] indicates whether or not species is exposed to sampling

    a0[i] ~ dnorm(a0.mean, tau.a0)
    a1[i] ~ dnorm(a1.mean, tau.a1)

    b0[i] ~ dnorm(b0.mean, tau.b0)
    b1[i] ~ dnorm(b1.mean, tau.b1)
    
    #Loop within a loop to estimate abund of spec i at site j
    for(j in 1:J){
    log(lambda[j,i]) <- a0[i] + a1[i]*cov1[j]
    mu.lambda[j,i] <- lambda[j,i]*w[i]
    Z[j,i] ~ dpois(mu.lambda[j,i])
    
    #Loop within loops for estimating det of spec i at site j at time k
    for(k in 1:K[j]){
    p[j,k,i] <- b0[i] + b1[i]*detcov[j,k]
    logit.p[j,k,i] <- 1 / (1 + exp(-p[j,k,i]))
    obsdata[j,k,i] ~ dbin(logit.p[j,k,i], Z[j,i])
    }
    }
    }

    #Estimate total richness N by adding observed (n) and unobserved (n0) species
    n0<-sum(w[(n+1):(n+n.aug)])
    N<-n+n0

    }
    ", file = "augmitcovs.txt")

#Compile data into list
datalist<-list(J=J, K=Ks, obsdata = augdata, n = n, n.aug = n.aug, 
               cov1 = cov1scale, detcov = detcov)

#Specify parameters
params<-list('Z','lambda','a0','b0', 'a0.mean', 'b0.mean', 'a1', 'b1', 'N')

#Specify initial values
init.values<-function(){
  omega.guess <- runif(1,0,1)
  lambda.guess <- runif(1,0,5)
  list(omega = omega.guess, 
       w=c(rep(1,n), rbinom(n = n.aug, size=1, prob=omega.guess)),
       a0 = rnorm(n = (n+n.aug), mean = mean(alpha0)),
       a1 = rnorm(n = (n+n.aug), mean = mean(alpha1)),
       b0 = rnorm(n = (n+n.aug), mean = runif(1,0,1)),
       b1 = rnorm(n = (n+n.aug), mean = runif(1,0,1)),
       Z = maxobs
  )
}

model <- bugs(model.file = "augmitcovs.txt", data = datalist, n.chains = 3,
              parameters.to.save = params, inits = init.values, n.burnin = 7500,
              n.iter = 15000, debug = T)

saveRDS(model, file = "modaugcovs.rds")

model <- readRDS(file = "modaugcovs.rds")

# Model evaluation: abundance/richness/diversity -----------------------------
#Regional richness
Ns <- model$sims.list$N
mean(Ns); quantile(Ns, c(0.025, 0.25, 0.75, 0.975))

ggplot()+
  geom_histogram(aes(x = Ns), binwidth = 1)

#Site-level richness
#Extract abundance matrices
Zs <- model$sims.list$Z

#Convert to occurrence matrices
to.occ <- function(x){
  y <- x
  y[which(y > 1)] <- 1
  return(y)
}

Zocc <- to.occ(x = Zs)
n.occ <- to.occ(x = ns)
obs.occ <- to.occ(x = maxobs)

#Get mean site-level richness
Zrich <- apply(Zocc, c(1,2), sum)
mean.rich <- apply(Zrich, 2, mean)

n.rich <- apply(n.occ, 2, sum)

obs.rich <- apply(obs.occ, 1, sum)

#Create data frame and plot
richness <- data.frame(True = n.rich, Estimated = mean.rich, Observed = obs.rich,
                       Sites = rank(n.rich))

ggplot(data = richness, aes(x = Sites))+
  geom_point(aes(y = True, color = "True"))+
  geom_smooth(aes(y = True, color = "True", fill = "True"))+
  geom_point(aes(y = Estimated, color = "Estimated"))+
  geom_smooth(aes(y = Estimated, color = "Estimated", fill = "Estimated"))+
  geom_point(aes(y = Observed, color = "Observed"))+
  geom_smooth(aes(y = Observed, color = "Observed", fill = "Observed"))+
  scale_color_manual(breaks = c("True", "Estimated", "Observed"), 
                    values = c("black", "blue", "red"))+
  scale_fill_manual(breaks = c("True", "Estimated", "Observed"), 
                    values = c("black", "blue", "red"), guide = F)+
  theme_bw()

#Are there differences in richness?
richness.melt <- gather(richness, True:Observed, key = Source, 
                        value = Richness)

rich.aov <- aov(data = richness.melt, Richness ~ Source + Sites)
summary(rich.aov) #Yep, there are differences
richness.groups <- TukeyHSD(rich.aov, which = "Source")
plot(richness.groups)
#Observed is different than estimated and true richness values

#Abundance
est.abund <- apply(Zs, c(1,2), sum)
mean.abund <- apply(est.abund, 2, mean)

true.abund <- colSums(ns)

obs.abund <- rowSums(maxobs)

abund <- data.frame(Site = rank(true.abund), True = true.abund, 
                    Estimated = mean.abund, Observed = obs.abund)

ggplot(data = abund, aes(x = Site))+
  geom_point(aes(y = True, color = "True"))+
  geom_smooth(aes(y = True, color = "True", fill = "True"))+
  geom_point(aes(y = Estimated, color = "Estimated"))+
  geom_smooth(aes(y = Estimated, color = "Estimated", fill = "Estimated"))+
  geom_point(aes(y = Observed, color = "Observed"))+
  geom_smooth(aes(y = Observed, color = "Observed", fill = "Observed"))+
  scale_color_manual(breaks = c("True", "Estimated", "Observed"), 
                     values = c("black", "blue", "red"))+
  scale_fill_manual(breaks = c("True", "Estimated", "Observed"), 
                    values = c("black", "blue", "red"), guide = F)+
  theme_bw()

#Differences in abundance?
abundance.melt <- gather(abund, True:Observed, key = Source, 
                        value = Abundance)

abund.aov <- aov(data = abundance.melt, Abundance ~ Source + Site)
summary(abund.aov) #Yep, there are differences
abund.groups <- TukeyHSD(abund.aov, which = "Source")
plot(abund.groups)

#Diversity
divers <- function(x){
  y <- diversity(x, index = "shannon", MARGIN = 1)
  return(y)
}

Zmat <- apply(Zs, c(2,3), mean)

Zshan <- divers(x = Zmat)
trueshan <- divers(x = t(ns))
obsshan <- divers(x = maxobs)

divers.frame <- data.frame(Site = rank(trueshan), True = trueshan, 
                           Estimated = Zshan, Observed = obsshan)

ggplot(data = divers.frame, aes(x = Site))+
  geom_point(aes(y = True, color = "True"))+
  geom_smooth(aes(y = True, color = "True", fill = "True"))+
  geom_point(aes(y = Estimated, color = "Estimated"))+
  geom_smooth(aes(y = Estimated, color = "Estimated", fill = "Estimated"))+
  geom_point(aes(y = Observed, color = "Observed"))+
  geom_smooth(aes(y = Observed, color = "Observed", fill = "Observed"))+
  scale_color_manual(breaks = c("True", "Estimated", "Observed"), 
                     values = c("black", "blue", "red"))+
  scale_fill_manual(breaks = c("True", "Estimated", "Observed"), 
                    values = c("black", "blue", "red"), guide = F)+
  theme_bw()

# Differences in diversity?
divers.melt <- gather(divers.frame, True:Observed, key = Source, 
                         value = Shannon)

divers.aov <- aov(data = divers.melt, Shannon ~ Source + Site)
summary(divers.aov) #Yep, there are differences
divers.groups <- TukeyHSD(divers.aov, which = "Source")
plot(divers.groups)

# Model evaluation: covariate estimates ----------------------------------------
#Abundance covariate
alpha1.est <- model$sims.list$a1

alpha1.mean <- apply(alpha1.est, 2, mean)
alpha1.lo <- apply(alpha1.est, 2, quantile, 0.025)
alpha1.hi <- apply(alpha1.est, 2, quantile, 0.975)
 
alpha1s <- data.frame(Mean = alpha1.mean, Lo = alpha1.lo, Hi = alpha1.hi)

ggplot(data = alpha1s, aes(x = seq(1:specs)))+
  geom_point(aes(y = Mean))+
  geom_point(aes(y = alpha1), color = "red")+
  geom_errorbar(aes(ymin = Lo, ymax = Hi))+
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  theme_bw()

#Detection covariate
beta1.est <- model$sims.list$b1

beta1.mean <- apply(beta1.est, 2, mean)
beta1.lo <- apply(beta1.est, 2, quantile, 0.025)
beta1.hi <- apply(beta1.est, 2, quantile, 0.975)

beta1s <- data.frame(Mean = beta1.mean, Lo = beta1.lo, Hi = beta1.hi)

ggplot(data = beta1s, aes(x = seq(1:specs)))+
  geom_point(aes(y = Mean))+
  geom_point(aes(y = beta1), color = "red")+
  geom_errorbar(aes(ymin = Lo, ymax = Hi))+
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  theme_bw()
