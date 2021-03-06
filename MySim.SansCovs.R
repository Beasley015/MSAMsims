##########################################
# Simulated abundance data for evaluating multi-species abundance model.
# No covariates, no data augmentation.
##########################################

#install and load packages ####
library(vcdExtra)
library(vegan)

#Prelim data: sites, survey, seed -----------------------------------------------
J <- 30 #sites
K <- 3 #surveys per site
specs<-10 #Number of species

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
mean.det <- runif(n = specs, min = 0.4, max = 0.8) #simulate mean detection probs
#These are mid to high detection probabilities
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
obsdata<-array(as.numeric(unlist(L)), dim=c(J, K, specs))

#Hells yeah it worked!!!

#Create an observation matrix using maximum abundance per species per site
maxobs<-apply(obsdata, c(1,3), max)
maxobs<-t(maxobs)
#This will be used as an initial value for the gibbs sampler

#Write model script -----------------------------------------------------
#install/load packages
# install.packages("R2OpenBUGS")
library(R2OpenBUGS)
#I loaded a package mid-script. Deal with it.

#Write model in BUGS language
setwd("c:/users/beasley/dropbox/MSAMsims")
cat("
    model{

    #Define hyperprior distributions

    a0.mean ~ dunif(0,1)
    mu.a0 <- log(a0.mean)-log(1-a0.mean)
    tau.a0 ~ dgamma(0.1, 0.1)

    b0.mean ~ dunif(0,1)
    mu.b0 <- log(b0.mean)-log(1-b0.mean)
    tau.b0 ~ dgamma(0.1, 0.1)

    for(i in 1:specs){
        #create priors from distributions above
            a0[i] ~ dnorm(mu.a0, tau.a0)

            b0[i] ~ dnorm(mu.b0, tau.b0)

        #Loop within a loop to estimate abund of spec i at site j
        for(j in 1:J){
            lambda[j,i] <- exp(a0[i])
            Z[j,i] ~ dpois(lambda[j,i])
            #Z is the estimated abundance matrix
            
            #Loop within loops for estimating det of spec i at site j at time k
            for(k in 1:K[j]){
                logit(p[j,k,i]) <- b0[i]
                obsdata[j,k,i] ~ dbin(p[j,k,i], Z[j,i])
            }
        }
    }
}
", file = "abundsanscovs.txt")

#Send model to Gibbs sampler

#Compile data into list
datalist<-list(specs=specs, J=J, K=Ks, obsdata=obsdata)

#Specify parameters you want to kick back to R
params<-list('Z','lambda','a0','b0', 'mu.a0', 'mu.b0', 'tau.a0', 'tau.b0')

#Specify initial values
init.values<-function(){
  list(a0 = rnorm(n = specs, mean = mean(alpha0)),
       b0 = rnorm(n = specs, mean = mean(beta0)),
       Z = t(maxobs)
  )
}

#Run the model in bugs and pray to every god created by humankind

# model3<-bugs(model.file = "abundsanscovs.txt", data = datalist, inits = init.values,
#              parameters.to.save = params, n.chains = 3, n.burnin = 2000,
#              n.iter = 15000, debug = T)

saveRDS(model3, file = "modsanscovs3.rds")

model3 <- readRDS(file = "modsanscovs3.rds")
#model 3 has best effective sample sizes and convergence

#Compare parameter estimates ---------------------------------------------------

all.lamb3<-model3$sims.list$lambda
mu.lambda3<-apply(all.lamb3, 3, mean)
lamb3low<-apply(all.lamb3, 3, quantile, probs = 0.025)
lamb3hi<-apply(all.lamb3, 3, quantile, probs = 0.975)
#1 and 3 are damn close. Will need to look at % error to check

#look at percent error
perc.err<-function(x,y){
  (x - y)/y
}

lamb.err3<-perc.err(x = mu.lambda3, y = mean.lambdas)

#After looking at means... 6-7% error.

#Compare Z estimates with true values and observation data 

#Extract estimated abundance matrix
Zs3<-model3$sims.list$Z

#Get mean estimated abundance matrices
mean.Z3<-apply(Zs3, c(2,3), mean)

obserr<-perc.err(x = maxobs, y = ns)
Z3err<-perc.err(x = mean.Z3, y = ns); Z3err[Z3err == Inf]<-NA

#Observed error is 24%, estimated is 16-17%

#Check to make sure estimated richness/abundance doesn't go below observed 
#Abundance: get credible intervals
spec.abunds <- apply(Zs3, c(1,3), sum)
spec.mean <- apply(spec.abunds, 2, mean)
spec.obs <- apply(maxobs, 1, sum)
spec.true <- rowSums(ns)

spec.comp <- data.frame(rank(spec.true), spec.mean, spec.obs, spec.true)
colnames(spec.comp) <- c("rank.true", "spec.mean", "spec.obs", "spec.true")

library(ggplot2)

ggplot(data = spec.comp, aes(x = rank.true))+
  geom_point(aes(y = spec.mean, fill = "Estimated", color = "Estimated"))+
  geom_smooth(aes(y = spec.mean, fill = "Estimated", color = "Estimated"))+
  geom_point(aes(y = spec.obs, fill = "Observed", color = "Observed"))+
  geom_smooth(aes(y = spec.obs, fill = "Observed", color = "Observed"))+
  geom_point(aes(y = spec.true, fill = "True", color = "True"))+
  geom_smooth(aes(y = spec.true, fill = "True", color = "True"))+
  scale_fill_manual(breaks = c("Estimated", "Observed", "True"),
                    values = c("blue", "black", "red"), guide = F)+
  scale_color_manual(breaks = c("Estimated", "Observed", "True"),
                     values = c("blue", "black", "red"))+
  theme_bw()+
  theme(legend.title = element_blank())


#Now check abundances per site
site.abunds <- apply(Zs3, c(1,2), sum)
site.mean <- apply(site.abunds, 2, mean)
site.true <- colSums(ns)
site.obs <- apply(maxobs, 2, sum)

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

#Check for 1:1 relationship between estimated and true abundances 
linears <- data.frame(as.vector(t(ns)),as.vector(apply(Zs3, c(2,3), mean)))
colnames(linears) <- c("Trues", "Estimated")

linears2 <- data.frame(rowSums(ns), colSums(apply(Zs3, c(2,3), mean)))
colnames(linears2) <- c("Trues", "Estimated")

ggplot(data = linears, aes(x = Trues, y = Estimated))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_bw()

ggplot(data = linears2, aes(x = Trues, y = Estimated))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_bw()

summary(lm(data = linears, Estimated~Trues))
summary(lm(data = linears2, Estimated~Trues))

#Check accuracy of diversity estimates ---------------------------------------
divers.true <- diversity(t(ns))
divers.est <- diversity(mean.Z3)
divers.obs <- diversity(t(maxobs))

divers.err <- perc.err(x = divers.est, y = divers.true)
summ.divers <- c(min(abs(divers.err)), max(abs(divers.err)), mean(divers.err))

hist(divers.err)

#Mean 2% error, most is less than 5%. 
#Community summary statistics are more accurate than site or species parameters,
#which is to be expected.

#What does the observed data look like?
divers.obs.err <- perc.err(x = divers.obs, y = divers.true)
summ.divers2 <- c(min(abs(divers.obs.err)), max(abs(divers.obs.err)), 
                  mean(divers.obs.err))

hist(divers.obs.err)

#So the maximum error is smaller, but there's more variation in error. 
#Estimated error rates appear to be skewed by an outlier

#Plot true vs estimated vs observed diversity
divers <- data.frame(Rank = rank(divers.true), divers.true, divers.est, divers.obs)

ggplot(data = divers, aes(x = Rank))+
  geom_point(aes(y = divers.true, color = "True"))+
  geom_smooth(aes(y = divers.true, color = "True"))+
  geom_point(aes(y = divers.est, color = "Estimates"))+
  geom_smooth(aes(y = divers.est, color = "Estimates"))+
  geom_point(aes(y = divers.obs, color = "Observed"))+
  geom_smooth(aes(y = divers.obs, color = "Observed"))+
  theme_bw()

#Diversity is overestimated by observed data and the model.
#Abundance of most common species isn't underestimated overall
#But maybe site-specific abundances are

true.hi.abund <- ns[6,]
est.hi.abund <- mean.Z3[,6]
obs.hi.abund <- maxobs[6,]

hi.abund <- data.frame(rank = rank(true.hi.abund), true.hi.abund, est.hi.abund,
                       obs.hi.abund)

ggplot(data = hi.abund, aes(x = rank))+
  geom_point(aes(y = true.hi.abund, color = "True"))+
  geom_point(aes(y = est.hi.abund, color = "Estimated"))+
  geom_point(aes(y = obs.hi.abund, color = "Observed"))
  
#Doesn't look like most abundant species is underestimated. 
#Don't know what's happening there
