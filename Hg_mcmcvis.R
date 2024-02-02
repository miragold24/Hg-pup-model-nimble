### Outputs/plots for MCMC runs / visualizations
###############################################
library(coda)
library(MCMCvis)
library(bayesplot)
###
###visualize priors
####

Nborn.prior<-rgamma(10000,shape=485, scale=1)
quantile(Nborn.prior)
plot(density(Nborn.prior))
mean(Nborn.prior)
Nborn.prior<-rgamma(10000,shape=4950, rate=10)
quantile(Nborn.prior)
mean(Nborn.prior)
lines(density(Nborn.prior))

start.prior<-runif(10000,25,45)
quantile(start.prior)
plot(density(start.prior))
mubday.prior<-runif(10000,35,50)
quantile(mubday.prior)
plot(density(mubday.prior))


sd.prior<-runif(10000,7,12)
quantile(sd.prior)
plot(density(sd.prior))
#########################################################
##mcmc output 
mcmc<-mcmc.out$samples
coda::gelman.diag(mcmc[,"sd.bday"])
# # # # # # # ## look at trace plots
plot(mcmc[,"Nborn"])
#trace and posterior density of params:
MCMCtrace(object = mcmc,
          pdf = FALSE, 
          ind = TRUE, # separate density lines per chain
          Rhat=F,
          n.eff=F,
          params = "all")



#95% credible intervals:
quantile(mcmc$chain1[,'Nborn'],probs=c(2.5,97.5)/100)


# check means:
Nborn.samples<-  c(mcmc$chain1[,'Nborn'], 
                   mcmc$chain2[,'Nborn'])
mean(Nborn.samples)

 start.samples<-  c(mcmc$chain1[,'start'], 
                     mcmc$chain2[,'start'])
 mean(start.samples)
mu.bday.samples<-  c(mcmc$chain1[,'mu.bday'], 
                    mcmc$chain2[,'mu.bday'])
 mean(mu.bday.samples)
 
sd.bday.samples<-  c(mcmc$chain1[,'sd.bday'], 
                      mcmc$chain2[,'sd.bday'])
mean(sd.bday.samples)

###########################################

#Pair plots to check correlation:

mcmc_pairs(mcmc)

