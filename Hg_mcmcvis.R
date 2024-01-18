### Outputs/plots for MCMC runs / visualizations
###############################################
library(coda)
library(MCMCvis)
###
###visualize priors
####

#Nborn

#test git update
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
sd(mcmc$chain1[,"sd.bday"])

head(mcmc)

# isolate with param samples:
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

 #############
 ## MCMC configuration: what nodes samplers operate on, any parameters needed, etc
 #HgMCMCconf <- configureMCMC(model)
 # HgMCMCconfiguration$printSamplers()
 #build MCMC object
 HgMCMC <- buildMCMC(model,monitors=params,enableWAIC = TRUE)
 #runMCMC(HgMCMC,niter=1)#set =1 bc with one calculation inits are used
 #& can diagnose there
 # # #compile the model object first for the MCMC to compile 
 Cmodel<- compileNimble(model,showCompilerOutput = F)
 CHgMCMC <- compileNimble(HgMCMC,project=Cmodel)
 samples<-runMCMC(CHgMCMC,niter=12000)
 #print(system.time(CHgMCMC$run(10000)))
 samples<-as.matrix(CHgMCMC$mvSamples)
 mean(samples[,'Nborn'])
 ## plot MCMC chains
 #par(mfrow = c(1, 3))
 plot(samples[ , 'Nborn'], type = 'l', xlab = 'iteration',
      ylab = expression(Nborn))
 
 # plot(samples[ , 'start'], type = 'l', xlab = 'iteration',
 #      ylab = expression(start))
 # plot(samples[ , 'mu.bday'], type = 'l', xlab = 'iteration',
 #      ylab = expression(mu.bday))
 # plot(samples[ , 'sd.bday'], type = 'l', xlab = 'iteration',
 #      ylab = expression(sd.bday))
 ## continue the same run of the MCMC for 10000 more iterations: 
 CHgMCMC$run(10000,reset=FALSE)
samples_update<- as.matrix(CHgMCMC$mvSamples)
plot(samples_update[ , 'Nborn'], type = 'l', xlab = 'iteration',  
        ylab = expression(Nborn))
 mean(samples_update[,'Nborn'])
 mean(samples_update[,'start'])
 mean(samples_update[,'mu.bday'])
 mean(samples_update[,'sd.bday'])

 ### plot priors with posteriors 
