### Outputs/plots for MCMC runs

#########################################################
##mcmc output 
mcmc<-mcmc.out$samples
coda::gelman.diag(mcmc[,"sd.bday"])
# # # # # # # ## look at trace plots
plot(mcmc[,"sd.bday"])
#trace and posterior density of params:
MCMCtrace(object = mcmc,
          pdf = FALSE, 
          ind = TRUE, # separate density lines per chain
          Rhat=F,
          n.eff=F,
          params = "mu.bday")

#95% credible intervals:
quantile(mcmc$chain1[,'Nborn'],probs=c(2.5,97.5)/100)
sd(mcmc$chain1[,"Nborn"])



# isolate with param samples:
 Nborn.samples<-  c(mcmc$chain1[,'Nborn'], 
                     mcmc$chain2[,'Nborn'])
 mean(Nborn.samples)


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
 # plot(samples_update[ , 'Nborn'], type = 'l', xlab = 'iteration',  
 #       ylab = expression(Nborn))
 mean(samples_update[,'Nborn'])
