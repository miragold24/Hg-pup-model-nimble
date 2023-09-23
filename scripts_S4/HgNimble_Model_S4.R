#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
library(nimble)
#library(coda)
#getwd() 
#source("input_HgNim.R")
source("HgSim_V3_MG.R")

code <- nimbleCode({
  
  #likelihood
  for (t in 1:ndays){
    born[t] <- Pb_sum[t]*Nborn
  }
  #specify priors
  Nborn ~ dpois(500)
  start ~ dnorm(32,5)
  mu.bday ~ dnorm(45,sd.bday)
  sd.bday ~ dnorm(5,5)
  

  lCV <- sd.bday/(mu.bday-start)
  lsd <- sqrt(log(lCV^2+1))
  lmu <- log(mu.bday-start)-(0.5*(lsd^2))
  #================================================================
  ### birth process model:
  # account for probabilities of birth across season and 
  # transitions to S4, moult, and leave stages 
  #================================================================

  for(t in 2:ndays){
    Pb[t] <- plnorm((t-start), lmu, lsd) - plnorm((t-start-1), lmu, lsd)

  }
  Pb_sum[1:ndays] <- Pb[1:ndays]/sum(Pb[1:ndays])
  
  #deterministic and continuous birth and transition processes

  
  for (i in 1:ndays) {
    cumsum_born[i] <- sum(born[1:i])
  }
      ###
      ###  PS4 probabilities 
      ###

  for (t in 2:amax){ #for 2:123
    PS4[t] <- pnorm(t, mu.s4day, sd.s4day) - pnorm(t-1, mu.s4day, sd.s4day)

  }

  PS4_sum[1:amax] <- PS4[1:amax]/sum(PS4[1:amax])

      ###
      ###  Pm probabilities
      ###

      for (t in 1:amax){
        Pm[t] <- pnorm(t, mu.mday, sd.mday) - pnorm(t-1, mu.mday, sd.mday)
      }

      Pm_sum[1:amax] <- Pm[1:amax]/sum(Pm[1:amax])
      
      ###
      ###  Pl probabilities 
      ###
      
      for (t in 1:amax){
        Pl[t] <- pnorm(t, mu.lday, sd.lday) - pnorm(t-1, mu.lday, sd.lday)
      }
      
      Pl_sum[1:amax] <- Pl[1:amax]/sum(Pl[1:amax])
    
      
    whitePups[1:ndays] <- whiteSum(born[1:ndays], PS4_sum[1:amax], Pm_sum[1:amax], Pl_sum[1:amax])
    S4Pups[1:ndays] <- S4Sum(born[1:ndays], PS4_sum[1:amax], Pm_sum[1:amax], Pl_sum[1:amax])
    moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax])
  #================================================================
  ### observation model :
  # all white and moulted counted on each survey day
  # account for imperfect detection and imperfect classification error
  #================================================================

     #
    for (obs in 1:length(obsdays)){
    wObs[obs] ~ dbinom(POWhite, whiteobs[obs])
    mObs[obs] ~ dbinom(POMoult, moultobs[obs])

}
#for exact observations
# for (obs in 1:length(obsdays)){
# 
#   wObs[obs]<- (whitePups[obsdays[obs]]*POWhite)+(moultPups[obsdays[obs]]* POMoult*(1-PCMoult))+
#     (S4Pups[obsdays[obs]]*POWhite*PCWhiteS4)
# 
#   mObs[obs]<- (moultPups[obsdays[obs]]*POMoult*PCMoult)+(whitePups[obsdays[obs]]*
#     POWhite*(1-PCWhiteS4))
# }


})
    
## supply data, constants and inits
data_list=list(POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult,PCWhiteS4=PCWhiteS4,
               whiteobs=round(simS4$whiteObs),moultobs=round(simS4$moultObs))

constants_list=list(obsdays = c(35, 40, 50, 65, 70, 80, 95),ndays=ndays,amax=amax,
                    mu.s4day=mu.s4day,sd.s4day=sd.s4day,mu.mday=mu.mday,sd.mday=sd.mday,
                    mu.lday=mu.lday,sd.lday=sd.lday)
#NOTE: bc truncated lognorm distribution of births, start needs to be less than bday for log calcs above
inits_list <- list(Nborn=100,start = runif(1,1,40),mu.bday = runif(1,41,100),
                   sd.bday = runif(1, 5, 10),
                   Pb=rep(0,ndays),PS4=rep(0,amax),Pm=rep(0,amax),Pl=rep(0,amax),
                   whitePups=rep(1,ndays),S4Pups=rep(1,ndays),moultPups=rep(1,ndays))

# MCMC settings
nc <- 2   ;   ni <- 12000   ;   nb <- 1000   ;   nt <- 1
# Parameters monitored
#params = c('cumsum_born','Pb','PS4','Pm','Pl','whitePups','S4Pups','moultPups','wObs','mObs')


##debug 
params = c('Nborn','born','whitePups','S4Pups','moultPups',
           'wObs','mObs','start','mu.bday','sd.bday')

#build model object from NimbleCode
model <- nimbleModel(code=code,inits=inits_list,data=data_list,constants=constants_list,
                     dimensions=list(born=ndays,PS4_sum=ndays,Pm_sum=ndays,Pl_sum=ndays))

# #build MCMC object
HgMCMC <- buildMCMC(model,monitors=params)#,monitors = parentNodes)
# samples<-
runMCMC(HgMCMC,niter=1)# can run this now uncompiled, but will be slow. 
                                #set =1 bc with one calculation inits are used
                                #and something wrong with initial values 

## MCMC configuration: what nodes samplers operate on, any parameters needed, etc
# HgMCMCconfiguration <- configureMCMC(model)
# HgMCMCconfiguration$printSamplers()
# # #compile the model object first for the MCMC to compile 
# Cmodel<- compileNimble(model,showCompilerOutput = TRUE)
# CHgMCMC <- compileNimble(HgMCMC,project=Cmodel)
# # now runMCMC on compiled model object/plot results: 
# 
# runMCMC_samples <- runMCMC(CHgMCMC,nchains=2,nburnin=500,niter=1000,samplesAsCodaMCMC = TRUE,summary=TRUE)
# mcmc.samples<-as.mcmc(runMCMC_samples)
# plot(mcmc.samples$chain1[ , 'Nborn'], type = 'l', xlab = 'iteration',  
# ylab = expression(start))
# 
mcmc.out <-nimbleMCMC(model=model,nchains =nc,niter = ni, nburnin = nb,
                      monitors = params, samplesAsCodaMCMC = TRUE,summary=TRUE)
# ,WAIC=TRUE

mcmc.out$summary
# ### check convergence and look at plots 
mcmc<-mcmc.out$samples
coda::gelman.diag(mcmc[,"mu.bday"])
# ## look at trace plots
plot(mcmc[,"sd.bday"])
