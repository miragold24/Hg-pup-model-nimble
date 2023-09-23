#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
library(nimble)
library(coda)
library(MCMCvis)

#getwd() 

#source("HgSim_noS4_MG.R")  #simulate pups born over season () 
#source("HgNimble_input_noS4.R") # Nimble functions to count whitecoat and moulted pups present on any given day
#source("HgSim_S3.R")
code <- nimbleCode({
  #priors
  Nborn ~ dgamma(shape=10000,rate=20)
  #Nborn ~ dgamma(shape=50000, rate=100) #still mean of 500, but more diffuse
  start ~ dunif(1,55) # start is 32 in simulation
  mu.bday ~ dnorm(45,sd=sd.bday)
  sd.bday ~ dunif(9,10)
  #constraints
  one ~ dconstraint(start < mu.bday)
  # for (obs in 1:length(obsdays)){
  # two[obs] ~ dconstraint(wObs[obs] < whiteCount[obs])
  # three[obs] ~ dconstraint(mObs[obs] < moultCount[obs])
  # }
  
  # deterministic and continuous births over breeding season (ndays)
  for (t in 1:ndays){ 
    born[t] <- Pb_sum[t]*Nborn
  }
  
  
  ###
  ### derived quantities for log scale calculations 
  ### 
  lCV <- sd.bday/(mu.bday-start)
  lsd <- sqrt(log(lCV^2+1))
  lmu <- log(mu.bday-start)-(0.5*(lsd^2))
  
  # Calculate probability of being born on any given day
  for(t in 2:ndays){
    Pb[t] <- plnorm((t-start), lmu, lsd) - plnorm((t-start-1), lmu, lsd)
    
  }
  Pb_sum[1:ndays] <- Pb[1:ndays]/sum(Pb[1:ndays])
  
  
  # cumulative sum of pups born on given day 
  for (i in 1:ndays) {
    cumsum_born[i] <- sum(born[1:i])
  }
  # account for probabilities of birth across season and 
  # transitions to moult and leave stages 
  
  ###
  ###  Probability of moulting across season
  ###
  
  for (t in 1:amax){
    Pm[t] <- pnorm(t, mu.mday, sd.mday) - pnorm(t-1, mu.mday, sd.mday)
  }
  
  Pm_sum[1:amax] <- Pm[1:amax]/sum(Pm[1:amax])
  
  ###
  ###  Probability of leaving colony across season
  ###
  
  for (t in 1:amax){
    Pl[t] <- pnorm(t, mu.lday, sd.lday) - pnorm(t-1, mu.lday, sd.lday)
  }
  
  Pl_sum[1:amax] <- Pl[1:amax]/sum(Pl[1:amax])
  
  # calculate how many white and moulted pups present on any given day, given age
  # Nimble functions in HGNimble_input_noS4.R
  whitePups[1:ndays] <- whiteSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax],whiteObs[])
  moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax],moultObs[])
  
  
  #N[1:ndays] <-whitePups[1:ndays]+moultPups[1:ndays]
  #================================================================
  ### observation model :
  # account for imperfect observation/classification error 
  # #==============================================================
  ### Observation error: 
  ### POWhite = probability of observing a whitecoat
  ### POMoult = probability of observing a moulted pup 
  ### Classification error:
  ### PCMoult = 0.82 Probability of correctly classifying a moulted as moulted
  ### 1-PCMoult is the probability of incorrectly classifying moulted as whitecoat 
  ### (add to whitecoat obs count)
  ### PCWhite = 1.0
  # you cannot have more pups observed than were there to be observed on survey day
  #but classification error allows for more observed whitepups than are actually 
  #there bc of misclassification 
  # =================================================================
  ## observation error only:
  for (obs in 1:length(obsdays)) {
    wObs[obs] ~ dbinom(POWhite, round(whitePups[obsdays[obs]]))
    mObs[obs] ~ dbinom(POMoult, round(moultPups[obsdays[obs]]))
    }

  ## debug to check if issue is with binomial process
  # for (obs in 1:length(obsdays)){
  #   wObs[obs] ~ dnorm(whitePups[obsdays[obs]],sd=0.01)
  #   mObs[obs] ~ dnorm(moultPups[obsdays[obs]],sd=0.03)
  #   
  # }
  
  
  #################################################################
  ###OBS + CLASSIFICATION ERROR 
  ### account for obs error
  #  whiteCount[1:ndays] <- round((whitePups[1:ndays]*POWhite))
  #  moultCount[1:ndays] <- round((moultPups[1:ndays]*POMoult))
  #  #NCount[1:ndays] <- whiteCount[1:ndays]+moultCount[1:ndays]
  # 
  # #add classification error
  # for (obs in 1:length(obsdays)){
  #  #WhiteCount * prob of incorrectly classifying moulted as white 
  #  wObs[obs] ~ dbinom((1-PCMoult), whiteCount[obsdays[obs]])
  #  mObs[obs] ~ dbinom( PCMoult, moultCount[obsdays[obs]])
  # 
  #  }
  
   })


## supply data, constants and inits
data_list=list(wObs=sim_S3$whiteObs,mObs=sim_S3$moultObs,
               whiteObs=sim_S3$whiteObs,moultObs=sim_S3$moultObs)

constants_list=list(ndays=ndays,amax=amax,obsdays=obsdays,
                    #start=start,
                    #mu.bday=mu.bday,
                    #sd.bday=sd.bday, #params to be estimated
                    mu.mday=mu.mday,sd.mday=sd.mday,mu.lday=mu.lday,sd.lday=sd.lday,
                    POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult)

inits_list <- list(Nborn=runif(1,570,600),
                   Pb=rep(0,ndays),Pm=rep(0,amax),Pl=rep(0,amax),
                   start = runif(1,31,33),
                   mu.bday=runif(1,44,46),
                   sd.bday= 10,
                   one=1)
################
model <- nimbleModel(code=code,inits=inits_list,data=data_list,
                    constants=constants_list,
                    dimensions=list(born=ndays,
                                    #whitePups=ndays,moultPups=ndays,
                                   # moultCount=length(obsdays),
                    Pb_sum=ndays,Pm_sum=ndays,Pl_sum=ndays))
                   
nc <- 2   ;   ni <- 11000   ;   nb <- 1000   ;   nt <- 0
 # Parameters monitored
                   
 params = c('Nborn','whitePups','moultPups',
            'start','mu.bday', 'sd.bday')
 HgMCMC <- buildMCMC(model,monitors=params)
 runMCMC(HgMCMC,niter=1)#set =1 bc with one calculation inits are used
 #& can diagnose there
 # # #compile the model object first for the MCMC to compile 
#  Cmodel<- compileNimble(model,showCompilerOutput = F)
#  CHgMCMC <- compileNimble(HgMCMC,project=Cmodel)
# #samples<-runMCMC(CHgMCMC,niter=10000)
#  print(system.time(CHgMCMC$run(10000)))
#  samples<-as.matrix(CHgMCMC$mvSamples)
#  mean(samples[,'Nborn']) 
#############
mcmc.out <-nimbleMCMC(model=model,nchains =nc,niter = ni, nburnin = nb,
                  monitors = params, samplesAsCodaMCMC = TRUE,summary=TRUE)
#                    
#                    # check convergence and look at plots 
 mcmc<-mcmc.out$samples
 coda::gelman.diag(mcmc[,"Nborn"])
# #look at trace plots
  plot(mcmc[,"Nborn"])
