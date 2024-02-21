#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
#library(nimble)
#library(coda)
#getwd()

#simulation code
#source("HgSim_S4.R")
# Nimble functions to count whitecoat and moulted pups present on any given day
#source("HgNimble_input_S4.R") 

# EXAMPLE

code <- nimbleCode({
   #priors
   #Nborn ~ dgamma(shape=500,scale=1) #Nborn (N)=500 in simulation
   Nborn ~ dgamma(shape=485, scale=1) 
   mu.bday ~ dunif(33,60)
   skew.bday ~ dunif(2,6)
   #shifted scale beta for sd.bday
   standard_beta ~ dbeta(alpha, beta)  
   # Shifting and scaling the distribution
   sd.bday <- shift + (scale * standard_beta)
   
   
  #constraints
  # one ~ dconstraint(start < mu.bday)
  # deterministic and continuous births over breeding season (ndays)
  for (t in 1:ndays){ 
    born[t] <- Pb_sum[t]*Nborn
  }
  
  #=============================================================================
  ### Calculate birth probabilities 

   ###
   ### Skew normal
   ###
  
   # convert paramaters for skew normal
   # delta <- skew.bday/sqrt(1+(skew.bday^2)) #location
   # omega<-sqrt(sd.bday^2/(1-(2*delta^2)/pi)) #scale
   # xi <- mu.bday - omega*delta*sqrt(2/pi) #shape

   # #calculate PDF of skew normal
   # for (x in 1:ndays){
   #   # next, calc the PDF of the skew normal (*scaling factor omega)
   #   #at (x - xi) / omega
   #   Pb[x]<- (2/omega)*dnorm((x-xi)/omega) *
   #     #calc CDF up to the point skew.bday * ((x - xi)/omega)
   #     pnorm(skew.bday*((x-xi)/omega))
   # }

  #calculate PDF of skew normal then integrate over intervals for Pb on given day
     Pb[1:ndays]<- slicePb(ndays,nslice,skew.bday,mu.bday,sd.bday) # using EJ method
     #Pb[1:ndays] <- intPb(ndays,int,skew.bday,mu.bday,sd.bday)

 #==============================================================================
  
  Pb_sum[1:ndays] <- Pb[1:ndays]/sum(Pb[1:ndays])
  

  # cumulative sum of pups born on given day 
  for (i in 1:ndays) {
    cumsum_born[i] <- sum(born[1:i])
  }
   
   
  # account for probabilities of birth across season and 
  # transitions to moult and leave stages 

  ###  PS4 as function of age (Probability of being S4)
  PS4[1:amax] <- pnorm(1:amax, mu.s4day, sd.s4day) - pnorm(0:(amax - 1), mu.s4day, sd.s4day)
    
  PS4_sum[1:amax] <- PS4[1:amax]/sum(PS4[1:amax])
  

  # Probability of moulting as function of age (Probability of finishing moult)
  Pm[1:amax] <- pnorm(1:amax, mu.mday, sd.mday) - pnorm(0:(amax - 1), mu.mday, sd.mday)
  #Ensure probabilities sum to 1
  Pm_sum[1:amax] <- Pm[1:amax]/sum(Pm[1:amax])
  
  # Note Pl_sum is given in data from "leaveFunction.R"
  
  # calculate how many white and moulted pups present on any given day, given age
  # Nimble functions in HGNimble_input_noS4.R
  whitePups[1:ndays] <- whiteSum(born[1:ndays],PS4_sum[1:amax], Pm_sum[1:amax])
  S4Pups[1:ndays] <- S4Sum(born[1:ndays], PS4_sum[1:amax], Pm_sum[1:amax])
  moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax])
 
  


  #################################################################
  # OBSERVATION PROCESS

  #####################################################################
 # binomial variance N *p*(1-p)
  # n (number of trials)* p(prob of success) * (1 - p) (prob of failure)
  for (obs in 1:length(obsdays)) {
    whiteCount[obs] <- whitePups[obsdays[obs]] * POWhite + 
      moultPups[obsdays[obs]] * POMoult * (1 - PCMoult) + 
      S4Pups[obsdays[obs]]*POWhite*PCWhiteS4

    sdwhite[obs] <- sqrt(whitePups[obsdays[obs]] * POWhite * (1 - POWhite) + 
      moultPups[obsdays[obs]] * POMoult*(1 - PCMoult)* (1 - POMoult*(1 - PCMoult))+
      S4Pups[obsdays[obs]] * POWhite * PCWhiteS4 * (1 - POWhite*PCWhiteS4))

    moultCount[obs] <- moultPups[obsdays[obs]] * POMoult * PCMoult + 
      S4Pups[obsdays[obs]]*POWhite*(1-PCWhiteS4)
   
     sdmoult[obs] <- sqrt(moultPups[obsdays[obs]] * POMoult * PCMoult * (1 - POMoult * PCMoult)+ 
      S4Pups[obsdays[obs]]*POWhite*(1-PCWhiteS4)*(1 - POWhite*(1-PCWhiteS4)))

    
    # normal approximation using binomial variance
    # wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs])
    # mObs[obs] ~ dnorm(moultCount[obs], sd = sdmoult[obs])
    # normal distribution with small variance (should get 500 pups when 
    # exact obs process in simulation):
    wObs[obs] ~ dnorm(whiteCount[obs], sd = 0.0316)
    mObs[obs] ~ dnorm(moultCount[obs], sd = 0.0316)
    
      
    }
  

})

## supply data, constants and inits
data_list=list(wObs=sim_S4$wObs,mObs=sim_S4$mObs,Pl_sum=PmPl2$PlGivenm)

constants_list=list(ndays=ndays,amax=amax,obsdays=obsdays,nslice=1000,
                    #int=100,
                    # Parameters for the standard beta distribution
                    alpha=2 ,  # Shape parameter (alpha > 0)
                    beta=5,    # Shape parameter (beta > 0)
                    shift=7,   # Shift the distribution to start at 7
                    scale=5,   # Scale the distribution to end at 12 (12-7)
                    #skew.bday=skew.bday,
                    #mu.bday=mu.bday,
                    #sd.bday=sd.bday,
                    pi=pi,
                    #sd.bday=sd.bday, #params to be estimated
                    mu.mday=mu.mday,sd.mday=sd.mday,mu.lday=mu.lday,sd.lday=sd.lday,
                    mu.s4day=mu.s4day,sd.s4day=sd.s4day,
                    POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult,PCWhiteS4=PCWhiteS4)



inits_list <- list(Nborn=runif(1,500,510),Pb=rep(0,ndays),
                   Pb_sum=rep(0,ndays),PS4=rep(0,amax),Pm=rep(0,amax),
                   #start = runif(1,25,45),
                    mu.bday=runif(1,33,60),
                    skew.bday=runif(1,2,6),
                    # #sd.bday= runif(1,7,12),
                    # #one=1,
                    standard_beta=rbeta(1,2,5) 
                   )
################
model <- nimbleModel(code=code,inits=inits_list,data=data_list,
                    constants=constants_list,
                    dimensions=list(born=ndays,Pb_sum=ndays,PS4_sum=ndays,
                                    Pm_sum=ndays))
                   
nc <- 2   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 0
 
# Parameters monitored
params = c('Nborn','mu.bday', 'sd.bday','skew.bday')
 
 HgMCMC <- buildMCMC(model,monitors=params)
 runMCMC(HgMCMC,niter=1)#set =1 bc with one calculation inits are used
#  #& can diagnose there
#  # # #compile the model object first for the MCMC to compile 
# Cmodel<- compileNimble(model,showCompilerOutput = F)
# CHgMCMC <- compileNimble(HgMCMC,project=Cmodel)
# CHgMCMC$run(10000,reset=FALSE)
# samples_update<- as.matrix(CHgMCMC$mvSamples)
# plot(samples_update[ , 'sd.bday'], type = 'l', xlab = 'iteration',
#      ylab = expression(Nborn))
# mean(samples_update[,'Nborn'])
# Nborn_samples <- samples[, "Nborn"]
# plot(density(Nborn_samples))
# start_samples <- samples[, "start"]
# plot(density(start_samples))
# mu.bday_samples <- samples[, "mu.bday"]
# plot(density(mu.bday_samples))
# # Plot the density
# plot(density, xlab = "Nborn", ylab = "Density")
# #samples<-runMCMC(CHgMCMC,niter=10000)
#  print(system.time(CHgMCMC$run(10000)))
#  samples<-as.matrix(CHgMCMC$mvSamples)
#  mean(samples[,'Nborn'])

#############
 #Run MCMC chains
mcmc.out <-nimbleMCMC(model=model,nchains =nc,niter = ni, nburnin = nb,
                  monitors = params, samplesAsCodaMCMC = TRUE,summary=TRUE)
# #                    
# #                    # check convergence and look at plots 
 mcmcS4<-mcmc.out$samples
#  coda::gelman.diag(mcmcS4[,"Nborn"])
# # #look at trace plots
#   plot(mcmcS4[,"Nborn"])
#   coda::gelman.diag(mcmc[,"mu.bday"])
#   # #look at trace plots
#   plot(mcmcS4[,"mu.bday"])
  library(MCMCvis)
 coda::gelman.diag(mcmcS4[,"mu.bday"])
  #trace and posterior density of params:
  MCMCtrace(object = mcmcS4,
            pdf = FALSE, 
            ind = TRUE, # separate density lines per chain
            Rhat=F,
            n.eff=F,
            params = "all")
  