#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
#library(nimble)
#library(MASS) 

#getwd() 
#source("HgSim_S4.R")
#source("HgNimble_input_S4.R") # Nimble functions to count whitecoat and moulted pups present on any given day

# EXAMPLE

code <- nimbleCode({
   #priors
   #Nborn ~ dgamma(shape=500,scale=1) #Nborn (N)=500 in simulation
   Nborn ~ dgamma(shape=485, scale=1) 
   start ~ dunif(25,45)
   mu.bday ~ dunif(20,60)
   #sd.bday ~ dunif(7,12) 
   #shifted scale beta for sd.bday
   # beta distribution
   standard_beta ~ dbeta(alpha, beta)  
   # Shifting and scaling the distribution
   sd.bday <- shift + (scale * standard_beta)
   
   
  #constraints
     one ~ dconstraint(start < mu.bday)
  # alt constraint - truncated normal for start
     #start ~ T(dnorm(0, 1), 0, mu.bday)
     #T(dnorm(0, 1), -15, mu.bday) #start is within 15 days of mu.bday
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
    Pb[t] <- plnorm((t-start+1), lmu, lsd) - plnorm((t-start), lmu, lsd)
    
  }
  Pb_sum[1:ndays] <- Pb[1:ndays]/sum(Pb[1:ndays])
  

  # cumulative sum of pups born on given day 
  for (i in 1:ndays) {
    cumsum_born[i] <- sum(born[1:i])
  }
  # account for probabilities of birth across season and 
  # transitions to moult and leave stages 
  
  
  ###
  ###  PS4 probabilities 
  ###
  
  for (t in 1:amax){ #for 2:123
    PS4[t] <- pnorm(t, mu.s4day, sd.s4day) - pnorm(t-1, mu.s4day, sd.s4day)
    
  }
  
  PS4_sum[1:amax] <- PS4[1:amax]/sum(PS4[1:amax])
  
  ###
  ###  Probability of moulting across season
  ##
  
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
  whitePups[1:ndays] <- whiteSum(born[1:ndays],PS4_sum[1:amax], Pm_sum[1:amax])
  S4Pups[1:ndays] <- S4Sum(born[1:ndays], PS4_sum[1:amax], Pm_sum[1:amax])
  moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax],moultObs[])
  #totPups[1:ndays] <- whitePups[1:ndays]+moultPups[1:ndays]
  


  #################################################################
  # ONE CALL TO DBINOM FOR ALL OBSERVATIONS 
    # NOTE: CONTINUOUS QUANTITIES whitePups & moultPups rounded for binom
  #===============

  #####################################################################
 # binomial variance formula
  #N *p*(1-p)
  # n (number of trials)* p(prob of success) * (1 - p) (prob of failure)
  for (obs in 1:length(obsdays)) {
    whiteCount[obs] <- whitePups[obsdays[obs]] * POWhite*PCWhite + moultPups[obsdays[obs]] * POMoult * (1 - PCMoult) + S4Pups[obsdays[obs]]*POWhite*PCWhiteS4

    sdwhite[obs] <- sqrt(whitePups[obsdays[obs]] * POWhite*PCWhite * (1 - POWhite*PCWhite) + moultPups[obsdays[obs]] * POMoult*(1 - PCMoult)* (1 - POMoult*(1 - PCMoult))+
      S4Pups[obsdays[obs]] * POWhite * PCWhiteS4 * (1 - POWhite*PCWhiteS4))

    moultCount[obs] <- moultPups[obsdays[obs]] * POMoult * PCMoult + S4Pups[obsdays[obs]]*POWhite*(1-PCWhiteS4)
    sdmoult[obs] <- sqrt(moultPups[obsdays[obs]] * POMoult * PCMoult * (1 - POMoult * PCMoult)+ 
      S4Pups[obsdays[obs]]*POWhite*(1-PCWhiteS4)*(1 - POWhite*(1-PCWhiteS4)))

    
    wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs])
    mObs[obs] ~ dnorm(moultCount[obs], sd = sdmoult[obs])
  }
  
  


})

## supply data, constants and inits
data_list=list(wObs=sim_S4$whiteObs,mObs=sim_S4$moultObs,moultObs=sim_S4$moultObs)

constants_list=list(ndays=ndays,amax=amax,obsdays=obsdays,
                    # Parameters for the standard beta distribution
                    alpha=2 ,  # Shape parameter (alpha > 0)
                    beta=5,    # Shape parameter (beta > 0)
                    shift=7,   # Shift the distribution to start at 7
                    scale=5,   # Scale the distribution to end at 12 (12-7)
                     #start=start,
                     #mu.bday=mu.bday,
                     #sd.bday=sd.bday, #params to be estimated
                    mu.mday=mu.mday,sd.mday=sd.mday,mu.lday=mu.lday,sd.lday=sd.lday,
                    mu.s4day=mu.s4day,sd.s4day=sd.s4day,
                    POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult,PCWhiteS4=PCWhiteS4)



inits_list <- list(Nborn=runif(1,500,510),
                   Pb=rep(0,ndays),PS4=rep(0,amax),Pm=rep(0,amax),Pl=rep(0,amax),
                   start = runif(1,25,45),
                   mu.bday=runif(1,20,60),
                    sd.bday= runif(1,7,12),
                    one=1,
                    standard_beta=rbeta(1,2,5) 
                   )
################
model <- nimbleModel(code=code,inits=inits_list,data=data_list,
                    constants=constants_list,
                    dimensions=list(born=ndays,
                                    #whitePups=ndays,moultPups=ndays,
                                #moultCount=length(obsdays), whiteCount=length(obsdays),
                    Pb_sum=ndays,PS4_sum=ndays,Pm_sum=ndays,Pl_sum=ndays))
                   
nc <- 2   ;   ni <- 11000   ;   nb <- 1000   ;   nt <- 0
 # Parameters monitored
                   
 #' params = c('cumsum_born','Nborn',
 #'            'whitePups','moultPups',#'totPups',
 #'            #'whiteCount','moultCount',# 'vwhite','vmoult',
 #'            'start','mu.bday', 'sd.bday'
 #'            #,'mm','mw',)
 #' )

params = c('Nborn','start','mu.bday', 'sd.bday')
 
 HgMCMC <- buildMCMC(model,monitors=params)
 runMCMC(HgMCMC,niter=1)#set =1 bc with one calculation inits are used
 #& can diagnose there
 # # #compile the model object first for the MCMC to compile 
Cmodel<- compileNimble(model,showCompilerOutput = F)
CHgMCMC <- compileNimble(HgMCMC,project=Cmodel)
CHgMCMC$run(10000,reset=FALSE)
samples_update<- as.matrix(CHgMCMC$mvSamples)
plot(samples_update[ , 'sd.bday'], type = 'l', xlab = 'iteration',
     ylab = expression(Nborn))
mean(samples_update[,'Nborn'])
Nborn_samples <- samples[, "Nborn"]
plot(density(Nborn_samples))
start_samples <- samples[, "start"]
plot(density(start_samples))
mu.bday_samples <- samples[, "mu.bday"]
plot(density(mu.bday_samples))
# Plot the density
plot(density, xlab = "Nborn", ylab = "Density")
#samples<-runMCMC(CHgMCMC,niter=10000)
 print(system.time(CHgMCMC$run(10000)))
 samples<-as.matrix(CHgMCMC$mvSamples)
 mean(samples[,'Nborn'])

#############
 #Run MCMC chains
mcmc.out <-nimbleMCMC(model=model,nchains =nc,niter = ni, nburnin = nb,
                  monitors = params, samplesAsCodaMCMC = TRUE,summary=TRUE)
# #                    
# #                    # check convergence and look at plots 
 mcmc<-mcmc.out$samples
 coda::gelman.diag(mcmc[,"Nborn"])
# #look at trace plots
  plot(mcmc[,"Nborn"])
  coda::gelman.diag(mcmc[,"mu.bday"])
  # #look at trace plots
  plot(mcmc[,"mu.bday"])
  library(MCMCvis)
  #trace and posterior density of params:
  MCMCtrace(object = mcmc,
            pdf = FALSE, 
            ind = TRUE, # separate density lines per chain
            Rhat=F,
            n.eff=F,
            params = "all")
  