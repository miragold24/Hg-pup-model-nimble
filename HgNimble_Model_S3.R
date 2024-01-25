#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
#library(nimble)
#library(MASS) 

#getwd() 
#source("HgSim_S3.R")
#source("HgNimble_input_S3.R") # Nimble functions to count whitecoat and moulted pups present on any given day

code <- nimbleCode({
  #priors
   #Nborn ~ dgamma(shape=500,scale=1) #Nborn (N)=500 in simulation
   Nborn ~ dgamma(shape=485, scale=1) 
   start ~ dunif(1,45) # start is 32 in simulation
   mu.bday ~ dunif(20,60)
   #sd.bday ~ dunif(7,12) 
   #shifted scale beta for sd.bday
   # beta distribution
   standard_beta ~ dbeta(alpha, beta)  
   # Shifting and scaling the distribution
   sd.bday <- shift + (scale * standard_beta)
   
   
  # # # #constraints
     one ~ dconstraint(start < mu.bday)

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
  whitePups[1:ndays] <- whiteSum(born[1:ndays], Pm_sum[1:amax])#,Pl_sum[1:amax],whiteObs[])
  moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax])#,moultObs[])
  #totPups[1:ndays] <- whitePups[1:ndays]+moultPups[1:ndays]
  
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
  # ============================================================================

# for (obs in 1:length(obsdays)){
#   wObs[obs] ~ dbinom(POWhite, round(whitePups[obsdays[obs]]))
#   mObs[obs] ~ dbinom(POMoult, round(moultPups[obsdays[obs]]))
# 
# }

  

  #################################################################
  # ONE CALL TO DBINOM FOR ALL OBSERVATIONS 
    # NOTE: CONTINUOUS QUANTITIES whitePups & moultPups rounded for binom
  #===============

  #####################################################################
  #model$whitePups[obsdays]*POWhite+model$moultPups[obsdays]*POMoult*(1-PCMoult)
  
  ## approximation of normal distribution with binomial variance
  # for (obs in 1:length(obsdays)){
  #   whiteCount[obs]<- whitePups[obsdays[obs]]*POWhite+moultPups[obsdays[obs]]*POMoult*(1-PCMoult)
  #   vwhite[obs] <- pow(whitePups[obsdays[obs]]*POWhite*(1-POWhite),2) + pow(moultPups[obsdays[obs]]*(1-PCMoult)*POMoult*(1-(1-PCMoult)*POMoult),2)
  #   sdwhite[obs] <- sqrt(vwhite[obs])
  # 
  #   moultCount[obs] <- moultPups[obsdays[obs]]*POMoult*PCMoult#+whitePups[obsdays[obs]]*POWhite*(1-PCWhite)
  #   vmoult[obs] <- pow(moultPups[obsdays[obs]]*POMoult*PCMoult*(1-POMoult*PCMoult),2)#+pow(whitePups[obsdays[obs]]*POWhite*(1-PCWhite)*1-(POWhite*(1-PCWhite)))
  #   sdmoult[obs] <- sqrt(vmoult[obs])
  # 
  #   wObs[obs] ~ dnorm(whiteCount[obs],sd=sdwhite[obs])
  #   mObs[obs] ~ dnorm(moultCount[obs],sd=sdmoult[obs])
  # }
  for (obs in 1:length(obsdays)) {
    whiteCount[obs] <- whitePups[obsdays[obs]] * POWhite + moultPups[obsdays[obs]] * POMoult * (1 - PCMoult)
    sdwhite[obs] <- sqrt(whitePups[obsdays[obs]] * POWhite * (1 - POWhite) + moultPups[obsdays[obs]] * (1 - PCMoult) * POMoult * (1 - (1 - PCMoult) * POMoult))

    moultCount[obs] <- moultPups[obsdays[obs]] * POMoult * PCMoult
    sdmoult[obs] <- sqrt(moultPups[obsdays[obs]] * POMoult * PCMoult * (1 - POMoult * PCMoult))

    wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs])
    mObs[obs] ~ dnorm(moultCount[obs], sd = sdmoult[obs])
  }


  ## debug to check if issue is with binomial process
  # for (obs in 1:length(obsdays)){
  #   wObs[obs] ~ dnorm(whitePups[obsdays[obs]],sd=0.01) 
  #   mObs[obs] ~ dnorm(moultPups[obsdays[obs]],sd=0.03)
  # 
  # }
  

  
    })

## supply data, constants and inits
data_list=list(wObs=sim_S3$whiteObs,mObs=sim_S3$moultObs,
               whiteObs=sim_S3$whiteObs,moultObs=sim_S3$moultObs)

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
                    POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult)


inits_list <- list(Nborn=runif(1,500,510),
                   # MM=rep(0,length(obsdays)),MW=rep(0,length(obsdays)),WW=rep(0,
                   #length(obsdays)),WM=rep(0,length(obsdays)),
                   #whiteCount=rep(1,length(obsdays)),moultCount=rep(1,length(obsdays)),
                   Pb=rep(0,ndays),Pm=rep(0,amax),Pl=rep(0,amax),
                    start = runif(1,1,45),
                    mu.bday=runif(1,20,60),
                    #sd.bday= runif(1,7,12),
                    one=1,
                    standard_beta=rbeta(1,2, 5) 
                   )
################
model <- nimbleModel(code=code,inits=inits_list,data=data_list,
                    constants=constants_list,
                    dimensions=list(born=ndays,
                                    #whitePups=ndays,moultPups=ndays,
                                # moultCount=length(obsdays), whiteCount=length(obsdays),
                    Pb_sum=ndays,Pm_sum=ndays,Pl_sum=ndays))
                   
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
# Cmodel<- compileNimble(model,showCompilerOutput = F)
# CHgMCMC <- compileNimble(HgMCMC,project=Cmodel)
# CHgMCMC$run(10000,reset=FALSE)
# samples_update<- as.matrix(CHgMCMC$mvSamples)
# plot(samples_update[ , 'Nborn'], type = 'l', xlab = 'iteration',  
#      ylab = expression(Nborn))
# mean(samples_update[,'Nborn'])
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
 mcmc<-mcmc.out$samples
 coda::gelman.diag(mcmc[,"Nborn"])
# #look at trace plots
  plot(mcmc[,"Nborn"])
