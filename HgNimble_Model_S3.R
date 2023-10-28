#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
library(nimble)
library(coda)
library(MCMCvis)

#getwd() 

#source("HgNimble_input_S3.R") # Nimble functions to count whitecoat and moulted pups present on any given day
#source("HgSim_S3.R")
code <- nimbleCode({
  #priors
  #Nborn ~ dgamma(shape=550,scale=1) #Nborn (N)=500 in simulation
  Nborn ~ dgamma(shape=50000, rate=100) #still mean of 500
  start ~ dunif(10,55) # start is 32 in simulation
  mu.bday ~ dunif(35,50)
  sd.bday ~ dunif(9,11) #shifted scale beta(?)
  # #constraints
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
  whitePups[1:ndays] <- whiteSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax],whiteObs[])
  moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax],moultObs[])
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

    
  ################################################################
  #CLASSIFICATION BEFORE OBSERVATION:
# 
#   for (obs in 1:length(obsdays)){
#  # do the classifications add up to moultPups+whitePups?
#   mm[obs] <- round(moultPups[obsdays[obs]] * PCMoult)
#   mw[obs] <- round(moultPups[obsdays[obs]]) - mm[obs]
#   ww[obs] <- round(whitePups[obsdays[obs]] * PCWhite)
#   wm[obs] <- round(whitePups[obsdays[obs]]) - ww[obs]
#   #whiteCount[obs] <- ww[obs]+mw[obs]
#   #moultCount[obs] <- mm[obs] +wm[obs]
#   }
# 
#   #function that constrains whiteCount to be > by 1 than wObs & mObs
#   whiteCount[1:length(obsdays)] <- wCount(whiteObs[],ww[],mw[])
#   moultCount[1:length(obsdays)] <- mCount(moultObs[],mm[],wm[])
# 
#   for (obs in 1:length(obsdays)){
#   wObs[obs] ~ dbinom(POWhite,whiteCount[obs])
#   mObs[obs] ~ dbinom(POMoult,moultCount[obs])
#   }
  

  #################################################################
  # ONE CALL TO DBINOM FOR ALL OBSERVATIONS 
    # NOTE: CONTINUOUS QUANTITIES whitePups & moultPups rounded for binom
  #===============
  #debug code
  # for (obs in 1:length(obsdays)){
  #   wmObs<- rep(0,length(obsdays))
  #   p<- rep(0,length(obsdays))
  #   totPups <- model$whitePups[obsdays[obs]] + model$moultPups[obsdays[obs]]
  # wmObs[obs] <-  (model$whitePups[obsdays[obs]]*Pwhite) + (model$moultPups[obsdays[obs]]*Pmoult)
  # p[obs]<-wmObs[obsdays[obs]]/totPups[obsdays[obs]]
  # 
  # totObs[obs] ~ dbinom(p[obs],round(totPups[obsdays[obs]]))
  # #
  # }
  #================
  # 
  # for (obs in 1:length(obsdays)){
  #   wmObs<- rep(0,length(obsdays))
  #      p<- rep(0,length(obsdays))
  # whiteCount[obs]<- whitePups[obsdays[obs]]*POWhite+moultPups[obsdays[obs]]*POMoult*(1-PCMoult)
  # moultCount[obs] <- moultPups[obsdays[obs]]*POMoult*PCMoult#+whitePups[obsdays[obs]]*POWhite*(1-PCWhite)
  # wmObs[obs] <- whiteCount[obs]+moultCount[obs]
  # p[obs]<-wmObs[obs]/totPups[obsdays[obs]]
  # 
  # totObs[obs] ~ dbinom(p[obs],round(totPups[obsdays[obs]]))
  # 
  # }
  #####################################################################
  #model$whitePups[obsdays]*POWhite+model$moultPups[obsdays]*POMoult*(1-PCMoult)
  
  ## approximation of normal distribution with binomial variance
  for (obs in 1:length(obsdays)){
    whiteCount[obs]<- whitePups[obsdays[obs]]*POWhite+moultPups[obsdays[obs]]*POMoult*(1-PCMoult)
    vwhite[obs] <- pow(whitePups[obsdays[obs]]*POWhite*(1-POWhite),2) + pow(moultPups[obsdays[obs]]*(1-PCMoult)*POMoult*(1-(1-PCMoult)*POMoult),2)
    sdwhite[obs] <- sqrt(vwhite[obs])

    moultCount[obs] <- moultPups[obsdays[obs]]*POMoult*PCMoult#+whitePups[obsdays[obs]]*POWhite*(1-PCWhite)
    vmoult[obs] <- pow(moultPups[obsdays[obs]]*POMoult*PCMoult*(1-POMoult*PCMoult),2)#+pow(whitePups[obsdays[obs]]*POWhite*(1-PCWhite)*1-(POWhite*(1-PCWhite)))
    sdmoult[obs] <- sqrt(vmoult[obs])

    wObs[obs] ~ dnorm(whiteCount[obs],sd=sdwhite[obs])
    mObs[obs] ~ dnorm(moultCount[obs],sd=sdmoult[obs])
  }


  ### NOTE: wObs + mObs should equal whiteCount+moultCount bc same total observed 
  # ##################################################################
  # PWhite[1:ndays]<- Pb_sum[]*(1-Pm_sum[]*(1-Pl_sum[])*POWhite*PCWhite)+(Pm_sum[]*
  #                    (1-Pl_sum[])*POMoult*(1-PCMoult))
  # PMoult[1:ndays]<- Pb_sum[]*(Pm_sum[]*(1-Pl_sum[])*POMoult*PCMoult)+(1-Pm_sum[])*
  #                   (1-Pl_sum[])*POWhite*(1-PCWhite)
  # 
  # 
  #  for (obs in 1:length(obsdays)){
  #   wObs[obs] ~ dbinom(PWhite[obsdays[obs]],round(totPups[obsdays[obs]]))
  #   mObs[obs] ~ dbinom(PMoult[obsdays[obs]],round(totPups[obsdays[obs]]))
  #  }

  ## debug to check if issue is with binomial process
  # for (obs in 1:length(obsdays)){
  #   wObs[obs] ~ dnorm(whitePups[obsdays[obs]],sd=0.01) dflkjs
  #   mObs[obs] ~ dnorm(moultPups[obsdays[obs]],sd=0.03)
  # 
  # }
  
  # for (obs in 1:length(obsdays)){
  #   wObs[obs] ~ dbinom(POWhite, round(whitePups[obsdays[obs]]))
  #   mObs[obs] ~ dbinom(POMoult, round(moultPups[obsdays[obs]]))
  #   
  # }
  
    })

## supply data, constants and inits
data_list=list(wObs=sim_S3$whiteObs,mObs=sim_S3$moultObs,
               whiteObs=sim_S3$whiteObs,moultObs=sim_S3$moultObs)

constants_list=list(ndays=ndays,amax=amax,obsdays=obsdays,
                    # start=start,
                    # mu.bday=mu.bday,
                    # sd.bday=sd.bday, #params to be estimated
                    mu.mday=mu.mday,sd.mday=sd.mday,mu.lday=mu.lday,sd.lday=sd.lday,
                    POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult,PCWhite=PCWhite)
                    #Pwhite=POWhite*PCWhite*(POMoult*(1-PCMoult)),
                    #Pmoult=POMoult*PCMoult)

inits_list <- list(Nborn=runif(1,500,510),
                   # MM=rep(0,length(obsdays)),MW=rep(0,length(obsdays)),WW=rep(0,
                   #  length(obsdays)),WM=rep(0,length(obsdays)),
                   #whiteCount=rep(1,length(obsdays)),moultCount=rep(1,length(obsdays)),
                   Pb=rep(0,ndays),Pm=rep(0,amax),Pl=rep(0,amax),
                   # start = runif(1,31,33),
                   # mu.bday=runif(1,44,46),
                   # sd.bday= 10,
                   # one=1)
                   start = runif(1,31,33),
                   mu.bday=runif(1,44,46),
                   sd.bday= runif(1,9,10),
                   one=1
                   )
################
model <- nimbleModel(code=code,inits=inits_list,data=data_list,
                    constants=constants_list,
                    dimensions=list(born=ndays,
                                    #wmObs=length(obsdays),
                                    #whitePups=ndays,moultPups=ndays,
                                # moultCount=length(obsdays), whiteCount=length(obsdays),
                    Pb_sum=ndays,Pm_sum=ndays,Pl_sum=ndays))
                   
nc <- 2   ;   ni <- 21000   ;   nb <- 1000   ;   nt <- 0
 # Parameters monitored
                   
 params = c('Nborn','cumsum_born',
            #'whitePups','moultPups',,#'totPups',
            #'whiteCount','moultCount')#),# 'vwhite','vmoult',
            'start','mu.bday', 'sd.bday')
            #,'mm','mw',)
 
 
 HgMCMC <- buildMCMC(model,monitors=params)
 runMCMC(HgMCMC,niter=1)#set =1 bc with one calculation inits are used
 #& can diagnose there
 # # #compile the model object first for the MCMC to compile 
#Cmodel<- compileNimble(model,showCompilerOutput = F)
#CHgMCMC <- compileNimble(HgMCMC,project=Cmodel)
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
