#####################
### NIMBLE Hg Sim
#####################
#getwd()
#StA drive file path
#work_dir<-("C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/06 - Scripts")
#setwd(work_dir)
#========================================================================
# add real data or use simulated data here 
#datafile <- "C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/04 - Raw Data/Hg_pup_counts_2018_temp.csv"
#Load the data and get the unique colony-year combinations.
# 1) Set up the loop to go through each combination.
# 2) Extract the data for the current colony-year.
# 4) Define the model inputs (constants, data, initial values).
# 5) Run the NIMBLE model for each iteration.
# 6) Store the results.
# #========================================================================
## load libraries, input data and functions
library(nimble)
library(coda)
#getwd()

nim_dskew <- nimbleFunction(
  run = function(x = double(1), param = double(1)) {
    skew.bday <- param[1]
    mu.bday   <- param[2]
    sd.bday   <- param[3]
    delta <- skew.bday/sqrt(1+(skew.bday^2))
    omega <- sqrt(sd.bday^2/(1-(2*delta^2)/pi))
    xi <- mu.bday - omega*delta*sqrt(2/pi)
    pb <- (2/omega)*dnorm((x-xi)/omega) * pnorm(skew.bday*((x-xi)/omega))
    returnType(double(1))
    return(pb)
  })



#generic function to calculate cumulative prob of being born ('days' are whatever you want to integrate over)
nim_pskew <- nimbleFunction(
  run = function(days = double(1), skew.bday = double(), mu.bday = double(), sd.bday = double()) {
    pb <- rep(0, length(days))
    pars <- c(skew.bday, mu.bday, sd.bday)
    for(t in 1:length(days)) {
      pb[t] <- nimIntegrate(nim_dskew, -Inf, days[t], pars)[1]
    }
    returnType(double(1))
    return(pb)
  })

intPb <- nimbleFunction(
  run = function(days = double(1), skew.bday = double(), mu.bday = double(), sd.bday = double()) {
    Pb <- rep(0, length(days))
    pars <- c(skew.bday, mu.bday, sd.bday)
    for(t in 1:length(days)) {
      Pb[t] <- nimIntegrate(nim_dskew, t-1, t, pars)[1]
    }
    returnType(double(1))
    return(Pb)
  })

code <- nimbleCode({
  
  # set prior on daily survival intercept and slope
  #surv.int.real ~ dbeta(175,1)
  surv.int <-  4    # log(surv.int.real/(1-surv.int.real))       # logit scale
  surv.slope <- 0.1 # ~ dunif(0,.5)                               # daily increase in survival (logit scale) 
  
  # set prior
  
  # set prior on detection and misclassification processes
  Pobs <- 0.95   # ~ dunif(0.9,1)
  PCMoult <-0.84 #~ dunif(0.8,0.9)
  
  # BIRTH DISTRIBUTION
  mu.bday ~ dunif(10,70) #dnorm(45,0.01)T(20,60)      ## Prior for mean of birth distribution (1 would mean born just before day 0 [t=1], 10 would mean born just before day 9 [t=10])
  sd.bday ~ dunif(1,20)                 ## Prior for standard deviation of birth distribution
  skew.bday ~ dunif(1,10)
  
  
  nborn ~ dunif(100,50000)       ## prior for mean number of births (total across season)
  
  mu.days.on.colony ~ dunif(min.age, max.age)       ## Prior for mean of the number of days on the colony
  sd.days.on.colony  ~ dunif(1, 15)        ## Prior for standard deviation of days on colony
  
  
  
  #=================================
  # set the probability of being moulted as a function of age
  #===================================
  
  pmoult.int <- -12   # ~ dnorm(5,1)               # probability of being moult at birth
  pmoult.slope <- 0.6 # ~ dunif(-0.5,0)          # daily decrease in the probability of being moult (increase in probability of moulting)
  for(a in 1:ndays){   # loop past the last possible age
    logit(pmoult[a]) <- pmoult.int + pmoult.slope*(a-1)    # deterministic for now
  }
  #=================================
  # set the probability of survival as a function of age
  #===================================
  
  for(a in 1:(max.age+1)){
    logit(surv[a]) <- surv.int + surv.slope*(a-1)
  }
  
  #=======================================
  # set the scheduling distributions 
  #=======================================
  
  # BIRTH DISTRIBUTION
  #normal 
  # nconst_Pb <- 1 - ( pnorm(0,mu.bday, sd.bday) + (1 - pnorm((max.day.birth),mu.bday, sd.bday) ) )  # constant for truncating the births distribution
  # Pb[1] <- 0  # can't be born after 0 time has passed in the season
  # Pb[2:(max.day.birth+1)] <- (pnorm(1:(max.day.birth), mu.bday, sd.bday) - pnorm(0:(max.day.birth-1), mu.bday, sd.bday))/nconst_Pb   ## Birth schedule

  #skew normal using generic function
  pb[1:(max.day.birth+1)] <- intPb(1:(max.day.birth+1), skew.bday, mu.bday, sd.bday)
  Pb[1:(max.day.birth+1)] <- pb[1:(max.day.birth+1)]/sum(pb[1:(max.day.birth+1)])
  
  # snconst_Pb <- 1 - (nim_pskew(0,skew.bday,mu.bday,sd.bday) + (1 - nim_pskew((max.day.birth+1),skew.bday,mu.bday,sd.bday) ) ) # constant for truncating the births distribution
  # Pb[1] <- 0  # can't enter after 0 time has passed in the season
  # Pb[2:(max.day.birth+1)] <- (nim_pskew(1:(max.day.birth),skew.bday,mu.bday,sd.bday) - nim_pskew(0:(max.day.birth-1),skew.bday,mu.bday,sd.bday))/snconst_Pb ## Birth schedule

  

  # TIME ON COLONY
  # define the probability mass of remaining on the colony for a certain number of days (having left the colony between the previous day and the current day)
  nconst_birth.to.leave <- 1 - ( pnorm(min.age,mu.days.on.colony, sd.days.on.colony) + (1 - pnorm(max.age,mu.days.on.colony, sd.days.on.colony) ) )
  for(day in 1:(min.age+1)){
    pmass.birth.to.leave[day] <- 0   # can not leave before minimum number of days
    cprob.remain[day] <- 1
  }
  for(day in (min.age+2):(max.age+1)){   # 1 is day of birth, 2 is 1 day after birth, so on
    pmass.birth.to.leave[day] <- (pnorm(day-1, mu.days.on.colony, sd.days.on.colony) -  pnorm(day-2, mu.days.on.colony, sd.days.on.colony))/nconst_birth.to.leave      # define probability of leaving on any given age in days
    cprob.remain[day] <- 1-sum(pmass.birth.to.leave[1:day])  # cumulative probability of remaining white for a specific number of days 
  }
  
  #=======================
  # place cohorts in the proper places over time
  #=======================
  
  # cohort.size ~ dmulti(Pb,trunc(nborn))  # might be able to put this back?
  
  for(c in 1:(max.day.birth+1)){
    cohort.size[c] ~ dpois(nborn*Pb[c])
  }
  
  for(c in 1:(max.day.birth+1)){    # loop through cohorts
    # SET ALL PRE-BIRTH ABUNDANCES TO ZERO
    for(a in 1:(c-1)){    # loop through all days prior to birth
      for(s in 1:3){      # loop through stages
        stageN_c[c,a,s] <- 0      # set to zero. the second dimension is converting to absolute time (vs cohort time)
      }
    }
    # SET ABUNDANCES WHILE POTENTIALLY ON COLONY
    stageN_c[c,c,1] <- cohort.size[c] * 1 * 1 * 1     # probability of being white after a days
    stageN_c[c,c,2] <- cohort.size[c] * 1 * 1 * 0   # moulted
    for(a in 2:min( (ndays-c),(max.age+1) ) ){ 
      stageN_c[c,c+a-1,1] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:a]) * (1-pmoult[a])     # probability of being white after a days
      stageN_c[c,c+a-1,2] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:a])* pmoult[a]          # moulted
    }
    
    for(a in (max.age+2):(ndays-c)){
      for(s in 1:2){      # loop through stages
        stageN_c[c,c+a-1,s] <- 0      #  set to zero after the maximum age
      }
    }
    
    for(a in 1:(ndays-c)){
      stageN_c[c,c+a-1,3] <- (nborn*Pb[c])-(stageN_c[c,c+a-1,1]+stageN_c[c,c+a-1,2])  # dead or left
    }
  }
  
  # compute the expected number in each stage for each day of the season (final state process step)
  
  for(t in 1:(ndays-1)){ # loop through days of the season
    stageN[t,1] <- sum(stageN_c[1:(max.day.birth+1),t,1])  # tot white across all cohorts
    stageN[t,2] <- sum(stageN_c[1:(max.day.birth+1),t,2])  # tot grey across all cohorts
    stageN[t,3] <- sum(stageN_c[1:(max.day.birth+1),t,3])  # tot leave across all cohorts
  }
  
  # ========================
  # Observation process
  # ========================
  
  # normal approximation
  
  # for(obs in 1:(length(obsdays))) {
  #   white.mu[obs] <- stageN[obsdays[obs],1] * Pobs
  #   white.sd[obs] <- pow(stageN[obsdays[obs],1] * Pobs,.5)*(1-Pobs)
  #   white.sd[obs] <- pow(white.sd[obs],-2)
  #   moult.mu[obs] <- stageN[obsdays[obs],2] * Pobs
  #   moult.sd[obs] <- pow(stageN[obsdays[obs],2] * Pobs,.5)*(1-Pobs)
  #   moult.sd[obs] <- pow(moult.sd[obs],-2)
  #   wObs[obs] ~ dnorm(white.mu[obs],white.sd[obs])
  #   mObs[obs] ~ dnorm(moult.mu[obs],moult.sd[obs])
  # }
  
  # for(obs in 1:(length(obsdays))) {
  #  wObs[obs] ~ dbinom(Pobs, trunc(stageN[obsdays[obs],1]))   # data node
  #  mObs[obs] ~ dbinom(Pobs, trunc(stageN[obsdays[obs],2]))   # data node
  # }
  
  
  
  WW_MW[1:(ndays-1)] <- stageN[1:(ndays-1),1]+stageN[1:(ndays-1),2]*(1-PCMoult)
  MM[1:(ndays-1)] <- stageN[1:(ndays-1),2]*PCMoult  #add misclassification to mObs
  
  for(obs in 1:(length(obsdays))) {
    wObs[obs] ~ dbinom(Pobs, trunc(WW_MW[obsdays[obs]]))   # data node
    mObs[obs] ~ dbinom(Pobs, trunc(MM[obsdays[obs]]))   # data node
  }
  
})

## supply data, constants and inits
## supply data, constants and inits
data_list=list(wObs=sim_wm$wObs,
               mObs=sim_wm$mObs)
               #Pb=sim_wm$Pb)

constants_list=list(     ndays=ndays,
                         max.age = max.age,
                         min.age = min.age,
                         max.day.birth = max.day.birth,
                         obsdays=seq(40,90,12)
)

inits_list <- list(     mu.bday = 45,
                        sd.bday = 15,
                        skew.bday=5,
                        mu.days.on.colony = 30, 
                        sd.days.on.colony= 10,  
                        nborn = 2000,
                        #Pobs = 0.95
                        cohort.size  = round(10000*(sim_wm$Pb)) 
)
################
model <- nimbleModel(code=code,inits=inits_list,data=data_list,
                     constants=constants_list
                     #dimensions=list(Pb=max.day.birth+1)
                     )
# ## make sure this gives you a real, negative number (not NaN or -Inf):
model$calculate()

#nc <- 1   ;   ni <- 10000   ;   nb <- 1000  ;   nt <- 5

# Parameters monitored
params <- c( "mu.bday", "sd.bday","nborn",'skew.bday',
             "mu.days.on.colony","sd.days.on.colony",
             #"surv.int","surv.slope","pwhite.int","pwhite.slope",
             #"pmass.birth.to.leave","cprob.remain","Pb",
             "stageN","cohort.size")


DEBUG=T

if(DEBUG){

  nc <- 1   ;   ni <- 10000   ;   nb <- 5000   ;   nt <- 10 ; na <-1000


  mcmc.out <-nimbleMCMC(model=model,nchains =nc,niter = ni, nburnin = nb,
                                  monitors = params, samplesAsCodaMCMC = TRUE,summary=TRUE)
}else{

  nc <- 3   ;   ni <- 80000   ;   nb <- 10000   ;   nt <- 5 ; na <-1000

  system.time(mcmc.out <-nimbleMCMC(model=model,nchains =nc,niter = ni, nburnin = nb,
                                    monitors = params, samplesAsCodaMCMC = TRUE,summary=TRUE)
              )
}
# save results 
mcmc<-mcmc.out$samples

#update without restarting MCMC

# Cmcmc$run(10000,reset=FALSE)
# samples_update<- as.matrix(HgMCMC$mvSamples)
# plot(samples_update[ , 'nborn'], type = 'l', xlab = 'iteration',
#      ylab = expression(Nborn))
# mean(samples_update[,'Nborn'])

library(MCMCvis)
# par(mfrow=c(3,3))
MCMCtrace(object = mcmc,
            pdf = FALSE,
            ind = TRUE, # separate density lines per chain
            Rhat=F,
            n.eff=F,
            exact = TRUE,
            #gvals = GV,
            params = c("mu.bday", "sd.bday","nborn","skew.bday",
             "mu.days.on.colony","sd.days.on.colony"))

# GV <- c(500, 45, 10,5)
# #GV <- c(500)
#  MCMCtrace(object = samples,
#            pdf = FALSE, 
#            type='density',
#            ind = TRUE, # separate density lines per chain
#            Rhat=F,
#            n.eff=F,
#            exact = TRUE,
#           # priors = PR,
#            gvals = GV,
#            post_zm = T,
#            params = c("nborn","mu.bday","sd.bday","mu.days.on.colony"),
#           main_den = c("Total pups born over season", "Mean birth date", "Standard deviation of birth date", "Skew birth date"))
# 