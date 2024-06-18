############################################
###
### The problem with this code is every new data set (so each simulation) recompiles.
### trouble updating the data outside the 
###
############################################

# setwd(work_dir)
#library(tidyverse)
library(nimble)
library(coda)
#source("HgSim_recode9_wm.R")
#num_simulations <- 10  # Define the number of simulations


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
#dskew_pb <- compileNimble(nim_dskew)
#nim_dskew(1:100,c(skew.bday,mu.bday,sd.bday))
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
  skew.bday <- 5 #dunif(1,10)
  
  
  nborn ~ dunif(100,50000)       ## prior for mean number of births (total across season)
  
  mu.days.on.colony ~ dunif(min.age, max.age)       ## Prior for mean of the number of days on the colony
  sd.days.on.colony <- 2 #~ dunif(4, 15)        ## Prior for standard deviation of days on colony
  
  
  
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
  #Pb<-rep(0,(max.day.birth+1))
  pb[1:(max.day.birth+1)] <- intPb(1:(max.day.birth+1), skew.bday, mu.bday, sd.bday)
  Pb[1:(max.day.birth+1)] <- pb[1:(max.day.birth+1)]/sum(pb[1:(max.day.birth+1)])
  
  # snconst_Pb <- 1 - (nim_pskew(0,skew.bday,mu.bday,sd.bday) + (1 - nim_pskew((max.day.birth),skew.bday,mu.bday,sd.bday) ) ) # constant for truncating the births distribution
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
      stageN_c[c,c+a-1,3] <- (nborn*Pb[c])-(stageN_c[c,c+a-1,1]+stageN_c[c,c+a-1,2])  # left
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
  # 
  # for (obs in 1:length(obsdays)) {
  #   whitemu[obs] <-  stageN[obsdays[obs],1]*Pobs + stageN[obsdays[obs],2]*Pobs*(1-PCMoult)
  #   
  #   sdwhite[obs] <- sqrt(stageN[obsdays[obs],1]* Pobs * (1 - Pobs) +
  #                          stageN[obsdays[obs],2]*Pobs*(1-PCMoult)* (1 - Pobs*(1 - PCMoult)))
  #                        
  #   moultmu[obs] <- stageN[obsdays[obs],2]*Pobs*PCMoult
  #                        
  #   sdmoult[obs] <- sqrt(stageN[obsdays[obs],2]*Pobs*PCMoult * (1 - Pobs * PCMoult))
  #                        
  #                        
  #    # normal approximation using binomial variance
  #    wObs[obs] ~ dnorm(whitemu[obs], sd = sdwhite[obs])
  #    mObs[obs] ~ dnorm(moultmu[obs], sd = sdmoult[obs])
  #    }
  
  
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


# Define the list of different obsdays combinations
obsdays_combinations <- list(
  seq(40, 80, by = 10),
  seq(60, 100, by = 10),
  seq(40, 90, by = 12),
  seq(30, 80, by = 12)
)


ndays <- 123
max.age <- 60
min.age <- 20
max.day.birth <- 100

# results_list <- list()
#dir.create("results")

# Loop through each combination and run the generate_obs function and NIMBLE model
for (obsdays in obsdays_combinations) {
  for (sim in 1:10) {
    cat("Running simulation", sim, "for obsdays", obsdays, "\n")
    
    # Generate observations using the generate_obs function
    sim_wm <- generate_obs(
      ndays = 123,
      min.age = 20,
      max.age = 60,
      max.day.birth = 100,
      obsdays = obsdays,
      mu.bday = 45,
      sd.bday = 10,
      skew.bday = 5,
      nborn = 700,
      mu.days.on.colony = 31.5,
      sd.days.on.colony = 2,
      surv.int = 4,
      surv.slope = 0.1,
      pmoult.int = -12,
      pmoult.slope = 0.6,
      Pobs = 0.95,
      PCMoult = 0.84
    )
    
    wObs <- sim_wm$wObs
    mObs <- sim_wm$mObs
    
    # Set data
    data <- list(
      wObs = wObs,
      mObs = mObs
    )
    
    # Define initial values
    inits_list <- list(     mu.bday = 45,
                            sd.bday = 15,
                            skew.bday=5,
                            mu.days.on.colony = 30, 
                            sd.days.on.colony= 10,  
                            nborn = 4000,
                            #Pobs = 0.95
                            cohort.size  = round(10000*(sim_wm$Pb)) 
    )
    
    # Constants will include current observation day scenario
    constants_list <- list(ndays = ndays,
                           max.age = max.age,
                           min.age = min.age,
                           max.day.birth = max.day.birth,
                           obsdays = obsdays)
    
    # Define parameters to monitor
    params <- c(
      "mu.bday", "sd.bday", "nborn",#,skew.bday,
      "mu.days.on.colony",# "sd.days.on.colony",
      "stageN", "cohort.size"
    )
    

    
    DEBUG <- T  
    if (DEBUG) {
      nc <- 1; ni <- 5000; nb <- 1000; nt <- 10
    } else {
      nc <- 3; ni <- 500000; nb <- 100000; nt <- 10
    }
    
    # Initialize and compile the NIMBLE model
    model <- nimbleModel(   # speed up compilation: can use calcualte=FALSE to disable calc of deterministics nodes and log-likelihood
      code = code,
      inits = inits_list,
      data = data,
      constants = constants_list
    )
    
    Cmodel <- compileNimble(model)
    
    # Configure and compile MCMC
    mcmc_conf <- configureMCMC(model, monitors = params)
    Rmcmc <- buildMCMC(mcmc_conf)
    Cmcmc <- compileNimble(Rmcmc, project = model)
    
    # Run the MCMC and measure execution time
    system.time(
      mcmc_samples <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, thin = nt, setSeed = TRUE)
    )
    
    # Save each result to an RDS file
    saveRDS(mcmc_samples, file = paste0("results/obsdays_", paste(obsdays, collapse = "_"), "_sim_", sim, "_", Sys.Date(), ".rds"))
  }
}



# 
# 
# results_dir <- "results/"
# # List all RDS files in the directory
# rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
# 
# # Initialize a list to store all loaded data frames
# results_list <- lapply(rds_files, readRDS)
# # Optionally save all results to a file
