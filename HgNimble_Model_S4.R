#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
#library(nimble)
#library(MASS) 
library(brms)
library(coda)
#getwd() 
#source("HgSim_S4.R")
#source("HgNimble_input_S4.R") # Nimble functions to count whitecoat and moulted pups present on any given day

# EXAMPLE

code <- nimbleCode({
   #priors
   #Nborn ~ dgamma(shape=500,scale=1) #Nborn (N)=500 in simulation
   Nborn ~ dgamma(shape=485, scale=1) 
   #start ~ dunif(25,45)
   mu.bday ~ dunif(33,60)
   skew.bday ~ dunif(3,7)
   #sd.bday ~ dunif(7,12) 
   #shifted scale beta for sd.bday
   # beta distribution
   standard_beta ~ dbeta(alpha, beta)  
   # Shifting and scaling the distribution
   sd.bday <- shift + (scale * standard_beta)
   
   
  #constraints
    # one ~ dconstraint(start < mu.bday)
  # alt constraint - truncated normal for start
     #start ~ T(dnorm(0, 1), 0, mu.bday)
     #T(dnorm(0, 1), -15, mu.bday) #start is within 15 days of mu.bday
  # deterministic and continuous births over breeding season (ndays)
  for (t in 1:ndays){ 
    born[t] <- Pb_sum[t]*Nborn
  }
  
  #=============================================================================
  ### Calculate birth probabilities 
   
  ###
  ### Log normal
  ###
  ### derived quantities for log scale calculations 
  ### 
  # lCV <- sd.bday/(mu.bday-start)
  # lsd <- sqrt(log(lCV^2+1))
  # lmu <- log(mu.bday-start)-(0.5*(lsd^2))
  # 
  # # Calculate probability of being born on any given day
  # for(t in 2:ndays){
  #   Pb[t] <- plnorm((t-start+1), lmu, lsd) - plnorm((t-start), lmu, lsd)
  #   
  # }
  # 

   ###
   ### Skew normal
   ###
   #using BRMS package (not estimating mu.bday and sd.bday)
   # for(t in 1:ndays){
   #   Pb[t] <- pskew_normal(q = t+1, mu = mu.bday, sigma = sd.bday, alpha = skew.bday) -
   #     pskew_normal(q = t, mu = mu.bday, sigma = sd.bday, alpha = skew.bday)
   # }
   # 
   # 
   # 
   # 
  # convert paramaters for skew normal
   delta <- skew.bday/sqrt(1+(skew.bday^2)) #location
   omega<-sqrt(sd.bday^2/(1-(2*delta^2)/pi)) #scaling factor
   xi <- mu.bday - omega*delta*sqrt(2/pi) #skewness
   # #x<- 1:123
   
   # #calculates Pb occurring at a certain point and
   # #being less than or equal to that point
   # #joint probability of a value being at (x - xi)/omega and
   # #being less than or equal to skew.bday * ((x - xi)/omega)
   #Pb_sum<- rep(0,ndays)
   for (x in 1:ndays){
     # next, calc the PDF of the skew normal (*scaling factor omega)
     #at the point (x - xi) / omega
     Pb[x]<- (2/omega)*dnorm((x-xi)/omega) *
       #calc CDF up to the point skew.bday * ((x - xi)/omega)
       pnorm(skew.bday*((x-xi)/omega))
   }
   
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
  moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax],moultObs[])
  #totPups[1:ndays] <- whitePups[1:ndays]+moultPups[1:ndays]
  


  #################################################################
  # OBSERVATION PROCESS

  #####################################################################
 # binomial variance formula
  #N *p*(1-p)
  # n (number of trials)* p(prob of success) * (1 - p) (prob of failure)
  for (obs in 1:length(obsdays)) {
    whiteCount[obs] <- whitePups[obsdays[obs]] * POWhite + moultPups[obsdays[obs]] * POMoult * (1 - PCMoult) + S4Pups[obsdays[obs]]*POWhite*PCWhiteS4

    sdwhite[obs] <- sqrt(whitePups[obsdays[obs]] * POWhite * (1 - POWhite) + moultPups[obsdays[obs]] * POMoult*(1 - PCMoult)* (1 - POMoult*(1 - PCMoult))+
      S4Pups[obsdays[obs]] * POWhite * PCWhiteS4 * (1 - POWhite*PCWhiteS4))

    moultCount[obs] <- moultPups[obsdays[obs]] * POMoult * PCMoult + S4Pups[obsdays[obs]]*POWhite*(1-PCWhiteS4)
    sdmoult[obs] <- sqrt(moultPups[obsdays[obs]] * POMoult * PCMoult * (1 - POMoult * PCMoult)+ 
      S4Pups[obsdays[obs]]*POWhite*(1-PCWhiteS4)*(1 - POWhite*(1-PCWhiteS4)))

    
    # normal approximation using binomial variance
    wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs])
    mObs[obs] ~ dnorm(moultCount[obs], sd = sdmoult[obs])
    # obs ~ normal distribution using tiny variance (should get 500 pups when exact obs process):
    # wObs[obs] ~ dnorm(whiteCount[obs], sd = 0.001)
    # mObs[obs] ~ dnorm(moultCount[obs], sd = 0.001)

     #==================================
     #LOG SUM EXP 
     #==================================
     #  # Integration range
      # a[obs] <- ceiling(whiteCount[obs]) 
      # b[obs] <- floor(whiteCount[obs])
      # 
      # # Logsumexp trick for white pups
      # c[obs] <- ld(mean = whiteCount[obs], sd = sdwhite[obs], x = (b[obs] + a[obs]) / 2)
      # 
      # 
      # for (x in seq(a[obs], b[obs])) {
      #   sum_ld[obs] <- c[obs] + log(exp(ld(x = x, mean = whiteCount[obs], sd = sdwhite[obs]) - c[obs]))+log(0.5)
      # }
      # # Use the integrated density in the likelihood
      # #wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs]) # Use the original distribution here
      # #wObsLogLik[obs] <- wObs[obs] +ld_w # Add the integrated density to the log-likelihood
      # wObs[obs] ~ dnorm(whiteCount[obs], sdwhite[obs])
      # wObsLogLik[obs] <- wObs[obs] + sum_ld[obs]
      # 
       ###########################
      # # # Integration range
      # e[obs] <- ceiling(moultCount[obs])
      # f[obs] <- floor(moultCount[obs])
      # 
      # # Logsumexp trick for moult pups
      # cm[obs] <- ld(mean = moultCount[obs], sd = sdmoult[obs], x = (e[obs] + f[obs]) / 2)

      # for (x in seq(e[obs], f[obs])) {
      #   sum_ldm[obs] <- cm[obs] + log(exp(ld(x = x, mean = moultCount[obs], sd = sdmoult[obs]) - cm[obs]))+log(0.5)
      # }
      # # Use the integrated density in the likelihood
      # #wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs]) # Use the original distribution here
      # #wObsLogLik[obs] <- wObs[obs] +ld_w # Add the integrated density to the log-likelihood
      # mObs[obs] ~ dnorm(moultCount[obs], sdmoult[obs])
      # mObsLogLik[obs] <- mObs[obs] + sum_ldm[obs]
      
    }
  
  
  


})

## supply data, constants and inits
data_list=list(wObs=sim_S4$whiteObs,mObs=sim_S4$moultObs,moultObs=sim_S4$moultObs,Pl_sum=PmPl2$PlGivenm)

constants_list=list(ndays=ndays,amax=amax,obsdays=obsdays,
                    # Parameters for the standard beta distribution
                    alpha=2 ,  # Shape parameter (alpha > 0)
                    beta=5,    # Shape parameter (beta > 0)
                    shift=7,   # Shift the distribution to start at 7
                    scale=5,   # Scale the distribution to end at 12 (12-7)
                    #start=start,
                    #skew.bday=skew.bday,
                    #mu.bday=mu.bday,
                    sigma=sd.bday,
                    pi=pi,
                    #sd.bday=sd.bday, #params to be estimated
                    mu.mday=mu.mday,sd.mday=sd.mday,mu.lday=mu.lday,sd.lday=sd.lday,
                    mu.s4day=mu.s4day,sd.s4day=sd.s4day,
                    POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult,PCWhiteS4=PCWhiteS4)



inits_list <- list(Nborn=runif(1,500,510),Pb=rep(0,ndays),
                   Pb_sum=rep(0,ndays),PS4=rep(0,amax),Pm=rep(0,amax),
                   #start = runif(1,25,45),
                    mu.bday=runif(1,33,60),
                    #sd.bday= runif(1,7,12),
                    #one=1,
                    standard_beta=rbeta(1,2,5) 
                   )
################
model <- nimbleModel(code=code,inits=inits_list,data=data_list,
                    constants=constants_list,
                    dimensions=list(born=ndays,
                                #whitePups=ndays,moultPups=ndays,
                                #moultCount=length(obsdays), 
                                #whiteCount=length(obsdays),
                    Pb_sum=ndays,PS4_sum=ndays,Pm_sum=ndays))
                   
nc <- 2   ;   ni <- 16000   ;   nb <- 5000   ;   nt <- 0
 # Parameters monitored
                   
 #' params = c('cumsum_born','Nborn',
 #'            'whitePups','moultPups',#'totPups',
 #'            #'whiteCount','moultCount',# 'vwhite','vmoult',
 #'            'start','mu.bday', 'sd.bday'
 #'            #,'mm','mw',)
 #' )

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
  