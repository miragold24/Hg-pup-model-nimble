#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
library(nimble)
#library(coda)
#getwd()

#simulation code
#source("HgSim_S4.R")
# Nimble functions to count whitecoat and moulted pups present on any given day
source("HgNimble_input_S4.R") 
source('leaveFunction.R')
# EXAMPLE

dskew <- nimbleFunction(
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

int_pb <- nimbleFunction(
  run = function(ndays = double(), skew.bday = double(), mu.bday = double(), sd.bday = double()) {
    pb <- numeric(ndays)
    pars <- c(skew.bday, mu.bday, sd.bday)
    for(t in 1:ndays) {
      pb[t] <- nimIntegrate(dskew, t-1, t, pars)[1]
    }
    returnType(double(1))
    return(pb)
  })

code <- nimbleCode({
   #priors
   #gamma prior chosen for Nborn in discrete model because gamma is continuous
   Nborn ~ dgamma(shape=485, scale=1) 
   #mu.bday ~ T(dnorm(45,10),33,70)
    # mu.bday ~ dunif(33,60)
    # skew.bday ~ dunif(2,6)
   #  mu.bday ~ dunif(0,ndays)
   #  skew.bday ~ dunif(0,ndays)
   #  sd.bday ~ dunif(0, ndays)   
   # #shifted scale beta for mu.bday,sd.bday,skew.bday parameters
   standard_beta.m ~ dbeta(alpha.m, beta.m)
   standard_beta.sd ~ dbeta(alpha.sd, beta.sd)
   standard_beta.s ~ dbeta(alpha.s, beta.s)
   # # Shifting and scaling the distribution
   mu.bday <- shift.m + (scale.m * standard_beta.m)
   sd.bday <- shift.sd + (scale.sd * standard_beta.sd)
   skew.bday <- shift.s + (scale.s * standard_beta.s)
   # 
   
  #constraints
  # one ~ dconstraint(start < mu.bday)
  # deterministic and continuous births over breeding season (ndays)
   
  ## for adding a mortality rate, could multiply Pb_sum[t]*mortality rate
  for (t in 1:ndays){ 
    born[t] <- Pb_sum[t]*Nborn
  }
  
  #=============================================================================
  ### Calculate birth probabilities 


  #calculate PDF of skew normal then integrate over intervals for Pb on given day
     Pb[1:ndays] <- int_pb(ndays, skew.bday, mu.bday, sd.bday)

 #==============================================================================
  
  Pb_sum[1:ndays] <- Pb[1:ndays]/sum(Pb[1:ndays])
  

  # cumulative sum of pups born on given day 
  # for (i in 1:ndays) {
  #   cumsum_born[i] <- sum(born[1:i])
  # }
   
   
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
 
  
  #=================================================================
  # OBSERVATION PROCESS
  #=================================================================
  # binomial with obs error only 
  # for (obs in 1:length(obsdays)) {
  #   wObs[obs] <- dbinom(round(whitePups[obsdays[obs]] + S4Pups[obsdays[obs]]), prob = POWhite)
  #   mObs[obs] <- dbinom(round(whitePups[obsdays[obs]] + S4Pups[obsdays[obs]]), prob = POWhite) 
  # }
  


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
    wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs])
    mObs[obs] ~ dnorm(moultCount[obs], sd = sdmoult[obs])
  #   # normal distribution with small variance (should get 500 pups when
  #   # exact obs process in simulation):
    # wObs[obs] ~ dnorm(whiteCount[obs], sd = 0.001)
    # mObs[obs] ~ dnorm(moultCount[obs], sd = 0.001)
  # 
  #     
   }
  
  

})

## supply data, constants and inits
data_list=list(wObs=sim_S4$wObs,mObs=sim_S4$mObs,Pl_sum=PmPl2$PlGivenm)

constants_list=list(ndays=ndays,amax=amax,obsdays=obsdays,#nslice=1000,
                    # Parameters for the standard beta distribution
                    #mu bday
                    alpha.m=2,
                    beta.m=2,
                    shift.m=20, # Shift the distribution to start at 33
                    # The "scale" parameter in a shifted beta distribution
                    # multiplies all values in the original beta distribution by 
                    # a constant factor. 
                    scale.m=30, 
                    # sd.bday
                    alpha.sd=2,  # Shape parameter (alpha > 0)
                    beta.sd=2,    # Shape parameter (beta > 0)
                    shift.sd=6,   # Shift the distribution to start at 6
                    scale.sd=5,   # 'stretch' to 11
                    # skew.bday
                    # alpha.s=5,
                    # beta.s=5,
                    # shift.s=2,  
                    # scale.s=9,  
                    alpha.s = 2,
                    beta.s = 2,
                    shift.s = 1,
                    scale.s = 5,
                    mu.mday=mu.mday,sd.mday=sd.mday,mu.lday=mu.lday,sd.lday=sd.lday,
                    mu.s4day=mu.s4day,sd.s4day=sd.s4day,
                    POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult,PCWhiteS4=PCWhiteS4)



inits_list <- list(Nborn=runif(1,500,520),Pb=rep(0,ndays),
                   Pb_sum=rep(0,ndays),PS4=rep(0,amax),Pm=rep(0,amax),
                   #start = runif(1,25,45),
                    #mu.bday=ndays/5,
                    # mu.bday=runif(1,33,60),
                    # skew.bday=runif(1,2,6),
                     #skew.bday=1,
                   standard_beta.m=rbeta(1,2,2),
                   standard_beta.s=rbeta(1,2,2),
                   # sd.bday= 5,
                    standard_beta.sd=rbeta(1,2,2)
                   # #one=1
                   )
################
model <- nimbleModel(code=code,inits=inits_list,data=data_list,
                    constants=constants_list,
                    dimensions=list(born=ndays,Pb_sum=ndays,PS4_sum=ndays,
                                    Pm_sum=ndays))
                   
nc <- 2   ;   ni <- 5000   ;   nb <- 1000  ;   nt <- 0
 
# Parameters monitored
params = c('Nborn','mu.bday', 'sd.bday','skew.bday')
 
 HgMCMC <- buildMCMC(model,monitors=params)
 ## make sure this gives you a real, negative number (not NaN or -Inf):
 model$calculate()
 
# Rmcmc <- buildMCMC(conf)
 
model <- compileNimble(model)
Cmcmc <- compileNimble(HgMCMC, project =model)
 
 set.seed(0)
 samples <- runMCMC(Cmcmc, nchains=3, niter = 20000, nburnin = 10000, thin = 5)
 
 ## adding priors and simulated values to posterior plots

 # Number of samples
n_samples <- 6000

alpha.sd <- 2
beta.sd <- 2
shift.sd <- 6
scale.sd <- 5

alpha.m <- 2
beta.m <- 2
shift.m <- 10
scale.m <- 40

alpha.s=2
beta.s=2
shift.s=1
scale.s=5
Nborn.prior<-rgamma(n_samples,shape=485, scale=1)
mubday.prior <- rbeta(n_samples, alpha.m, beta.m) * scale.m + shift.m
sdbday.prior <- rbeta(n_samples, alpha.sd, beta.sd) * scale.sd + shift.sd
skew.prior <- rbeta(n_samples, alpha.s, beta.s) * scale.s + shift.s
# mubday.prior<-runif(10000,33,60)
# sd.prior<-runif(10000,7,12)
# skew.prior<-runif(10000,4,8)
PR=cbind(Nborn.prior,mubday.prior,sdbday.prior,skew.prior)
GV <- c(500, 45, 10,5)
#GV <- c(500)
 MCMCtrace(object = samples,
           pdf = FALSE, 
           type='density',
           ind = TRUE, # separate density lines per chain
           Rhat=F,
           n.eff=F,
           exact = TRUE,
          # priors = PR,
           gvals = GV,
           post_zm = T,
           params = "all",
          main_den = c("Total pups born over season", "Mean birth date", "Standard deviation of birth date", "Skew birth date"))

 MCMCsummary(samples,params="all")
#runMCMC(HgMCMC,niter=500)#set =1 bc with one calculation inits are used
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
# mcmc.out <-nimbleMCMC(model=model,nchains =nc,niter = ni, nburnin = nb,
#                   monitors = params, samplesAsCodaMCMC = TRUE,summary=TRUE)
 
## check convergence and look at plots 
# mcmcS4<-mcmc.out$samples
# mcmcs4allobs<-mcmc.out$samples
#  coda::gelman.diag(mcmcS4[,"Nborn"])
# # #look at trace plots
#   plot(mcmcS4[,"Nborn"])
#   coda::gelman.diag(mcmc[,"mu.bday"])
##  look at trace plots
#   plot(mcmcS4[,"mu.bday"])
#   library(MCMCvis)
# coda::gelman.diag(mcmcS4[,"mu.bday"])

 #trace and posterior density of params:
 
 # generating values for each parameter used to simulate data
# GV <- c(500, 45, 10,5)
#   MCMCtrace(object = mcmcs4allobs,
#             pdf = FALSE, 
#             ind = TRUE, # separate density lines per chain
#             Rhat=F,
#             n.eff=F,
#             exact = TRUE,
#             #gvals = GV,
#             params = "all")
#   
# MCMCsummary(mcmcs4allobs,params="all")
# par(mfrow=c(1,1))
# Define colorblind-friendly palette
# cbPalette <- c("#E69F00", "#cc8400", "#56B4E9","#3c7da3")
# # Plot the lines for whitePups and S4Pups
# plot(model$whitePups + S4Pups, typ = 'l', col = cbPalette[1],lwd=2, ylab = '# pups', xlab = 'day of season', main = 'Present on colony vs observed pups')
# lines(model$moultPups, col = cbPalette[3],lwd=2)
# 
# # Add points for whiteCount and moultCount on observation days
# #obsdays <- c(35, 50, 75, 90, 105, 120)
# points(x = obsdays, y = model$whiteCount, col = "#cc8400", pch = 19)
# points(x = obsdays, y = model$moultCount, col = "#3c7da3", pch = 19)
# 
## Create legend with points
# legend("topright", legend = c("White present", "White observed", "Moult present", "Moult observed"), col = cbPalette, lty = c(1,NA,1,NA), lwd=2, pch = c(NA, 19, NA, 19))

