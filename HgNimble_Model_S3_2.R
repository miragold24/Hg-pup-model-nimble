#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
library(nimble)
library(dplyr)
#StA drive file path
work_dir<-("C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/07 - GitHub Repositories/Hg_pup_model")
setwd(work_dir)
#getwd() 

###################################
###
### Choose data 
###
#simulated
#source("HgSim_S3.R")
source("HgNimble_input_S3.R") # Nimble functions to count whitecoat and moulted pups present on any given day

#real data 
#"C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/04 - Raw Data/
datafile <- "HgPupCountsForProdModel_2018_NEEfs.csv"
data <- read.csv(datafile)
colyears <- data %>% select(Colony, Year) %>% distinct()
# Convert Date column to Date type
data$Date <- as.Date(data$Date, format = "%d/%m/%Y")
# Calculate the day number from October 1
start_date <- as.Date("01/10/2018", format = "%d/%m/%Y")
data$Day <- as.numeric(difftime(data$Date, start_date, units = "days")) + 1
###################################


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
  Nborn ~ dunif(100,50000) 
  #mean.nborn ~  dgamma(shape=485, scale=1) 
  #Nborn ~ dpois(mean.nborn)  
  #Nborn ~ dgamma(shape=485, scale=1) 
  mu.bday ~ dunif(20,60)
  skew.bday ~ dunif(2,6)
  sd.bday ~ dunif(5,15)
  #shifted scale beta 
  #standard_beta.m ~ dbeta(alpha.m, beta.m)
  # standard_beta ~ dbeta(alpha, beta)  
  # #standard_beta.s ~ dbeta(alpha.s, beta.s)
  # #Shifting and scaling the distribution
  # mu.bday <- shift.m + (scale.m * standard_beta)
  # sd.bday <- shift + (scale * standard_beta)
  # skew.bday <- shift.s + (scale.s * standard_beta)
  
  
  # constraints
  #one ~ dconstraint(start < mu.bday)
  
  # deterministic and continuous births over breeding season (ndays)
  for (t in 1:ndays){ 
    born[t] <- Pb_sum[t]*Nborn
  }

  #calculate PDF of skew normal then integrate over intervals for Pb on given day
  Pb[1:ndays] <- int_pb(ndays,skew.bday,mu.bday,sd.bday)
  
   Pb_sum[1:ndays] <- Pb[1:ndays]/sum(Pb[1:ndays])
   
  ###
  ###  Probability of moulting across season
  ##
  
  # Probability of moulting as function of age (Probability of finishing moult)
  Pm[1:amax] <- pnorm(1:amax, mu.mday, sd.mday) - pnorm(0:(amax - 1), mu.mday, sd.mday)
  #Ensure probabilities sum to 1
  Pm_sum[1:amax] <- Pm[1:amax]/sum(Pm[1:amax])
  
  # # *Use PmPl2.PlGivenm (prob of leaving given moulted) in leaveFunction.R
  Pl[1:amax] <- pnorm(1:amax, mu.lday, sd.lday) - pnorm(0:(amax - 1), mu.lday, sd.lday)
  #Ensure probabilities sum to 1
  Pl_sum[1:amax] <- Pl[1:amax] / sum(Pl[1:amax])
  
  
  # calculate how many white and moulted pups present on any given day, given age
  # Nimble functions in HGNimble_input_noS4.R
  whitePups[1:ndays] <- whiteSum(born[1:ndays], Pm_sum[1:amax],whiteObs[])#,Pl_sum[1:amax],
  moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax],moultObs[])#,
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
  
  # #   # add misclassification
  #   # whiteMC[obs] <- whitePups[obsdays[obs]] + moultPups[obsdays[obs]] * (1-PCMoult)
  #   # moultMC[obs] <- moultPups[obsdays[obs]] * PCMoult
  #   # wObs[obs] ~ dbinom(POWhite,round(whiteMC[obs]))
  #   # mObs[obs] ~ dbinom(POMoult,round(moultMC[obs]))
  # 
  
  #####################################################################
  ###
  ### Normal approximation using binomial variance
  ###
  #####################################################################
  for (obs in 1:length(obsdays)) {
    whiteCount[obs] <- whitePups[obsdays[obs]] * POWhite + moultPups[obsdays[obs]] * POMoult * (1 - PCMoult)
    sdwhite[obs] <- sqrt(whitePups[obsdays[obs]] * POWhite * (1 - POWhite) + moultPups[obsdays[obs]] * (1 - PCMoult) * POMoult * (1 - POMoult * (1 - PCMoult)))

    moultCount[obs] <- moultPups[obsdays[obs]] * POMoult * PCMoult
    sdmoult[obs] <- sqrt(moultPups[obsdays[obs]] * POMoult * PCMoult * (1 - POMoult * PCMoult))

     wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs])
     mObs[obs] ~ dnorm(moultCount[obs], sd = sdmoult[obs])

    ## small variance and exact obs process :
    # wObs[obs] ~ dnorm(whiteCount[obs],sd=0.1)
    # mObs[obs] ~ dnorm(moultCount[obs],sd=0.1)
    #


  }
  
})

## supply data, constants and inits

## single colony run 
# Extract data for the specified colony and year
selected_colony="Farne Islands"
selected_year="2018"
realdata <- subset(data, Colony == selected_colony & Year == selected_year)
wObs <- realdata$Whitecoats + realdata$DeadPups
mObs <- realdata$Moulteds
obsdays <- unique(realdata$Day)

###################################################################
# Define constants
amax <- 123 
ndays<-123
mu.mday <- 23
sd.mday <- 5
mu.lday <- 31.5
sd.lday <- 7
POWhite <- 0.95  
POMoult <- 0.95  
PCMoult <- 0.82  
# Directory where you want to save the .rds files from laptop
#output_dir <- "C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/05 - Results/real data"

# Directory where you want to save the .rds files from remote desktop:
output_dir <- "D:/mg349/results"


#simulated ddata
# data_list=list(wObs=sim_S3$wObs,mObs=sim_S3$mObs,whiteObs=sim_S3$wObs,moultObs=sim_S3$mObs,
#                Pl_sum=sim_S3$Pl)
#real data 
data_list=list(wObs=wObs,mObs=mObs,whiteObs=wObs,moultObs=mObs)

constants_list=list(ndays=ndays,amax=amax,obsdays=obsdays,
                    mu.mday=mu.mday,sd.mday=sd.mday,
                    mu.lday=mu.lday,sd.lday=sd.lday,
                    POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult
                      )
                    # constants for the standard beta distribution
                    #mu bday
                    # alpha.m=4,
                    # beta.m=4,
                    # shift.m=10, # Shift the distribution to start at 33
                    # # The "scale" parameter in a shifted beta distribution
                    # # multiplies all values in the original beta distribution by
                    # # a constant factor.
                    # scale.m=45,
                    # # sd.bday
                    # alpha=5,  # Shape parameter (alpha > 0)
                    # beta=5,    # Shape parameter (beta > 0)
                    # shift=5,
                    # scale=7,
                    # # skew.bday
                    # # alpha.s=5,
                    # # beta.s=5,
                    # shift.s=1,  # Shift the distribution to start at 2
                    # scale.s=6,  # Scale the distribution to end at 7 (7-2)
                


inits_list <- list(Nborn=runif(1,100,50000),
                   # Nborn =500,
                   # mean.nborn = ndays*100,
                   # Nborn = ndays*100,
                   mu.bday= runif(1,20,60),
                   skew.bday=runif(1,2,6),
                   sd.bday=runif(1,5,15)
                   #skew.bday=runif(1,2,6),
                   #mu.bday=rnorm(1,45,60),S
                   #sd.bday= runif(1,7,12),
                   # standard_beta.m=rbeta(1,4,4),
                   # standard_beta.s=rbeta(1,5,5),
                   # #sd.bday= runif(1,7,12),
                   # #one=1,
                   #standard_beta=rbeta(1,5,5)
)

################
model <- nimbleModel(code=code,
                     inits=inits_list,
                     data=data_list,
                     constants=constants_list,
                    dimensions=list(born=ndays,
                    Pb_sum=amax,Pm_sum=amax,Pl_sum=amax))
model$calculate()


#nc <- 1   ;   ni <- 10000   ;   nb <- 5000   ;   nt <- 10


# # Parameters monitored
params = c('Nborn','mu.bday', 'sd.bday',"skew.bday")
# 
# #Build and compile MCMC, update without reset
HgMCMC <- buildMCMC(model,monitors=params)

#Cmodel <- compileNimble(model)
#Cmcmc <- compileNimble(HgMCMC, project = Cmodel)
set.seed(0)
conf <- configureMCMC(model, useConjugacy = FALSE)
Rmcmc <- buildMCMC(conf)  # Create MCMC configuration based on the compiled model
Cmodel <- compileNimble(model)  # Compile the NIMBLE model, translating it to C++
Cmcmc <- compileNimble(Rmcmc, project = model)  # Compile the MCMC configuration and link it to the already compiled model
runMCMC(HgMCMC,niter=1)#set =1 bc with one calculation inits are used

time_S3 <- system.time({
  samples <- runMCMC(Cmcmc, nchains = 1, niter = 10000, nburnin = 5000, thin = 10)
})
print(time_S3)

sample_summary<-as.matrix(Cmcmc$mvSamples)
samplesSummary(sample_summary)
# The summary function returns a matrix, so we can convert it directly
samps_summary_df <- as.data.frame(samps_summary)
library(MCMCvis)
MCMCtrace(object = samples,
          pdf = FALSE,
          #type='density',
          #ind = TRUE, # separate density lines per chain
          Rhat=F,
          n.eff=F,
          exact = TRUE,
          # priors = PR,
          #gvals = GV,
          #post_zm = T,
          params = "all")
# ## traceplots of parameters
library(tidyverse)
samples %>%
  as.data.frame() %>%
  mutate(iter = 1:nrow(samples)) %>%
  pivot_longer(-iter, names_to = 'param', values_to = 'value') %>%
  ggplot(aes(x = iter, y = value)) +
  geom_line() +
  facet_wrap(~param, scales = 'free')


# update MCMC chain to improve convergence: 
niter_ad <-10000
Cmcmc$run(niter_ad, reset=FALSE)
# extract matrix of previous MCMC samples augmented with new ones: 
more_samples <- as.matrix(Cmcmc$mvSamples)
samplesSummary(more_samples)


