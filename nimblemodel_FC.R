

library(nimble)
###
### functions to calculate probability of birth on given day of season
##  as skew normal distribution
###
dskew <- nimbleFunction(
    run = function(x = double(1), param = double(1)) {
        skew.bday <- param[1]
        mu.bday   <- param[2]
        sd.bday   <- param[3]
        delta <- skew.bday/sqrt(1+(skew.bday^2))
        omega <- sqrt(sd.bday^2/(1-(2*delta^2)/pi))
        xi <- mu.bday - omega*delta*sqrt(2/pi)
        Pb <- (2/omega)*dnorm((x-xi)/omega) * pnorm(skew.bday*((x-xi)/omega))
        returnType(double(1))
        return(Pb)
    })

intPb <- nimbleFunction(
    run = function(ndays = double(), skew.bday = double(), mu.bday = double(), sd.bday = double()) {
        Pb <- numeric(ndays)
        nday<-123
        pars <- c(skew.bday, mu.bday, sd.bday)
        for(t in 1:ndays) {
            Pb[t] <- nimIntegrate(dskew, t-1, t, pars)[1]
        }
        returnType(double(1))
        return(Pb)
    })


## =====================================================================
## ============== insert simulated / real data here ===============
## =====================================================================

## option 1
##  random data:
ndays <- 30
amax <- ndays
obsdays2 <- 11:27
wObs2 <- c(1, 1, 3, 4, 4, 5, 5, 3, 2, 1, 1, 1, 1, 0, 0, 0, 0)
mObs2 <- c(0, 0, 0, 0, 0, 1, 2, 4, 5, 5, 5, 4, 2, 2, 2, 2, 1)
obsdays <- 1:ndays
wObs <- mObs <- rep(0, ndays)
wObs[obsdays2] <- wObs2
mObs[obsdays2] <- mObs2

## option 2
# simulated data:
# sim_S4 <- simSealsS4()
# wObs <- sim_S4$wObs
# mObs <- sim_S4$mObs
# obsdays <- 1:length(wObs)

## =====================================================================
## =====================================================================

code <- nimbleCode({
  # Define priors for parameters
  # mu.bday ~ dunif(0, ndays)       ## Prior for mean of birth date 
  # skew.bday ~ dunif(0, ndays)     ## Prior for skewness of birth date
  # sd.bday ~ dunif(0, ndays)       ## Prior for standard deviation of birth date
  # scale ~ dunif(0, 200)           ## Prior for scale parameter
  
  #shifted scale beta for mu.bday,sd.bday,skew.bday parameters
  standard_beta.m ~ dbeta(alpha.m, beta.m)
  standard_beta.sd ~ dbeta(alpha.sd, beta.sd)
  standard_beta.s ~ dbeta(alpha.s, beta.s)
  # # Shifting and scaling the distribution
  mu.bday <- shift.m + (scale.m * standard_beta.m)
  sd.bday <- shift.sd + (scale.sd * standard_beta.sd)
  skew.bday <- shift.s + (scale.s * standard_beta.s)
  # Calculate probabilities for birth, moult, and leave events
  Pb[1:ndays]<- intPb(ndays, skew.bday, mu.bday, sd.bday)
  # Pb[1:ndays] <- dskew(x=1:ndays,param=c(skew.bday,mu.bday,sd.bday)) ##trouble shooting nimIntegrate issue
  # mu.sday ~ dunif(0, ndays)       ## Prior for mean of stage 4 tranisition date
  # sd.sday ~ dunif(0, ndays)       ## Prior for standard deviation of stage 4 transition date
  
  ###
  ### Transition probabilities follow CDF of normal distribution and are a function of age (amax)
  ###  mean stage transition day, and sd of transition day
  ###
  ## Probability of transitioning to stage 4 (moulting) 
  Ps[1:amax] <- pnorm(1:amax, mu.sday, sd.sday) - pnorm(0:(amax-1), mu.sday, sd.sday) 
  # mu.mday ~ dunif(0, ndays)       ## Prior for mean of moult date
  # sd.mday ~ dunif(0, ndays)       ## Prior for standard deviation of moult date
  
  ## Probability of completing moult (moulted) 
  Pm[1:amax] <- pnorm(1:amax, mu.mday, sd.mday) - pnorm(0:(amax-1), mu.mday, sd.mday) 
  # mu.lday ~ dunif(0, ndays)       ## Prior for mean of leave date
  # sd.lday ~ dunif(0, ndays)       ## Prior for standard deviation of leave date
  
  ## Probability of leaving 
  Pl[1:amax] <- pnorm(1:amax, mu.lday, sd.lday) - pnorm(0:(amax-1), mu.lday, sd.lday) 
  # Prior distributions for observation parameters
  # POWhite ~ dunif(0.7, 1)
  # POMoult ~ dunif(0.7, 1)
  # PCMoult ~ dunif(0.7, 1)
  # PCStge4 ~ dunif(0.7, 1)
  # # Set constraints
  # cons[1] ~ dconstraint(mu.sday < mu.mday) ## Ensure stage 4 comes before moult
  # cons[2] ~ dconstraint(mu.mday < mu.lday) ## and moult comes before leaving
  
  ##############################################################################
  # Define the dynamic process model
  # The number of pups born on day t follows a Poisson distribution.
  # The lambda parameter of the Poisson distribution is determined by the 
  # scale parameter multiplied by the birth probability for that day.
 ###############################################################################
  for(t in 1:ndays) {
    # in notation write up 'white' is same as 'b' or 'born'. All pups born into 'white' class
    white[t,1] ~ dpois(scale * Pb[t])
    
    # Initialize pups for each state  to 0
    stge4[t,1] <- 0                     
    moult[t,1] <- 0                        
    leave[t,1] <- 0                        
    
    # Calculate daily aggregates
    bornDaily[t] <- white[t,1]             ## Store the number of white born on the day
    whiteDaily[t] <- sum(white[t,1:t])     ## Calculate the total white count up to this day
    stge4Daily[t] <- sum(stge4[t,1:t])     ## and so on for all age classes
    moultDaily[t] <- sum(moult[t,1:t])     
    leaveDaily[t] <- sum(leave[t,1:t])    
  }
  
  
  #############################################################################
  # STORE NUMBER OF PUPS TRANSITIONING FROM ONE STATE TO ANOTHER -> STORE TOTAL 
  # NUMBER OF PUPS PRESENT IN COLONY IN EACH CLASS AT AGE A ON DAY T
  #############################################################################
  # Loop through each day (except the last one) and each possible age from 1 
  # to the minimum value of t and amax (except the last)
  for(t in 1:(ndays-1)) {
    for(a in 1:min(t,amax-1)) {
      
      # Transition probabilities from one state to another
      # For each combination of day and age, simulate transitions between 
      # different states (white->stage 4->moult->leave) based on transition probabilities above.
      white_s[t,a] ~ dbinom(Ps[a], white[t,a])   # the number of white pups (of age a) which transition to stage4 on day t
      stge4_m[t,a] ~ dbinom(Pm[a], stge4[t,a])   # the number of stage4 pups (of age a) which transition to moulted on day t
      moult_l[t,a] ~ dbinom(Pl[a], moult[t,a])   # the number of moulted pups (of age a) which leave the colony on day t
      
      # Update counts based on transitions
      white[t+1,a+1] <- white[t,a] - white_s[t,a]   # total number of white pups present in the colony, of age a on day t 
      stge4[t+1,a+1] <- stge4[t,a] - stge4_m[t,a] + white_s[t,a]  # total number of stage4 pups present in the colony, of age a on day t 
      moult[t+1,a+1] <- moult[t,a] - moult_l[t,a] + stge4_m[t,a]  # total number of moulted pups present in the colony, of age a on day t 
      leave[t+1,a+1] <- moult_l[t,a]                # total number of pups that leave the colony, of age a on day t
    }
  }
  
  # Calculate total births and leaves
   bornTotal  <- sum(bornDaily[1:ndays])    ## Total births calculated by summing daily births
  leaveTotal <- sum(leaveDaily[1:ndays])   ## Total leaves calculated by summing daily leaves
  
  
  
  # Observation process
  for(i in 1:length(obsdays)) {
  # Calculate mean and binomial variance for white and moulted pups
    whiteM[i] <- whiteDaily[obsdays[i]] * POWhite + stge4Daily[obsdays[i]] * POWhite * PCStge4 + moultDaily[obsdays[i]] * POMoult * (1 - PCMoult)
    whiteV[i] <- whiteDaily[obsdays[i]] * POWhite * (1 - POWhite) + stge4Daily[obsdays[i]] * POWhite * PCStge4 * (1 - POWhite * PCStge4) + moultDaily[obsdays[i]] * POMoult * (1 - PCMoult) * (1 - POMoult * (1 - PCMoult))
    whiteS[i] <- max(sqrt(whiteV[i]), 0.00000001)
    moultM[i] <- stge4Daily[obsdays[i]] * POWhite * (1 - PCStge4) + moultDaily[obsdays[i]] * POMoult * PCMoult
    moultV[i] <- stge4Daily[obsdays[i]] * POWhite * (1 - PCStge4) * (1 - POWhite * (1 - PCStge4)) + moultDaily[obsdays[i]] * POMoult * PCMoult * (1 - POMoult * PCMoult)
    moultS[i] <- max(sqrt(moultV[i]), 0.00000001)
   
    wObs[i] ~ dnorm(whiteM[i], sd = whiteS[i])
    mObs[i] ~ dnorm(moultM[i], sd = moultS[i])
  }
})

# Define constants
constants <- list(  
  ndays = ndays,
  amax = amax,
  obsdays = obsdays,
  mu.sday = ndays/4,
  sd.sday = 20,
  mu.mday = ndays/3,
  sd.mday = 20,
  mu.lday = ndays/2,
  sd.lday = 20,
  POWhite=1,POMoult=0.94,PCMoult=0.82,PCStge4=0.75
 )


# Define data
data <- list(
  wObs = wObs,
  mObs = mObs,
  cons = c(1,1)
)


useInits <- 4

# Initialize variables if needed
if(useInits == 1) {
  zeroArray <- array(0, c(ndays,amax))
  # assign zero-filled array to variables 'white_init', 'white_s_init', 
  # 'stge4_m_init', and 'moult_l_init', effectively initializing them with zeros.
  white_init <- white_s_init <- stge4_m_init <- moult_l_init <- zeroArray
  white_init
}

if(useInits == 2) {
    zeroArray <- array(0, c(ndays,amax))
    white_init <- white_s_init <- stge4_m_init <- moult_l_init <- zeroArray
    # 'white_init[1:10,1]' sets the initial value to 2 for the first 10 days 
    # in the 'white_init' array, column 1.
    white_init[1:10,1] <- 2
    white_s_init[10,1:10] <- 2
    stge4_m_init[20,11:20] <- 2
    moult_l_init[30,21:30] <- 2
}

if(useInits == 3) {
    zeroArray <- array(0, c(ndays,amax))
    white_init <- white_s_init <- stge4_m_init <- moult_l_init <- zeroArray
    white_init[1:10,1] <- 3
    white_s_init[10,1:10] <- 3
    stge4_m_init[20,11:20] <- 3
    moult_l_init[30,21:30] <- 3
}

if(useInits == 4) {
  # Calculate minimum required initial populations based on observations
    wMin <- ceiling(max(wObs/obsdays))
    mMin <- ceiling(max(mObs/obsdays))
    bMin <- ndays * max(wMin, mMin)
    zeroArray <- array(0, c(ndays,amax))
    white_init <- white_s_init <- stge4_m_init <- moult_l_init <- zeroArray
    white_init[1:ndays,1] <- bMin
    for(t in 1:(ndays-1))   white_s_init[t,1:min(t,amax-1)] <- wMin
    for(t in 2:(ndays-1))   stge4_m_init[t,2:min(t,amax-1)] <- mMin
    for(t in 3:(ndays-1))   moult_l_init[t,3:min(t,amax-1)] <- 1
}

inits <- list(
  mu.bday = ndays/5,
  skew.bday = 1,
  sd.bday = 20,
  scale = 10,
  alpha.m=2,
  beta.m=2,
  shift.m=10, # Shift the distribution to start at 33
  # The "scale" parameter in a shifted beta distribution
  # multiplies all values in the original beta distribution by 
  # a constant factor. 
  scale.m=40, 
  # sd.bday
  alpha.sd=2,  # Shape parameter (alpha > 0)
  beta.sd=2,    # Shape parameter (beta > 0)
  shift.sd=6,   # Shift the distribution to start at 6
  scale.sd=5,   # 'stretch' to 11
  alpha.s = 2,
  beta.s = 2,
  shift.s = 1,
  scale.s = 5,
    # mu.sday = ndays/4,
    # sd.sday = 20,
    # mu.mday = ndays/3,
    # sd.mday = 20,
    # mu.lday = ndays/2,
    # sd.lday = 20,
    # POWhite = 0.9,
    # POMoult = 0.9,
    # PCMoult = 0.9,
    # PCStge4 = 0.9,
    white   = white_init,
    white_s = white_s_init,
    stge4_m = stge4_m_init,
    moult_l = moult_l_init,
    moult = zeroArray,
    stge4 = zeroArray,
    leave = zeroArray
)

Rmodel <- nimbleModel(code, constants, data, inits)

## check this is real, negative number:
Rmodel$calculate()

conf <- configureMCMC(Rmodel, useConjugacy = FALSE)
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples <- runMCMC(Cmcmc, 500)

samplesSummary(samples)




  