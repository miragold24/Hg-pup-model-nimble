

library(nimble)

## choose dataset
## 1: my random data
## 2: your simulated (sim_S4) data
useData <- 2

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

dIndividual <- nimbleFunction(
    run = function(x = double(1), pb = double(1), ps = double(1), pm = double(1), pl = double(1), log = integer(0, default = 0)) {
        ## x: (born day, stage4 day, moult day, leave day)
        ndays <- length(pb)
        lp <- 0
        if(0 < x[1] & x[1] < x[2] & x[2] < x[3] & x[3] < x[4] & x[4] <= ndays) {
            age_s <- x[2] - x[1] + 1
            age_m <- x[3] - x[1] + 1
            age_l <- x[4] - x[1] + 1
            lp <- log(pb[x[1]]) + log(ps[age_s]) + log(pm[age_m]) + log(pl[age_l])
        }
        returnType(double())
        if(log) return(lp) else return(exp(lp))
    })

rIndividual <- nimbleFunction(
    run = function(n = integer(), pb = double(1), ps = double(1), pm = double(1), pl = double(1)) {
        stop('not implemented')
        returnType(double(1))
        return(c(1,2,3,4))
    })

registerDistributions(
    list(dIndividual = list(
             BUGSdist = 'dIndividual(pb, ps, pm, pl)',
             types = c('value = double(1)', 'pb = double(1)', 'ps = double(1)', 'pm = double(1)', 'pl = double(1)'),
             discrete = TRUE,
             mixedSizes = TRUE,
             pqAvail = FALSE
         ))
)

calcDaily_boundary <- nimbleFunction(
    run = function(ndays = double(), z = double(1)) {
        ans <- numeric(ndays)
        for(t in 1:ndays)   ans[t] <- sum(z == t)
        returnType(double(1))
        return(ans)
    }
)

calcDaily_intermediate <- nimbleFunction(
    run = function(ndays = double(), z = double(2)) {
        ans <- numeric(ndays)
        for(t in 1:ndays)   ans[t] <- sum(  (z[,1] <= t)  &  (z[,2] > t)  )
        returnType(double(1))
        return(ans)
    }
)

sampler_add_remove <- nimbleFunction(
    name = 'sampler_add_remove',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        tries <- if(!is.null(control$tries)) control$tries else 1
        calcNodes <- model$getDependencies(target)
        if(target != 'z')  stop('target should only be z node')
        ndays <- length(model$pb)
        M <- dim(model$z)[1]
        totalTries <- 0
        totalAccept <- 0
        totalAdd <- 0
        totalRemove <- 0
    },
    run = function() {
        for(i in 1:tries) {
            lp0 <- model$getLogProb(calcNodes)
            nInd <- model$nInd[1]
            a <- runif(1)
            if(a < 0.5) {  ## propose adding an individual
                if(nInd == M) { print('Error: trying to add an individual, already M individuals in model, MCMC invalid'); return() }
                z_new <- generate_z()
                model$z[nInd+1,] <<- z_new
                lpFwd <- dIndividual(z_new, model$pb, model$ps, model$pm, model$pl, 1)
                lpRev <- -log(nInd+1)
            } else {              ## propose removing an individual
                ind <- ceiling(nInd*runif(1))
                z_old <- model$z[ind,]
                model$z[ind,] <<- model$z[nInd,]
                model$z[nInd,] <<- c(0, 0, 0, 0)
                lpFwd <- -log(nInd)
                lpRev <- dIndividual(z_old, model$pb, model$ps, model$pm, model$pl, 1)
            }
            lp1 <- model$calculate(calcNodes)
            logMHR <- lp1 + lpRev - lp0 - lpFwd
            jump <- decide(logMHR)
            if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
            ## diagnostics:
            totalTries <<- totalTries + 1
            if(jump) {
                totalAccept <<- totalAccept + 1
                if(a < 0.5)   totalAdd <<- totalAdd + 1   else   totalRemove <<- totalRemove + 1
            }
        }
    },
    methods = list(
        generate_z = function() {
            z_new <- numeric(4)
            z_new[1] <- ceiling(ndays*runif(1))
            z_new[2] <- ceiling(ndays*runif(1))
            while(z_new[2] == z_new[1])   z_new[2] <- ceiling(ndays*runif(1))
            z_new[3] <- ceiling(ndays*runif(1))
            while(z_new[3] == z_new[1] | z_new[3] == z_new[2])   z_new[3] <- ceiling(ndays*runif(1))
            z_new[4] <- ceiling(ndays*runif(1))
            while(z_new[4] == z_new[1] | z_new[4] == z_new[2] | z_new[4] == z_new[3])   z_new[4] <- ceiling(ndays*runif(1))
            ## sort:
            for(i in 1:3)   if(z_new[i] > z_new[i+1])   z_new[i:(i+1)] <- c(z_new[i+1], z_new[i])
            for(i in 1:3)   if(z_new[i] > z_new[i+1])   z_new[i:(i+1)] <- c(z_new[i+1], z_new[i])
            for(i in 1:3)   if(z_new[i] > z_new[i+1])   z_new[i:(i+1)] <- c(z_new[i+1], z_new[i])
            returnType(double(1))
            return(z_new)
        },
        reset = function() {
            totalTries <<- 0
            totalAccept <<- 0
            totalAdd <<- 0
            totalRemove <<- 0
        }
    )
)

sampler_modify <- nimbleFunction(
    name = 'sampler_modify',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        tries <- if(!is.null(control$tries)) control$tries else 1
        calcNodes <- model$getDependencies(target)
        if(target != 'z')  stop('target should only be z node')
        ndays <- length(model$pb)
        z_new <- numeric(4)
        totalTries <- 0
        totalAccept <- 0
    },
    run = function() {
        for(i in 1:tries) {
            nInd <- model$nInd
            ind <- ceiling(nInd*runif(1))
            index <- ceiling(4*runif(1))
            delta <- 2*floor(2*runif(1))-1
            for(j in 1:4)   z_new[j] <<- model$z[ind,j][1]
            z_new[index] <<- z_new[index] + delta
            valid <- 0 < z_new[1] & z_new[1] < z_new[2] & z_new[2] < z_new[3] & z_new[3] < z_new[4] & z_new[4] <= ndays
            if(valid) {
                lp0 <- model$getLogProb(calcNodes)
                model$z[ind,] <<- z_new
                lp1 <- model$calculate(calcNodes)
                logMHR <- lp1 - lp0
                jump <- decide(logMHR)
                if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
                else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
            } else jump <- FALSE
            ## diagnostics:
            totalTries <<- totalTries + 1
            if(jump)   totalAccept <<- totalAccept + 1
        }
    },
    methods = list(
        reset = function() {
            totalTries <<- 0
            totalAccept <<- 0
        }
    )
)

########################################################################
################# insert your simulated / real data here ###############
########################################################################

## option 1
## my random data:
if(useData == 1) {
    ndays <- 30
    wObs <- c(1, 1, 3, 4, 4, 5, 5, 3, 2, 1, 1, 1, 1, 0, 0, 0, 0)
    mObs <- c(0, 0, 0, 0, 0, 1, 2, 4, 5, 5, 5, 4, 2, 2, 5, 2, 1)
    obsdays <- 5:21
}

## option 2
## your simulated data:
if(useData == 2) {
   # sim_S4 <- simSealsS4()
    wObs <- round(sim_S4$wObs)
    mObs <- round(sim_S4$mObs)
    ndays <- length(wObs)
    obsdays <- 1:length(wObs)
}

########################################################################
########################################################################

code <- nimbleCode({
  
  ###
  ### Priors 
  ###
  
  #shifted beta priors
  # standard_beta.m ~ dbeta(alpha.m, beta.m)
  # standard_beta.sd ~ dbeta(alpha.sd, beta.sd)
  # standard_beta.s ~ dbeta(alpha.s, beta.s)
  # # Shifting and scaling the distribution
  # mu.bday <- shift.m + (scale.m * standard_beta.m)
  # sd.bday <- shift.sd + (scale.sd * standard_beta.sd)
  # skew.bday <- shift.s + (scale.s * standard_beta.s)

    
    #uninformative priors 
    mu.bday ~ dunif(0, ndays)
    skew.bday ~ dunif(0, ndays)
    sd.bday ~ dunif(0, ndays)
    mu.sday ~ dunif(0, ndays)
    sd.sday ~ dunif(0, ndays)
    mu.mday ~ dunif(0, ndays)
    sd.mday ~ dunif(0, ndays)
    mu.lday ~ dunif(0, ndays)
    sd.lday ~ dunif(0, ndays)
    # # 
    cons[1] ~ dconstraint(mu.sday < mu.mday)
    cons[2] ~ dconstraint(mu.mday < mu.lday)
    #observation priors 
    # POWhite ~ dunif(.5, 1)
    # POMoult ~ dunif(.5, 1)
    # PCMoult ~ dunif(.5, 1)
    # PCStge4 ~ dunif(.5, 1)
    
    #calculate birth and transition probabilities 
    pb[1:ndays] <- int_pb(ndays, skew.bday, mu.bday, sd.bday)

    ps1[1:ndays] <- pnorm(1:ndays, mu.sday, sd.sday) - pnorm(0:(ndays-1), mu.sday, sd.sday)
    ps[1:ndays] <- ps1[1:ndays]/sum(ps1[1:ndays])
  
    pm1[1:ndays] <- pnorm(1:ndays, mu.mday, sd.mday) - pnorm(0:(ndays-1), mu.mday, sd.mday)
    pm[1:ndays] <- pm1[1:ndays]/sum(pm1[1:ndays])

   
    pl1[1:ndays] <- pnorm(1:ndays, mu.lday, sd.lday) - pnorm(0:(ndays-1), mu.lday, sd.lday)
    pl[1:ndays] <- pl1[1:ndays]/sum(pl1[1:ndays])
   

    ##
    for(i in 1:M)   z[i, 1:4] ~ dIndividual(pb[1:ndays], ps[1:ndays], pm[1:ndays], pl[1:ndays])
    ##
    bornDaily [1:ndays] <- calcDaily_boundary    (ndays, z[1:M,1  ])
    whiteDaily[1:ndays] <- calcDaily_intermediate(ndays, z[1:M,1:2])
    stge4Daily[1:ndays] <- calcDaily_intermediate(ndays, z[1:M,2:3])
    moultDaily[1:ndays] <- calcDaily_intermediate(ndays, z[1:M,3:4])
    leaveDaily[1:ndays] <- calcDaily_boundary    (ndays, z[1:M,  4])
    nInd <- sum(bornDaily[1:ndays])
    #whiteDaily1_4[1:ndays]<-whiteDaily[1:ndays]+stge4Daily[1:ndays]

    ###
    ### OBSERVATION MODEL
    ### 
          
          # summing binomials (?):
          # for (i in 1:length(obsdays)){
       
          #   wObs[i] ~ dbinom(POWhite, whiteDaily[i])+ dbinom(PCStge4*POWhite, stge4Daily[i]) + dbinom((1-PCMoult)*POMoult,moultDaily[i]) 
          #   mObs[i] ~ dbinom(POMoult*PCMoult, moultDaily[i]) +dbinom((1-PCStge4)*POWhite,whiteDaily[i])
          # 
          # }
    
    
        # Poisson
         for (i in 1:length(obsdays)){
              whiteM[i] <-  whiteDaily[obsdays[i]] * POWhite               + stge4Daily[obsdays[i]] * POWhite*PCStge4                       + moultDaily[obsdays[i]] * POMoult*(1-PCMoult)
              moultM[i] <- stge4Daily[obsdays[i]] * POWhite*(1-PCStge4)                         + moultDaily[obsdays[i]] * POMoult*PCMoult
              wObs[i] ~ dpois(whiteM[i])
              mObs[i] ~ dpois(moultM[i])
         }


    ## norm apprx
        # for(i in 1:length(obsdays)) {
        #     whiteM[i] <- whiteDaily[obsdays[i]] * POWhite               + stge4Daily[obsdays[i]] * POWhite*PCStge4                       + moultDaily[obsdays[i]] * POMoult*(1-PCMoult)
        #     whiteV[i] <- whiteDaily[obsdays[i]] * POWhite * (1-POWhite) + stge4Daily[obsdays[i]] * POWhite*PCStge4 * (1-POWhite*PCStge4) + moultDaily[obsdays[i]] * POMoult*(1-PCMoult) * (1-POMoult*(1-PCMoult))
        #     whiteS[i] <- max(sqrt(whiteV[i]), 0.00000001)
        #     moultM[i] <- stge4Daily[obsdays[i]] * POWhite*(1-PCStge4)                         + moultDaily[obsdays[i]] * POMoult*PCMoult
        #     moultV[i] <- stge4Daily[obsdays[i]] * POWhite*(1-PCStge4)*(1-POWhite*(1-PCStge4)) + moultDaily[obsdays[i]] * POMoult*PCMoult * (1-POMoult*PCMoult)
        #     moultS[i] <- max(sqrt(moultV[i]), 0.00000001)
        # 
        #   # With small variance  model expects the observations to be very close to the means whiteM[i] and moultM[i].
        #   # If the observations are not close to the means, the log prob will be a large negative number
        #     wObs[i] ~ dnorm(whiteM[i], sd = 0.001)
        #     mObs[i] ~ dnorm(moultM[i], sd = 0.001)
        #     # wObs[i] ~ dnorm(whiteM[i], sd = whiteS[i])
        #     # mObs[i] ~ dnorm(moultM[i], sd = moultS[i])
        # }

    
    
    ## BINOMIAL WITH CLASSIFICATION FIRST -> POBS IN BINOMIAL

            # for(i in 1:length(obsdays)) {
            #   whiteM[i] <- whiteDaily[obsdays[i]]+ stge4Daily[obsdays[i]] * PCStge4                       + moultDaily[obsdays[i]] * (1-PCMoult)
            #   # whiteV[i] <- whiteDaily[obsdays[i]] + stge4Daily[obsdays[i]] * POWhite*PCStge4 * (1-PCStge4) + moultDaily[obsdays[i]] *(1-PCMoult) * (1-(1-PCMoult))
            #   # whiteS[i] <- max(sqrt(whiteV[i]), 0.00000001)
            #   moultM[i] <- stge4Daily[obsdays[i]] * (1-PCStge4)                         + moultDaily[obsdays[i]] * PCMoult
            #   # moultV[i] <- stge4Daily[obsdays[i]] * (1-PCStge4)*(1-(1-PCStge4)) + moultDaily[obsdays[i]] * PCMoult * (1-PCMoult)
            #   # moultS[i] <- max(sqrt(moultV[i]), 0.00000001)
            # 
            #   wObs[i] ~ dbinom(POWhite,whiteM[i])
            #   mObs[i] ~ dbinom(POMoult,moultM[i])
            # 
            # }
         
    
    
    
 })

M <- max(ceiling(max(wObs+mObs) * 1.2), 20)       ## maximum possible number of individuals this model will allow
constants <- list(ndays = ndays, M = M, obsdays = obsdays,
                   # ,sd.sday=5,sd.mday=5,sd.lday = 7,mu.sday=17,mu.mday=23,mu.lday = 31.5,
                  POWhite=0.98,POMoult=0.94,PCMoult=0.82,PCStge4=0.75)

data <- list(wObs = wObs, mObs = mObs, cons = c(1,1))
z_init <- array(0, c(M,4))

#for (1 in 1:x) x must be less than the total generated nInd
if(useData == 1) {    for(i in 1:15 )   z_init[i,1:4] <- c( 5, 10, 13,  22)    }
if(useData == 2) {    for(i in 1:400)   z_init[i,1:4] <- c(25, 40, 55, 120)    }
#if(useData == 2) {    for(i in 1:400)   z_init[i,1:4] <- c(45, 62, 68, 76)    }


inits <- list(
  # shifted beta prior inits 
   # standard_beta.m=rbeta(1,2,2),
   # standard_beta.s=rbeta(1,2,2),
   # standard_beta.sd=rbeta(1,2,2),
   # alpha.m=2,
   # beta.m=2,
   # shift.m=10,
   # scale.m=40,
   # # sd.bday
   # alpha.sd=2,
   # beta.sd=2,
   # shift.sd=6,
   # scale.sd=5,
   # alpha.s = 2,
   # beta.s = 2,
   # shift.s = 1,
   # scale.s = 5,
  
    
    # Dt inits for data 1 with uninformative priors
    # mu.bday =  45, #(1/5) * ndays,
    # skew.bday = 5, #2,
    # sd.bday = 10, #20,
    mu.sday = 17,#(2/5) * ndays,
    sd.sday = 10,#20,
    mu.mday = 23,#(3/5) * ndays,
    sd.mday = 7, #20,
    mu.lday = 31.5,#(4/5) * ndays,
    sd.lday = 10, #20,
#     POWhite = 0.9,
#     POMoult = 0.9,
#     PCMoult = 0.9,
#     PCStge4 = 0.9,
     z = z_init
 )

Rmodel <- nimbleModel(code, constants, data, inits,dimensions=list(whiteM=ndays,moultM=ndays))
# calculate() provides the log likelihood of the model 
# the closer to zero the log-likelihood, the better the model fits the data.
Rmodel$calculate()   ## make sure this gives you a real number (not NaN or -Inf)

conf <- configureMCMC(Rmodel, nodes = Rmodel$getNodeNames(topOnly = TRUE))
conf$addMonitors('nInd')
conf$addSampler('z', 'add_remove', tries = 10)
conf$addSampler('z', 'modify',     tries = 10)

Rmcmc <- buildMCMC(conf)
runMCMC(Rmcmc,niter=1)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
samples1 <- runMCMC(Cmcmc, niter = 10000, nburnin = 5000, thin = 2)

#check initial values

#samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 5000, thin = 10, nchains = 2)
# Run the code and capture time information
#run_time <- system.time(samples <- runMCMC(Cmcmc, niter = 10000, nburnin = 5000, thin = 10, nchains = 2))
# Access the elapsed time from the list
#cat("Elapsed time for runMCMC:", run_time[["elapsed"]], "seconds\n")
############################################################################# 
####################### plots and diagnostics ############################### 
############################################################################# 

## final value of individual histories:
Cmodel$z

wObs2 <- mObs2 <- rep(NA, ndays)
wObs2[obsdays] <- wObs
mObs2[obsdays] <- mObs

## all observation data:
data.frame(
    day  = 1:ndays,
    wObs = wObs2,
    mObs = mObs2
)

## observation data and (final) daily counts from MCMC:
data.frame(
    born  = Cmodel$bornDaily,
    wObs = wObs2,
    white = Cmodel$whiteDaily,
    stge4 = Cmodel$stge4Daily,
    mObs = mObs2,
    moult = Cmodel$moultDaily,
    leave = Cmodel$leaveDaily
)

library(coda)
library(MCMCvis)
#par(mfrow=c(2,4))
# 
# MCMCtrace(object = samples1,
#                    pdf =F, 
#                     type='density',
#                    ind = TRUE, # separate density lines per chain
#                    Rhat=F,
#                    n.eff=F,
#                     params = 'all') #c('nInd','mu.bday','sd.bday','skew.bday'))

MCMCsummary(samples1)
## posterior statistics of parameters:
# samplesSummary(samples)
# 
# ## traceplots of parameters
library(tidyverse)
samples1 %>%
    as.data.frame() %>%
    mutate(iter = 1:nrow(samples1)) %>%
    pivot_longer(-iter, names_to = 'param', values_to = 'value') %>%
    ggplot(aes(x = iter, y = value)) +
    geom_line() +
    facet_wrap(~param, scales = 'free')



