
library(nimble)

code <- nimbleCode({
    n ~ dpois(10)
    p ~ dunif(0, 1)

    y_mean <- n*p
    #binomial variance calculation
    y_sd <- sqrt(n*p*(1-p))
    
    for(i in 1:5) {
        y[i] ~ dnorm(y_mean, sd = y_sd)
    }
    # for (i in 1:5) {
    #   y[i] ~ dpois(y_mean)  # Sample counts from Poisson distribution
    # }
    # for (i in 1:5) {
    #   wi[i] <- dbinom(round(whitePups[idays[i]] + S4Pups[idays[i]]), prob = POWhite) * PCWhite + 
    #     dbinom(round(moultPups[idays[i]]),  prob = POMoult) * (1 - PCMoult)
    #   mObs[i] <- dbinom(round(whitePups[idays[i]] + S4Pups[idays[i]]), prob = POWhite) * (1 - PCWhite) + 
    #     dbinom(round(moultPups[idays[i]]),  prob = POMoult) * PCMoult
    # }
    
    
    
})

constants <- list()
data <- list(y = c(10, 10, 10, 10, 15))
inits <- list(n = 15, p = 0.5)

Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
samples <- runMCMC(Cmcmc, 10000)

samplesSummary(samples)
table(samples[, 'n'])
mean(samples[, 'n'] < 15)    ## nearly half of MCMC chain represents n < 15

