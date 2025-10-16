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
#Functions 

source("HgNimble_input_S3.R") # Nimble functions to count whitecoat and moulted pups present on any given day

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

#####################################

code <- nimbleCode({
  #priors
  Nborn ~ dunif(100,10000) 
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
    sdwhite[obs] <- sqrt(whitePups[obsdays[obs]] * POWhite * (1 - POWhite) + moultPups[obsdays[obs]] * (1 - PCMoult) * POMoult * (1 - (1 - PCMoult) * POMoult))

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

######################################################################

library(nimble)
library(dplyr)
library(MCMCvis)
library(ggplot2)

# directory to save results tto remote desktop
#output_dir <- "D:/mg349/results"
# directory local comp
output_dir <- "C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/07 - GitHub Repositories/Hg_pup_model/results/HgNimble_Model_S3_simulation"
# Create a function to run the model and save results for a given dataset
run_model_simulation <- function(simulation_number) {
  # Generate a new dataset
  # Generate observations using the generate_obs function
  sim_S3 <- simSealsS3(
    N = 2000, 
    ndays = 123, amax = 123, mu.bday = 45, sd.bday = 10,
    birthcurve = "skewnormal",
    skew.bday = 5, 
    mu.mday = 23, sd.mday = 5,
    mu.lday = 31.5, sd.lday = 7,
    POWhite = 0.98,
    POMoult = 0.95, 
    PCWhite = 1, 
    PCMoult = 0.82,  
    obsdays =  c(33,49,61,75,96),
    #popproc: biological process can be generated as: 
    #         "detcont"=deterministic and continuous
    #         "stocdis"= stochastic and discrete
    popproc = "detcont", 
    # obsproc: observation process can be generated as: 
    #          "exact"= obs generated with no variation, 
    #          "normappx"= obs are normal approximated using binomial variance
    #          "binomial" = obs are binomially distributed 
    obsproc = "normappx" 
  )
  
  wObs <- sim_S3$wObs
  mObs <- sim_S3$mObs
  obsdays <- c(33,49,61,75,96)
  
  # Print the generated data
  cat("Sim:", simulation_number, "\n")
  cat("Obs days:", obsdays, "\n")
  cat("wObs:", wObs, "\n")
  cat("mobs:", mObs, "\n")
  cat("========================================\n")
  
  # Define the data, constants, and initial values
  data_list <- list(wObs = wObs, mObs = mObs, whiteObs = wObs, moultObs = mObs)
  
  constants_list <- list(
    ndays = 123,
    amax = 123,
    obsdays = obsdays,
    mu.mday = 23,
    sd.mday = 5,
    mu.lday = 31.5,
    sd.lday = 7,
    POWhite = 0.95,
    POMoult = 0.95,
    PCMoult = 0.82
  )
  
  inits_list <- list(
    Nborn = runif(1, 100, 10000),
    mu.bday = runif(1, 20, 60),
    skew.bday = runif(1, 2, 6),
    sd.bday = runif(1, 5, 10)
  )
  
  # Define the model
  code <- nimbleCode({
    # Priors
    Nborn ~ dunif(100, 10000)
    mu.bday ~ dunif(20, 60)
    skew.bday ~ dunif(2, 6)
    sd.bday ~ dunif(5, 15)
    
    # Model specifics
    for (t in 1:ndays) { 
      born[t] <- Pb_sum[t] * Nborn
    }
    
    
    Pb[1:ndays] <- int_pb(ndays, skew.bday, mu.bday, sd.bday)
    Pb_sum[1:ndays] <- Pb[1:ndays] / sum(Pb[1:ndays])
    
    Pm[1:amax] <- pnorm(1:amax, mu.mday, sd.mday) - pnorm(0:(amax - 1), mu.mday, sd.mday)
    Pm_sum[1:amax] <- Pm[1:amax] / sum(Pm[1:amax])
    
    Pl[1:amax] <- pnorm(1:amax, mu.lday, sd.lday) - pnorm(0:(amax - 1), mu.lday, sd.lday)
    Pl_sum[1:amax] <- Pl[1:amax] / sum(Pl[1:amax])
    
    whitePups[1:ndays] <- whiteSum(born[1:ndays], Pm_sum[1:amax], whiteObs[])
    moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax], moultObs[])
    
    
    
    for (obs in 1:length(obsdays)) {
      whiteCount[obs] <- whitePups[obsdays[obs]] * POWhite + moultPups[obsdays[obs]] * POMoult * (1 - PCMoult)
      sdwhite[obs] <- sqrt(whitePups[obsdays[obs]] * POWhite * (1 - POWhite) + moultPups[obsdays[obs]] * (1 - PCMoult) * POMoult * (1 - POMoult * (1 - PCMoult)))
      moultCount[obs] <- moultPups[obsdays[obs]] * POMoult * PCMoult
      sdmoult[obs] <- sqrt(moultPups[obsdays[obs]] * POMoult * PCMoult * (1 - POMoult * PCMoult))
      
      wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs])
      mObs[obs] ~ dnorm(moultCount[obs], sd = sdmoult[obs])
    }
  })
  
  ## Parameters monitored
  params = c('Nborn','mu.bday', 'sd.bday',"skew.bday")
  # Create and compile the model and MCMC configuration
  model <- nimbleModel(code = code, inits = inits_list, data = data_list, constants = constants_list)
  Cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = params, useConjugacy = FALSE)
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # Run MCMC
  set.seed(0)
  samples <- runMCMC(Cmcmc, nchains = 1, niter = 10000, nburnin = 5000, thin = 10)
  
  # Save results
  saveRDS(samples, file = file.path(output_dir, paste0("HgNimble_Model_S3_","sim_", simulation_number,"_", Sys.Date(), ".rds")))
}

# Run simulations
n_simulations <- 100
for (i in 1:n_simulations) {
  run_model_simulation(i)
}


current_date <- Sys.Date()

# Initialize a list to store all loaded data frames
all_mcmc_results <- list()

# Function to process MCMC results
process_mcmc_results <- function(filename) {
  # Read the MCMC results from the RDS file
  mcmc_data <- readRDS(filename)
  
  return(mcmc_data)
}

#library(coda)
#mcmc<-as.mcmc(mcmc_results)
# 
# MCMCtrace(object=mcmc, 
#           pdf = FALSE,
#           #type='density',
#           #ind = TRUE, # separate density lines per chain
#           Rhat=F,
#           n.eff=F,
#           exact = TRUE,
#           # priors = PR,
#           #gvals = GV,
#           #post_zm = T,
#           params = "all")

# Plotting using ggplot2
#Define the true values for each parameter
true_values <- list(
  "mu.bday" = 45, # Example true value for mu.bday
  "sd.bday" = 10, # Example true value for sd.bday
  "nborn" = 500, # Example true value for nborn
  "mu.days.on.colony" = 31.5 # Example true value for mu.days.on.colony
)
###############################################
#this plot plots all simulations for each scenario for each parameter
# Plotting using ggplot2

library(HDInterval)
# Function to calculate HPD interval
calculate_hpd <- function(plot_data, param, prob = 0.95) {
  hpd_interval <-hdi(plot_data[[param]], prob = prob)
  return(hpd_interval)
}

for (param in params_of_interest) {
  mean_value <- mean(plot_data[[param]])
  hpd_interval <- calculate_hpd(plot_data, param)
  hpd_intervals[param, c("mean_value", "lower_hpd", "upper_hpd")] <- c(mean_value, hpd_interval[1], hpd_interval[2])
}


# lots <- list()
for (param in params_of_interest) {
  p <- ggplot(plot_data, aes_string(x = param, color = "simulation")) +
    geom_density() +
    facet_wrap(~ scenario, scales = "free") +
    labs(title = paste("Posterior Density of", param),
         x = param,
         y = "Density") +
    scale_color_discrete(guide = "none") +
    theme_minimal() +
    theme(legend.position = "topright") +
    guides(color = guide_legend(title.position = "top", label.position = "top"))
  
  # Add vertical line for the true value
  true_value <- true_values[[param]]
  p <- p + geom_vline(xintercept = true_value, linetype = "dashed", color = "red", size = 1)
  
  # Add confidence intervals
  # p <- p + geom_ribbon(aes_string(ymin = "lower_bound", ymax = "upper_bound"), alpha = 0.2, fill = "blue")
  
  plots[[param]] <- p
}

# Display the plots
for (param in params_of_interest) {
  print(plots[[param]])
}

