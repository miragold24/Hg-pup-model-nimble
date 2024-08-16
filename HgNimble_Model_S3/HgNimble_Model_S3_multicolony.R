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
# Extract data for the specified colony and year

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


# Define the NIMBLE model code
code <- nimbleCode({
  #priors
  Nborn ~ dunif(100, 10000)
  mu.bday ~ dunif(20,60)
  skew.bday ~ dunif(2,6)
  sd.bday ~ dunif(5,15)
  
  # deterministic and continuous births over breeding season (ndays)
  for (t in 1:ndays) { 
    born[t] <- Pb_sum[t] * Nborn
  }
  
  # calculate PDF of skew normal then integrate over intervals for Pb on given day
  Pb[1:ndays] <- int_pb(ndays, skew.bday, mu.bday, sd.bday)
  Pb_sum[1:ndays] <- Pb[1:ndays] / sum(Pb[1:ndays])
  
  # Probability of moulting as function of age (Probability of finishing moult)
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

# Output directory
output_dir <- "D:/mg349/results"

# Parameters monitored
params <- c('Nborn', 'mu.bday', 'sd.bday', 'skew.bday')

# Measure overall time for running the models
overall_time <- system.time({
# Loop through each colony and year
for (i in 1:nrow(colyears)) {
  selected_colony <- colyears$Colony[i]
  selected_year <- colyears$Year[i]
  
  print(paste("Running model for colony:", selected_colony, "year:", selected_year))
  
  realdata <- subset(data, Colony == selected_colony & Year == selected_year)
  wObs <- realdata$Whitecoats + realdata$DeadPups
  mObs <- realdata$Moulteds
  obsdays <- unique(realdata$Day)
  ndays <- max(obsdays)
  
  # Define data, constants, and initial values
  #real data 
  data_list=list(wObs=wObs,mObs=mObs,whiteObs=wObs,moultObs=mObs)
  constants_list=list(ndays=ndays,amax=amax,obsdays=obsdays,
                      mu.mday=mu.mday,sd.mday=sd.mday,
                      mu.lday=mu.lday,sd.lday=sd.lday,
                      POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult)
  
  
  inits_list <- list(Nborn=runif(1,100,10000),
                     mu.bday= runif(1,20,60),
                     skew.bday=runif(1,2,6),
                     sd.bday=runif(1,5,15)
  )
  # Compile the model and MCMC configuration
  model <- nimbleModel(code = code, 
                       inits = inits_list, 
                       data = data_list, 
                       constants = constants_list)
  
  Cmodel <- compileNimble(model)
  conf <- configureMCMC(model, useConjugacy = FALSE)
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # Run the model
  set.seed(0)
 # time_S3 <- system.time({
    samples <- runMCMC(Cmcmc, nchains = 3, niter = 10000, nburnin = 5000, thin = 10)
 # })
  
  # Save the results
  saveRDS(samples, file =  file.path(output_dir, paste0(selected_colony, "_samples_HgNimble_Model_S3_", Sys.Date(), ".rds")))
  }
})
# Print overall time
print(paste("Overall time for all models:", overall_time))
###############################################################################
# Load the saved samples and generate MCMC trace plots for each colony

for (i in 1:nrow(colyears)) {
  selected_colony <- colyears$Colony[i]
  selected_year <- colyears$Year[i]
  
  # Load the samples
  sample_file <- file.path(output_dir, paste0(selected_colony, "_samples_HgNimble_Model_S3_", Sys.Date(), ".rds"))
  if (file.exists(sample_file)) {
    loaded_samples <- readRDS(sample_file)
    
    # Generate MCMC trace plots
    pdf(file.path(output_dir, paste0(selected_colony, "_MCMC_trace_plots_", Sys.Date(), ".pdf")))
    MCMCtrace(object = loaded_samples,
              pdf = FALSE,
              ind = TRUE,
              Rhat = TRUE,
              n.eff = TRUE,
              exact = TRUE,
              params = "all")
    dev.off()
  } else {
    print(paste("File not found:", sample_file))
  }
}

samps<-as.matrix(Cmcmc$mvSamples)
samplesSummary(samps)
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
########################################################################################################
# Parallel execution of NIMBLE model and MCMC sampling



library(parallel)
nbcores <-detectCores() -75
my_cluster <-makeCluster(nbcores)
# Placeholder for filename dataframe
filenames_df <- tibble()
results <- list()
# Export necessary data to parallel workers
clusterExport(my_cluster, c("data", "colyears", "dskew", "int_pb","whiteSum","moultSum"))

results <- parLapply(my_cluster, 1:nrow(colyears), function(i) {
  library(nimble)
  library(coda)
  library(dplyr)
  
  
  code <- nimbleCode({
    #priors
    Nborn ~ dunif(100,5000) 
    #mean.nborn ~  dgamma(shape=485, scale=1) 
    #Nborn ~ dpois(mean.nborn)  
    #Nborn ~ dgamma(shape=485, scale=1) 
    mu.bday ~ dunif(10,60)
    skew.bday ~ dunif(2,6)
    sd.bday ~ dunif(1,15)
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
    
    
    # cumulative sum of pups born on given day 
    # for (i in 1:ndays) {
    #   cumsum_born[i] <- sum(born[1:i])
    # }
    # account for probabilities of birth across season and 
    # transitions to moult and leave stages 
    
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
    #   #wObs[obs] ~ dbinom(POWhite, round(whitePups[obsdays[obs]]))
    #   #mObs[obs] ~ dbinom(POMoult, round(moultPups[obsdays[obs]]))
    # 
    # 
    # 
    # # #   # add misclassification
    #   whiteMC[obs] <- whitePups[obsdays[obs]] + moultPups[obsdays[obs]] * (1-PCMoult)
    #   moultMC[obs] <- moultPups[obsdays[obs]] * PCMoult
    #   wObs[obs] ~ dbinom(POWhite,round(whiteMC[obs]))
    #   mObs[obs] ~ dbinom(POMoult,round(moultMC[obs]))
    # #
    # }
    #####################################################################
    ###
    ### Normal approximation using binomial variance
    ###
    #####################################################################
    for (obs in 1:length(obsdays)) {
      whiteCount[obs] <- whitePups[obsdays[obs]] * POWhite + moultPups[obsdays[obs]] * POMoult * (1 - PCMoult)
      sdwhite[obs] <- sqrt(whitePups[obsdays[obs]] * POWhite * (1 - POWhite) + moultPups[obsdays[obs]] * (1 - PCMoult) * POMoult * (1-POMoult*(1-PCMoult)))

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
  
# Extract data for the current colony
realdata <- subset(data, Colony == colyears$Colony[i] & Year == colyears$Year[i])
wObs <- realdata$Whitecoats + realdata$DeadPups
mObs <- realdata$Moulteds
obsdays <- unique(realdata$Day)

data_list <- list(wObs = wObs, mObs = mObs)

constants_list <- list(
  whiteObs = wObs,
  moultObs = mObs,
  obsdays = obsdays,
  ndays = max(data$Day),
  mu.mday=23,
  sd.mday=5,
  mu.lday=31.5,
  sd.lday=7,
  amax = 123,
  POWhite = 0.94,
  POMoult = 0.98,  
  PCMoult = 0.84
)

inits_list <- list(
  mu.bday= 45,#runif(1,1,60),
  skew.bday=runif(1,2,6),
  sd.bday=runif(1,1,15),
  Nborn = 5000
)


set.seed(0)  # Set seed for reproducibility

model <- nimbleModel(
  code = code,
  inits = inits_list,
  data = data_list,
  constants = constants_list
)
# Calculate model results
model$calculate()

conf <- configureMCMC(model, useConjugacy = FALSE)
Rmcmc <- buildMCMC(conf)  # Create MCMC configuration based on the compiled model
Cmodel <- compileNimble(model)  # Compile the NIMBLE model, translating it to C++
Cmcmc <- compileNimble(Rmcmc, project = model)  # Compile the MCMC configuration and link it to the already compiled model

# Update model data
Cmodel$setData(data_list)

# Run the MCMC
mcmc.out <- runMCMC(Cmcmc,
                    niter = 20000,
                    nburnin = 10000,
                    nchains = 3,
                    setSeed = FALSE,
                    samples = TRUE,
                    samplesAsCodaMCMC = FALSE)



#cat("MCMC time for", colyears$Colony[i], ":", total_time, "\n")

# Store the results
result <- list(
  mcmc = mcmc.out
  #total_time = total_time
)

return(result)
})



# Close the parallel cluster
stopCluster(my_cluster)

names(results) <- paste(colyears$Colony, colyears$Year, sep = "_")

# Directory where you want to save the .rds files from remote desktop:
output_dir <- "D:/mg349/results"
# Iterate over each colony's results

for (colony_name in names(results)) {
  mcmc_object <- results[[colony_name]]$mcmc  # Extract mcmc object
  
  # Construct the file name for saving
  file_name <- file.path(output_dir, paste0(colony_name, "_mcmc_results_HgNimble_Model_S3_", Sys.Date(), ".rds"))
  
  # Save the mcmc results for the current colony as .rds file
  saveRDS(mcmc_object, file = file_name)
}




summary_stats <- list()
library(coda)
# Iterate over each colony's MCMC results
for (colony_name in names(results)) {
  mcmc_object <- as.mcmc(do.call(rbind, results[[colony_name]]$mcmc))  # Combine chains and convert to mcmc object
  
  # Plot trace plots for each parameter
  par(mfrow = c(4, 2))  # Set layout for 2x2 grid
  for (param in c("mu.bday", "sd.bday", "Nborn", "skew.bday")) {
    traceplot(mcmc_object[, param], main = paste("Trace Plot of", param, "for", colony_name))
    
    # Density plot with mean value
    mean_value <- mean(mcmc_object[, param])
    dens <- density(mcmc_object[, param])
    plot(dens, main = paste("Density Plot of", param, "for", colony_name))
    #abline(v = mean_value, col = "red", lwd = 2)
  }
  
  # Calculate summary statistics for each parameter
  param_summary <- list()
  for (param in c("mu.bday", "sd.bday", "Nborn", "skew.bday")) {
    param_summary[[param]] <- summary(mcmc_object[, param])
  }
  
  # Store the summary statistics for the colony
  summary_stats[[colony_name]] <- param_summary
}
par(mfrow = c(1, 1))  # Reset layout

# # Print summary statistics for each colony
for (colony_name in names(summary_stats)) {
  cat("Summary statistics for", colony_name, ":\n")
  print(summary_stats[[colony_name]])
  cat("\n")
}

# Prepare to write all summary statistics to CSV file
summary_df <- data.frame()  # Initialize an empty data frame

# Iterate over each colony's summary statistics
for (colony_name in names(summary_stats)) {
  param_summary <- summary_stats[[colony_name]]
  
  # Create a data frame for the colony's summary statistics
  colony_summary_df <- data.frame(
    colony = colony_name,
    parameter = c("mu.bday", "sd.bday", "Nborn", "skew.bday"),
    mean = sapply(param_summary, function(x) x$statistics["Mean"]),
    sd = sapply(param_summary, function(x) x$statistics["Std. Dev."]),
    min = sapply(param_summary, function(x) x$statistics["Min."]),
    max = sapply(param_summary, function(x) x$statistics["Max."])
  )
  
  # Append the colony's summary statistics to the overall summary data frame
  summary_df <- rbind(summary_df, colony_summary_df)
}


# Save the summary data frame as a CSV file
output_csv_file <- "D:/mg349/results/summary_stats_S3_Farnes2018.csv"
write.csv(summary_df, file = output_csv_file, row.names = FALSE)


