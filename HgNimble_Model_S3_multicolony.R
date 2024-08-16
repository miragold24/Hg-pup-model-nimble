# Load necessary libraries
library(nimble)
library(coda)
library(dplyr)
library(parallel)
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


# Placeholder for filename dataframe
filenames_df <- tibble()
results <- list()

# Iterate over each row in colyears to process each colony and year sequentially
for (i in 1:nrow(colyears)) {
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
    mu.mday = 23,
    sd.mday = 5,
    mu.lday = 31.5,
    sd.lday = 7,
    amax = 123,
    POWhite = 0.94,
    POMoult = 0.98,
    PCMoult = 0.84
  )
  
  inits_list <- list(
    mu.bday = 45,
    sd.bday = 10,
    skew.bday = 4,
    Nborn = 5000
  )
  
  code <- nimbleCode({
    #priors
    Nborn ~ dunif(100, 10000)
    mu.bday ~ dunif(33, 60)
    skew.bday ~ dunif(2, 6)
    sd.bday ~ dunif(7, 15)
    
    # deterministic and continuous births over breeding season (ndays)
    for (t in 1:ndays){ 
      born[t] <- Pb_sum[t] * Nborn
    }
    
    #calculate PDF of skew normal then integrate over intervals for Pb on given day
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
      sdwhite[obs] <- sqrt(whitePups[obsdays[obs]] * POWhite * (1 - POWhite) + moultPups[obsdays[obs]] * (1 - PCMoult) * POMoult * (1 - (1 - PCMoult) * POMoult))
      
      moultCount[obs] <- moultPups[obsdays[obs]] * POMoult * PCMoult
      sdmoult[obs] <- sqrt(moultPups[obsdays[obs]] * POMoult * PCMoult * (1 - POMoult * PCMoult))
      
      wObs[obs] ~ dnorm(whiteCount[obs], sd = sdwhite[obs])
      mObs[obs] ~ dnorm(moultCount[obs], sd = sdmoult[obs])
    }
  })
  
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
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(model)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
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
  
  # Store the results
  result <- list(
    mcmc = mcmc.out
  )
  
  # Save the results
  colony_name <- paste(colyears$Colony[i], colyears$Year[i], sep = "_")
  results[[colony_name]] <- result
  
  # Directory where you want to save the .rds files from laptop
  #output_dir <- "C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/05 - Results/real data"
  
  # Directory on remote desktop (blue whale)
  #output_dir <- "D:/mg349/results"
  file_name <- file.path(output_dir, paste0(colony_name, "_mcmc_results_HgNimble_Model_S3_", Sys.Date(), ".rds"))
  
  # Save the mcmc results for the current colony as .rds file
  saveRDS(result$mcmc, file = file_name)
}