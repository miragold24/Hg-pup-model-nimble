
library(nimble)
library(coda)
library(dplyr)


######################################################
#add data here
datafile <- "C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/04 - Raw Data/HgPupCountsForProdModel_2018_NEEfs.csv"
data <- read.csv(datafile)
colyears <- data %>% select(Colony, Year) %>% distinct()
# Convert Date column to Date type
data$Date <- as.Date(data$Date, format = "%d/%m/%Y")
# Calculate the day number from October 1
start_date <- as.Date("01/10/2018", format = "%d/%m/%Y")
data$Day <- as.numeric(difftime(data$Date, start_date, units = "days")) + 1
##############################################################################


# Define custom NIMBLE functions
nim_dskew <- nimbleFunction(
  run = function(x = double(1), param = double(1)) {
    skew.bday <- param[1]
    mu.bday   <- param[2]
    sd.bday   <- param[3]
    delta <- skew.bday / sqrt(1 + (skew.bday^2))
    omega <- sqrt(sd.bday^2 / (1 - (2 * delta^2) / pi))
    xi <- mu.bday - omega * delta * sqrt(2 / pi)
    pb <- (2 / omega) * dnorm((x - xi) / omega) * pnorm(skew.bday * ((x - xi) / omega))
    returnType(double(1))
    return(pb)
  })

intPb <- nimbleFunction(
  run = function(days = double(1), skew.bday = double(), mu.bday = double(), sd.bday = double()) {
    Pb <- rep(0, length(days))
    pars <- c(skew.bday, mu.bday, sd.bday)
    for (t in 1:length(days)) {
      Pb[t] <- nimIntegrate(nim_dskew, t-1, t, pars)[1]
    }
    returnType(double(1))
    return(Pb)
  })
###
### Parallel
###
library(parallel)
nbcores <-detectCores() -3
my_cluster <-makeCluster(nbcores)
# Placeholder for filename dataframe
filenames_df <- tibble()
results <- list()
# Export necessary data to parallel workers
clusterExport(my_cluster, c("data", "colyears", "nim_dskew", "intPb"))

results <- parLapply(my_cluster, 1:nrow(colyears), function(i) {
  library(nimble)
  library(coda)
  library(dplyr)

code <- nimbleCode({
  
  # set prior on daily survival intercept and slope
  #surv.int.real ~ dbeta(175,1)
  #surv.int <-  4    # log(surv.int.real/(1-surv.int.real))       # logit scale
  #surv.slope <- 0.1 # ~ dunif(0,.5)                               # daily increase in survival (logit scale) 
  
  # set prior
  
  # set prior on detection and misclassification processes
  #Pobs<- 0.95 # ~ dunif(0.8,1) #   #
  Pobs_m <- 0.95   # ~ dbeta(95,5)
  Pobs_w <- 0.95   # ~ dbeta(95,5)
  PCMoult<-0.91# ~ dunif(0.7,0.9) # #~ 
  
  # BIRTH DISTRIBUTION
  mu.bday ~ dunif(10,70) #dnorm(45,0.01)T(20,60)      ## Prior for mean of birth distribution (1 would mean born just before day 0 [t=1], 10 would mean born just before day 9 [t=10])
  sd.bday ~ dunif(1,20)                 ## Prior for standard deviation of birth distribution
  skew.bday ~ dunif(2,10) #<- 5
  nborn ~ dunif(100,50000)       ## prior for mean number of births (total across season)
  mu.days.on.colony <- 31.5#~ dunif(min.age, max.age)       ## Prior for mean of the number of days on the colony
  sd.days.on.colony <- 7#~ dunif(4, 15)        ## Prior for standard deviation of days on colony
  
  
  
  #=================================
  # set the probability of being moulted as a function of age
  #===================================
  
  pmoult.int <- -12   # ~ dnorm(5,1)               # probability of being moult at birth
  pmoult.slope <- 0.6 # ~ dunif(-0.5,0)          # daily decrease in the probability of being moult (increase in probability of moulting)
  for(a in 1:ndays){   # loop past the last possible age
    logit(pmoult[a]) <- pmoult.int + pmoult.slope*(a-1)    # deterministic for now
  }
  #=================================
  # set the probability of survival as a function of age - in this version all pups survive with p=1
  #===================================
  
  # for(a in 1:(max.age+1)){
  #   logit(surv[a]) <- surv.int + surv.slope*(a-1)
  # }
  # 
  #=======================================
  # set the scheduling distributions 
  #=======================================
  
  # BIRTH DISTRIBUTION
  #normal 
  # nconst_Pb <- 1 - ( pnorm(0,mu.bday, sd.bday) + (1 - pnorm((max.day.birth),mu.bday, sd.bday) ) )  # constant for truncating the births distribution
  # Pb[1] <- 0  # can't be born after 0 time has passed in the season
  # Pb[2:(max.day.birth+1)] <- (pnorm(1:(max.day.birth), mu.bday, sd.bday) - pnorm(0:(max.day.birth-1), mu.bday, sd.bday))/nconst_Pb   ## Birth schedule
  
  #skew normal using generic function
  pb[1:(max.day.birth+1)] <- intPb(1:(max.day.birth+1), skew.bday, mu.bday, sd.bday)
  Pb[1:(max.day.birth+1)] <- pb[1:(max.day.birth+1)]/sum(pb[1:(max.day.birth+1)])
  
  # snconst_Pb <- 1 - (nim_pskew(0,skew.bday,mu.bday,sd.bday) + (1 - nim_pskew((max.day.birth+1),skew.bday,mu.bday,sd.bday) ) ) # constant for truncating the births distribution
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
      stageN_c[c,c+a-1,1] <- cohort.size[c] * cprob.remain[a]  * (1-pmoult[a])   #* prod(surv[1:a])  # probability of being white after a days
      stageN_c[c,c+a-1,2] <- cohort.size[c] * cprob.remain[a] * pmoult[a]        #* prod(surv[1:a])  # moulted
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
  
  # for(obs in 1:(length(obsdays))) {
  #   white.mu[obs] <- stageN[obsdays[obs],1] * Pobs
  #   white.sd[obs] <- pow(stageN[obsdays[obs],1] * Pobs,.5)*(1-Pobs)
  #   white.sd[obs] <- pow(white.sd[obs],-2)
  #   moult.mu[obs] <- stageN[obsdays[obs],2] * Pobs
  #   moult.sd[obs] <- pow(stageN[obsdays[obs],2] * Pobs,.5)*(1-Pobs)
  #   moult.sd[obs] <- pow(moult.sd[obs],-2)
  #   wObs[obs] ~ dnorm(white.mu[obs],white.sd[obs])
  #   mObs[obs] ~ dnorm(moult.mu[obs],moult.sd[obs])
  # }
  
  # for(obs in 1:(length(obsdays))) {
  #  wObs[obs] ~ dbinom(Pobs, trunc(stageN[obsdays[obs],1]))   # data node
  #  mObs[obs] ~ dbinom(Pobs, trunc(stageN[obsdays[obs],2]))   # data node
  # }
  
  
  
  WW_MW[1:(ndays-1)] <- stageN[1:(ndays-1),1]+stageN[1:(ndays-1),2]*(1-PCMoult)
  MM[1:(ndays-1)] <- stageN[1:(ndays-1),2]*PCMoult  #add misclassification to mObs
  
  for(obs in 1:(length(obsdays))) {
    wObs[obs] ~ dbinom(Pobs_w, trunc(WW_MW[obsdays[obs]]))   # data node
    mObs[obs] ~ dbinom(Pobs_m, trunc(MM[obsdays[obs]]))   # data node
  }
  
})

# Parameters monitored
params <- c( "mu.bday", "sd.bday","nborn",
              'skew.bday'
             #"sd.days.on.colony",
             #"mu.days.on.colony",
             #"surv.int","surv.slope","pwhite.int","pwhite.slope",
             #"pmass.birth.to.leave","cprob.remain","Pb",
             #"stageN")#,"Pobs","cohort.size")
)
# ndays <- 123
# max.age <- 60
# min.age <- 20
# max.day.birth <- 100


# Parallel execution of NIMBLE model and MCMC sampling
#total_time <- system.time({

    # Extract data for the current colony
    realdata <- subset(data, Colony == colyears$Colony[i] & Year == colyears$Year[i])
    wObs <- realdata$Whitecoats + realdata$DeadPups
    mObs <- realdata$Moulteds
    obsdays <- unique(realdata$Day)
    
    data_list <- list(wObs = wObs, mObs = mObs)
    
    constants_list <- list(
      ndays = 123, 
      max.age = 60, 
      min.age = 20,  
      max.day.birth = 100,
      #mu.days.on.colony=31.5,
      obsdays = unique(realdata$Day)
    )
    
    inits_list <- list(
      mu.bday = 45,
      sd.bday = 15,
      skew.bday= 5,
      #mu.days.on.colony = 31.5,
      #sd.days.on.colony = 10,
      nborn = 5000
      #cohort.size= round(10000*(sim.wm$Pb)) 
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
                        niter = 15000,
                        nburnin = 5000,
                        nchains = 2,
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

########################
##
## Save results
##
########################

names(results) <- paste(colyears$Colony, colyears$Year, sep = "_")
# Directory where you want to save the .rds files
output_dir <- "C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/07 - GitHub Repositories/Hg_pup_model/S3_LR/real data results"
# Iterate over each colony's results
for (colony_name in names(results)) {
  mcmc_object <- results[[colony_name]]$mcmc  # Extract mcmc object
  
  # Construct the file name for saving
  file_name <- file.path(output_dir, paste0(colony_name, "_mcmc_results3chain_PCMoult.91_", Sys.Date(), ".rds"))
  
  # Save the mcmc results for the current colony as .rds file
  saveRDS(mcmc_object, file = file_name)
}

####
# Save model info 

###########################
## RELOAD MCMC OBJECT
###########################

# Get a list of all .rds files in the directory
rds_files <- list.files(output_dir, pattern = "_mcmc_results3chain_PCMoult.5_.*.rds", full.names = TRUE)

# Load each .rds file into a list
results <- list()
for (file_path in rds_files) {
  # Extract colony name from the file name, replacing spaces with underscores
  colony_name <- sub("_mcmc_results3chain_PCMoult.5_.*", "", basename(file_path))
  #colony_name <- gsub(" ", "_", colony_name)  # Replace spaces with underscores
  
  # Load the MCMC object from the .rds file
  results[[colony_name]] <- list(mcmc = readRDS(file_path))
}

# Verify the list has entries for all five colonies
print(names(results))  # Should show all five colony names with underscores
###########################
## Look at trace/density
###########################
summary_stats <- list()
param_list<- c("mu.bday", "nborn", "sd.bday","skew.bday") #,"mu.days.on.colony")
# Iterate over each colony's MCMC results
for (colony_name in names(results)) {
  mcmc_object <- as.mcmc(do.call(rbind, results[[colony_name]]$mcmc))  # Combine chains and convert to mcmc object
  
  # Plot trace plots for each parameter
  par(mfrow = c(4, 2))  # Set layout for 2x2 grid
  for (param in param_list) {
    traceplot(mcmc_object[, param], main = paste("Trace Plot of", param, "for", colony_name))
    
    # Density plot with mean value
    mean_value <- mean(mcmc_object[, param])
    dens <- density(mcmc_object[, param])
    plot(dens, main = paste("Density Plot of", param, "for", colony_name))
    abline(v = mean_value, col = "red", lwd = 2)
  }
  
  # Calculate summary statistics for each parameter
  param_summary <- list()
  for (param in param_list) {
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

# # Prepare to write all summary statistics to CSV file
# summary_df <- data.frame()  # Initialize an empty data frame
# 
# # Iterate over each colony's summary statistics
# for (colony_name in names(summary_stats)) {
#   param_summary <- summary_stats[[colony_name]]
#   
#   # Create a data frame for the colony's summary statistics
#   colony_summary_df <- data.frame(
#     colony = colony_name,
#     parameter = c("mu.bday", "nborn", "sd.bday"),
#     mean = sapply(param_summary, function(x) x$statistics["Mean"]),
#     sd = sapply(param_summary, function(x) x$statistics["Std. Dev."]),
#     min = sapply(param_summary, function(x) x$statistics["Min."]),
#     max = sapply(param_summary, function(x) x$statistics["Max."])
#   )
#   
#   # Append the colony's summary statistics to the overall summary data frame
#   summary_df <- rbind(summary_df, colony_summary_df)
# }

# Write all summary statistics to the CSV file
#write.csv(summary_df, csv_file, row.names = FALSE)


# Iterate over each colony's MCMC results
for (colony_name in names(results)) {
  mcmc_object <- as.mcmc(do.call(rbind, results[[colony_name]]$mcmc))  # Combine chains and convert to mcmc object
  
  # Set up PDF file path
  pdf_file <- file.path(output_dir, paste(date,"HgNimble9_multicolony2_parallel_rm_surv_TTL31.5_plots_", colony_name, ".pdf", sep = ""))
  
  # Start PDF device
  pdf(pdf_file)
  
  # Plot trace plots and density plots for each parameter
  par(mfrow = c(4, 2))  # Set layout for 4x2 grid
  
  for (param in param_list) {
    # Trace plot
    traceplot(mcmc_object[, param], main = paste("Trace Plot of", param, "for", colony_name))
    
    # Density plot with mean value
    mean_value <- mean(mcmc_object[, param])
    dens <- density(mcmc_object[, param])
    plot(dens, main = paste("Density Plot of", param, "for", colony_name))
    abline(v = mean_value, col = "red", lwd = 2)
  }
  
  # Close PDF device
  dev.off()
}





# Load necessary libraries
library(coda)
library(HDInterval)

# Create an empty list to store summary statistics
summary_stats2 <- list()

# Iterate over each colony's MCMC results
for (colony_name in names(results)) {
  mcmc_object <- as.mcmc(do.call(rbind, results[[colony_name]]$mcmc))  # Combine chains and convert to mcmc object
  
  # Initialize a list to store parameter summaries
  param_summary <- list()
  
  for (param in param_list) {
    # Extract parameter samples
    param_samples <- mcmc_object[, param]
    
    # Compute summary statistics
    mean_value <- mean(param_samples)
    median_value <- median(param_samples)
    sd_value <- sd(param_samples)
    
    # Compute HPDI
    hpdi_interval <- hdi(param_samples, credMass = 0.95)  # 95% HDI
    hpdi_low <- hpdi_interval[1]
    hpdi_high <- hpdi_interval[2]
    
    # Store statistics in list
    param_summary[[param]] <- list(
      Mean = mean_value,
      Median = median_value,
      SD = sd_value,
      HPDI_Low = hpdi_low,
      HPDI_High = hpdi_high
    )
  }
  
  # Store the summary statistics for the colony
  summary_stats2[[colony_name]] <- param_summary
}

# Convert summary statistics list to a data frame
# Prepare data frame to store all results
results_df <- data.frame()

# Iterate over summary statistics to build the data frame
for (colony_name in names(summary_stats2)) {
  colony_summary <- summary_stats2[[colony_name]]
  
  for (param in names(colony_summary)) {
    param_stats <- colony_summary[[param]]
    results_df <- rbind(results_df, data.frame(
      Colony = colony_name,
      Parameter = param,
      Mean = param_stats$Mean,
      Median = param_stats$Median,
      SD = param_stats$SD,
      lower_hpd = param_stats$HPDI_Low,
      higher_hpd = param_stats$HPDI_High
    ))
  }
}

# Write the results to a CSV file
#write.csv(results_df, file = "HgNimble9_multicolony2_parallel_rm_surv_summarystats.csv", row.names = FALSE)

# Metadata
script_filename <- "HgNimble9_multicolony2_parallel_rm_surv.R"
rdsfile_dir <- "C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/07 - GitHub Repositories/Hg_pup_model/S3_LR/real data results"
priors <- "mu.bday ~ dunif(10,70); sd.bday ~ dunif(1,20); skew.bday ~ (2,10); nborn ~ dunif(100,50000); "
# # skew.bday ~ dunif(2,10) ;
initial_values <- paste0("mu.bday = ", inits_list$mu.bday, "; sd.bday = ", inits_list$sd.bday, "; skew.bday = ", inits_list$skew.bday, "; Nborn = ", inits_list$nborn) #, "; mu.days.on.colony = ", inits_list$mu.days.on.colony, ";,"; sd.days.on.colony = ", inits_list$sd.days.on.colony)
mu.days.on.colony <- 31.5
pcmoult <- 0.5
chains <- 2
iterations <- 10000
burnin <- 5000
output_dir <-"C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/07 - GitHub Repositories/Hg_pup_model/S3_LR/real data results"
date <- Sys.Date()
# Create a data frame with metadata
metadata <- data.frame(
  ScriptFilename = script_filename,
  Priors = priors,
  InitialValues = initial_values,
  mu.days.on.colony = mu.days.on.colony,
  pcmoult = pcmoult,
  Chains = chains,
  Iterations = iterations,
  Burnin = burnin,
  date=date,
  OutputDir = output_dir,
  rdsfiles = rdsfile_dir,
  stringsAsFactors = FALSE
)


# Create a filename with the date
filename <- paste0(date,"_HgNimble9_multicolony2_parallel_rm_surv_summarystats_pcmoult0.5", ".csv")

# Combine the output directory with the filename
filepath <- file.path(output_dir, filename)

# Write the metadata to the CSV file in the output directory
write.csv(metadata, file = filepath, row.names = FALSE)

# Append the results data frame to the CSV file in the output directory
write.table(results_df, file = filepath, append = TRUE, sep = ",", row.names = FALSE, col.names = TRUE)
