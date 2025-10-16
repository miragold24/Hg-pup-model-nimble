#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
library(nimble)
library(dplyr)
library(tidyr)
library(MCMCvis)
library(coda)
library(ggplot2)
library(HDInterval)
#StA drive file path
#work_dir<-("C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/07 - GitHub Repositories/Hg_pup_model")
#setwd(work_dir)
#getwd() 

###################################
#Functions 

source("HgNimble_input_S3.R") # Nimble functions to count whitecoat and moulted pups present on any given day
source("HgSim_S3.R")
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
  whitePups[1:ndays] <- whiteSum(born[1:ndays], Pm_sum[1:amax])#whiteObs[]),Pl_sum[1:amax],
  moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax])#moultObs[]),
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


# directory to save results tto remote desktop
output_dir <- "D:/mg349/results/HgNimble_Model_S3"
# directory local comp
#output_dir <- "C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/07 - GitHub Repositories/Hg_pup_model/results/HgNimble_Model_S3_simulation"
# Create a function to run the model and save results for a given dataset
run_model_simulation <- function(simulation_number) {
  # Generate a new dataset
  # Generate observations using the generate_obs function
  sim_S3 <- simSealsS3(
    N = 3000, 
    ndays = 123, amax = 123, mu.bday = 45, sd.bday = 10,
    birthcurve = "skewnormal",
    skew.bday = 5, 
    mu.mday = 23, sd.mday = 5,
    mu.lday = 31.5, sd.lday = 7,
    POWhite = 0.95,
    POMoult = 0.95, 
    PCWhite = 1, 
    PCMoult = 0.91,  
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
  cat("mObs:", mObs, "\n")
  cat("========================================\n")
  
  # Define the data, constants, and initial values
  data_list <- list(wObs = wObs, mObs = mObs) #whiteObs = wObs, moultObs = mObs)
  
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
    PCMoult = 0.91
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
    
    whitePups[1:ndays] <- whiteSum(born[1:ndays], Pm_sum[1:amax])#, whiteObs[])
    moultPups[1:ndays] <- moultSum(born[1:ndays], Pm_sum[1:amax], Pl_sum[1:amax])#, moultObs[])
    
    
    
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
  
  
  samples <- runMCMC(Cmcmc, nchains = 3, niter = 10000, nburnin = 5000, thin = 5)
  
  # Save results
  saveRDS(samples, file = file.path(output_dir, paste0("HgNimble_Model_S3_","sim_", simulation_number,"_", Sys.Date(), ".rds")))
}

# Run simulations
n_simulations <- 100
for (i in 1:n_simulations) {
  run_model_simulation(i)
}

# Convert each chain to an mcmc object
mcmc_list <- lapply(samples, as.mcmc)
# Combine into a single mcmc.list object
mcmc_list <- as.mcmc.list(mcmc_list)
# Gelman-Rubin diagnostic to check convergence
gelman.diag(mcmc_list)
################################################################################
# read in simulation, visualize plots
# directory to save results tto remote desktop
output_dir <- "D:/mg349/results/HgNimble_Model_S3"
# directory local comp
date <- Sys.Date()-1

# List all RDS files for the current date
rds_files <- list.files(output_dir, pattern = paste0("HgNimble_Model_S3_.*_", date, "\\.rds$"), full.names = TRUE)

# Initialize a list to store all MCMC chains
all_chains <- list()

# Load and combine MCMC chains
for (file in rds_files) {
  samples <- readRDS(file)
  
  # Convert each chain to mcmc object
  chains <- lapply(samples, function(chain) as.mcmc(chain))
  
  # Combine chains into a single mcmc.list object
  combined_chains <- mcmc.list(chains)
  
  all_chains <- append(all_chains, list(combined_chains))
}

# Initialize a list to store all samples
all_samples <- list()

# Load and combine MCMC chains
for (file in rds_files) {
  # Read the MCMC samples
  samples <- readRDS(file)
  
  # Convert each sample to a data frame
  samples_df <- lapply(samples, as.data.frame)
  
  # Add a column for simulation number
  samples_df <- lapply(seq_along(samples_df), function(i) {
    samples_df[[i]] %>%
      mutate(simulation = basename(file)) %>%
      mutate(chain = i)
  })
  
  # Combine the data frames for this file
  combined_samples <- bind_rows(samples_df)
  
  # Append to the list of all samples
  all_samples <- append(all_samples, list(combined_samples))
}

# Combine all data frames into one
combined_samples <- bind_rows(all_samples)

# Convert to long format for easier summarization
long_samples <- combined_samples %>%
  pivot_longer(-c(simulation, chain), names_to = "parameter", values_to = "value")

# Define the true values for each parameter
true_values <- list(
  "Nborn" = 3000,
  "sd.bday" = 10, 
  "mu.bday" = 45, 
  "skew.bday" = 5
)

# Calculate summary statistics across all simulations
summary_stats <- long_samples %>%
  group_by(parameter) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    .groups = 'drop'
  )

# Function to calculate HPD intervals
calculate_hpd <- function(param_values) {
  hpd_interval <- hdi(param_values, prob = 0.95)
  return(c(lower = hpd_interval[1], upper = hpd_interval[2]))
}

# Calculate Relative Bias and HPD intervals
result_summary <- summary_stats %>%
  rowwise() %>%
  mutate(
    true_value = ifelse(parameter %in% names(true_values), true_values[[parameter]], NA),
    relative_bias = ifelse(!is.na(true_value), (mean - true_value) / true_value* 100, NA),
    hpd_interval = list(calculate_hpd(long_samples$value[long_samples$parameter == parameter])),
    lower_hpd = hpd_interval[[1]],
    upper_hpd = hpd_interval[[2]]
  ) %>%
  ungroup() %>%
  # Round numeric columns to 3 decimal places
  mutate(
    mean = round(mean, 3),
    sd = round(sd, 3),
    median = round(median, 3),
    relative_bias = round(relative_bias, 3),
    lower_hpd = round(lower_hpd, 3),
    upper_hpd = round(upper_hpd, 3)
  )
# Compute precision 1/ (sd^2)
result_summary <- result_summary %>%
  mutate(precision = 1 / (sd^2))
# Print the result summary
print(result_summary)

# Define the path to save the CSV file
#output_csv_path <- "D:/mg349/results/HgNimble_Model_S3/summary_statistics_HgNimble_S3_sim.csv"
output_csv_path <-file.path(output_dir, paste0(Sys.Date(),"_HgNimble_Model_S3_","summarystats_sim", ".csv"))
# Remove the 'hpd_interval' column
result_summary <- result_summary %>%
  select(-hpd_interval)

# Write the modified result summary to a CSV file
readr::write_csv(result_summary, output_csv_path)

#############################################################
params_of_interest <- c("Nborn","sd.bday","mu.bday","skew.bday")
plots <- list()
# Plot results for each parameter
for (param in params_of_interest) {
  # Filter long_samples for the current parameter
  param_samples <- long_samples %>% filter(parameter == param)
  
  # Check if HPD interval columns are available
  param_summary <- summary_stats %>% filter(parameter == param)
  
  # Debug: Print param_summary to check column names
  print(param_summary)
  
  # Create plot
  p <- ggplot(param_samples, aes(x = value, color = simulation)) +
    geom_density() +
    facet_wrap(~ parameter, scales = "free") +
    labs(title = paste("Posterior Density of", param),
         x = param,
         y = "Density") +
    scale_color_discrete(guide = "none") +
    theme_minimal() +
    theme(legend.position = "none") +      
    guides(color = guide_none())        
  
  # Add vertical line for the true value
  true_value <- true_values[[param]]
  p <- p + geom_vline(xintercept = true_value, linetype = "dashed", color = "red", size = 1)
  
  # Add HPD intervals if available
  if ("lower_hpd" %in% colnames(param_summary) && "upper_hpd" %in% colnames(param_summary)) {
    p <- p + geom_ribbon(data = param_summary,
                         aes(xmin = lower_hpd, xmax = upper_hpd, ymin = 0, ymax = Inf),
                         alpha = 0.2, fill = "blue", color = NA)
  } else {
    warning(paste("HPD interval columns not found for parameter:", param))
  }
  
  plots[[param]] <- p
}

# Display the plots
for (param in params_of_interest) {
  print(plots[[param]])
}
######################


## plot histogram and posterior density
# Define the descriptive names for the parameters
param_labels <- c(
  "Nborn" = "Total Pups Born",
  "sd.bday" = "SD Birthday", 
  "mu.bday" = "Mean Birthday", 
  "skew.bday" = "Skew Birthday"
)

# Function to calculate HPD intervals
calculate_hpd <- function(values) {
  hdi(values, prob = 0.95)
}

# Calculate HPD intervals for all parameters
hpd_intervals <- long_samples %>%
  group_by(parameter) %>%
  summarize(lower_hpd = calculate_hpd(value)[1], upper_hpd = calculate_hpd(value)[2], .groups = 'drop')

  # Plot results for each parameter
  for (param in params_of_interest) {
    # Filter long_samples for the current parameter
    param_samples <- long_samples %>% filter(parameter == param)
    
    # Calculate density for the current parameter
    density_data <- density(param_samples$value)
    
    # Get HPD intervals for the current parameter
    param_hpd <- hpd_intervals %>% filter(parameter == param)
    
    # Create a data frame for the density and HPD intervals
    density_df <- data.frame(x = density_data$x, y = density_data$y)
    density_df$within_hpd <- ifelse(density_data$x >= param_hpd$lower_hpd & density_data$x <= param_hpd$upper_hpd, density_data$y, NA)
    
    # Create plot with density and histogram
    p <- ggplot(param_samples, aes(x = value)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "lightgray", color = "black", alpha = 0.5) +
      geom_line(data = density_df, aes(x = x, y = y), color = "blue") +
      geom_ribbon(data = density_df, aes(x = x, ymin = 0, ymax = within_hpd), fill = "blue", alpha = 0.2) +
      labs(#title = ,  # Use descriptive name for the title
        x = param_labels[[param]],
        y = "Frequency") +
      theme_minimal() +
      theme(
        axis.title.x = element_text( size = 14),  # face = "bold", Make x-axis label bold and increase font size
        axis.title.y = element_text( size = 14),  # face = "bold",
        axis.text.x = element_text(size = 12),  # Increase size of x-axis tick labels
        axis.text.y = element_text(size = 12),  
        #plot.title = element_text(hjust = 0.5),  # Center the title
        theme(legend.position = "none")
      )
    # Add vertical line for the true value
    true_value <- true_values[[param]]
    p <- p + geom_vline(xintercept = true_value, linetype = "dashed", color = "red", size = 1)
    
  # Save the plot as PDF
  pdf_file <- paste0("D:/mg349/results/HgNimble_Model_S3/", Sys.Date(), "_", param, "_sim_100_hist_95CI.pdf")
  ggsave(pdf_file, plot = p, device = "pdf", width = 10, height = 6)
}

