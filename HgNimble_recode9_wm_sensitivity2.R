# Set up model and MCMC configuration outside the loop
library(nimble)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

nim_dskew <- nimbleFunction(
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
#dskew_pb <- compileNimble(nim_dskew)
#nim_dskew(1:100,c(skew.bday,mu.bday,sd.bday))
intPb <- nimbleFunction(
  run = function(days = double(1), skew.bday = double(), mu.bday = double(), sd.bday = double()) {
    Pb <- rep(0, length(days))
    pars <- c(skew.bday, mu.bday, sd.bday)
    for(t in 1:length(days)) {
      Pb[t] <- nimIntegrate(nim_dskew, t-1, t, pars)[1]
    }
    returnType(double(1))
    return(Pb)
  })


#######################################################################

run_simulation <- function(scenario_index, obsdays_combinations) {
code <- nimbleCode({
  
  # set prior on daily survival intercept and slope
  #surv.int.real ~ dbeta(175,1)
  surv.int <-  4    # log(surv.int.real/(1-surv.int.real))       # logit scale
  surv.slope <- 0.1 # ~ dunif(0,.5)                               # daily increase in survival (logit scale) 
  
  # set prior
  
  # set prior on detection and misclassification processes
  Pobs <- 0.95   # ~ dunif(0.9,1)
  PCMoult <-0.84 #~ dunif(0.8,0.9)
  

  # BIRTH DISTRIBUTION
  mu.bday ~ dunif(10,70) #dnorm(45,0.01)T(20,60)      ## Prior for mean of birth distribution (1 would mean born just before day 0 [t=1], 10 would mean born just before day 9 [t=10])
  sd.bday ~ dunif(1,20)                 ## Prior for standard deviation of birth distribution
  skew.bday ~ dunif(1,10)
  
  
  nborn ~ dunif(100,50000)       ## prior for mean number of births (total across season)
  
  mu.days.on.colony ~ dunif(min.age, max.age)       ## Prior for mean of the number of days on the colony
  sd.days.on.colony  ~ dunif(1, 15)        ## Prior for standard deviation of days on colony
  
  
  
  #=================================
  # set the probability of being moulted as a function of age
  #===================================
  
  pmoult.int <- -12   # ~ dnorm(5,1)               # probability of being moult at birth
  pmoult.slope <- 0.6 # ~ dunif(-0.5,0)          # daily decrease in the probability of being moult (increase in probability of moulting)
  for(a in 1:ndays){   # loop past the last possible age
    logit(pmoult[a]) <- pmoult.int + pmoult.slope*(a-1)    # deterministic for now
  }
  #=================================
  # set the probability of survival as a function of age
  #===================================
  
  for(a in 1:(max.age+1)){
    logit(surv[a]) <- surv.int + surv.slope*(a-1)
  }
  
  #=======================================
  # set the scheduling distributions 
  #=======================================
  
  # BIRTH DISTRIBUTION
  #normal 
  # nconst_Pb <- 1 - ( pnorm(0,mu.bday, sd.bday) + (1 - pnorm((max.day.birth),mu.bday, sd.bday) ) )  # constant for truncating the births distribution
  # Pb[1] <- 0  # can't be born after 0 time has passed in the season
  # Pb[2:(max.day.birth+1)] <- (pnorm(1:(max.day.birth), mu.bday, sd.bday) - pnorm(0:(max.day.birth-1), mu.bday, sd.bday))/nconst_Pb   ## Birth schedule
  
  #skew normal using generic function
  #Pb<-rep(0,(max.day.birth+1))
  pb[1:(max.day.birth+1)] <- intPb(1:(max.day.birth+1), skew.bday, mu.bday, sd.bday)
  Pb[1:(max.day.birth+1)] <- pb[1:(max.day.birth+1)]/sum(pb[1:(max.day.birth+1)])
  
  # snconst_Pb <- 1 - (nim_pskew(0,skew.bday,mu.bday,sd.bday) + (1 - nim_pskew((max.day.birth),skew.bday,mu.bday,sd.bday) ) ) # constant for truncating the births distribution
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
      stageN_c[c,c+a-1,1] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:a]) * (1-pmoult[a])     # probability of being white after a days
      stageN_c[c,c+a-1,2] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:a])* pmoult[a]          # moulted
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
  # 
  # for (obs in 1:length(obsdays)) {
  #   whitemu[obs] <-  stageN[obsdays[obs],1]*Pobs + stageN[obsdays[obs],2]*Pobs*(1-PCMoult)
  #   
  #   sdwhite[obs] <- sqrt(stageN[obsdays[obs],1]* Pobs * (1 - Pobs) +
  #                          stageN[obsdays[obs],2]*Pobs*(1-PCMoult)* (1 - Pobs*(1 - PCMoult)))
  #                        
  #   moultmu[obs] <- stageN[obsdays[obs],2]*Pobs*PCMoult
  #                        
  #   sdmoult[obs] <- sqrt(stageN[obsdays[obs],2]*Pobs*PCMoult * (1 - Pobs * PCMoult))
  #                        
  #                        
  #    # normal approximation using binomial variance
  #    wObs[obs] ~ dnorm(whitemu[obs], sd = sdwhite[obs])
  #    mObs[obs] ~ dnorm(moultmu[obs], sd = sdmoult[obs])
  #    }
  
  
  # for(obs in 1:(length(obsdays))) {
  #  wObs[obs] ~ dbinom(Pobs, trunc(stageN[obsdays[obs],1]))   # data node
  #  mObs[obs] ~ dbinom(Pobs, trunc(stageN[obsdays[obs],2]))   # data node
  # }
  
  
  
  WW_MW[1:(ndays-1)] <- stageN[1:(ndays-1),1]+stageN[1:(ndays-1),2]*(1-PCMoult)
  MM[1:(ndays-1)] <- stageN[1:(ndays-1),2]*PCMoult  #add misclassification to mObs

  for(obs in 1:(length(obsdays))) {
    wObs[obs] ~ dbinom(Pobs, trunc(WW_MW[obsdays[obs]]))   # data node
    mObs[obs] ~ dbinom(Pobs, trunc(MM[obsdays[obs]]))   # data node
  }
  
})

# obsdays_combinations <- list(
#   #c(35,45, 55, 65, 75,85)
#   seq(30, 80, by = 12),  # Observation days starting at 10, going to 50, spaced 10 days apart
#   # seq(30, 70, by = 10),
#   #seq(40, 90, by = 12),
#   seq(50, 90, by = 10)
#   # seq(30, 100, by = 15),
#   # seq(40, 100, by = 15)
# )

# Define the list of different obsdays combinations
obsdays_combinations <- list(
  seq(10, 70, by = 10),
  seq(60, 100, by = 12),
  seq(40, 90, by = 20),
  seq(30, 90, by = 15)
)

# Loop over each scenario (obsdays combination)

# Parameters
# Parameters monitored
params <- c( "mu.bday", "sd.bday","nborn",'skew.bday',
             "mu.days.on.colony","sd.days.on.colony",
             #"surv.int","surv.slope","pwhite.int","pwhite.slope",
             #"pmass.birth.to.leave","cprob.remain","Pb",
             "stageN","cohort.size")


num_simulations <- 10  # Define the number of simulations
ndays <- 123
max.age <- 60
min.age <- 20
max.day.birth <- 100



# Define initial values for each simulation
inits_list <- list(     mu.bday = 45,
                        sd.bday = 15,
                        skew.bday=5,
                        mu.days.on.colony = 30, 
                        sd.days.on.colony= 10,  
                        nborn = 4000
                        #Pobs = 0.95
                        #cohort.size  = round(10000*(Pb)) 
)
# Placeholder for filename dataframe
filenames_df <- tibble()

# Loop through each scenario
for (scenario in 1:length(obsdays_combinations)) {
  this_obsdays <- obsdays_combinations[[scenario]]
  
  # Update constants to include the current observation days
  constants_list <- list(ndays = ndays,
                         max.age = max.age,
                         min.age = min.age,
                         max.day.birth = max.day.birth,
                         obsdays = this_obsdays)
  
  # Generate initial data for model setup
  initial_data <- generate_obs(ndays = 123,
                               min.age = 20,
                               max.age = 60,
                               max.day.birth = 100,
                               obsdays = this_obsdays,
                               mu.bday = 45,
                               sd.bday = 10,
                               skew.bday = 5,
                               nborn = 700,
                               mu.days.on.colony = 31.5,
                               sd.days.on.colony = 2,
                               surv.int = 4,
                               surv.slope = 0.1,
                               pmoult.int = -12,
                               pmoult.slope = 0.6,
                               Pobs = 0.95,
                               PCMoult = 0.84)
  
  set.seed(0)
  model <- nimbleModel(   # speed up compilation: can use calcualte=FALSE to disable calc of deterministics nodes and log-likelihood
    code = code,
    inits = inits_list,
    data = initial_data,
    constants = constants_list
  )
  model$calculate()
  conf <- configureMCMC(model, useConjugacy = FALSE) #conjugacy=false disables search for conjugate samplers
  Rmcmc <- buildMCMC(conf) # Create MCMC configuration based on the compiled model
  Cmodel <- compileNimble(model) # Compile the NIMBLE model, translating it to C++
  Cmcmc <- compileNimble(Rmcmc, project = model) # Compile the MCMC configuration and link it to the already compiled model
  
  # Loop to run simulations
  for (sim_index in 1:num_simulations) {
    data <- generate_obs(ndays = 123,
                         min.age = 20,
                         max.age = 60,
                         max.day.birth = 100,
                         obsdays = this_obsdays,
                         mu.bday = 45,
                         sd.bday = 10,
                         skew.bday = 5,
                         nborn = 700,
                         mu.days.on.colony = 31.5,
                         sd.days.on.colony = 2,
                         surv.int = 4,
                         surv.slope = 0.1,
                         pmoult.int = -12,
                         pmoult.slope = 0.6,
                         Pobs = 0.95,
                         PCMoult = 0.84)
    
    cat("Running simulation", sim_index, "for obsdays", this_obsdays, "\n")
    
    print(data)
    # Update model data
    Cmodel$setData(data)
    
    # Run the MCMC
    mcmc.out <- runMCMC(Cmcmc,
                        niter = 5000,
                        nburnin = 1000,
                        nchains = 1,
                        setSeed = FALSE,
                        samples = TRUE,
                        samplesAsCodaMCMC = FALSE)
  
    # Store the MCMC results
    # tempdf <- tibble(
    #   filename <- sprintf("C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/07 - GitHub Repositories/results/test_obsdays_%s_sim_%s.RDS", this_obsdays, sim_index, Sys.Date()),
    #   simulation = sim_index,
    #   scenario = scenario,
    #   obsdays = list(this_obsdays)
    # )
    # saveRDS(mcmc.out, file = tempdf$filename)
    # filenames_df <- rbind(filenames_df, tempdf)
    }
  }
} #end run simulation function


################################################################################
# Set up parallel 
cl <- makeCluster(detectCores()-3)  # Create a cluster using all available cores
clusterExport(cl, c("run_simulation", "obsdays_combinations","intPb","nim_dskew", "num_simulations", "code","generate_obs"))  # Export necessary functions and variables to the cluster
clusterEvalQ(cl, { library(nimble); library(tibble); })  # Load libraries 

# Execute simulations in parallel
results<-parLapply(cl, seq_along(obsdays_combinations), function(idx) {
  run_simulation(idx, obsdays_combinations)
}
)

# Clean up
stopCluster(cl)
################################################################################



# Example process_mcmc_results function
process_mcmc_results <- function(filename) {
  # Read the MCMC results from the CSV file
  mcmc_data <- readRDS(filename)
  
  # Perform any additional processing if necessary
  # For example, converting columns to appropriate types
  # or filtering out certain rows/columns
  
  return(mcmc_data)
}
# #for ()
# mcmc_areas(mcmc.out$samples,
#            pars = 'nborn',
#            prob = 0.9) + plot_title

# Placeholder for plotting data
plot_data <- tibble()

# Loop through each scenario in filenames_df
for (s in 1:nrow(filenames_df)) {
  filename <- filenames_df$filename[s]
  sim_index <- filenames_df$simulation[s]
  scenario <- filenames_df$scenario[s]
  
  # Read MCMC results from file
  mcmc_results <- process_mcmc_results(filename)
  
  # Extract parameters of interest
  params_of_interest <- c("mu.bday", "sd.bday", "nborn", "mu.days.on.colony")
  
  # Combine scenario and simulation info with parameter samples
  param_samples <- as.data.frame(mcmc_results) %>%
    select(all_of(params_of_interest)) %>%
    mutate(simulation = factor(sim_index), scenario = scenario)
  
  # Append to plot_data
  plot_data <- bind_rows(plot_data, param_samples)
}


# Convert scenario and simulation to factor for better plotting
plot_data <- plot_data %>%
  mutate(scenario = factor(scenario),
         simulation = factor(simulation))

# Calculate means for ablines
mean_values <- plot_data %>%
  summarise(across(all_of(params_of_interest), mean))

# Plotting using ggplot2
#Define the true values for each parameter
true_values <- list(
  "mu.bday" = 45, # Example true value for mu.bday
  "sd.bday" = 10, # Example true value for sd.bday
  "nborn" = 700, # Example true value for nborn
  "mu.days.on.colony" = 31.5 # Example true value for mu.days.on.colony
)
###############################################
#this plot plots all simulations for each scenario for each parameter
# Plotting using ggplot2
plots <- list()
library(ggplot2)
for (param in params_of_interest) {
  p <- ggplot(plot_data, aes_string(x = param, color = "simulation")) +
    geom_density() +
    facet_wrap(~ factor(scenario, levels = 1:2, labels = sapply(obsdays_combinations, paste, collapse = "-")), scales = "free") +
    labs(title = paste("Posterior Density of", param),
         x = param,
         y = "Density") +
    scale_color_discrete(guide = "none") +
    theme_minimal() +
    #theme_light() +  # Change to theme_light() for a white background
    theme(legend.position = "topright") +  # Place legend at the top right
    guides(color = guide_legend(title.position = "top", label.position = "top"))  # Set legend title
  # Add vertical line for the true value
  true_value <- true_values[[param]]
  p <- p + geom_vline(xintercept = true_value, linetype = "dashed", color = "red", size = 1)
  
  plots[[param]] <- p
  
}
# Save all plots
# for (plot_name in names(plots)) {
#   plot_filename <- sprintf("C:/Users/mg349/OneDrive - University of St Andrews/01 - Grey Seal Project/07 - GitHub Repositories/Hg_pup_model/sim_results/%s.png", plot_name)
#   ggsave(plot_filename, plot = plots[[plot_name]], width = 8, height = 6, units = "in", dpi = 300)
# }

# Display one of the plots as an example
print(plots[["mu.bday"]])
plots[["sd.bday"]]
plots[["nborn"]]
plots[["mu.days.on.colony"]]

###############################
#calculate bias: 

# Calculate the bias for each parameter
#bias_values <- mean_values - true_values
bias <- data.frame(parameter = names(mean_values),
                   mean_value = unlist(mean_values),
                   true_value = unlist(true_values),
                   bias = unlist(mean_values) - unlist(true_values))

library(HDInterval)
# Function to calculate HPD interval
calculate_hpd <- function(plot_data, param, prob = 0.95) {
  hpd_interval <-hdi(plot_data[[param]], prob = prob)
  return(hpd_interval)
}

# Print HPD intervals
print(hpd_interval)
hpd_intervals <- data.frame(
  mean_value = numeric(),
  lower_hpd = numeric(),
  upper_hpd = numeric(),
  stringsAsFactors = FALSE
)


for (param in params_of_interest) {
  mean_value <- mean(plot_data[[param]])
  hpd_interval <- calculate_hpd(plot_data, param)
  hpd_intervals[param, c("mean_value", "lower_hpd", "upper_hpd")] <- c(mean_value, hpd_interval[1], hpd_interval[2])
}

print(hpd_intervals)

###########################################################################
#BASE R PLOT OF DIFFFERENT SCENARIOS AND ALL SIMULATIONS
# Observation days mapping
# obsdays_mapping <- list(
#   `Scenario A` = c(30, 40, 50, 60, 70),
#   `Scenario B` = c(40, 50, 60, 70, 80)
# )
# observation_days <- unlist(lapply(obsdays_mapping, function(x) paste("(", paste(x, collapse = ", "), ")")))
# 
# # Legend labels with observation days
# legend_labels <- paste(names(obsdays_mapping), observation_days)
# # Define plot colors for scenarios
# color_blind_friendly <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
# custom_colors <- setNames(color_blind_friendly[seq_along(observation_days)], observation_days)

# Reset graphics parameters
#par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# for (param in params_of_interest) {
#   # Plot initialization
#   plot(NULL, xlim = range(plot_data[[param]]), ylim = c(0, 2), 
#        main = paste("Posterior Density of", param), 
#        xlab = param, ylab = "Density")
#   
#   # Plot density for each scenario and simulation
#   for (scenario in levels(plot_data$scenario)) {
#     scenario_index <- which(levels(plot_data$scenario) == scenario)
#     
#     for (simulation in levels(plot_data$simulation)) {
#       sim_data <- as.numeric(unlist(plot_data[plot_data$scenario == scenario & plot_data$simulation == simulation, param]))
#       
#       # Ensure sim_data is numeric and not empty
#       if (is.numeric(sim_data) && length(sim_data) > 0) {
#         density_values <- density(sim_data)
#         lines(density_values, col = custom_colors[scenario_index], lwd = 2)
#        }
#      }
#    }
#   
#   # Add true value line
#   abline(v = true_values[[param]], lty = 2, col = "black", lwd = 1.5)
#   
#   # Add legend
#   legend("topright", legend = legend_labels, col = custom_colors, lty = 1, lwd = 2, cex = 0.5, bty = "n")
# }
# 
