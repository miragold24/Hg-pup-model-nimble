###
### Load in MCMC results from any of the scripts that compiles each time (HgNimble_recode9_sensitivity3.R)
###

###
library(tibble)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(MCMCvis)
obsdays_combinations <- list(
  seq(10, 70, by = 10),
  seq(60, 100, by = 12),
  seq(40, 90, by = 20),
  seq(30, 90, by = 15)
)

results_dir <- "results/"
current_date <- Sys.Date()
yesterday <- Sys.Date() - 1

num_simulations<- 50
# List all RDS files in the directory that include either today's or yesterday's date
rds_files_today <- list.files(results_dir, pattern = paste0("_", current_date, "\\.rds$"), full.names = TRUE)
rds_files_yesterday <- list.files(results_dir, pattern = paste0("_", yesterday, "\\.rds$"), full.names = TRUE)
rds_files <- c(rds_files_today, rds_files_yesterday)

# Generate custom labels for the scenarios
scenario_labels <- sapply(obsdays_combinations, function(x) paste(x, collapse = "-"))
# # Initialize a list to store all loaded data frames



# Placeholder for filenames_df dataframe
filenames_df <- tibble(
  filename = rds_files,
  simulation = rep(1:50, each = length(obsdays_combinations)),
  scenario = rep(scenario_labels, times = 50)
)

##############################################################################
###
###  BAYESPLOT/MCMCVIS STUFFS. 
###
##############################################################################
# Function to load an RDS file and convert to MCMC object
load_mcmc_object <- function(filename) {
  # Read the RDS file
  mcmc_samples <- readRDS(filename)
  
 
  #if (!inherits(mcmc_samples, "mcmc")) {
    mcmc_samples <- coda::mcmc(mcmc_samples)
 # }
  
  return(mcmc_samples)
}


# Initialize an empty list to store MCMC objects
mcmc_objects <- list()
summaries<-list()
# Loop through each scenario in filenames_df
for (s in 1:nrow(filenames_df)) {
  filename <- filenames_df$filename[s]
  sim_index <- filenames_df$simulation[s]
  scenario <- filenames_df$scenario[s]
  
  # Load the MCMC object
  mcmc_sample <- load_mcmc_object(filename)
  
  # Create a nested list structure to store the MCMC objects
  if (!is.list(mcmc_objects[[scenario]])) {
    mcmc_objects[[scenario]] <- list()
  }
  mcmc_objects[[scenario]][[paste0("sim_", sim_index)]] <- mcmc_sample
  

}

# Initialize a list to store summaries
summaries <- list()

# Extract parameters of interest
params_of_interest <- c("mu.bday", "sd.bday", "skew.bday", "nborn", "mu.days.on.colony", "sd.days.on.colony")

# Loop through each scenario in mcmc_objects
for (scenario in names(mcmc_objects)) {
  # Initialize a list to store filtered MCMC samples
  filtered_samples <- vector("list", length = length(mcmc_objects[[scenario]]))
  
  # Loop through each simulation in the current scenario
  for (sim in seq_along(mcmc_objects[[scenario]])) {
    # Filter MCMC samples for params_of_interest
    filtered_samples[[sim]] <- mcmc_objects[[scenario]][[sim]][, params_of_interest]
  }
  
  # Combine all simulations for the current scenario
  combined_mcmc_scenario <- mcmc(do.call(rbind, filtered_samples))
  
  # Compute the summary for the combined MCMC sample
  summaries[[scenario]] <- MCMCsummary(combined_mcmc_scenario,   
                                                HPD = TRUE )
  
  # # Generate density plot using MCMCvis
  # plots[[scenario]] <- MCMCtrace(filtered_samples, 
  #                                params = params_of_interest,
  #                                type="both",
  #                                gvals="true_values")

}

# Print summaries for each scenario
for (scenario in names(summaries)) {
  cat("Scenario:", scenario, "\n")
  print(summaries[[scenario]])
  
  print(plots[[scenario]])
}



library(bayesplot)
#mcmc_pairs(combined_mcmc$`60-72-84-96`[, params_of_interest])


##############################################################################
###
###  GGPLOT STUFFS. (Plot simulations)
###
##############################################################################
# process_mcmc_results function
process_mcmc_results <- function(filename) {
  # Read the MCMC results from the CSV file
  mcmc_data <- readRDS(filename)
  
  return(mcmc_data)
}
# Placeholder for plotting data
plot_data <- tibble()

# Loop through each scenario in filenames_df
for (s in 1:nrow(filenames_df)) {
  filename <- filenames_df$filename[s]
  sim_index <- filenames_df$simulation[s]
  scenario <- filenames_df$scenario[s]
  
  # Read MCMC results from file
  mcmc_results <- process_mcmc_results(filename)
  mcmc_objects <- lapply(filename, load_mcmc_object)
  # Optionally, you can name the list elements for easier reference
  names(mcmc_objects) <- filename
  # Extract parameters of interest
  params_of_interest <- c("mu.bday", "sd.bday","skew.bday", "nborn", "mu.days.on.colony","sd.days.on.colony")
  
  # Combine scenario and simulation info with parameter samples
  param_samples <- as.data.frame(mcmc_results) %>%
    select(all_of(params_of_interest)) %>%
    mutate(simulation = factor(sim_index), scenario)
  
  # Append to plot_data
  plot_data <- bind_rows(plot_data, param_samples)
}


# Convert scenario and simulation to factor for better plotting
plot_data <- plot_data %>%
  mutate(scenario = factor(scenario),
         simulation = factor(simulation))



###
### CALCULATE BIAS AND PRECISION
###

library(HDInterval)
# Function to calculate HPD interval
calculate_hpd <- function(plot_data, param, prob = 0.95) {
  hpd_interval <-hdi(plot_data[[param]], prob = prob)
  return(hpd_interval)
}


summary_stats <- list()
hpd_intervals <- list()

# Iterate over scenarios
for (scenario in levels(plot_data$scenario)) {
  # Subset data for the current scenario
  scenario_data <- plot_data[plot_data$scenario == scenario, ]
  
  # Iterate over parameters of interest
  for (param in params_of_interest) {
    # Ensure param exists as a column in scenario_data
    if (param %in% colnames(scenario_data)) {
      # Calculate mean value
      mean_value <- mean(scenario_data[[param]])
      
      # Calculate HPD interval
      hpd_interval <- calculate_hpd(scenario_data, param)
      
      # Calculate bias (if true_values are available)
      # Example: Replace true_values[[param]] with your actual structure for true values
      bias <- mean_value - true_values[[param]]  # Adjust as per your true_values structure
      
      # Store results in lists
      summary_stats[[paste(scenario, param, sep = "_")]] <- c(mean_value, hpd_interval[1], hpd_interval[2], bias)
      hpd_intervals[[paste(scenario, param, sep = "_")]] <- hpd_interval
    } else {
      # Handle case where param does not exist in scenario_data
      warning(paste("Parameter", param, "not found in scenario", scenario, ". Skipping..."))
    }
  }
}
# Convert summary_stats to data frame
summary_stats_df <- tibble(
  parameter = gsub(".*_", "", names(summary_stats)),
  scenario = gsub("_.*", "", names(summary_stats)),
  mean_value = sapply(summary_stats, "[[", 1),
  lower_hpd = sapply(summary_stats, "[[", 2),
  upper_hpd = sapply(summary_stats, "[[", 3),
  bias = sapply(summary_stats, "[[", 4)
)



# Calculate summary measures across scenarios
summary_measures <- summary_stats_df %>%
  group_by(scenario) %>%
  summarize(
    mean_bias = mean(bias, na.rm = TRUE),
    min_hpd_width = min(upper_hpd - lower_hpd, na.rm = TRUE)
  ) %>%
  arrange(mean_bias, min_hpd_width)

# Print summary measures for comparison
print(summary_measures)



###########################################################################

#Function to create a point plot with error bars for mean values for a given parameter
for (param in params_of_interest) {
  param_data <- summary_stats_df[summary_stats_df$parameter == param, ]
  
  p<- ggplot(param_data, aes(x = scenario, y = mean_value, color = scenario)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower_hpd, ymax = upper_hpd), width = 0.2) +
    labs(x = "Scenario", y = "Mean Value", title = paste("Mean Values with HPD Intervals for", param), color = "Scenario") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # Save the plot
  pdf_file1 <- paste0(results_dir, param, "_", Sys.Date(), "_sensitivity3_50_mean_hpd_plot.pdf")
  ggsave(pdf_file1, plot = p, device = "pdf", width = 10, height = 6)
  
  }


# Create point plots for each parameter
for (param in params_of_interest) {
  print(plot_mean_values(param))
}
# List of parameters to plot
params_of_interest <- unique(summary_stats_df$parameter)

# # Create bar plots for each parameter
# for (param in params_of_interest) {
#   print(plot_hpd_intervals(param))
# }
###
### PLOTTING
###
# Plotting using ggplot2
#Define the true values for each parameter
true_values <- list(
  "mu.bday" = 45, # Example true value for mu.bday
  "sd.bday" = 10, # Example true value for sd.bday
  "skew.bday"=5,
  "nborn" = 700, # Example true value for nborn
  "mu.days.on.colony" = 31.5, # Example true value for mu.days.on.colony
  "sd.days.on.colony"=2
)
###############################################
#this plot plots all simulations for each scenario for each parameter
# Plotting using ggplot2
library(tidyverse)
library(bayesplot)
library(gridExtra)

# # Function to generate and display posterior density plots using ggplot2
# density_plots <- function(params, data) {
#   plots_list <- list()
# 
#   for (scenario in unique(data$scenario)) {
#     plots <- lapply(params, function(param) {
#       scenario_data <- data %>% filter(scenario == scenario, parameter == param)
# 
#       # Create density plot
#       density_plot <- ggplot(scenario_data, aes(x = value, fill = simulation)) +
#         geom_density(alpha = 0.5) +
#         labs(title = paste("Posterior Density for", param, "in Scenario", scenario),
#              x = param, y = "Density") +
#         theme_minimal() +
#         theme(legend.position = "none")
# 
#       return(density_plot)
#     })
# 
#     plots_list[[scenario]] <- do.call(grid.arrange, c(plots, ncol = 3)) #length(params)))
#   }
# 
#   return(plots_list)
# }
# 
# # Specify parameters of interest
# params_of_interest <- c("mu.bday", "sd.bday", "skew.bday", "nborn", "mu.days.on.colony", "sd.days.on.colony")
# 
# # Generate and display density plots for each scenario
# #scenario_plots <- density_plots(params_of_interest, plot_data_long)
# 
# # Display plots for each scenario
# for (i in seq_along(scenario_plots)) {
#   print(scenario_plots[[i]])
# }

# Combine scenario and simulation info with parameter samples

library(ggplot2)
#plots <- list()
for (param in params_of_interest) {
  p <- ggplot(plot_data, aes_string(x = param , color = "simulation")) +
    geom_density() +
    facet_wrap(~ scenario, scales = "free") +
    labs(title = paste("Posterior Density of", param),
         x = param,
         y = "Density") +
    scale_color_discrete(guide = "none") +
    theme_minimal()
  
  # Add vertical line for the true value
  true_value <- true_values[[param]]
  p <- p + geom_vline(xintercept = true_value, linetype = "dashed", color = "red", linewidth = 1)
  
  # Save the plot as PDF
  pdf_file <- paste0("results/", param,"_", Sys.Date(),"_sensitivity3","_50","_density_plots.pdf")
  ggsave(pdf_file, plot = p, device = "pdf", width = 10, height = 6)
  #plots[[param]] <- p
  }
# # Display the plots
# for (param in params_of_interest) {
#   print(plots[[param]])
#  }

# # Display one of the plots as an example
 print(plots[["mu.bday"]])
# plots[["sd.bday"]]
# plots[["nborn"]]
# plots[["mu.days.on.colony"]]

