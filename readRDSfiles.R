###
### Load in MCMC results from any of the scripts that compiles each time (HgNimble_recode9_sensitivity3.R)
###

###
library(tibble)
library(dplyr)
library(coda)
library(MCMCvis)
library(ggplot2)
library(purrr)
# obsdays_combinations <- list(
#   seq(10, 70, by = 10),
#   seq(60, 100, by = 12),
#   seq(40, 90, by = 20),
#   seq(30, 90, by = 15)
# )

results_dir <- "results/"
date <-  Sys.Date() - 1
current_date <- Sys.Date()
yesterday <- Sys.Date() - 1

num_simulations<- 2
# List all RDS files in the directory that include either today's or yesterday's date
#rds_files_today <- list.files(results_dir, pattern = paste0("_", current_date, "\\.rds$"), full.names = TRUE)
# rds_files_yesterday <- list.files(results_dir, pattern = paste0("_", yesterday, "\\.rds$"), full.names = TRUE)
# rds_files <- c(rds_files_today, rds_files_yesterday)
rds_files<- list.files(results_dir, pattern = paste0("_", date, "\\.rds$"), full.names = TRUE)
#extract associated data for each simulation/scenario
data_files<- list.files(results_dir, pattern = paste0(".*_", date, ".*_data_\\.rds$"), full.names = TRUE)

# Function to extract scenario and simulation from the filename
extract_scenario_simulation <- function(filename) {
  parts <- unlist(strsplit(basename(filename), "_"))
  scenario <- paste(parts[2:(length(parts)-4)], collapse = "_")
  simulation <- as.numeric(parts[length(parts)-2])
  return(list(scenario = scenario, simulation = simulation))
}

#Function to read data and add scenario and simulation info
process_file <- function(file) {
  info <- extract_scenario_simulation(file)
  data <- readRDS(file)
  tibble(
    scenario = info$scenario,
    simulation = info$simulation,
    data.wObs = data$wObs,
    data.mObs = data$mObs
  )
}

# Process all files and combine into a single data frame
data_df <- map_dfr(data_files, process_file)
# this next line of code is because im dumb and didn't name things properly for this test simulation!
# Remove "22" and everything before it in scenario
data_df$scenario <- sub(".*?22", "", data_df$scenario)

# Aggregate data.wObs and data.mObs by scenario and simulation
aggregated_data <- data_df %>%
  group_by(scenario, simulation) %>%
  summarise(wObs = list(data.wObs), mObs = list(data.mObs), .groups = 'drop')

# Unique scenarios from data_df
unique_scenarios <- unique(data_df$scenario)
# Initialize an empty list to store obsdays combinations
obsdays_combinations <- vector("list", length(unique_scenarios))

# Loop through unique scenarios and construct obsdays for each scenario
for (i in seq_along(unique_scenarios)) {
  # Split the scenario string by underscores and convert to numeric
  obsdays <- as.numeric(strsplit(unique_scenarios[i], "_")[[1]])
  # Store obsdays in the corresponding position in obsdays_combinations
  obsdays_combinations[[i]] <- obsdays
}

# Generate custom labels for the scenarios
scenario_labels <- sapply(obsdays_combinations, function(x) paste(x, collapse = "-"))
# # Initialize a list to store all loaded data frames

# Placeholder for filenames_df dataframe
filenames_df <- tibble(
  filename = rds_files,
  simulation = rep(1:num_simulations, each = length(obsdays_combinations)),
  scenario = rep(scenario_labels, times = num_simulations)
)


##############################################################################
###
###  BAYESPLOT STUFFS. 
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
  
  # Extract stageN samples and store them in a nested list for use to do ppd later
  if (!is.list(stageN_samples[[scenario]])) {
    stageN_samples[[scenario]] <- list()
  }
  # 
  # # Convert mcmc_sample to matrix
  mcmc_matrix <- as.matrix(mcmc_sample)
  # # Extract stageN samples
  stageN_samples[[scenario]] <- mcmc_matrix[, grepl("stageN", colnames(mcmc_matrix))]
  
}



# Initialize a list to store mean results
mean_results <- list()

# Loop through each combination of obsdays_combinations
for (i in seq_along(obsdays_combinations)) {
  # Get the current combination of days
  days <- obsdays_combinations[[i]]
  
  # Generate column names for stageN[, 1] and stageN[, 2] from the current days
  column_names_1 <- paste0("stageN[", days, ", 1]")
  column_names_2 <- paste0("stageN[", days, ", 2]")
  
  # Find column indices corresponding to column_names in stageN_samples
  col_indices_1 <- which(colnames(stageN_samples) %in% column_names_1)
  col_indices_2 <- which(colnames(stageN_samples) %in% column_names_2)
  
  # Subset stageN_samples based on column indices for stageN[, 1]
  stageN_1_mean <- colMeans(stageN_samples[, col_indices_1])
  
  # Subset stageN_samples based on column indices for stageN[, 2]
  stageN_2_mean <- colMeans(stageN_samples[, col_indices_2])

  # Store results in mean_results
  mean_results[[i]] <- list(
    scenario = days,
    mean_stageN_1 = stageN_1_mean,
    mean_stageN_2 = stageN_2_mean 
  )
}

# Print mean_results in the specified format
for (i in seq_along(mean_results)) {
  cat("Scenario:", scenario, "\n")
  cat("mean white\n")
  print(mean_results[[i]]$mean_stageN_1)
  cat("mean moult\n")
  print(mean_results[[i]]$mean_stageN_2)
}


########################################
###
### posterior predictive check
###
########################################


## use mean_stageN_1 and mean_stageN_2 as simulated data to test against observed





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
  summaries[[scenario]] <- MCMCvis::MCMCsummary(combined_mcmc_scenario, HPD = TRUE)
}

# Print summaries for each scenario
for (scenario in names(summaries)) {
  cat("Scenario:", scenario, "\n")
  print(summaries[[scenario]])
}















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
#### 
#Plotting using ggplot2
#Define the true values for each parameter
true_values <- list(
  "mu.bday" = 45, # Example true value for mu.bday
  "sd.bday" = 10, # Example true value for sd.bday
  "skew.bday"=5,
  "nborn" = 700, # Example true value for nborn
  "mu.days.on.colony" = 31.5, # Example true value for mu.days.on.colony
  "sd.days.on.colony"=2
)

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
  #pdf_file1 <- paste0(results_dir, param, "_", Sys.Date(), "_sensitivity3_50_mean_hpd_plot.pdf")
  #ggsave(pdf_file1, plot = p, device = "pdf", width = 10, height = 6)
  
  }
#plot(p)

# Create point plots for each parameter
# for (param in params_of_interest) {
#   print(plot_mean_values(param))
# }
# List of parameters to plot
params_of_interest <- unique(summary_stats_df$parameter)

# # Create bar plots for each parameter
# for (param in params_of_interest) {
#   print(plot_hpd_intervals(param))
# }
###
### PLOTTING

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
plots <- list()
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
 # pdf_file <- paste0("results/", param,"_", Sys.Date(),"_sensitivity3","_50","_density_plots.pdf")
 # ggsave(pdf_file, plot = p, device = "pdf", width = 10, height = 6)
  plots[[param]] <- p
  }
# # Display the plots
for (param in params_of_interest) {
  print(plots[[param]])
 }

# # Display one of the plots as an example
 print(plots[["mu.bday"]])
# plots[["sd.bday"]]
# plots[["nborn"]]
# plots[["mu.days.on.colony"]]

