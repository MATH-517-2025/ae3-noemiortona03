library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Set seed
set.seed(953)

# Define the true regression function m(x) 
true_m <- function(x) {
  sin(1 / (x / 3 + 0.1))
}

#Function to generate data
# Generates n samples from Y = m(X) + epsilon, with X ~ Beta(alpha, beta) and epsilon ~ N(0, sigma_sq)
generate_data <- function(n, alpha, beta, sigma_sq = 1) {
  # Generate the covariate X from a Beta distribution
  x <- rbeta(n, shape1 = alpha, shape2 = beta)
  
  # Generate the error epsilon from a Normal distribution
  epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma_sq))
  
  # Calculate the response Y
  y <- true_m(x) + epsilon
  
  # Return the data as a dataframe
  return(data.frame(x = x, y = y))
}


#Function to find the optimal number of blocks 
find_optimal_N <- function(sim_data, min_block_size = 20) {
  n <- nrow(sim_data)
  N_max <- max(min(floor(n / min_block_size), 5), 1) # N_max respects min_block_size
  
  if (N_max < 1) return(1) # Edge case for very small n
  
  rss_values <- numeric(N_max) #create a vector to save rrs_values
  
  for (N in 1:N_max) {
    if (n / N < min_block_size && N > 1) {
      rss_values[N] <- NA
      next
    }
    
    # Create the 'block' column. 
    # If N=1, assign 1 to all rows --> every row has the exact same identifier 
    sim_data_N <- sim_data %>% 
      mutate(block = if (N > 1) cut(x, breaks = N, labels = FALSE, include.lowest = TRUE) else 1)
    
    # Calculate RSS for each block
    #For the current block (.), it attempts to fit a linear model (lm) with a fourth-degree polynomial 
    #If the fit fails (too few data points in the block), the program will store an error object (no stop)
    #The sum of the RRS of the blocks at the end of each iteration is saved in the vector 'rss_values'
    
    rss_by_block <- sim_data_N %>%
      group_by(block) %>%
      do({
        fit <- try(lm(y ~ poly(x, 4, raw = TRUE), data = .), silent = TRUE)
        if (inherits(fit, "try-error")) { data.frame(rss = NA) } else { data.frame(rss = sum(residuals(fit)^2)) }
      })
    rss_values[N] <- sum(rss_by_block$rss, na.rm = TRUE)
  }
  
  if(all(is.na(rss_values))) return(1)
  
  valid_indices <- which(!is.na(rss_values))
  last_valid_N_max <- max(valid_indices)
  
  if ((n - 5 * last_valid_N_max) <= 0) return(1) 
  rss_Nmax <- rss_values[last_valid_N_max]
  
  #Define a vector in which we save the Mallow's Cp for each N
  cp_values <- (rss_values[valid_indices] / (rss_Nmax/ (n - 5 * last_valid_N_max))) - (n - 10 * valid_indices)
  
  # Find the index of the minimum Cp among the valid ones
  best_index_in_valid_set <- which.min(cp_values)
  
  # Return the correct value of N
  N_star <- valid_indices[best_index_in_valid_set]
  
  return(N_star)
}

estimate_parameters <- function(sim_data, N) {
  n <- nrow(sim_data)
  
  #Estimate theta_22
  #Divide the data into N blocks based on the value of x
  sim_data_N <- sim_data %>% 
    mutate(block = if (N > 1) cut(x, breaks = N, labels = FALSE, include.lowest = TRUE) else 1)
  
  #Function to calculate the squared second derivative:
  get_m_double_prime_sq <- function(data_subset) {
    # Fit the 4th-degree polynomial model.
    fit <- try(lm(y ~ poly(x, 4, raw = TRUE), data = data_subset), silent = TRUE)
    
    # Error handling: if the fit fails, return NA
    if (inherits(fit, "try-error") || length(coef(fit)) < 5) {
      return(rep(NA, nrow(data_subset)))
    }
    
    # Extract the five estimated coefficients (beta_0, beta_1, ...)
    coeffs <- coef(fit)
    
    # Apply the formula for the second derivative of a quartic polynomial:
    # m''(x) = 2*beta_2 + 6*beta_3*x + 12*beta_4*x^2
    m_double_prime <- 2 * coeffs[3] + 6 * coeffs[4] * data_subset$x + 12 * coeffs[5] * data_subset$x^2
    
    # Return the square of this value (as required by the theta_22 formula)
    return(m_double_prime^2)
  }
  
  # Apply the helper function to each block of data
  # 'group_map' creates a list where each element is the result for one block
  # 'unlist' combines this list of vectors into a single large vector (inner sum of teh given formula)
  m_double_prime_sq_values <- unlist(
    sim_data_N %>% group_by(block) %>% group_map(~ get_m_double_prime_sq(.x))
  )
  
  # Calculate theta_22 by taking the mean of all the squared second derivative values
  # Theta_22 = (1/n) * sum(m''(xi)^2)
  theta_22_hat <- mean(m_double_prime_sq_values, na.rm = TRUE)
  
  #Estimate sigma_sq
  # We re-calculate the RSS for the chosen number of blocks, N.
  rss_by_block <- sim_data_N %>%
    group_by(block) %>%
    do({
      fit <- try(lm(y ~ poly(x, 4, raw = TRUE), data = .), silent = TRUE)
      if (inherits(fit, "try-error")) { data.frame(rss = NA) } else { data.frame(rss = sum(residuals(fit)^2)) }
    })
  rss_N <- sum(rss_by_block$rss, na.rm = TRUE)
  
  # Calculate sigma^2 using the RSS (as in the formula)
  if ((n - 5 * N) <= 0) {
    sigma_sq_hat <- NA
  } else {
    sigma_sq_hat <- rss_N / (n - 5 * N)
  }
  
  #Return all estimated parameters in a list
  return(list(
    sigma_sq_hat = sigma_sq_hat,
    theta_22_hat = theta_22_hat,
    RSS = rss_N
  ))
}
#-------------------------------------------------------
#How does h_AMISE behave when N grows?

# Function to calculate h_AMISE for each potential N
analyze_h_vs_N <- function(n, alpha, beta, sigma2) {
  sim_data <- generate_data(n = n, alpha = alpha, beta = beta, sigma_sq = sigma2) %>% arrange(x)
  N_max <- max(min(floor(n / 20), 5), 1)
  h_results <- vector("list", N_max)
  
  for (N in 1:N_max) {
    params <- estimate_parameters(sim_data, N) # Estimate parameters for the current N
    
    if(is.na(params$theta_22_hat) || params$theta_22_hat == 0 || is.na(params$sigma_sq_hat)) {
      h_amise_N <- NA
    } else {
      h_amise_N <- (n^(-1/5)) * ((35 * params$sigma_sq_hat) / abs(params$theta_22_hat))^(1/5)
    }
    h_results[[N]] <- data.frame(N = N, h_amise_at_N = h_amise_N)
  }
  return(bind_rows(h_results))
}

# Run the analysis for a large dataset (n=2000)
h_vs_N_data <- analyze_h_vs_N(n = 2000, alpha = 2, beta = 5, sigma2 = 1)

# Plot the results
plot_h_vs_N <- ggplot(h_vs_N_data, aes(x = N, y = h_amise_at_N)) +
  geom_line(linewidth = 1, color = "purple") +
  geom_point(size = 4, color = "purple") +
  labs(
    title = "Behavior of h_AMISE as N Grows",
    x = "Number of Blocks (N)",
    y = "Estimated h_AMISE"
  ) +
  scale_x_continuous(breaks = 1:max(h_vs_N_data$N, na.rm=TRUE)) +
  theme_minimal()

print(plot_h_vs_N)

#---------------------------------------------------------------
#Should N depend on n? Why?
        
# Vector of sample sizes to test
n_values <- seq(200, 3000, by = 200)
# Number of replications for each n to average out randomness
N_REPLICATIONS <- 200

# Fixed parameters for the Beta distribution
ALPHA <- 2
BETA <- 5

results_list <- list()

# OUTER LOOP: iterates over each value of n
for (n_i in n_values) {
  
  print(paste("Running", N_REPLICATIONS, "replications for n =", n_i))
  
  # Vector to store the results of each replication for the current n
  n_opt_replicates <- numeric(N_REPLICATIONS)
  
  # INNER LOOP: runs the simulation multiple times for the same n
  for (rep in 1:N_REPLICATIONS) {
    # 1. Generate data
    sim_data <- generate_data(n = n_i, alpha = ALPHA, beta = BETA) %>% arrange(x)
    
    # 2. Find the optimal N for this specific dataset
    N_opt <- find_optimal_N(sim_data)
    
    # 3. Store the result of this replication
    n_opt_replicates[rep] <- N_opt
  }
  
  # Calculate the average of the results for the current n
  average_n_opt <- mean(n_opt_replicates, na.rm = TRUE)
  
  # Add the averaged result for this n_i to the list
  results_list[[length(results_list) + 1]] <- data.frame(
    n = n_i,
    avg_N_opt = average_n_opt
  )
}

#Combine all results into a dataframe
results_averaged <- bind_rows(results_list)

# Plot: Average Optimal N vs. n
plot_N_vs_n <- ggplot(results_averaged, aes(x = n, y = avg_N_opt)) +
  geom_line(color="darkgreen", linewidth = 1.2) +
  geom_point(color="darkgreen", size = 3) +
  labs(title = "Average Optimal N vs. Sample Size (n)",
       subtitle = paste("Based on", N_REPLICATIONS, "replications per point"),
       x = "Sample Size (n)",
       y = "Average Optimal Number of Blocks (N_opt)") +
  theme_minimal()

print(plot_N_vs_n) 

#-------------------------------------------------------------------
#What happens when the number of observations varies a lot between different regions in the support of X? 
#How is this linked to the parameters of the Beta distribution?
  
beta_params <- list(
  "Uniform (a=1, b=1)" = c(1, 1), 
  "Unimodal (a=5, b=5)" = c(5, 5), 
  "U-shaped (a=0.5, b=0.5)" = c(0.5, 0.5), 
  "Asymmetric (a=5, b=2)" = c(5, 2)
)
n_fixed <- 500
results_dist_list <- list()

for (name in names(beta_params)) {
  print(paste("  Simulating for distribution:", name))
  params_beta <- beta_params[[name]]
  for (rep in 1:N_REPLICATIONS) {
    sim_data <- generate_data(n = n_fixed, alpha = params_beta[1], beta = params_beta[2], sigma_sq = 1) %>% arrange(x)
    N_star <- find_optimal_N(sim_data)
    if (is.na(N_star)) {
      h_amise <- NA
    } else {
      params <- estimate_parameters(sim_data, N_star)
      if(is.na(params$theta_22_hat) || params$theta_22_hat == 0 || is.na(params$sigma_sq_hat)) {
        h_amise <- NA
      } else {
        h_amise <- (n_fixed^(-1/5)) * ((35 * params$sigma_sq_hat) / abs(params$theta_22_hat))^(1/5)
      }
    }
    results_dist_list[[length(results_dist_list) + 1]] <- data.frame(distribution = name, h_amise = h_amise)
  }
}

results_dist_df <- bind_rows(results_dist_list) %>% filter(!is.na(h_amise))

results_dist_df$distribution <- factor(results_dist_df$distribution, levels = names(beta_params))
plot_h_vs_dist <- ggplot(results_dist_df, aes(x = distribution, y = h_amise, fill = distribution)) +
  geom_boxplot() +
  labs(title = paste("Impact of n° of observations in different regions on h_AMISE"), 
       x = "Distribution of X", 
       y = "Estimated Optimal Bandwidth") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 10, hjust = 1))
print(plot_h_vs_dist)
