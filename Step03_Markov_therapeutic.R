# =============================================================================
# Step03_Markov_therapeutic.R
# Purpose : (1) Same two-state Markov model for THERAPEUTIC AMU with 3 cut points;
#           (2) MCMC diagnostics; 
#           (3) daily initiation probability (p4);
#           (4) posterior-predictive proportion of days under treatment per phase (p5)
# Input   : thera_logistic_data, thera_poisson_data (Step00);
#           model/Therapeutic_AMU_three_CPs.txt
# Output  : output/Supp_Figure_S12_Trace_plot_therapeutic.tiff
#           output/Supp_Figure_S14_Pairs_plot_therapeutic.tiff
#           JAGS summary -> log; p4, p5 are saved in Step04 (Figure 5)
# Packages: loaded in Master_Run_All.R | Authors/license/dates: see README.md
# =============================================================================

################################Data preparation########################################
# Import data (objects created in Step00_data_process.R)
logistic_data <- thera_logistic_data
poisson_data <- thera_poisson_data

# variable description
usage_today <- as.numeric(logistic_data$Drug_usage_final)
farm_id <- as.numeric(logistic_data$farm_id)
T <- as.numeric(logistic_data$Admin_days_cycle)              # see README issue 3 (column name)
duration <- as.numeric(poisson_data$duration)
farm_event <- as.numeric(poisson_data$farm_id)     # see README issue 3 (column name)

##########Therapeutic AMU Bayesian Three state Markov JAGS Model Fitting ####################
#### List of data - function ###
data_jags <- list(
  # Logistic
  N_no_use_days = length(usage_today),
  usage_today = usage_today,
  farm_id = farm_id, ##farm level
  T = T,
  CP1=27,   # HARD-CODED cut point 1 (day of age): estimated from the GAMM
  CP2=53,   # HARD-CODED cut point 2 (day of age): estimated from the GAMM
  CP3=83,   # HARD-CODED cut point 3 (day of age): estimated from the GAMM
  
  # Poisson
  N_usage_events = length(duration),
  duration = duration,
  farm_event = farm_event, ##farm level

  #farm and cycle
  N_farms = 60,   # HARD-CODED: number of study farms
  SD = 2          
)

### Parameters to save - Function ###
parameters <- c(
  "beta0",
  "beta_interval0",
  "beta_interval1",
  "beta_interval2",
  "beta_interval3",
  "logistic_farm_sigma",
  "alpha0",
  "poisson_farm_sigma"
)

### Defining initial value of rjags function ###
inits <- function(){
  list(
    beta0 = rnorm(1,0,0.1), 
    beta_interval0 = rnorm(1,0,0.1),
    beta_interval1 = rnorm(1,0,0.1),
    beta_interval2 = rnorm(1,0,0.1),
    beta_interval3 = rnorm(1,0,0.1),
    logistic_farm_tau = runif(1, min = 0.03, max = 4),
    alpha0 = rnorm(1,0,0.1),
    poisson_farm_tau = runif(1, min = 0.03, max = 4)
  )
}

model.file <- "model/Therapeutic_AMU_three_CPs.txt"

### Run JAGS and save outputs ###
registerDoParallel(cores = 3)
runjags.options(silent.jags = FALSE, modules = "glm")

jags_fit <- jags.parallel(
  data = data_jags,
  inits = inits,
  parameters.to.save = parameters,
  model.file = model.file,
  n.chains = 4,
  n.iter = 20000,
  n.burnin = 5000,
  n.thin = 2,
  n.cluster = 3
)

##############################outputs##########################################################
# Extract MCMC samples
sims_array <- jags_fit$BUGSoutput$sims.array  # Iterations x Chains x Parameters
n_chains <- dim(sims_array)[2]  # number of chains

selected_params <- c("beta0", 
                     "beta_interval0","beta_interval1","beta_interval2","beta_interval3",
                     "logistic_farm_sigma", 
                     "alpha0", 
                     "poisson_farm_sigma")

# Convert to mcmc.list for trace and pairs plots
mcmc_combined <- as.mcmc.list(lapply(1:n_chains, function(chain) {
  m <- sims_array[, chain, selected_params, drop = FALSE]
  as.mcmc(matrix(m, nrow = dim(m)[1], ncol = length(selected_params),
                 dimnames = list(NULL, selected_params)))
}))

### Trace Plot ###
custom_colors <- c("#3D9F3C", "#9ED17B", "#367DB0", "#9DC7DD")
tiff("output/Supp_Figure_S12_Trace_plot_therapeutic.tiff", width = 15, height = 10, units = "in", res = 300, compression = "lzw")
trace_plot <- mcmc_trace(mcmc_combined, pars = selected_params) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 16)
print(trace_plot)
dev.off()

### Pairs Plot ###
if (length(selected_params) > 1) {
  tiff("output/Supp_Figure_S14_Pairs_plot_therapeutic.tiff", width = 18, height = 14, units = "in", res = 300, compression = "lzw")
  pairs_plot <- mcmc_pairs(
    mcmc_combined,
    pars = selected_params,
    off_diag_args = list(size = 1, alpha = 0.5)
  )
  print(pairs_plot)
  dev.off()
} else {
  warning("Not enough parameters for pairs plot.")
}

### Daily probability of starting a therapeutic treatment (mode, 95% HDI) ###
CP1 <- 27
CP2 <- 53
CP3 <- 83
max_days <- 161

# Extract posterior samples from JAGS model (all samples)
jags_samples <- jags_fit$BUGSoutput$sims.list
beta0_samples <- jags_samples$beta0
beta_interval0_samples <- jags_samples$beta_interval0
beta_interval1_samples <- jags_samples$beta_interval1
beta_interval2_samples <- jags_samples$beta_interval2
beta_interval3_samples <- jags_samples$beta_interval3

# Function to calculate mode
calculate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

# Initialize storage
results_list <- list()

# Calculate probability distribution for each day
for (day in 0:max_days) {
  # Calculate probability for each posterior iteration
  logit_p <- beta0_samples + 
    beta_interval0_samples * day + 
    beta_interval1_samples * pmax(day - CP1, 0) + 
    beta_interval2_samples * pmax(day - CP2, 0) +
    beta_interval3_samples * pmax(day - CP3, 0) 
  
  probs <- plogis(logit_p)
  
  # Calculate mode and 95% HDI
  mode_val <- calculate_mode(probs)
  hdi <- hdi(probs, credMass = 0.95)
  
  # Store results
  results_list[[day + 1]] <- data.frame(
    day = day,
    mode = mode_val,
    hdi_lower = hdi[1],
    hdi_upper = hdi[2]
  )
}

# Combine results
results_df <- do.call(rbind, results_list)

p4 <- ggplot(results_df, aes(x = day)) +
  geom_ribbon(aes(ymin = hdi_lower, ymax = hdi_upper), 
              fill = "#A2D94D", alpha = 0.2) +
  geom_line(aes(y = mode), color = "#0073B1", linewidth = 1.2) +
  geom_line(aes(y = hdi_lower), color = "#A2D94D", 
            linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = hdi_upper), color = "#A2D94D", 
            linetype = "dashed", linewidth = 0.8) +
  # vertical lines at the three cut points
  geom_vline(aes(xintercept = CP1, linetype = paste0("Day ", CP1)), 
             color = "darkred", linewidth = 0.8) +
  geom_vline(aes(xintercept = CP2, linetype = paste0("Day ", CP2)), 
             color = "darkblue", linewidth = 0.8) +
  geom_vline(aes(xintercept = CP3, linetype = paste0("Day ", CP3)), 
             color = "black", linewidth = 0.8) +
  # three linetype values
  scale_linetype_manual(
    name = "Cutpoints",
    values = c("dashed", "dotted", "longdash")
  ) +
  labs(
    x = "Days of age",
    y = "Probability of farm starting a therapeutic treatment"
  ) +
  theme_classic() +
  theme(
    plot.title        = element_text(size = 16, hjust = 0.5),
    axis.title        = element_text(size = 16),
    axis.text         = element_text(size = 16),
    legend.title      = element_text(size = 16),   # "Cutpoints"
    legend.text       = element_text(size = 16),   # "Day 27" / "Day 53" / "Day 83"
    legend.position   = c(0.85, 0.78),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size   = unit(1.2, "lines")         # longer legend keys
  ) +
  scale_y_continuous(limits = c(0, 0.2)) +
  scale_x_continuous(limits = c(0, 161)) +
  annotate("text", x = 137, y = 0.12, 
           label = "Blue line: mode\nShaded area: 95% HDI", 
           hjust = 0.5, size = 16 / .pt, color = "black")

p4

### Proportion of days under therapeutic treatment per phase based on the posterior-predictive simulation ###
### Posterior-predictive simulation ###
# read the production cycles for different farms 
farm_cycles <- cycle$Days_cycle   # one value per farm (created in Step00)
n_days_max  <- max(farm_cycles)   # longest production cycle across farms
chosen_iter <- sample(seq_len(length(jags_samples$beta0)), 5000) # HARD-CODED: 5000 posterior draws for the posterior-predictive simulation (this can be modified if needed)
n_iter <- 5000
n_farms <- 60

# Initialize data structures for posterior predictions
post_pred_durations_t <- array(list(), dim = c(n_iter, n_farms))
post_pred_usage_t <- array(NA, dim = c(n_iter, n_farms, n_days_max + 1))

# Loop through each posterior iteration
for (iter in 1:n_iter) {
  beta0 <- jags_samples$beta0[chosen_iter[iter]]
  beta_interval0 <- jags_samples$beta_interval0[chosen_iter[iter]]
  beta_interval1 <- jags_samples$beta_interval1[chosen_iter[iter]]
  beta_interval2 <- jags_samples$beta_interval2[chosen_iter[iter]]
  beta_interval3 <- jags_samples$beta_interval3[chosen_iter[iter]]
  alpha0 <- jags_samples$alpha0[chosen_iter[iter]]
  
  for (f in 1:n_farms) {
    n_days_farm <- farm_cycles[f]
    remaining_duration <- 0
    durations_farm <- integer(0)
    
    for (day in 0:n_days_farm) {
      array_index <- day + 1
      
      if (remaining_duration > 0) {
        post_pred_usage_t[iter, f, array_index] <- 1
        remaining_duration <- remaining_duration - 1
        
      } else {
        logit_p <- beta0 + beta_interval0 * day + 
          beta_interval1 * pmax(day - CP1, 0) + 
          beta_interval2 * pmax(day - CP2, 0) +
          beta_interval3 * pmax(day - CP3, 0) 
        
        start_prob <- plogis(logit_p)
        start_usage <- rbinom(n = 1, size = 1, prob = start_prob)
        post_pred_usage_t[iter, f, array_index] <- start_usage
        
        if (start_usage == 1) {
          log_lambda <- alpha0 
          lambda <- exp(log_lambda)
          duration <- rtpois(1, lambda = lambda, a = 0)
          duration <- min(duration, n_days_farm - day + 1)
          remaining_duration <- duration - 1
          durations_farm <- c(durations_farm, duration)
        }
      }
    }  
    if (n_days_farm < n_days_max) {
      post_pred_usage_t[iter, f, (n_days_farm + 2):(n_days_max + 1)] <- NA
    }
    post_pred_durations_t[[iter, f]] <- durations_farm
  }
}

### Calculate proportion of days with treatment in each phase based on the posterior-predictive simulation ###
# Initialize storage for proportions per iteration
prop_phase1_iter <- numeric(n_iter)
prop_phase2_iter <- numeric(n_iter)
prop_phase3_iter <- numeric(n_iter)
prop_phase4_iter <- numeric(n_iter)

# Calculate proportions for each iteration
for (iter in 1:n_iter) {
  # Storage for farm-level proportions in this iteration
  prop_phase1_farms <- numeric(n_farms)
  prop_phase2_farms <- numeric(n_farms)
  prop_phase3_farms <- numeric(n_farms)
  prop_phase4_farms <- numeric(n_farms)
  
  for (f in 1:n_farms) {
    n_days_farm <- farm_cycles[f]
    
    # Phase 1
    phase1_end <- min(CP1, n_days_farm)
    if (phase1_end >= 0) {
      phase1_days <- phase1_end + 1  # number of days in phase 1
      phase1_treatment <- sum(post_pred_usage_t[iter, f, 1:(phase1_end + 1)], na.rm = TRUE)
      prop_phase1_farms[f] <- phase1_treatment / phase1_days
    } else {
      prop_phase1_farms[f] <- NA
    }
    # Phase 2
    if (n_days_farm > CP1) {
      phase2_start <- CP1 + 1
      phase2_end <- min(CP2, n_days_farm)
      phase2_days <- phase2_end - CP1  # number of days in phase 2
      phase2_treatment <- sum(post_pred_usage_t[iter, f, (phase2_start + 1):(phase2_end + 1)], na.rm = TRUE)
      prop_phase2_farms[f] <- phase2_treatment / phase2_days
    } else {
      prop_phase2_farms[f] <- NA
    }
    # Phase 3
    if (n_days_farm > CP2) {
      phase3_start <- CP2 + 1
      phase3_end <- min(CP3, n_days_farm)
      phase3_days <- phase3_end - CP2  # number of days in phase 3
      phase3_treatment <- sum(post_pred_usage_t[iter, f, (phase3_start + 1):(phase3_end + 1)], na.rm = TRUE)
      prop_phase3_farms[f] <- phase3_treatment / phase3_days
    } else {
      prop_phase3_farms[f] <- NA
    }
    
    # Phase 4
    if (n_days_farm > CP3) {
      phase4_start <- CP3 + 1
      phase4_end <- n_days_farm
      phase4_days <- phase4_end - CP3  # number of days in phase 4
      phase4_treatment <- sum(post_pred_usage_t[iter, f, (phase4_start + 1):(phase4_end + 1)], na.rm = TRUE)
      prop_phase4_farms[f] <- phase4_treatment / phase4_days
    } else {
      prop_phase4_farms[f] <- NA
    }
  }
  
  # Average across farms for this iteration
  prop_phase1_iter[iter] <- mean(prop_phase1_farms, na.rm = TRUE)
  prop_phase2_iter[iter] <- mean(prop_phase2_farms, na.rm = TRUE)
  prop_phase3_iter[iter] <- mean(prop_phase3_farms, na.rm = TRUE)
  prop_phase4_iter[iter] <- mean(prop_phase4_farms, na.rm = TRUE)
}

# Calculate statistics for plotting
calc_stats <- function(data) {
  data.frame(
    median = median(data, na.rm = TRUE),
    mean = mean(data, na.rm = TRUE),
    q25 = quantile(data, 0.25, na.rm = TRUE),
    q75 = quantile(data, 0.75, na.rm = TRUE),
    hdi_low = hdi(data)[1],
    hdi_high = hdi(data)[2]
  )
}
phase_stats <- data.frame(
  Phase = c("P1\n(0-27)", "P2\n(28-53)", "P3\n(54-83)", "P4\n(84+)"),
  rbind(
    calc_stats(prop_phase1_iter),
    calc_stats(prop_phase2_iter),
    calc_stats(prop_phase3_iter),
    calc_stats(prop_phase4_iter)
  )
)

# Prepare data for plotting
proportion_data <- data.frame(
  Phase = rep(c("P1\n(0-27)", "P2\n(28-53)", "P3\n(54-83)", "P4\n(84+)"), 
              each = n_iter),
  Proportion = c(prop_phase1_iter, prop_phase2_iter, prop_phase3_iter, prop_phase4_iter)
)

# Create visualization
p5 <- ggplot(proportion_data, aes(x = Phase, y = Proportion, fill = Phase)) +
  # Violin plot
  geom_violin(alpha = 0.5) +
  # IQR box (25%-75% quantiles) with median line
  geom_crossbar(
    data = phase_stats,
    aes(x = Phase, y = median, ymin = q25, ymax = q75),
    width = 0.2,
    fatten = 2,
    inherit.aes = FALSE
  ) +
  # 95% HDI whiskers
  geom_errorbar(
    data = phase_stats,
    aes(x = Phase, ymin = hdi_low, ymax = hdi_high),
    width = 0.15,
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  # Mean diamond marker
  geom_point(
    data = phase_stats,
    aes(x = Phase, y = mean),
    shape = 23,
    size = 1.5,
    fill = "white",
    color = "black",
    stroke = 1,
    inherit.aes = FALSE
  ) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Phase", 
    y = "Proportion of days under therapeutic treatments"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.4, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 0.25), labels = scales::percent)

print(p5)
# p4 + p5 are combined and saved in Step04 (Figure 5)