# =============================================================================
# Step04_All_AMU_simulation.R
# Purpose : (1) risk ratio of initiating preventive vs therapeutic AMU by broiler age;
#           (2) final combined figure 5 (RR + p1/p2 + p4/p5);
#           (3) joint posterior-predictive simulation of P-AMU, T-AMU and any AMU (PPC)
# Input   : jags_fit_p, p1, p2 (Step02); jags_fit, p4, p5 (Step03); data (Step01)
#          
# Output  : output/Figure_5_RR_initiation_and_phase_proportions.tiff
#           output/Figure_S15_PPC_daily_number_of_farms.tiff
# Packages: loaded in Master_Run_All.R | Authors/license/dates: see README.md
# =============================================================================
# Step02 (P-AMU) and Step03 (T-AMU) must be run first!

################################ Joint simulation of P-AMU and T-AMU ################################
# production cycle length of farms
farm_cycles <- cycle$Days_cycle   # one value per farm (created in Step00)
n_days_max  <- max(farm_cycles)   # longest production cycle across farms
n_farms <- 60      # HARD-CODED: number of study farms

# posterior samples
jags_samples_p <- jags_fit_p$BUGSoutput$sims.list  # P-AMU
jags_samples_t <- jags_fit$BUGSoutput$sims.list    # T-AMU

# HARD-CODED cut points (same as in Step02 / Step03)
CP1_p <- 57
CP2_p <- 79

CP1_t <- 27
CP2_t <- 53
CP3_t <- 83

################ risk ratio of initiating preventive vs therapeutic AMU by broiler age ##################
max_days <- n_days_max
days     <- 0:max_days
n_days   <- length(days)

# P-AMU posterior samples
beta0_P          <- jags_samples_p$beta0
beta_interval0_P <- jags_samples_p$beta_interval0
beta_interval1_P <- jags_samples_p$beta_interval1
beta_interval2_P <- jags_samples_p$beta_interval2

# T-AMU posterior samples
beta0_T          <- jags_samples_t$beta0
beta_interval0_T <- jags_samples_t$beta_interval0
beta_interval1_T <- jags_samples_t$beta_interval1
beta_interval2_T <- jags_samples_t$beta_interval2
beta_interval3_T <- jags_samples_t$beta_interval3

n_iter_p <- length(beta0_P)
n_iter_t <- length(beta0_T)
stopifnot(n_iter_p == n_iter_t)
n_iter <- n_iter_p

# daily p_P, p_T, RR posterior (n_iter x n_days)
p_init_P <- matrix(NA_real_, nrow = n_iter, ncol = n_days)
p_init_T <- matrix(NA_real_, nrow = n_iter, ncol = n_days)
RR       <- matrix(NA_real_, nrow = n_iter, ncol = n_days)

for (day in 0:max_days) {
  # P-AMU
  logit_p_P <- beta0_P +
    beta_interval0_P * day +
    beta_interval1_P * pmax(day - CP1_p, 0) +
    beta_interval2_P * pmax(day - CP2_p, 0)
  probs_P <- plogis(logit_p_P)
  
  # T-AMU
  logit_p_T <- beta0_T +
    beta_interval0_T * day +
    beta_interval1_T * pmax(day - CP1_t, 0) +
    beta_interval2_T * pmax(day - CP2_t, 0) +
    beta_interval3_T * pmax(day - CP3_t, 0)
  probs_T <- plogis(logit_p_T)
  
  # Store results
  p_init_P[, day + 1] <- probs_P
  p_init_T[, day + 1] <- probs_T
  RR[, day + 1]       <- probs_P / probs_T   # ratio of daily initiation probabilities
}

# summary: median + 95% HDI
summarize_post <- function(mat, days) {
  q_med <- apply(mat, 2, median, na.rm = TRUE)
  q_lo  <- apply(mat, 2, function(z) {
    z <- z[is.finite(z)]; if (length(z) < 2) return(NA); hdi(z, 0.95)[1]
  })
  q_hi  <- apply(mat, 2, function(z) {
    z <- z[is.finite(z)]; if (length(z) < 2) return(NA); hdi(z, 0.95)[2]
  })
  data.frame(day = days, median = q_med, lo = q_lo, hi = q_hi)
}

sum_pP <- summarize_post(p_init_P, days)
sum_pT <- summarize_post(p_init_T, days)
sum_RR <- summarize_post(RR,       days)

daily_table <- data.frame(
  day       = days,
  pP_median = sum_pP$median, pP_lo = sum_pP$lo, pP_hi = sum_pP$hi,
  pT_median = sum_pT$median, pT_lo = sum_pT$lo, pT_hi = sum_pT$hi,
  RR_median = sum_RR$median, RR_lo = sum_RR$lo, RR_hi = sum_RR$hi
)
# plot
# last day (days of age) of the production cycle of each farm
farm_cycles <- cycle$Days_cycle

# 1) number of farms still within their cycle on each day
# farms active on day d = farms with cycle length >= d
n_active_per_day <- sapply(days, function(d) sum(farm_cycles >= d))

# 2) last day with >= 10 active farms
min_active  <- 10   # HARD-CODED: minimum number of active farms to display RR: [AUTHOR: justify]
max_valid_d <- max(days[n_active_per_day >= min_active])

cat(sprintf("Total farms loaded: %d\n", length(farm_cycles)))
cat(sprintf("Cutoff: day %d  (last day with >= %d active farms; n at cutoff = %d)\n",
            max_valid_d, min_active,
            n_active_per_day[which(days == max_valid_d)]))

sum_pP_f <- subset(sum_pP, day <= max_valid_d)
sum_pT_f <- subset(sum_pT, day <= max_valid_d)
sum_RR_f <- subset(sum_RR, day <= max_valid_d)

# Risk Ratio over age (truncated)
p_RR <- ggplot(sum_RR_f, aes(x = day)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#e07b39", alpha = 0.25) +
  geom_line(aes(y = median), color = "#c0531a", linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  scale_y_continuous(trans = "log10",
                     breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10, 25, 100, 500),
                     labels = c("0.1","0.25","0.5","1","2","4",
                                "10","25","100","500")) +
  scale_x_continuous(breaks = seq(0, max_valid_d, 20),
                     limits = c(0, max_valid_d),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Days of age",
       y = expression("RR = " * p[P](d) / p[T](d) ~ "  (log scale)"),
       subtitle = sprintf(
         "Posterior median and 95%% HDI; restricted to days with \u2265 %d active farms (day 0\u2013%d)",
         min_active, max_valid_d)) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold"))

n_df <- data.frame(day = days, n = n_active_per_day) |>
  subset(day <= max_valid_d)

p_n <- ggplot(n_df, aes(x = day, y = n)) +
  geom_area(fill = "grey80") +
  geom_line(color = "grey40") +
  geom_hline(yintercept = min_active, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(0, max_valid_d, 20),
                     limits = c(0, max_valid_d),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Days of age", y = "n farms\nat cycle") +
  theme_classic(base_size = 16)

### plot combine (creates Figure 5 (RR + preventive p1/p2 + therapeutic p4/p5) ###
p3 <- p1 + p2 # combine plots (side by side)      
p6 <- p4 + p5

# final combine
# 1) common font size (including legends)
common_theme <- theme(
  text          = element_text(size = 16),
  axis.title    = element_text(size = 16),
  axis.text     = element_text(size = 16),
  legend.title  = element_text(size = 16),
  legend.text   = element_text(size = 16),
  plot.title    = element_text(size = 16),
  plot.subtitle = element_text(size = 16),
  plot.caption  = element_text(size = 16),
  strip.text    = element_text(size = 16)
)

# 2) apply common font size inside p_rr, keep 4:1 height ratio
p_rr <- ((p_RR & common_theme) / (p_n & common_theme)) +
  plot_layout(heights = c(4, 1))

# 3) wrap_elements() treats p_rr as one element (prevents patchwork from flattening the nesting)
p_final <- wrap_elements(p_rr) /
  (p3 & common_theme) /
  (p6 & common_theme)

ggsave(
  filename = "output/Figure_5_RR_initiation_and_phase_proportions.tiff",
  plot     = p_final,
  device   = "tiff",
  width    = 15,
  height   = 18,
  dpi      = 300,
  units    = "in"
)

################### posterior predictive check: daily number of farms using AMU) ###########################
n_iter <- 5000
chosen_iter_p <- sample(seq_len(length(jags_samples_p$beta0)), n_iter)
chosen_iter_t <- sample(seq_len(length(jags_samples_t$beta0)), n_iter)

# initialise storage arrays
post_pred_usage_p <- array(NA, dim = c(n_iter, n_farms, n_days_max + 1))
post_pred_usage_t <- array(NA, dim = c(n_iter, n_farms, n_days_max + 1))
post_pred_usage_combined <- array(NA, dim = c(n_iter, n_farms, n_days_max + 1))

post_pred_durations_p <- array(list(), dim = c(n_iter, n_farms))
post_pred_durations_t <- array(list(), dim = c(n_iter, n_farms))

# main loop
for (iter in 1:n_iter) {
  
  # P-AMU parameters
  beta0_p <- jags_samples_p$beta0[chosen_iter_p[iter]]
  beta_interval0_p <- jags_samples_p$beta_interval0[chosen_iter_p[iter]]
  beta_interval1_p <- jags_samples_p$beta_interval1[chosen_iter_p[iter]]
  beta_interval2_p <- jags_samples_p$beta_interval2[chosen_iter_p[iter]]
  alpha0_p <- jags_samples_p$alpha0[chosen_iter_p[iter]]
  
  # T-AMU parameters
  beta0_t <- jags_samples_t$beta0[chosen_iter_t[iter]]
  beta_interval0_t <- jags_samples_t$beta_interval0[chosen_iter_t[iter]]
  beta_interval1_t <- jags_samples_t$beta_interval1[chosen_iter_t[iter]]
  beta_interval2_t <- jags_samples_t$beta_interval2[chosen_iter_t[iter]]
  beta_interval3_t <- jags_samples_t$beta_interval3[chosen_iter_t[iter]]
  alpha0_t <- jags_samples_t$alpha0[chosen_iter_t[iter]]
  
  for (f in 1:n_farms) {
    n_days_farm <- farm_cycles[f]

    ### P-AMU simulation ###
    remaining_duration_p <- 0
    durations_farm_p <- integer(0)
    
    for (day in 0:n_days_farm) {
      array_index <- day + 1
      
      if (remaining_duration_p > 0) {
        post_pred_usage_p[iter, f, array_index] <- 1
        remaining_duration_p <- remaining_duration_p - 1
      } else {
        logit_p <- beta0_p + 
          beta_interval0_p * day + 
          beta_interval1_p * pmax(day - CP1_p, 0) + 
          beta_interval2_p * pmax(day - CP2_p, 0)
        
        start_prob <- plogis(logit_p)
        start_usage <- rbinom(1, 1, start_prob)
        post_pred_usage_p[iter, f, array_index] <- start_usage
        
        if (start_usage == 1) {
          lambda <- exp(alpha0_p)
          duration <- rtpois(1, lambda = lambda, a = 0)
          duration <- min(duration, n_days_farm - day + 1)
          remaining_duration_p <- duration - 1
          durations_farm_p <- c(durations_farm_p, duration)
        }
      }
    }
    
    ### T-AMU simulation ###
    remaining_duration_t <- 0
    durations_farm_t <- integer(0)
    
    for (day in 0:n_days_farm) {
      array_index <- day + 1
      
      if (remaining_duration_t > 0) {
        post_pred_usage_t[iter, f, array_index] <- 1
        remaining_duration_t <- remaining_duration_t - 1
      } else {
        logit_p <- beta0_t + 
          beta_interval0_t * day + 
          beta_interval1_t * pmax(day - CP1_t, 0) + 
          beta_interval2_t * pmax(day - CP2_t, 0) +
          beta_interval3_t * pmax(day - CP3_t, 0)
        
        start_prob <- plogis(logit_p)
        start_usage <- rbinom(1, 1, start_prob)
        post_pred_usage_t[iter, f, array_index] <- start_usage
        
        if (start_usage == 1) {
          lambda <- exp(alpha0_t)
          duration <- rtpois(1, lambda = lambda, a = 0)
          duration <- min(duration, n_days_farm - day + 1)
          remaining_duration_t <- duration - 1
          durations_farm_t <- c(durations_farm_t, duration)
        }
      }
    }
    
    ### Combined: either P or T use ###
    for (day in 0:n_days_farm) {
      array_index <- day + 1
      post_pred_usage_combined[iter, f, array_index] <- 
        as.integer(post_pred_usage_p[iter, f, array_index] == 1 | 
                     post_pred_usage_t[iter, f, array_index] == 1)
    }
    
    # days beyond the farm's cycle set to NA
    if (n_days_farm < n_days_max) {
      post_pred_usage_p[iter, f, (n_days_farm + 2):(n_days_max + 1)] <- NA
      post_pred_usage_t[iter, f, (n_days_farm + 2):(n_days_max + 1)] <- NA
      post_pred_usage_combined[iter, f, (n_days_farm + 2):(n_days_max + 1)] <- NA
    }
    
    post_pred_durations_p[[iter, f]] <- durations_farm_p
    post_pred_durations_t[[iter, f]] <- durations_farm_t
  }
  
  if (iter %% 1000 == 0) {
    cat("Completed iteration:", iter, "/", n_iter, "\n")
  }
}

### summary of the daily number of farms using AMU ####
daily_count_p <- apply(post_pred_usage_p, c(1, 3), sum, na.rm = TRUE)
daily_count_t <- apply(post_pred_usage_t, c(1, 3), sum, na.rm = TRUE)
daily_count_combined <- apply(post_pred_usage_combined, c(1, 3), sum, na.rm = TRUE)

# summary statistics
summary_p <- data.frame(
  day = 0:n_days_max,
  median = apply(daily_count_p, 2, median, na.rm = TRUE),
  q025 = apply(daily_count_p, 2, quantile, 0.025, na.rm = TRUE),
  q25 = apply(daily_count_p, 2, quantile, 0.25, na.rm = TRUE),
  q75 = apply(daily_count_p, 2, quantile, 0.75, na.rm = TRUE),
  q975 = apply(daily_count_p, 2, quantile, 0.975, na.rm = TRUE)
)

summary_t <- data.frame(
  day = 0:n_days_max,
  median = apply(daily_count_t, 2, median, na.rm = TRUE),
  q025 = apply(daily_count_t, 2, quantile, 0.025, na.rm = TRUE),
  q25 = apply(daily_count_t, 2, quantile, 0.25, na.rm = TRUE),
  q75 = apply(daily_count_t, 2, quantile, 0.75, na.rm = TRUE),
  q975 = apply(daily_count_t, 2, quantile, 0.975, na.rm = TRUE)
)

summary_combined <- data.frame(
  day = 0:n_days_max,
  median = apply(daily_count_combined, 2, median, na.rm = TRUE),
  q025 = apply(daily_count_combined, 2, quantile, 0.025, na.rm = TRUE),
  q25 = apply(daily_count_combined, 2, quantile, 0.25, na.rm = TRUE),
  q75 = apply(daily_count_combined, 2, quantile, 0.75, na.rm = TRUE),
  q975 = apply(daily_count_combined, 2, quantile, 0.975, na.rm = TRUE)
)

# function to process observed data
process_observed_data <- function(observed_data) {
  all_combinations <- observed_data %>% distinct(Admin_days_cycle)
  
  used_data <- observed_data %>%
    filter(Drug_usage_final == 1) %>%
    group_by(Admin_days_cycle) %>%
    summarise(Usage_Count = n(), .groups = "drop")
  
  total_daily_usage <- all_combinations %>%
    left_join(used_data, by = "Admin_days_cycle") %>%
    mutate(Usage_Count = replace_na(Usage_Count, 0))
  
  return(total_daily_usage)
}

### observed data ###

# Observed data: preventive use only
observed_data_p <- data %>%
  filter(preventive == 1)

# Observed data: therapeutic use only
observed_data_t <- data %>%
  filter(therapeutic == 1)

# Observed data: all use regardless of purpose
observed_data_all <- data %>%
  filter(Drug_usage_final == 1)

# Calculate total daily usage
total_daily_usage_p   <- process_observed_data(observed_data_p)
total_daily_usage_t   <- process_observed_data(observed_data_t)
total_daily_usage_all <- process_observed_data(observed_data_all)

# merge simulated and observed data
plot_data_p <- summary_p %>%
  left_join(total_daily_usage_p, by = c("day" = "Admin_days_cycle")) %>%
  rename(actual_usage = Usage_Count)

plot_data_t <- summary_t %>%
  left_join(total_daily_usage_t, by = c("day" = "Admin_days_cycle")) %>%
  rename(actual_usage = Usage_Count)

plot_data_all <- summary_combined %>%
  left_join(total_daily_usage_all, by = c("day" = "Admin_days_cycle")) %>%
  rename(actual_usage = Usage_Count)

#  plotting function
create_ppc_plot <- function(plot_data, point_color, panel_label, y_label, subtitle_text) {
  
  ggplot(plot_data, aes(x = day)) +
    # whiskers: 2.5%-97.5% quantiles (q025, q975) - see README issue 7
    geom_linerange(aes(ymin = q025, ymax = q975), 
                   color = "gray30", linewidth = 0.3) +
    # IQR (boxes): q25 and q75
    geom_crossbar(aes(y = median, ymin = q25, ymax = q75), 
                  fill = "gray75", color = "gray50", 
                  width = 0.8, fatten = 0, linewidth = 0.2) +
    # Observed data points
    geom_point(aes(y = actual_usage), 
               color = point_color, size = 1.2) +
    # axes
    scale_x_continuous(
      breaks = seq(0, 160, by = 5),
      limits = c(-1, 162),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      breaks = seq(0, 40, by = 10),
      limits = c(0, 40),
      expand = c(0.02, 0)
    ) +
    # labels
    labs(
      x = "Days of age",
      y = y_label,
      subtitle = subtitle_text,
      tag = panel_label
    ) +
    # theme
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 11, color = "black", margin = margin(t = 5)),
      axis.title.y = element_text(size = 11, color = "black", margin = margin(r = 5)),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0.01, 0.98),
      plot.subtitle = element_text(size = 8, color = "black", hjust = 0, margin = margin(b = 5)),
      plot.margin = margin(5, 10, 5, 5)
    )
}

# create the three panels
p_A <- create_ppc_plot(
  plot_data = plot_data_p,
  point_color = "#CC79A7",  # magenta/pink
  panel_label = "A",
  y_label = "Number of farms with preventive AMU",
  subtitle_text = "Boxes: IQR (25%-75%), Whiskers: 95% HDI\nRed points: Observed numbers"
)

p_B <- create_ppc_plot(
  plot_data = plot_data_t,
  point_color = "#0072B2",  # blue
  panel_label = "B",
  y_label = "Number of farms with therapeutic AMU",
  subtitle_text = "Boxes: IQR (25%-75%), Whiskers: 95% HDI\nBlue points: Observed numbers"
)

p_C <- create_ppc_plot(
  plot_data = plot_data_all,
  point_color = "#56ba8c",  # green
  panel_label = "C",
  y_label = "Number of farms with AMU",
  subtitle_text = "Boxes: IQR (25%-75%), Whiskers: 95% HDI\nBlue points: Observed numbers"
)

# combine and save
combined_plot <- p_A / p_B / p_C

ggsave("output/Figure_S15_PPC_daily_number_of_farms.tiff", combined_plot, 
       width = 10, height = 10, dpi = 300, bg = "white")