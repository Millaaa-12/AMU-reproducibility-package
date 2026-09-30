# =============================================================================
# Step00_data_process.R
# Purpose : Build logistic (start of use) and Poisson (duration of use)
#           datasets for preventive and therapeutic AMU (used in Step02/Step03)
# Input   : data/VN_Daily_AMU.csv
# Output  : in-memory objects prevent_logistic_data, prevent_poisson_data,
#           thera_logistic_data, thera_poisson_data (no files written)
# Packages: loaded in Master_Run_All.R | Authors/license/dates: see README.md
# =============================================================================

### import dataset, distribution by relative date ###
data <- read.csv("data/VN_Daily_farm_records.csv", stringsAsFactors = FALSE)

# Make sure data are sorted by farm and time
data <- data %>%
  arrange(farm_id, Admin_days_cycle)

################# Preventive ##############

# Logistic regression dataset (preventive use)
prevent_logistic_data <- data %>%
  mutate(Drug_usage_final = ifelse(Drug_usage_final == 1 & preventive == 1, 1, 0)) %>%  # usage = 1 only when preventive = 1
  group_by(farm_id) %>%
  mutate(usage_yesterday = lag(Drug_usage_final, default = 0)) %>%  # usage status on the previous day
  filter(usage_yesterday == 0) %>%                                  # keep days with no use on the previous day
  dplyr::select(farm_id, Admin_days_cycle, Drug_usage_final, preventive) %>%
  ungroup()

# Poisson regression dataset (preventive use)
prevent_poisson_data <- data %>%
  mutate(Drug_usage_final = ifelse(Drug_usage_final == 1 & preventive == 1, 1, 0)) %>%  # usage = 1 only when preventive = 1
  group_by(farm_id) %>%
  mutate(event_change = Drug_usage_final != lag(Drug_usage_final, default = 0),  # flag status change
         event_id = cumsum(event_change & Drug_usage_final == 1)) %>%            # cumulative usage event ID
  filter(Drug_usage_final == 1) %>%                                              # keep usage days only
  group_by(farm_id, event_id, preventive) %>%
  summarise(duration = n(), .groups = "drop")                                    # duration (days) of each event

############## Therapeutic ##############

# Logistic regression dataset (therapeutic use)
thera_logistic_data <- data %>%
  mutate(Drug_usage_final = ifelse(Drug_usage_final == 1 & therapeutic == 1, 1, 0)) %>%  # usage = 1 only when therapeutic = 1
  group_by(farm_id) %>%
  mutate(usage_yesterday = lag(Drug_usage_final, default = 0)) %>%  # usage status on the previous day
  filter(usage_yesterday == 0) %>%                                  # keep days with no use on the previous day
  dplyr::select(farm_id, Admin_days_cycle, Drug_usage_final, therapeutic) %>%
  ungroup()

# Poisson regression dataset (therapeutic use)
thera_poisson_data <- data %>%
  mutate(Drug_usage_final = ifelse(Drug_usage_final == 1 & therapeutic == 1, 1, 0)) %>%  # usage = 1 only when therapeutic = 1
  group_by(farm_id) %>%
  mutate(event_change = Drug_usage_final != lag(Drug_usage_final, default = 0),  # flag status change
         event_id = cumsum(event_change & Drug_usage_final == 1)) %>%            # cumulative usage event ID
  filter(Drug_usage_final == 1) %>%                                              # keep usage days only
  group_by(farm_id, event_id, therapeutic) %>%
  summarise(duration = n(), .groups = "drop")  
  
############### Farm-level production cycle (one row per farm) #########################
# Each farm has multiple daily rows but a single Days_cycle value.
cycle <- data %>%
  distinct(farm_id, Days_cycle) %>%
  arrange(farm_id)            # order must match the farm index f used in the models