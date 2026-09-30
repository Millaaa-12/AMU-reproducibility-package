# =============================================================================
# Step01_GAMM.R
# Purpose : (1) GAMMs of preventive / therapeutic / all AMU vs. broiler age (day) and
#           within-flock mortality risk; 
#           (2)partial effects, 
#           (3) interaction,
#           (4) 3-D surfaces
#           (5) rate-of-change (marginal effect) of day
# Input   : data/data/VN_Daily_AMU.csv.csv
# Output  : GAM summaries (Table 3)
#           output/Figure_S4_GAMM_partial_effects.tiff
#           output/Figure_S5_Interaction between broiler age and within-flock mortality risk on the probability of AMU.tiff
#           output/Figure_S6_Three-dimensional graphics.tiff
#           output/Figure_S10_cut points plot -> Rate of Change by Broiler Age
# Packages: loaded in Master_Run_All.R | Authors/license/dates: see README.md
# =============================================================================

### data preparation###
data_updated <- data %>%
  group_by(farm_id) %>%
  mutate(farm_id = factor(farm_id)) %>%
  arrange(farm_id, Admin_days_cycle) %>% 
  mutate(
    mortality_risk = ifelse(alive_start_of_day > 0, 
                            dead / alive_start_of_day, 0),
    
    day = as.integer(Admin_days_cycle),
    
    preventive = as.factor(preventive),
    therapeutic = as.factor(therapeutic),
    all_AMU = as.factor(Drug_usage_final)
  ) %>%
  ungroup()

####################GAM MODEL 1 - PREVENTIVE#########################
# ---- This block creates Table X / in-text results: GAMM for preventive AMU ----
model1 <- gam(
  preventive ~ 
    s(day, k = 10) +  
    s(mortality_risk, k = 10) +  
    ti(day,  mortality_risk, k = c(4, 4)) +    
    s(farm_id, bs = "re"),
  data = data_updated,
  family = binomial(),
  method = "REML"
)
summary(model1)
#par(mfrow = c(2, 2))
#gam.check(model1)
#k.check(model1, subsample=5000, n.rep=400)

############################GAM MODEL 2 - THERAPEUTIC###################
model2 <- gam(
  therapeutic ~ 
    s(day, k = 10) + 
    s(mortality_risk, k = 10) +  
    ti(day,  mortality_risk, k = c(8, 5)) +       
    s(farm_id, bs = "re"),
  data = data_updated,
  family = binomial(),
  method = "REML"
)
summary(model2)
#par(mfrow = c(2, 2))
#gam.check(model2)
#k.check(model2, subsample=5000, n.rep=400)

#########################GAM MODEL 3 - ALL AMU############################
# ---- This block creates Table X / in-text results: GAMM for all AMU ----
model3 <- gam(
  all_AMU ~ 
    s(day, k = 10) +  
    s(mortality_risk, k = 10) +  
    ti(day,  mortality_risk, k = c(10, 10)) +       
    s(farm_id, bs = "re"),
  data = data_updated,
  family = binomial(),
  method = "REML"
)
summary(model3)
#par(mfrow = c(2, 2))
#gam.check(model3)
#k.check(model3, subsample=5000, n.rep=400)

#############Partial effect plot#####################################
p1 <- draw(model1)
p2 <- draw(model2)
p3 <- draw(model3)

combined_plot <- wrap_plots(
  wrap_elements(p1),
  wrap_elements(p2),
  wrap_elements(p3),
  plot_spacer(),
  ncol = 2, 
  nrow = 2
) + plot_annotation(tag_levels = "a")

combined_plot
ggsave("output/Figure_S4_GAMM_partial_effects.tiff", 
       plot = combined_plot,
       width = 20, 
       height = 18, 
       dpi = 300, 
       device = "tiff")

###### significance of the day x mortality interaction on the probability of AMU ############
plot_ti_sig <- function(model,
                        smooth   = "ti(day,mortality_risk)",
                        n        = 100,      # number of grid points per dimension
                        level    = 0.95,
                        too_far  = 0.1,      # hide grid points too far from the data; NULL = show all
                        title    = NULL) {
  
  crit <- qnorm((1 + level) / 2)
  sm <- smooth_estimates(model, select = smooth, n = n)
  if (!".estimate" %in% names(sm)) {
    sm <- sm |> rename(.estimate = est, .se = se)
  }
  sm <- sm |>
    mutate(
      lower = .estimate - crit * .se,
      upper = .estimate + crit * .se,
      sig = case_when(
        lower > 0 ~ "Above 0",
        upper < 0 ~ "Below 0",
        TRUE      ~ "CI includes 0"
      ),
      sig = factor(sig, levels = c("Above 0", "Below 0", "CI includes 0"))
   )
  
  # interaction significance plot
  if (!is.null(too_far)) {
    dat <- model$model
    far <- mgcv::exclude.too.far(sm$day, sm$mortality_risk,
                                 dat$day, dat$mortality_risk,
                                 dist = too_far)
    sm <- sm[!far, ]
  }
  ggplot(sm, aes(x = day, y = mortality_risk)) +
    geom_raster(aes(fill = sig)) +
    geom_contour(aes(z = .estimate), colour = "black",
                 linewidth = 0.3, bins = 10) +
    scale_fill_manual(
      values = c("Above 0"       = "#D7301F",
                 "Below 0"       = "#2B8CBE",
                 "CI includes 0" = "grey92"),
      drop = FALSE,
      name = paste0(level * 100, "% CI")
    ) +
    labs(x = "Day", y = "Mortality risk",
    ) +
    theme_bw()
}

p1_int <- plot_ti_sig(model1)
p2_int <- plot_ti_sig(model2)
p3_int <- plot_ti_sig(model3)
p_all <- (p1_int | p2_int | p3_int) + plot_layout(guides = "collect")

ggsave("output/Figure_S5_GAMM_interaction_significance.tiff", plot = p_all,
       width = 12, height = 3, units = "in",
       dpi = 300,
       device = "tiff",
       compression = "lzw")     

##########surface with standard error##########################################
# Open TIFF graphics device
tiff("output/Figure_S6_GAMM_3D_surface_SE.tiff", width = 15, height = 5, units = "in", res = 300)

# Set up a 1x3 layout
par(mfrow = c(1, 3))
vis.gam(model1,
        se = 2,
        main = "(a)",
        xlab = "Age of broiler (day)", 
        ylab = "Within-flock mortality",
        zlab = "The logit of the probability of preventive AMU")

vis.gam(model2,
        se = 2,
        main = "(b)",
        xlab = "Age of broiler (day)", 
        ylab = "Within-flock mortality",
        zlab = "The logit of the probability of therapeutic AMU")

vis.gam(model3,
        se = 2,
        main = "(c)",
        xlab = "Age of broiler (day)", 
        ylab = "Within-flock mortality",
        zlab = "The logit of the probability of AMU")

# Reset to single plot layout
par(mfrow = c(1, 1))

# Close the graphics device
dev.off()

#############Definition of cutpoint (critical change point)#############################
# 1. Slopes of the predicted probability with respect to day
get_slopes <- function(model) {
  s <- slopes(model,
              variables = "day",
              newdata   = datagrid(model = model,
                                   day = seq(min(model$model$day), max(model$model$day), length.out = 200)),
              type      = "response",
              exclude   = "s(farm_id)")   # typical farm (random effect = 0)
  s <- as.data.frame(s)
  s <- s[order(s$day), ]                      # sort by day
  s$sig <- s$conf.low > 0 | s$conf.high < 0   # significant if the CI excludes 0
  s
}

# 2. Intervals of significant change (consecutive significant days)
get_intervals <- function(df) {
  dir <- ifelse(df$conf.low > 0, "Increase",
                ifelse(df$conf.high < 0, "Decrease", "ns"))
  r     <- rle(dir)
  end   <- cumsum(r$lengths)
  start <- end - r$lengths + 1
  out <- data.frame(direction = r$values,
                    day_start = round(df$day[start], 1),
                    day_end   = round(df$day[end], 1))
  out[out$direction != "ns", ]
}

# 3. Turning points (slope = 0, linear interpolation)
get_turning <- function(df) {
  i <- which(diff(sign(df$estimate)) != 0)
  data.frame(
    day  = round(df$day[i] - df$estimate[i] * (df$day[i + 1] - df$day[i]) /
                   (df$estimate[i + 1] - df$estimate[i]), 1),
    type = ifelse(df$estimate[i] > 0, "Peak (up -> down)", "Trough (down -> up)")
  )
}

# 4. Plot (significant parts in red)
plot_slopes_df <- function(df, ylab) {
  ggplot(df, aes(day, estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "grey80") +
    geom_line(aes(colour = sig, group = 1), linewidth = 1) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red"), guide = "none") +
    labs(x = "Days of age", y = ylab) +
    theme_minimal(base_size = 12)
}

# 5. Run
s1 <- get_slopes(model1)
s2 <- get_slopes(model2)
s3 <- get_slopes(model3)

p1 <- plot_slopes_df(s1, "Rate of change in probability\nof a farm with preventive AMU")
p2 <- plot_slopes_df(s2, "Rate of change in probability\nof a farm with therapeutic AMU")
p3 <- plot_slopes_df(s3, "Rate of change in probability\nof a farm with AMU")

# 6. Combine and save
combined <- p1 + p2 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
combined

ggsave("output/Figure_10_cutpoint_GAMM_marginal_effect_day.tiff", combined,
       width = 8, height = 4, dpi = 300, compression = "lzw", bg = "white")