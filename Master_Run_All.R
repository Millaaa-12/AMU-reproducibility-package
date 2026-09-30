# =============================================================================
# Master_Run_All.R  -  runs the complete analysis from beginning to end
# HOW TO RUN
#   1. Install R [4.4.3] and JAGS [4.3.2] (https://mcmc-jags.sourceforge.io/)
#   2. Set the working directory to the root folder of this package
#      (the folder containing this file), e.g. setwd("path/to/AMU_Reproducibility_Package")
#   3. source("Master_Run_All.R")
#
# EXECUTION ORDER: Step00 -> Step01 -> Step02 -> Step03 -> Step04
#   Steps 02-04 use R objects created by earlier steps, so all steps must be
#   run in this order in ONE R session (this file does that).
# =============================================================================

rm(list = ls())
start_time <- Sys.time()
setwd("path/to/AMU_Reproducibility_Package")

############ 0. Check working directory and create output folders #################
if (!file.exists("data/VN_Daily_farm_records.csv")) {
  stop("Working directory must be the package root folder (the folder that contains 'data/').")
}
dir.create("output", showWarnings = FALSE)
dir.create("log",    showWarnings = FALSE)

########### 1. Log file (console output, warnings and messages) ####################
log_con <- file("log/Master_log.txt", open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
cat("Analysis started:", format(start_time), "\n")

########## 2. Packages: installed if missing, loaded once for all steps ############
pkgs <- c("dplyr", "tidyverse", "ggplot2", "tidyr", "gratia", "mgcv",
          "marginaleffects", "patchwork",                        # Step00-01
          "R2jags", "bayesplot", "rjags", "HDInterval", "coda",
          "parallel", "runjags", "extraDistr", "doParallel",    # Step02-03
          "scales", "tibble","flextable", "officer")                                              # Step04
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

############################### 3. Run all steps in sequence #######################
run_step <- function(f) {
  cat("\n\n==================== Running", f, "====================\n")
  t0 <- Sys.time()
  source(f, echo = TRUE, max.deparse.length = Inf, encoding = "UTF-8")
  cat("\n---", f, "finished in",
      format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "---\n")
}

run_step("Step00_data_process.R")
run_step("Step01_GAMM.R")
run_step("Step02_Markov_preventive.R")
run_step("Step03_Markov_therapeutic.R")
run_step("Step04_All_AMU_simulation.R")
run_step("Step05_AMU_quantification.R")

########### 4. Computing environment (recorded AFTER all packages are loaded) ########
cat("\n\n==================== Computing environment ====================\n")
print(sessionInfo())
for (pkg in pkgs) cat(sprintf("%-16s %s\n", pkg, as.character(packageVersion(pkg))))
try(cat("JAGS version:", as.character(rjags::jags.version()), "\n"))
cat("Logical CPU cores:", parallel::detectCores(), "\n")

end_time <- Sys.time()
cat("Analysis finished:", format(end_time), "\n")
cat("Total run time:", format(round(difftime(end_time, start_time, units = "hours"), 2)), "\n")

sink(type = "message"); sink(); close(log_con)
