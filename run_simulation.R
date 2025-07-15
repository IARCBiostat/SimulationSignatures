rm(list=ls())
renv::activate()
# Scenario 1 ------------------------------------------------------------
# config_file_1 = "config/scenario_1.R"
# 
# message("Running scenario: ", config_file_1)
# system2("Rscript", args = c("R/run_scenario.R", config_file_1))

# Scenario 2 ------------------------------------------------------------
config_file_2 = "config/scenario_2.R"

message("Running scenario: ", config_file_2)
system2("Rscript", args = c("R/run_scenario.R", config_file_2))

