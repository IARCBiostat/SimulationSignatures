rm(list=ls())
renv::activate()
# Scenario 3 (supplementary) --------------------------------------------
config_file_3 = "config/supplementary_scenario_3.R"

message("Running scenario: ", config_file_3)
system2("Rscript", args = c("R/run_supplementary_scenario.R", config_file_3))