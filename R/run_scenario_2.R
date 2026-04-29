rm(list = ls())
library(glmnet)
library(foreach)
library(tibble)
source("R/analysis_fxns.R")
source("R/generate_data.R")
source("config/scenario_2.R")

set.seed(1+i)

if(i <= nrow(params)){
  message(
    paste0(
      "Running scenario 2 | iteration=", iteration,
      " | seed=", 1 + i,
      "\nParameters: n=", n,
      ", P=", P,
      ", prel=", prel,
      ", rho_g=", rho_g,
      ", rho_l=", rho_l,
      ", beta=", beta,
      ", Bhiddenmothers=", Bhiddenmothers,
      ", Bstepbrothers=", Bstepbrothers,
      "\n"
    )
  )
  
  gendat <- generate_data_scen2(n=n, P=P, prel=prel, rho_g=rho_g, rho_l=rho_l, beta=beta, Bhiddenmothers=Bhiddenmothers, Bstepbrothers=Bstepbrothers)
  gendat.test <- generate_data_scen2(n=10000, P=P, prel=prel, rho_g=rho_g, rho_l=rho_l, beta=beta, Bhiddenmothers=Bhiddenmothers, Bstepbrothers=Bstepbrothers)
  testres <- run_scenario(gendat, gendat.test, n=n, P=P, prel=prel, rho_g=rho_g, rho_l=rho_l, beta=beta, Bhiddenmothers=Bhiddenmothers, Bstepbrothers=Bstepbrothers, iteration=iteration)
  savedir <- "results/scenario_2/"
  savepath <- paste0(savedir, paste("scenario2",n,P,prel,rho_g,rho_l,beta,Bhiddenmothers,Bstepbrothers,iteration, sep="_"),".rds")
  
  if(!dir.exists(savedir)){dir.create(savedir, recursive=TRUE)}
  saveRDS(testres, savepath)
  
  message("Finished")
}

