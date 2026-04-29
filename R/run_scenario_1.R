args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Please provide a config file path.")

source(args[1])  # Load config

library(tidyverse)
library(mvtnorm)
library(glmnet)
library(parallel)
library(doParallel)
library(foreach)

source("R/generate_data.R")
source("R/analysis_fxns.R")

AllParams  <- expand.grid(NSIM, N, PMRELS, BETAS, rhoM)
colnames(AllParams) <- c("nsim", "n", "pMrel", "beta", "rhoM")

clusSim <- makeCluster(NCore)
registerDoParallel(clusSim, outfile='')

# Run simulations for all combinations of parameters
all_res <- foreach (i = 1:nrow(AllParams), .combine=rbind, .packages=c("glmnet", "mvtnorm", "tidyverse")) %dopar% {
  cat(i, "\n")
  set.seed(i)
  pMrel <- AllParams[i, "pMrel"]
  beta  <- AllParams[i, "beta"]
  n     <- AllParams[i, "n"]
  nsim  <- AllParams[i, "nsim"]
  rhoM  <- AllParams[i, "rhoM"]
  
  dat.train <- generate.data(n = n, pM = pM, pMrel = pMrel,
                             SecondLayerM = SecondLayerM,
                             beta = beta,
                             rhoM = rhoM)
  dat.test <- generate.data(n = 1e4, pM = pM, pMrel = pMrel,
                            SecondLayerM = SecondLayerM,
                            beta = beta,
                            rhoM = rhoM)
  
  tryCatch({
    all_coefs <- estimate_coef_sig_w_wo_prior_screening(dat.train)
    VecE1onM <- dat.train$VecE1onM
    
    temp_tibble <- tibble(
      Approach = names(all_coefs),
      CardS    = sapply(all_coefs, function(x) sum(x != 0)),
      Sensi    = sapply(all_coefs, function(x) sum(x != 0 & VecE1onM != 0) / sum(VecE1onM != 0)),
      Speci    = sapply(all_coefs, function(x) sum(x == 0 & VecE1onM == 0) / sum(VecE1onM == 0)),
      Corr     = sapply(all_coefs, function(x) cor_sig_exp_test(x, dat.test)),
      n, pM, pMrel, SecondLayerM, beta, rhoM, nsim
    )
    return(temp_tibble)
  }, error = function(e) {
    write(paste0(pE, "\t", SecondLayerM, "\t", i),
          file = paste0("results/Error_SecondLayerM_", SecondLayerM, "_rhoM_", rhoM, ".txt"), append = TRUE)
    return(NULL)
  })
}

saveRDS(all_res, file = paste0("results/Res_SecondLayerM_", SecondLayerM, ".rds"))
stopCluster(clusSim)
