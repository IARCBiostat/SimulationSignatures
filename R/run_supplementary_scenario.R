args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Please provide a config file path.")

source(args[1])  # Load config

library(tidyverse)
library(mvtnorm)
library(glmnet)
library(parallel)
library(doParallel)
library(foreach)

source("R/generate_latent_data.R")
source("R/analysis_fxns.R")

AllParams  <- expand.grid(NSIM, N, PLS, pM, PLRELS, BETAS, RHOLS, typeCorM)
colnames(AllParams) <- c("nsim", "n", "pL","pM", "pLrel", "beta", "rhoL", "typeCorM")

clusSim     <- makeCluster(NCore)
registerDoParallel(clusSim, outfile='')
all_res <- foreach (i = 1:nrow(AllParams), .combine=rbind, .packages=c("glmnet", "mvtnorm", "tidyverse")) %dopar%
  {
    cat(i, "\n")
    set.seed(i)
    pL    <- AllParams[i, "pL"]
    pM    <- AllParams[i, "pM"]
    pLrel <- AllParams[i, "pLrel"]
    beta  <- AllParams[i, "beta"]
    rhoL  <- AllParams[i, "rhoL"]
    n     <- AllParams[i, "n"]
    nsim  <- AllParams[i, "nsim"]
    typeCorM <- AllParams[i, "typeCorM"]
    
    dat.train <- generate.data.withlatent(n = n,pL=pL, pM = pM, pLrel = pLrel,
                                          SecondLayerL = SecondLayerL, typeCorM=typeCorM,
                                          beta = beta,
                                          rhoL = rhoL,
                                          rhoM = rhoM)
    
    dat.test <- generate.data.withlatent(n = 1e4,pL=pL, pM = pM, pLrel = pLrel,
                                         SecondLayerL = SecondLayerL, typeCorM = typeCorM,
                                         beta = beta,
                                         rhoL = rhoL,
                                         rhoM = rhoM)
    
    
    tryCatch(
      {
        all_coefs <- estimate_coef_sig_w_wo_prior_screening(dat.train)
        
        VecE1onL        <- dat.train$VecE1onL
        ind_Mrel        <- dat.train$ind_M_rel
        VecE1onM        <- rep(0, pM)
        VecE1onM[unlist(ind_Mrel)] <- 1
        sensitivityVV_L <- function(pred){
          sum(sapply(ind_Mrel, function(templll){1*sum(pred[templll]!=0)>0}))/length(ind_Mrel)
        }
        sensitivityVV   <- function(pred){sum(pred!=0 & VecE1onM!=0)/sum(VecE1onM!=0)}
        specificityVV   <- function(pred){sum(pred==0 & VecE1onM==0)/sum(VecE1onM==0)}
        Sensi           <- do.call(rbind, lapply(all_coefs, sensitivityVV))
        Sensi_L         <- do.call(rbind, lapply(all_coefs, sensitivityVV_L))
        Speci           <- do.call(rbind, lapply(all_coefs, specificityVV))
        Corrs <- do.call(rbind, lapply(all_coefs, function(x) cor_sig_exp_test(x, dat.test)))
        temp_tibble <- tibble(Approach = rownames(Sensi), Sensi=Sensi[, 1], Sensi_L=Sensi_L[, 1], Speci=Speci[, 1], Corr=Corrs[,1],
                              n, pE, pL, pM, pLrel, SecondLayerL, typeCorM, beta, rhoL, rhoM, nsim)# temp_tibble
        
        return(temp_tibble)
      }, error              = function(e){write(paste0(pE, "\t", SecondLayerL, "\t", as.character(i)),
                                                paste0("results/Error_withLatents_pE_", pE, "_typeCorM_", typeCorM, "_SecondLayerL_", SecondLayerL, "_typeCorrM_", typeCorM, "_withrelax.txt"), append=T)}
    )
  }
saveRDS(all_res, file= paste0("results/Res_withLatents_pE_", pE, "_SecondLayerL_", SecondLayerL, "_rhoM_", rhoM, "_typeCorrM_", typeCorM, "_sim", min(NSIM), "_", max(NSIM), ".rds"))
stopCluster(clusSim)