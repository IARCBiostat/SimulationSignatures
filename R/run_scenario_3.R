rm(list = ls())
library(tidyverse)
library(glmnet)
print(getwd())
source("config/scenario_3.R")
source("R/generate_data.R")
source("R/analysis_fxns.R")

AllParams  <- expand.grid(NREP, NSIM, N, BETAS, M17_latent)
colnames(AllParams) <- c("block", "nsim", "n", "beta", "M17_latent")


t1 <- Sys.time()
# Need to fix from here on----
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
task <- AllParams[task_id,]

cat("Beta:", task$beta, "\n",
    "n:", task$n, "\n",
    "N blocks:", task$block, "\n")

cat(sprintf(">>> Task %d: beta=%s, n=%s, block=%s\n",
            task_id, task$beta, task$n, task$block))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

sim_res <- run_simulation(blocks = task$block, 
                          B = task$nsim, 
                          beta=task$beta, 
                          M17_latent=task$M17_latent, 
                          n = task$n, 
                          q_screen = 0.05, 
                          use_lambda = "lambda.1se", 
                          seed0 = 19+task_id)

saveRDS(sim_res, file=paste0("results/Res_ToyEx_General_n_", task$n, 
                             "_beta_", task$beta,
                             "_M17_latent_" , task$M17_latent, 
                             "_block_", task$block ,".rds"))

t2 <- Sys.time()
cat(sprintf(">>> Completed in %s\n", format(t2 - t1)))

cat("\n=== Average metrics over replications ===\n")
print(sim_res$summary_metrics, digits = 3)

cat("\n=== Selection frequency by feature (first 16 rows) ===\n")
print(sim_res$freq[order(-sim_res$freq$sel_freq_no_screen), ], row.names = FALSE)

