# Config for toy example, expanded
NCore <- 2
NREP    <- c(55, 110)
NSIM  <- 100
N     <- c(500, 1000, 2500, 12500, 5*12500, 25*12500)
BETAS  <- seq(0.1, 0.5, 0.1)
M17_latent=TRUE

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