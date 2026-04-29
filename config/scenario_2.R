# Config for scenario 2

# ---------------------------------
# Parameters
# ---------------------------------
base <- expand.grid(
  n = 1000,
  P = 1000,
  rho_g = c(0,0.2),
  rho_l = 0.5,
  beta = c(1:5)/10,
  iteration = 1:100,
  stringsAsFactors = FALSE
)

# helper to add specific combos
add_params <- function(prel, nmother, nstepbro) {
  df <- expand.grid(
    prel = prel,
    Bhiddenmothers = nmother,
    Bstepbrothers = nstepbro
  )
  merge(base, df)
}

params <- rbind(
  add_params(5, 25, 7),
  add_params(5, 7, 25),
  add_params(25, 5, 7),
  add_params(25, 7, 5),
  add_params(125, 2, 3),
  add_params(125, 3, 2),
  add_params(25, 1:5, 1:7)
)

i <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

if (is.na(i) || i < 1 || i > nrow(params)) {
  stop("Invalid SLURM_ARRAY_TASK_ID: ", i)
}

n = params$n[i]
P = params$P[i]
prel = params$prel[i]
rho_g = params$rho_g[i]
rho_l = params$rho_l[i]
beta = params$beta[i]
Bhiddenmothers = params$Bhiddenmothers[i]
Bstepbrothers = params$Bstepbrothers[i]
iteration=params$iteration[i]
