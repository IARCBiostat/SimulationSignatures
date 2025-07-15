# Config for: No second layer M, rhoM = 0, 0.5

NCore <- 40
pM    <- 1000
NSIM  <- 1:100
N     <- c(1000)
PMRELS <- c(5, 25, 125)
BETAS  <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rhoM   <- c(0, 0.5)
SecondLayerM <- FALSE

# NCore <- 12
# pM    <- 100
# NSIM  <- 1:10
# N     <- c(100)
# PMRELS <- c(5, 25)
# BETAS  <- c(0.1, 0.3, 0.5)
# rhoM   <- c(0, 0.5)
# SecondLayerM <- FALSE