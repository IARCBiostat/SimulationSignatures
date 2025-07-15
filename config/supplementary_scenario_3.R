# Config for: Latent variable with no correlations among M

NCore <- 40
pE    <- 1
PLRELS <- 2
PLS    <- 10
SecondLayerL <- FALSE

rhoM  <- 0
pM    <- 1000
NSIM  <- 1:100
N     <- c(1000)

BETAS  <- c(0.1, 0.2, 0.3, 0.4, 0.5)
RHOLS   <- c(0, 0.5)

typeCorM <- "allM"

## Toy scenario for testing
# NCore <- 2
# pE    <- 1
# PLRELS <- 2
# PLS    <- 10
# SecondLayerL <- FALSE
# 
# rhoM  <- 0
# pM    <- 100
# NSIM  <- 1:5
# N     <- c(1000)
# 
# BETAS  <- c(0.1, 0.2)
# RHOLS   <- c(0, 0.5)
# 
# typeCorM <- "allM"
