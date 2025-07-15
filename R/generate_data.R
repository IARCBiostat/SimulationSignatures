generate.data <- function(n, pM, pMrel,
                          SecondLayerM,
                          beta,
                          rhoM)
{
  # Scenario 3-4
  if(SecondLayerM)
  {
    pFirstLayer  <- 1:pMrel
    pSecondLayer <- (pMrel+1):(2*pMrel)
    # betaFSLayer  <- beta/2
  }
  
  # Get indicator vector of related features
  VecE1onM = c(rep(1, pMrel), rep(0, pM-pMrel))
  # Indicator vector * beta
  MatBeta  <- beta * matrix(VecE1onM, nrow=1)
  
  # Exposure drawn from normal distribution
  ZE   <- matrix(rnorm(n), nrow=1)
  E    <- t(ZE)
  
  # p x p matrix of metabolites
  ZM   <- matrix(rnorm(n*pM), nrow=pM)
  # Matrix of correlations between metabolites (0 or 0.5)
  SigM <- matrix(rhoM, nrow = pM, ncol=pM) + diag(1-rhoM, pM, pM)
  # Get Cholesky decomposition of M correlation structure
  CM   <- t(chol(SigM))
  
  prmM <- CM[lower.tri(CM, diag = TRUE)]
  ltM <- ltMatrices(matrix(prmM, ncol = 1L),
                    diag = TRUE, ### has diagonal elements
                    byrow = FALSE) ### prm is column-major
  
  # Give M correlation structure ltM + exposure*beta
  M    <- t(Mult(ltM, ZM)) + E%*%MatBeta
  
  # Scenarios 3-4 with second layer of metabolites related to the first
  if(SecondLayerM){
    # Second layer is related to first layer w: unif(0.1,0.5)
    MatBetaFirstSecondLayer <- matrix(runif(length(pSecondLayer)*length(pFirstLayer), min=0.2, max=1), nrow = length(pFirstLayer))
    # Second layer has correlation structure as defined above
    M[, pSecondLayer]       <- M[, pSecondLayer] + M[, pFirstLayer]%*% MatBetaFirstSecondLayer
  }
  # Scale features
  M <- scale(M)
  
  # Give features names "M_" + number of the feature
  colnames(M)    <- paste0("M_", 1:pM)
  
  # Returns a list with the features, exposure, and which rows of M are associated with E
  gendat <- list(
      M        = M,
      E        = E,
      VecE1onM = VecE1onM)
  return(gendat)
}