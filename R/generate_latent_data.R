generate.data.withlatent <- function(n, pL, pM, pLrel,
                                     SecondLayerL, typeCorM,
                                     beta,
                                     rhoL,
                                     rhoM)
{
  # Generate indicator of E ~ L where E affects pLrel latent vars
  VecE1onL = c(rep(1, pLrel), rep(0, pL-pLrel))
  
  # Indicator with values
  MatBeta  <- beta * matrix(VecE1onL, nrow=1)
  
  # Random normal exposure
  ZE   <- matrix(rnorm(n), nrow=1)
  E    <- t(ZE)
  
  # Random normal latent variables
  ZL   <- matrix(rnorm(n*pL), nrow=pL)
  SigL <- matrix(rhoL, nrow = pL, ncol=pL) + diag(1-rhoL, pL, pL)
  CL   <- t(chol(SigL))
  prmL <- CL[lower.tri(CL, diag = TRUE)]
  ltL <- ltMatrices(matrix(prmL, ncol = 1L),
                    diag = TRUE, ### has diagonal elements
                    byrow = FALSE) ### prm is column-major
  # Latent variable with correlation structure ltL and relationship to exposure
  L    <- t(Mult(ltL, ZL)) + E%*%MatBeta
  
  if(SecondLayerL){
    betaFSLayer  <- beta/2
    for (k in 1:floor(pL/pLrel))
    {
      pFirstLayer  <- unique((k-1)*pLrel + c(1:pLrel))
      pSecondLayer <- (3*pLrel+1):(2*pLrel)
      MatBetaFirstSecondLayer <- betaFSLayer*matrix(runif(length(pSecondLayer)*length(pFirstLayer), min=0.1, max=0.5), nrow = length(pFirstLayer))
    }
    
    L[, pSecondLayer]       <- L[, pSecondLayer] + L[, pFirstLayer]%*% MatBetaFirstSecondLayer
  }
  
  #### Now that we have the Ls, we generate the Ms.
  # Create a grouping index
  groupsM <- cut(seq_along(1:pM), breaks=pL, labels=FALSE)
  # Split the vector into chunks
  chunksM <- split(1:pM, groupsM)
  ind_M_rel <- vector("list")
  M         <- matrix(NA, ncol=pM, nrow=n)
  
  # Metabolites within latent variables are correlated
  if (typeCorM =="amongL")
  {
    # For each latent variable
    for (ltmp in 1:pL) 
    {
      # Get metabolite indices related to latent variable "ltmp"
      ind_mtmp   <- chunksM[[ltmp]]
      # Get number of metabolites related
      pMtmp      <- length(ind_mtmp)
      # If it's a latent variable related to the exposure, save index of related M
      if (ltmp %in% 1:pLrel){ind_M_rel[[ltmp]] <- ind_mtmp}
      # Make vector of related metabolites
      MatBetatmp  <- matrix(2, nrow=1, ncol=pMtmp)
      
      # Multivariate normal residuals
      ZMtmp   <- matrix(rnorm(n*pMtmp), nrow=pMtmp)
      # For this set of metabolites, make the correlation matrix
      SigMtmp <- matrix(rhoM, nrow = pMtmp, ncol=pMtmp) + diag(1-rhoM, pMtmp, pMtmp)
      # Cholesky decomp of correlation matrix
      CMtmp   <- t(chol(SigMtmp))
      prmMtmp <- CMtmp[lower.tri(CMtmp, diag = TRUE)]
      ltMtmp  <- ltMatrices(matrix(prmMtmp, ncol = 1L),
                            diag = TRUE, ### has diagonal elements
                            byrow = FALSE) ### prm is column-major
      # Give residuals correlation structure of ltMtmp + latent variable effect
      Mtmp    <- t(Mult(ltMtmp, ZMtmp)) + L[, ltmp]%*%MatBetatmp
      M[, ind_mtmp] <- Mtmp
    }
  }else{
    # Correlation structure is across all metabolites
    if (typeCorM=="allM")
    {
      # Make M with correlation structure rhoM across all metabolites
      ZM   <- matrix(rnorm(n*pM), nrow=pM)
      SigM <- matrix(rhoM, nrow = pM, ncol=pM) + diag(1-rhoM, pM, pM)
      CM   <- t(chol(SigM))
      prmM <- CM[lower.tri(CM, diag = TRUE)]
      ltM  <- ltMatrices(matrix(prmM, ncol = 1L),
                         diag = TRUE, ### has diagonal elements
                         byrow = FALSE) ### prm is column-major
      M    <- t(Mult(ltM, ZM))
      for (ltmp in 1:pL) # Add back effect of latent variable for each metabolite
      {
        ind_mtmp   <- chunksM[[ltmp]]
        pMtmp      <- length(ind_mtmp)
        if (ltmp %in% 1:pLrel){ind_M_rel[[ltmp]] <- ind_mtmp}
        MatBetatmp  <- matrix(2, nrow=1, ncol=pMtmp)
        
        M[, ind_mtmp] <- M[, ind_mtmp]  + L[, ltmp]%*%MatBetatmp
      }
    }
  }
  M <- scale(M)
  colnames(M)    <- paste0("M_", 1:pM)

  gendat <- list(
      M         = M,
      E         = E,
      VecE1onL  = VecE1onL,
      ind_M_rel = ind_M_rel) # Which metabolites are related to the exposure through L
  
  return(gendat)
}
