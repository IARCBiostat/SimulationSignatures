generate.data <- function(n, pM, pMrel,
                          SecondLayerM,
                          beta,
                          rhoM)
{
  # Scenario 1
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
    # Second layer is related to first layer w: unif(0.005,0.01)
    MatBetaFirstSecondLayer <- matrix(runif(length(pSecondLayer)*length(pFirstLayer), min=0.005, max=0.01), nrow = length(pFirstLayer))
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

generate_exposure <- function(n = 1000, beta = 0.6, sigma_eps = 0.5){
  #--------------------------------------------------
  # Exposure
  # W1 -> E
  #--------------------------------------------------
  W1  <- rnorm(n, 0, 1)
  E <- beta * W1 + rnorm(n, 0, sigma_eps)
  return(list(E=E, W1=W1))
}

generate_DAG_data <- function(n = 1000, 
                              beta = 0.6,       # generic edge strength
                              seed = 123,
                              sigma_eps = 0.5, 
                              M17_latent = FALSE, # noise sd
                              E,
                              W1) { 
  
  set.seed(seed)
  
  #--------------------------------------------------
  # Latent variables
  #--------------------------------------------------
  U1 <- rnorm(n, 0, 1)
  U2 <- rnorm(n, 0, 1)
  
  #--------------------------------------------------
  # Exogenous observed variables
  #--------------------------------------------------
  M2  <- rnorm(n, 0, 1)
  M7  <- rnorm(n, 0, 1)
  M5  <- rnorm(n, 0, 1)
  M15 <- rnorm(n, 0, 1)
  
  #--------------------------------------------------
  # "Simple" Children of E
  #--------------------------------------------------
  M1  <- beta * E + beta * M2 + rnorm(n, 0, sigma_eps)
  M11 <- beta * E + beta * W1 + beta * U2 + rnorm(n, 0, sigma_eps)
  
  #--------------------------------------------------
  # Descendants
  #--------------------------------------------------
  M3  <- beta * M1 + rnorm(n, 0, sigma_eps)
  M4  <- beta * M3 + beta * M5 + rnorm(n, 0, sigma_eps)
  M12 <- beta * M11 + rnorm(n, 0, sigma_eps)
  M13 <- beta * M12 + rnorm(n, 0, sigma_eps)
  
  #--------------------------------------------------
  # More complex children of E, and their descendants
  #--------------------------------------------------
  M6  <- beta * E + beta * M3 + beta * M7 + beta * U1 + rnorm(n, 0, sigma_eps)
  M8  <- beta * M6 + rnorm(n, 0, sigma_eps)
  
  #--------------------------------------------------
  # Spouses via latent variables
  #--------------------------------------------------
  M9  <- beta * U1 + rnorm(n, 0, sigma_eps)
  M14 <- beta * U2 + beta * M15 + rnorm(n, 0, sigma_eps)
  
  #--------------------------------------------------
  # Other non-descendants
  #--------------------------------------------------
  M16 <- beta * W1 + rnorm(n, 0, sigma_eps)
  M10 <- beta * M9 + rnorm(n, 0, sigma_eps)
  
  
  #--------------------------------------------------
  # Features linked to M17
  #--------------------------------------------------
  M19 <- rnorm(n, 0, 1)
  if(M17_latent){M17 <- beta * E + beta * M19 + rnorm(n, 0, sigma_eps)}else{M17 <- rnorm(n, 0, 1)}
  M18 <- beta * M17 + rnorm(n, 0, sigma_eps)
  
  #--------------------------------------------------
  # Return observed data only (latents removed)
  #--------------------------------------------------
  datagen <- data.frame(
    M1, M2, M3, M4, M5,
    M6, M7, M8, M9, M10,
    M11, M12, M13, M14, M15, M16, 
    M18, M19
  )
  
  if (!M17_latent)
  {
    children          <- c("M1","M6","M11")
    descendants       <- c("M3","M4","M8","M12","M13")
    nondesc           <- c("M2","M5","M7","M9","M10","M14","M15","M16","M19")
    sign_rest         <- c("M1","M3","M6","M11")
    sign_full         <- c("M1","M3","M6","M11","M2","M7","M9","M14","M15")
  } else
  {
    children          <- c("M1","M6","M11","M18")
    descendants       <- c("M3","M4","M8","M12","M13")
    nondesc           <- c("M2","M5","M7","M9","M10","M14","M15","M16","M19")
    sign_rest         <- c("M1","M3","M6","M11","M18")
    sign_full         <- c("M1","M3","M6","M11","M2","M7","M9","M14","M15", "M18","M19")
  } 
  
  return(list(datagen = datagen, 
              children = children, 
              descendants = descendants, 
              nondesc     = nondesc, 
              sign_rest   = sign_rest, 
              sign_full   = sign_full))
}

generate_multiblock <- function(blocks = 50,
                                n = 1000, 
                                beta = 0.6,       # generic edge strength
                                seed = 123,
                                sigma_eps = 0.5, 
                                M17_latent = FALSE) {
  all_datagen   <- list()
  all_children  <- c()
  all_desc      <- c()
  all_nondesc   <- c()
  all_rest <- c()
  all_full <- c()
  
  # Generate initial exposure and confounder W1
  init_E_W <- generate_exposure(
    n = n, 
    beta = beta,
    sigma_eps = sigma_eps)
  
  old_Mnames <- c(1:16, 18:19)
  for (i in 1:blocks){
    # Generate one iteration of block
    tmp_blk <- generate_DAG_data(
      n = n,
      beta = beta,
      seed = seed + i,
      sigma_eps = sigma_eps,
      M17_latent = M17_latent,
      E = init_E_W$E,
      W1 = init_E_W$W1
    )
    
    # Calculate shift for metabolite numbers
    M_count <- 19
    shift <- (i-1)*M_count
    
    # New names for df
    new_Mnames <- paste0("M", old_Mnames +shift)
    colnames(tmp_blk$datagen) <- new_Mnames
    
    # Function to shift children names
    shift_vec <- function(v)
      paste0("M", as.numeric(sub("M","", v)) + shift)
    
    # Add identity for all nodes
    all_children <- c(all_children, shift_vec(tmp_blk$children))
    all_desc     <- c(all_desc,     shift_vec(tmp_blk$descendants))
    all_nondesc  <- c(all_nondesc,     shift_vec(tmp_blk$nondesc))
    all_rest     <- c(all_rest,     shift_vec(tmp_blk$sign_rest))
    all_full     <- c(all_full,     shift_vec(tmp_blk$sign_full))
    
    all_datagen[[i]] <- tmp_blk$datagen
  }
  all_datagen_fin <- cbind(E = init_E_W$E, W1 = init_E_W$W1, do.call(cbind, all_datagen))
  return(list(datagen     = all_datagen_fin, 
              children    = all_children, 
              descendants = all_desc, 
              nondesc     = all_nondesc, 
              sign_rest   = all_rest, 
              sign_full   = all_full))
}