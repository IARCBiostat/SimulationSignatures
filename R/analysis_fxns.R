# Run lasso
lasso_and_co <- function(x, y, n=1000, ...)
{
  pM    <- ncol(x)
  cv.lasso              <- cv.glmnet(x, y, alpha = 1, gamma = 1, relax=T, maxp=min(pM, n/5), ...)
  coef.lasso.min        <- coef(cv.lasso, s="lambda.min")[-1][1:pM]
  coef.lasso.1se        <- coef(cv.lasso, s="lambda.1se")[-1][1:pM]
  
  all_coefs <- list(lasso.min          = coef.lasso.min,
                    lasso.1se          = coef.lasso.1se)
  
  return(all_coefs)
}

# Estimate lasso results with and without screening
estimate_coef_sig_w_wo_prior_screening <- function(dat.train)
{
  M.train      <- dat.train$M
  E.train      <- dat.train$E
  y.train      <- E.train[, 1]
  pM           <- ncol(M.train)
  pE           <- ncol(E.train)
  
  # Screening
  p.val.screen.univar <- sapply(1:pM, function(j){
    summary(lm(M.train[,j]~E.train))$coef[2, "Pr(>|t|)"]
  })
  fdr.screen.univar <- p.adjust(p.val.screen.univar, method = "fdr")
  to.keep           <- order(fdr.screen.univar)[1:max(length(which(fdr.screen.univar <0.05)), 3)]
  M.train.screened  <- M.train[, to.keep]
  
  x.train                 <- M.train
  x.train.screened        <- M.train.screened
  
  ##### now run the lasso and other approaches
  all_coefs_noscreen      <- lasso_and_co(x = x.train, y = y.train)
  all_coefs_screen0       <- lasso_and_co(x = x.train.screened, y = y.train)
  
  all_coefs_screen        <- c(
    lapply(all_coefs_screen0[which(grepl("res", names(all_coefs_screen0)))], function(vec){
      vec1 <- rep(0, pM)
      vec1[to.keep.res]<- vec
      return(vec1)}),
    lapply(all_coefs_screen0[which(!grepl("res", names(all_coefs_screen0)))], function(vec){
      vec1 <- rep(0, pM)
      vec1[to.keep]<- vec
      return(vec1)}))
  
  names(all_coefs_noscreen) <- paste0(names(all_coefs_noscreen), "_noscreen")
  names(all_coefs_screen) <- paste0(names(all_coefs_screen), "_screen")
  all_coefs               <- c(all_coefs_screen, all_coefs_noscreen)
  return(all_coefs)
}

cor_sig_exp_test <- function(coef_sig, dat.test)
{
  pM       <- ncol(dat.test$M)
  Sig_test <- dat.test$M %*%coef_sig[1:pM]
  return(cor(Sig_test, dat.test$E[,1]))
}


## Returns the names of features passing BH FDR <= q
screen_univariate <- function(dat, q = 0.05, slow=F) {
  all_feats <- setdiff(colnames(dat), c("E", "W1"))
  pvals <- numeric(length(all_feats))
  names(pvals) <- all_feats
  if(slow){
    for (m in all_feats) {
      fit <- lm(E ~ W1 + dat[[m]], data = dat)
      # p-value for the feature coefficient
      pvals[m] <- summary(fit)$coefficients[3, 4]
    }
  } else{
    # Base variables
    W1 <- dat$W1
    Y  <- dat$E
    X  <- as.matrix(dat[all_feats])
    
    # Residualize Y on W1
    fit_y <- lm.fit(cbind(1, W1), Y)
    Y_res <- fit_y$residuals
    
    # Residualize ALL features on W1 at once
    fit_x <- lm.fit(cbind(1, W1), X)
    X_res <- fit_x$residuals   # matrix: n × p
    
    # Compute coefficients (vectorized)
    beta <- colSums(X_res * Y_res) / colSums(X_res^2)
    
    # Residual variance per feature
    n <- nrow(dat)
    df <- n - 3
    
    residuals_mat <- sweep(X_res, 2, beta, "*")
    residuals_mat <- Y_res - residuals_mat
    
    sigma2 <- colSums(residuals_mat^2) / df
    
    # Standard errors
    se <- sqrt(sigma2 / colSums(X_res^2))
    
    # t-stats and p-values
    tstat <- beta / se
    pvals <- 2 * pt(-abs(tstat), df = df)
  }
  qvals <- p.adjust(pvals, method = "BH")
  kept  <- names(qvals)[qvals <= q]
  return(list(kept = kept, p = pvals, q = qvals))
}
generate_data_scen2 <- function(n,P,prel,rho_g,rho_l,beta,Bhiddenmothers,Bstepbrothers=5){
  #Bstepbrothers   <- floor((P - prel) / (prel*Bhiddenmothers)) ### number of descendants per hiddenmother
  
  P_used <- prel + (prel*Bhiddenmothers)*Bstepbrothers
  VecE1onM = c(rep(1, prel), rep(0, P-prel))
  
  cat("Variables used:", P_used, " / ", P, "\n")
  
  # ---------------------------------
  # Latent variables
  # ---------------------------------
  E  <- rnorm(n)                 # exposure
  F0 <- rnorm(n)                 # global hidden mother
  F_local <- matrix(rnorm(n * Bhiddenmothers*prel), n, Bhiddenmothers*prel) # local hidden mothers
  
  # ---------------------------------
  # Initialize observed matrix
  # ---------------------------------
  M <- matrix(0, nrow = n, ncol = P)
  
  # ---------------------------------
  # Primary variables (1 ... prel)
  # ---------------------------------
  for (j in 1:prel) {
    if(Bhiddenmothers != 1){
      M[, j] <- sqrt(rho_g) * F0 +
        sqrt(rho_l) * rowSums(F_local[, ((j-1)*Bhiddenmothers + 1):(j*Bhiddenmothers)]) +
        beta * E
    } else{
      M[, j] <- sqrt(rho_g) * F0 +
        sqrt(rho_l) * F_local[, ((j-1)*Bhiddenmothers + 1):(j*Bhiddenmothers)] +
        beta * E
    }}
  
  # ---------------------------------
  # Secondary blocks
  # ---------------------------------
  for (j in 1:prel) for (b in 1:Bhiddenmothers){
    k     <- (j-1)*Bhiddenmothers + b
    start <- prel + (k - 1) * Bstepbrothers + 1
    end   <- prel + k * Bstepbrothers
    M[, start:end] <- sqrt(rho_l) * F_local[, k]
  }
  
  # ---------------------------------
  # Residual variances (ensure Var = 1)
  # ---------------------------------
  # contributions per variable:
  # primary: rho_g + rho_l + beta^2
  # secondary: rho_l
  # resid_var <- rep(NA, P)
  # resid_var[1:prel] <- 1 - (rho_g + rho_l + beta^2)
  # resid_var[(prel + 1):P_used] <- 1 - rho_l
  # 
  # # safety for numerical tolerance
  # resid_var <- pmax(resid_var, 0)
  
  # ---------------------------------
  # Add residual noise
  # ---------------------------------
  for (p in 1:P_used) {
    #if (resid_var[p] > 0) {
    M[, p] <- M[, p] + rnorm(n, sd = 0.5 ) #sqrt(resid_var[p])
    #}
  }
  
  # ---------------------------------
  # Optional: remaining variables (if P_used < P)
  # ---------------------------------
  if (P_used < P) {
    M[, (P_used + 1):P] <- matrix(rnorm(n * (P - P_used)), n)
  }
  
  
  # Returns a list with the features, exposure, and which rows of M are associated with E
  gendat <- list(
    M        = M,
    E        = matrix(E, n),
    VecE1onM = VecE1onM)
  return(gendat)
}

estimate_coef_sig <- function(dat.train, n=1000, ...)
{
  M.train      <- dat.train$M
  E.train      <- dat.train$E
  y.train      <- E.train[, 1]
  pM           <- ncol(M.train)
  pE           <- ncol(E.train)
  
  # Screening
  p.val.screen.univar <- fast_pvals(M.train, E.train)
  fdr.screen.univar <- p.adjust(p.val.screen.univar, method = "fdr")
  to.keep           <- order(fdr.screen.univar)[1:max(length(which(fdr.screen.univar <0.05)), 3)]
  M.train.screened  <- M.train[, to.keep]
  
  x.train                 <- M.train
  x.train.screened        <- M.train.screened
  
  ##### now run the lasso and other approaches
  all_coefs_noscreen      <- lasso_and_co(x = x.train, y = y.train, n = n, ...)
  all_coefs_screen0       <- lasso_and_co(x = x.train.screened, y = y.train, ...)
  
  all_coefs_screen        <- c(
    lapply(all_coefs_screen0[which(!grepl("res", names(all_coefs_screen0)))], function(vec){
      vec1 <- rep(0, pM)
      vec1[to.keep]<- vec
      return(vec1)}))
  
  names(all_coefs_noscreen) <- paste0(names(all_coefs_noscreen), "_noscreen")
  names(all_coefs_screen) <- paste0(names(all_coefs_screen), "_screen")
  all_coefs               <- c(all_coefs_screen, all_coefs_noscreen)
  return(all_coefs)
}

run_scenario <- function(dat, dat.test, n, P, prel, rho_g, rho_l, beta, Bhiddenmothers, Bstepbrothers, iteration){
  
  screened <- fast_pvals(dat$M, dat$E)
  fdr <- p.adjust(screened, method = "fdr")
  nfeat_screen <- sum(fdr < 0.05)
  
  all_coefs <- estimate_coef_sig(dat, n = 1000)
  all_coefs <- all_coefs[c("lasso.1se_screen", "lasso.1se_noscreen")] # Choose 1se only
  
  VecE1onM  <- dat$VecE1onM
  
  temp_tibble <- tibble(
    Approach = names(all_coefs),
    CardS    = sapply(all_coefs, function(x) sum(x != 0)),
    Sensi    = sapply(all_coefs, function(x) sum(x != 0 & VecE1onM != 0) / sum(VecE1onM != 0)),
    Speci    = sapply(all_coefs, function(x) sum(x == 0 & VecE1onM == 0) / sum(VecE1onM == 0)),
    Corr     = sapply(all_coefs, function(x) cor_sig_exp_test(x, dat.test)),
    Screen_child = sum((fdr < 0.05) & (VecE1onM == 1)),
    Screen_nonchild = sum((fdr < 0.05) & (VecE1onM == 0)),
    n, P, prel, rho_g, rho_l, beta, Bhiddenmothers, Bstepbrothers, iteration
  )
  return(temp_tibble)
}

fast_pvals <- function(M, E) {
  n <- nrow(M)
  
  # center
  M_centered <- scale(M, center = TRUE, scale = FALSE)
  E_centered <- as.numeric(scale(E, center = TRUE, scale = FALSE))
  
  # correlation
  cor_vals <- colSums(M_centered*E_centered) /
    sqrt(colSums(M_centered^2) * sum(E_centered^2))
  
  # t-stat
  t_stats <- cor_vals * sqrt((n - 2) / (1 - cor_vals^2))
  
  # p-values
  2 * pt(-abs(t_stats), df = n - 2)
}

## ---- 3) LASSO with W1 unpenalized
## x: matrix of features (columns = subset of all_feats + W1)
## y: E
## unpenalized: character vector of column names to set penalty = 0 (here "W1")
## returns vector of selected feature names (nonzero at lambda.1se), excluding intercept and W1
fit_lasso_unpenalized_W1 <- function(x, y, unpenalized = "W1",
                                     use_lambda = c("lambda.1se","lambda.min"), 
                                     x.test = NULL, y.test=NULL) {
  use_lambda <- match.arg(use_lambda)
  stopifnot(unpenalized %in% colnames(x))
  pen <- rep(1, ncol(x))
  pen[match(unpenalized, colnames(x))] <- 0
  cv <- cv.glmnet(x = x, y = y,
                  alpha = 1, family = "gaussian",
                  standardize = TRUE, intercept = TRUE,
                  penalty.factor = pen, nfolds = 10)
  lam <- switch(use_lambda, lambda.1se = cv$lambda.1se, lambda.min = cv$lambda.min)
  b   <- coef(cv, s = lam)
  # extract nonzero coefficients, drop intercept and W1
  nm  <- rownames(b)
  sel_idx <- which(b[,1] != 0)
  sel <- nm[sel_idx]
  sel <- setdiff(sel, c("(Intercept)", unpenalized))
  if(is.null(x.test) |is.null(y.test)){corr.test.all <- corr.test.noW <- NA}else
  {
    corr.test.all = as.numeric(cor(y.test, predict(cv, newx = x.test, s=lam)))
    x.testnoW     <- x.test[, -which(colnames(x.test)=="W1")]
    b.noW         <- b[-which(rownames(b)=="W1")]
    print(length(b.noW))
    print(length(y.test))
    print(ncol(x.testnoW))
    if(length(b.noW) == 0){
      corr.test.noW <- 0
    } else{corr.test.noW <- as.numeric(cor(y.test, b[1] + x.testnoW%*%b.noW[-1]))}
  }
  return(list(selected = sel, cv = cv, lambda = lam, corr.test.noW = corr.test.noW, corr.test.all = corr.test.all))
}

## ---- 4) One replication of both pipelines
## Returns selection sets and metrics
one_replication <- function(blocks = 50, n = 500, beta = 0.5, sigma_eps = 0.5, M17_latent=TRUE, seed = 1234, q_screen = 0.05, use_lambda = "lambda.1se") {
  
  # seed        <- seed + 1
  temp        <- generate_multiblock(blocks = blocks, n = n, beta = beta, seed = seed, sigma_eps = sigma_eps, M17_latent=M17_latent)
  temptest    <- generate_multiblock(blocks = blocks, n = 10000, beta = beta, seed = seed, sigma_eps = sigma_eps, M17_latent=M17_latent)
  
  dat         <- temp$datagen
  y           <- dat$E
  
  all_feats   <- setdiff(colnames(temptest$datagen), c("E", "W1"))
  x.test      <- as.matrix(temptest$datagen[, c(all_feats, "W1")])
  y.test      <- as.matrix(temptest$datagen$E)
  
  
  # oo <- lm(y~ dat$M19 + dat$W1); summary(oo)
  
  children    <- temp$children
  descendants <- temp$descendants
  nondesc     <- temp$nondesc
  sign_rest   <- temp$sign_rest
  sign_full   <- temp$sign_full
  
  ## LASSO without screening (all features + W1)
  X_all <- as.matrix(dat[, c(all_feats, "W1")])
  fitA  <- fit_lasso_unpenalized_W1(X_all, y, unpenalized = "W1", use_lambda = use_lambda, x.test = x.test, y.test = y.test)
  #selA  <- fitA$selected
  #corr.test.noW.A <- fitA$corr.test.noW
  #corr.test.all.A <- fitA$corr.test.all
  
  ## Screening (adjusted for W1), then LASSO on kept + W1
  sc    <- screen_univariate(dat, q = q_screen)
  kept  <- sc$kept
  if (length(kept) == 0) {
    selB <- character(0)
    corr.test.noW.B <- 0
    corr.test.all.B <- 0
  } else {
    X_scr <- as.matrix(dat[, c(kept, "W1")])
    fitB  <- fit_lasso_unpenalized_W1(x = X_scr, y = y, unpenalized = "W1", use_lambda = use_lambda, x.test = x.test[, c(kept, "W1")], y.test = y.test)
    #selB  <- fitB$selected
    #corr.test.noW.B <- fitB$corr.test.noW
    #corr.test.all.B <- fitB$corr.test.all
  }
  
  ## Metrics for each method
  eval_metrics <- function(fit) {
    sel           <- fit$selected
    corr.test.noW <- fit$corr.test.noW
    corr.test.all <- fit$corr.test.all
    tp_child     <- sum(sel %in% children)
    sel_indir    <- sum(sel %in% descendants)
    sel_nond     <- sum(sel %in% nondesc)
    sel_full_sgn <- sum(sel %in% sign_full)
    sel_rest_sgn <- sum(sel %in% sign_rest)
    
    n_selected = length(sel)
    tp_children = tp_child
    recall_children = tp_child / length(children)
    #      prec_children = ifelse(n_selected == 0, NA, tp_child / n_selected)
    selected_indirect = sel_indir
    selected_nondesc  = sel_nond
    tp_fullsign = sel_full_sgn
    recall_fullsign = sel_full_sgn / length(sign_full)
    prec_fullsign = ifelse(n_selected == 0, NA, sel_full_sgn / n_selected) 
    
    tp_restsign = sel_rest_sgn
    recall_restsign = sel_rest_sgn / length(sign_rest)
    prec_restsign = ifelse(n_selected == 0, NA, sel_rest_sgn / n_selected) 
    
    is_restsign   =ifelse(recall_restsign==1 & prec_restsign==1, 1, 0)
    is_fullsign   =ifelse(recall_fullsign==1 & prec_fullsign==1, 1, 0)
    list(
      n_selected = n_selected, 
      tp_children = tp_children, 
      recall_children = recall_children, 
      selected_indirect = selected_indirect, 
      selected_nondesc  = selected_nondesc, 
      tp_fullsign = tp_fullsign, 
      recall_fullsign = recall_fullsign, 
      prec_fullsign = prec_fullsign,
      tp_restsign = tp_restsign, 
      recall_restsign = recall_restsign, 
      prec_restsign = prec_restsign,
      is_restsign  = is_restsign, 
      is_fullsign = is_fullsign, 
      corr.test.noW = corr.test.noW, 
      corr.test.all = corr.test.all
    )
  }
  mA <- eval_metrics(fitA)
  mB <- eval_metrics(fitB)
  
  out <- list(sel_no_screen = fitA$sel, sel_with_screen = fitB$sel,
              metrics_no_screen = mA, metrics_with_screen = mB,
              kept_screen = kept)
  
  return(list(out = out, 
              children = children, 
              descendants = descendants, 
              nondesc     = nondesc))
}

## ---- 5) Run many replications and summarize
run_simulation <- function(blocks = 50, B = 5, beta=0.5, sigma_eps = 0.5, M17_latent=TRUE, n = 500, q_screen = 0.05, use_lambda = "lambda.1se", seed0 = 1234) {
  set.seed(seed0)
  
  tmp_dat <- generate_multiblock(blocks = blocks, n = n, beta = beta, sigma_eps = sigma_eps, M17_latent = M17_latent, seed = seed0)
  all_feats <- setdiff(colnames(tmp_dat$datagen), c("E", "W1"))
  # selection frequencies by feature
  freqA <- setNames(rep(0, length(all_feats)), all_feats)
  freqB <- setNames(rep(0, length(all_feats)), all_feats)
  
  # metrics accumulators
  acc <- data.frame(
    method = character(),
    n_selected = integer(),
    tp_children = integer(),
    recall_children = numeric(),
    prec_children = numeric(),
    selected_indirect = integer(),
    selected_nondesc = integer(),
    tp_fullsign = integer(),
    recall_fullsign = numeric(),
    prec_fullsign = numeric(), 
    is_children = integer(), 
    is_fullsign = integer(),
    corr.test.noW.full = numeric(), 
    corr.test.all.full = numeric(), 
    corr.test.noW.rest = numeric(), 
    corr.test.all.rest = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (b in 1:B) {
    seedy <- sample.int(1e7, 1)
    print(sprintf("...Starting iteration %s, seed=%s", b, seedy))
    temp <- one_replication(blocks = blocks, n = n, beta = beta, M17_latent = M17_latent, seed = sample.int(1e7, 1), q_screen = q_screen, use_lambda = use_lambda)
    res  <- temp$out
    children    <- temp$children
    descendants <- temp$descendants
    nondesc     <- temp$nondesc
    
    
    # frequencies
    if (length(res$sel_no_screen))  freqA[res$sel_no_screen]  <- freqA[res$sel_no_screen]  + 1
    if (length(res$sel_with_screen)) freqB[res$sel_with_screen] <- freqB[res$sel_with_screen] + 1
    
    # metrics
    acc <- rbind(acc,
                 data.frame(method = "no_screen",  t(unlist(res$metrics_no_screen)),  row.names = NULL),
                 data.frame(method = "with_screen", t(unlist(res$metrics_with_screen)), row.names = NULL))
    print(sprintf("Finished iteration %s", b))
  }
  
  # average metrics
  summary_metrics <- aggregate(. ~ method, data = acc, FUN = function(x) mean(x, na.rm = TRUE))
  rownames(summary_metrics) <- NULL
  
  # selection frequency tables
  freq_tab <- data.frame(
    feature = all_feats,
    sel_freq_no_screen  = freqA / B,
    sel_freq_with_screen = freqB / B,
    is_child     = all_feats %in% children,
    is_desc      = all_feats %in% descendants,
    is_nondesc   = all_feats %in% nondesc,
    stringsAsFactors = FALSE
  )
  
  list(summary_metrics = summary_metrics, freq = freq_tab, runs = acc)
}
