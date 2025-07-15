# Run lasso
lasso_and_co <- function(x, y)
{
  pM    <- ncol(x)
  cv.lasso              <- cv.glmnet(x, y, alpha = 1, gamma = 1, relax=T, maxp=min(pM, n/5))
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
