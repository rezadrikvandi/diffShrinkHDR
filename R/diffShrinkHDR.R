#' High Dimensional Regression with Different Shrinkage
#' @param X Design matrix of predictors (n x p)
#' @param Y Response vector of length n
#' @param Index_XI Optional vector of indices for parameters of interest. If NULL, they will be selected automatically
#' @param alpha Significance level for confidence intervals and hypothesis tests (e.g., 0.05)

diffShrinkHDR <- function(X,Y,Index_XI=NULL,alpha){
  n <- nrow(X)
  p <- ncol(X)
  if(is.null(Index_XI))
  {
   quiet <- function(x)
   {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
   }
   sis_res <- quiet(SIS(X, Y, family='gaussian', iter=TRUE))
   Index_XI <- sis_res$sis.ix0
   Index_XI <- sis_res$ix
  }
  XI <- X[,Index_XI]
  XN <- X[,-Index_XI]
  XI <- matrix(c(XI),nrow(X),length(Index_XI))
  XN <- matrix(c(XN),nrow(X),p-length(Index_XI))

  #CV for lambda1 and lambda2
  CV_prediction <- function(X,Y,X_test,Y_test,Index_XI,lambda_1,lambda_2)
  {
    n <- nrow(X)
    p <- ncol(X)
    XI <- X[,Index_XI]
    XN <- X[,-Index_XI]
    XI <- matrix(c(XI),n,length(Index_XI))
    XN <- matrix(c(XN),n,p-length(Index_XI))
    PI <- XI%*%ginv(t(XI)%*%XI+n*lambda_1*diag(ncol(XI)))%*%t(XI)
    Y_star <- sqrtm(diag(n)-PI)%*%Y
    XN_star <- sqrtm(diag(n)-PI)%*%XN
    fit_lasso_betaN <- glmnet(XN_star, Y_star, alpha = 1, family = "gaussian", intercept = FALSE, standardize = TRUE)
    betaN_hat <- as.vector(coef(fit_lasso_betaN, s=lambda_2)[-1])

    betaI_hat <- ginv(t(XI)%*%XI+n*lambda_1*diag(ncol(XI)))%*%t(XI)%*%(Y-XN%*%betaN_hat)
    betaI_hat <- as.vector(betaI_hat)

    X_test_I <- X_test[,Index_XI]
    X_test_N <- X_test[,-Index_XI]
    Y_test_hat <- X_test_I%*%betaI_hat+X_test_N%*%betaN_hat
    MSPE_CV <- mean((Y_test - Y_test_hat)^2)
    return(MSPE_CV)
  }

  grid_lambda1 <- 10^seq(-1.0,-1.4,by=-0.05)
  grid_lambda2 <- 10^seq(-0.7,-1.1,by=-0.05)
  # 10-fold cross validation
  f_MSPE <- function(i,lambda1,lambda2,rand_shuffl,folds)
  {
    X.CV <- X[rand_shuffl,]
    Y.CV <- Y[rand_shuffl]
    #Segment the data by fold using the which() function;
    testIndexes <- which(folds==i,arr.ind=TRUE)
    XCV_train <- X.CV[-testIndexes, ]
    XCV_valid <- X.CV[testIndexes, ]
    YCV_train <-  Y.CV[-testIndexes]
    YCV_valid <-  Y.CV[testIndexes]
    MSPE_CV <- CV_prediction(X=XCV_train,Y=YCV_train,X_test=XCV_valid,Y_test=YCV_valid,Index_XI=Index_XI,lambda_1=lambda1,lambda_2=lambda2)
    return(MSPE_CV)
  }
  CV_pred_results <- NA
  for(lambda1 in grid_lambda1)
  {
    for(lambda2 in grid_lambda2)
    {
      rand_shuffl <- sample(1:nrow(X))
      X.CV <- X[rand_shuffl,]
      Y.CV <- Y[rand_shuffl]
      #Create 10 equally-sized folds;
      folds <- cut(seq(1,nrow(X.CV)),breaks=10,labels=FALSE)
      #Perform 10 fold cross validation using sapply on the function "f_MSPE" above;
      MSPE_CV <- sapply(X=c(1:10),FUN=f_MSPE,lambda1=lambda1,lambda2=lambda2,rand_shuffl=rand_shuffl,folds=folds)
      pred_error_CV <- c(mean(MSPE_CV, na.rm=TRUE),lambda1,lambda2)
      CV_pred_results <- rbind(CV_pred_results,pred_error_CV)
    }
  }
  CV_pred_results <- CV_pred_results[-1,]
  CV_pred_results <- CV_pred_results[order(CV_pred_results[,1]),]
  lambda_1 <- as.numeric(CV_pred_results[1,2])
  lambda_2 <- as.numeric(CV_pred_results[1,3])

  #calculate the parameter estimates betaI_hat and betaN_hat with the new method in the paper
  PI <- XI%*%ginv(t(XI)%*%XI+n*lambda_1*diag(ncol(XI)))%*%t(XI)
  Y_star <- sqrtm(diag(n)-PI)%*%Y
  XN_star <- sqrtm(diag(n)-PI)%*%XN
  fit_lasso_betaN <- glmnet(XN_star, Y_star, alpha = 1, family = "gaussian", intercept = FALSE, standardize = TRUE)
  betaN_hat <- as.vector(coef(fit_lasso_betaN, s=lambda_2)[-1])

  betaI_hat <- ginv(t(XI)%*%XI+n*lambda_1*diag(ncol(XI)))%*%t(XI)%*%(Y-XN%*%betaN_hat)
  betaI_hat <- as.vector(betaI_hat)

  #debiasing the betaI_hat
  eigendecom <- eigen(t(XI)%*%XI+n*lambda_1*diag(ncol(XI)), symmetric=TRUE)
  Q <- eigendecom$vectors
  LAMBDA_eigendecom <- diag(eigendecom$values)
  betaI_hat_debiaedNEW <- Q%*%ginv(diag(ncol(XI))-n*lambda_1*solve(LAMBDA_eigendecom))%*%t(Q)%*%betaI_hat
  betaI_hat_debiaedNEW <- as.vector(betaI_hat_debiaedNEW)

  #Construct confidence intervals for parameters of interest using the new method with debiased estimate;
  alpha <- alpha
  adj_alpha <- alpha/length(betaI_hat_debiaedNEW[Index_XI]) #Bonferroni adjusted multiple CIs
  var_betaI_hat_debiaedNEW <- sigma^2*Q%*%ginv(LAMBDA_eigendecom-n*lambda_1*diag(ncol(XI)))%*%t(Q)
  CI <- cbind(betaI_hat_debiaedNEW-qnorm(1-adj_alpha/2,0,1)*sqrt(diag(var_betaI_hat_debiaedNEW)),betaI_hat_debiaedNEW+qnorm(1-adj_alpha/2,0,1)*sqrt(diag(var_betaI_hat_debiaedNEW)))

  #conduct hypothesis test for betaI;
  adj_alpha <- alpha/length(betaI_hat_debiaedNEW[Index_XI]) #Bonferroni adjusted multiple tests
  test_stat_beta <- numeric(length(betaI_hat_debiaedNEW[Index_XI]))
  testpower <- numeric(length(betaI_hat_debiaedNEW[Index_XI]))
  for(j in 1:length(betaI_hat_debiaedNEW[Index_XI]))
  {
    test_stat_beta[j] <- (betaI_hat_debiaedNEW[j])/sqrt(diag(var_betaI_hat_debiaedNEW)[j])
    testpower[j] <- ifelse(abs(test_stat_beta[j])>=qnorm(1-adj_alpha/2,0,1),1,0)
  }
  return(list(Index_XI=Index_XI,betaI_hat=betaI_hat,betaI_hat_debiaed=betaI_hat_debiaedNEW,betaN_hat=betaN_hat,CI=CI,testpower=testpower))
}





